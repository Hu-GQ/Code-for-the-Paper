function ctx = mwspm_run_init(root, expId, varargin)
%MWSPM_RUN_INIT  Create a unique run directory and freeze the run config.
%
%   ctx = mwspm_run_init(root, expId)
%   ctx = mwspm_run_init(root, expId, 'require_ready', false)
%
%   run_id = <YYYYmmddTHHMMSS local, UTC+8>_<sha256(config file) first 8 chars>
%   results/runs/<run_id>/{config.json, logs/, data/}
%
%   Guarantees required by the project record Sec.4:
%     * the directory is unique and never overwritten (a counter suffix is
%       added if a collision occurs within the same second);
%     * the config actually used is copied into the run directory;
%     * the experiment config file itself is NEVER modified (read-only input);
%     * a run whose config has status TODO/not enabled is refused unless the
%       caller explicitly opts in with 'require_ready',false.

p = inputParser;
p.addRequired('root', @(x) ischar(x) || isstring(x));
p.addRequired('expId', @(x) ischar(x) || isstring(x));
p.addParameter('require_ready', true, @(x) islogical(x) || isnumeric(x));
p.parse(root, expId, varargin{:});
opt = p.Results;

root = char(root);
expId = char(expId);

cfgPath = fullfile(root, 'config', 'experiments', [expId '.json']);
assert(exist(cfgPath, 'file') == 2, 'mwspm:run', 'experiment config not found: %s', cfgPath);

cfg = mwspm_readjson(cfgPath);
assert(isfield(cfg, 'experiment_id') && strcmp(cfg.experiment_id, expId), 'mwspm:run', ...
    'experiment_id in %s does not match its file name', cfgPath);

status = 'READY';
if isfield(cfg, 'status'), status = upper(char(cfg.status)); end
if logical(opt.require_ready) && ~strcmp(status, 'READY')
    blocked = '';
    if isfield(cfg, 'blocked_by') && ~isempty(cfg.blocked_by)
        blocked = [' blocked_by: ' char(string(cfg.blocked_by))];
    end
    error('mwspm:run:notReady', ...
        ['experiment "%s" has status %s and may not be executed yet.%s\n' ...
         'Pass ''require_ready'',false only for an explicitly authorised run.'], expId, status, blocked);
end

% Orthogonal to status: is the entry FUNCTION written at all?  Caught here so
% that an unimplemented experiment fails BEFORE a run directory is created
% (require_ready,false authorises the run, it does not implement it).
implementation = 'TODO';
if isfield(cfg, 'implementation'), implementation = upper(char(cfg.implementation)); end
if strcmp(implementation, 'TODO') && exist(char(string(getdef(cfg, 'entry_function', expId))), 'file') ~= 2
    error('mwspm:run:notImplemented', ...
        ['experiment "%s" is not implemented: implementation=%s and entry function "%s" ' ...
         'does not exist on the MATLAB path.\n' ...
         'Implement src/matlab/%s.m first (config/experiments/README.md, section 0).'], ...
        expId, implementation, char(string(getdef(cfg, 'entry_function', expId))), ...
        char(string(getdef(cfg, 'entry_function', expId))));
end

cfgHash = mwspm_sha256(cfgPath);

modelCardPath = fullfile(root, 'config', 'experiments', [char(cfg.model) '.json']);
assert(exist(modelCardPath, 'file') == 2, 'mwspm:run', 'model card not found: %s', modelCardPath);
modelCard = mwspm_readjson(modelCardPath);

assert(strcmp(modelCard.model_version, cfg.model_version), 'mwspm:run', ...
    'model version mismatch: card "%s" vs experiment "%s"', modelCard.model_version, cfg.model_version);

stamp = datestr(now, 'yyyymmddTHHMMSS');                     %#ok<TNOW1,DATST>
baseId = sprintf('%s_%s', stamp, cfgHash(1:8));

runId = baseId;
k = 2;
while exist(fullfile(root, 'results', 'runs', runId), 'dir') == 7
    runId = sprintf('%s_%d', baseId, k);
    k = k + 1;
end
runDir = fullfile(root, 'results', 'runs', runId);
mkdir(runDir);
mkdir(fullfile(runDir, 'logs'));
mkdir(fullfile(runDir, 'data'));

% resolved config actually used: experiment config + model card pointer
resolved = struct();
resolved.run_id            = runId;
resolved.experiment_id     = expId;
resolved.config_file       = sprintf('config/experiments/%s.json', expId);
resolved.config_sha256     = cfgHash;
resolved.model_card_file   = sprintf('config/experiments/%s.json', char(cfg.model));
resolved.model_card_sha256 = mwspm_sha256(modelCardPath);
resolved.project_root      = root;
resolved.experiment        = cfg;
resolved.model_card        = modelCard;
mwspm_writejson(fullfile(runDir, 'config.json'), resolved);

logPath = fullfile(runDir, 'logs', 'run_log.txt');
fid = fopen(logPath, 'w', 'n', 'UTF-8');
assert(fid > 0, 'mwspm:run', 'cannot open log file %s', logPath);

ctx = struct();
ctx.root      = root;
ctx.expId     = expId;
ctx.cfgPath   = cfgPath;
ctx.cfg       = cfg;
ctx.cfgHash   = cfgHash;
ctx.modelCard = modelCard;
ctx.modelCardPath = modelCardPath;
ctx.runId     = runId;
ctx.runDir    = runDir;
ctx.dataDir   = fullfile(runDir, 'data');
ctx.logPath   = logPath;
ctx.logFid    = fid;
ctx.startedAt = datestr(now, 'yyyy-mm-ddTHH:MM:SS');          %#ok<TNOW1,DATST>
ctx.tStart    = tic;

fprintf(fid, 'run_id      : %s\n', runId);
fprintf(fid, 'experiment  : %s\n', expId);
fprintf(fid, 'task        : %s\n', char(string(getdef(cfg, 'task', '-'))));
fprintf(fid, 'status      : %s\n', status);
fprintf(fid, 'implemented : %s\n', implementation);
fprintf(fid, 'config      : %s  (sha256 %s)\n', ctx.cfgPath, cfgHash);
fprintf(fid, 'model card  : %s  (sha256 %s)\n', modelCardPath, resolved.model_card_sha256);
fprintf(fid, 'claims      : %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));
fprintf(fid, '\n');
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
