function out = reproduce(expId, varargin)
%REPRODUCE  Unified MATLAB entry point for every numerical experiment.
%
%   reproduce()                          list the configured experiments
%   reproduce('list')                    same
%   reproduce(expId)                     run experiment <expId>
%   reproduce(expId, 'dry_run', true)    validate config/inputs only, no numerics
%   reproduce(expId, 'require_ready', false)
%                                        allow a TODO/not enabled config to run
%                                        (only for an explicitly authorised run)
%
%   Design constraints (the project record Sec.4):
%     * paths are resolved from the PROJECT LOCATION, never from pwd;
%     * no dependence on the MATLAB base workspace or on any GUI action;
%     * every run gets its own results/runs/<run_id>/ containing config.json,
%       manifest.json, logs/ and data/, plus input and code hashes;
%     * a failed run is recorded as FAILED and never leaves a success flag.
%
%   The experiment logic itself lives in src/matlab/exp_*.m, and the model
%   equations in src/matlab/mwspm_rhs.m - this file is orchestration only.

p = inputParser;
p.addOptional('expId', 'list', @(x) ischar(x) || isstring(x));
p.addParameter('dry_run', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('require_ready', true, @(x) islogical(x) || isnumeric(x));
p.parse(expId, varargin{:});
opt = p.Results;
expId = char(opt.expId);

root = reproduce_locate_root();
addpath(fullfile(root, 'src', 'matlab'), '-begin');
rehash; %#ok<REHASH>

cfgDir = fullfile(root, 'config', 'experiments');
assert(exist(cfgDir, 'dir') == 7, 'reproduce', 'missing config directory %s', cfgDir);

if any(strcmpi(expId, {'list', ''}))
    out = reproduce_list(cfgDir);
    return
end

cfgPath = fullfile(cfgDir, [expId '.json']);
assert(exist(cfgPath, 'file') == 2, 'reproduce', ...
    'no such experiment "%s" (config file %s does not exist)', expId, cfgPath);

fprintf('reproduce: project root = %s\n', root);
fprintf('reproduce: experiment  = %s\n', expId);

ctx = mwspm_run_init(root, expId, 'require_ready', opt.require_ready);
if logical(opt.dry_run)
    ctx.dryRun = true;
end
fprintf('reproduce: run_id      = %s\n', ctx.runId);
fprintf('reproduce: run_dir     = %s\n', ctx.runDir);

entryFcn = char(string(getdef(ctx.cfg, 'entry_function', expId)));
assert(exist(entryFcn, 'file') == 2, 'reproduce:noEntryFcn', ...
    ['entry function "%s" not found on the MATLAB path.\n' ...
     'Check config/experiments/%s.json: implementation = %s.\n' ...
     'A config with implementation=TODO has a frozen specification but no code yet - ' ...
     'the experiment must be implemented (see config/experiments/README.md) before it can run.'], ...
    entryFcn, expId, char(string(getdef(ctx.cfg, 'implementation', 'unknown'))));

payload = struct();
errInfo = [];
status = 'SUCCESS';
outputs = {};
try
    fn = str2func(entryFcn);
    payload = fn(ctx);
    if isfield(payload, 'outputs') && ~isempty(payload.outputs)
        outputs = cellstr(string(payload.outputs));
    end
catch ME
    status = 'FAILED';
    errInfo = ME;
end

man = mwspm_run_finalize(ctx, status, payload, errInfo, outputs);

fprintf('reproduce: status      = %s\n', man.status);
fprintf('reproduce: manifest    = %s\n', fullfile(ctx.runDir, 'manifest.json'));
fprintf('reproduce: log         = %s\n', ctx.logPath);

if strcmp(status, 'FAILED')
    rethrow(errInfo);
end

out = struct();
out.run_id = ctx.runId;
out.run_dir = ctx.runDir;
out.status = status;
out.manifest = man;
out.dry_run = logical(opt.dry_run);
end

% ---------------------------------------------------------------------------
function root = reproduce_locate_root()
% Locate the project root without depending on pwd: this file lives in
% <root>/scripts, so walk upwards until a directory containing
% config/project.json is found.
anchor = fileparts(mfilename('fullpath'));
cand = anchor;
for k = 1:10
    if exist(fullfile(cand, 'config', 'project.json'), 'file') == 2
        root = cand;
        return
    end
    parent = fileparts(cand);
    if strcmp(parent, cand)
        break
    end
    cand = parent;
end
error('reproduce:rootNotFound', ...
    'could not locate the project root above %s', anchor);
end

% ---------------------------------------------------------------------------
function out = reproduce_list(cfgDir)
d = dir(fullfile(cfgDir, '*.json'));
rows = {};
for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);
    c = mwspm_readjson(f);
    if ~isfield(c, 'experiment_id')
        continue
    end
    st = 'READY';
    if isfield(c, 'status'), st = upper(char(c.status)); end
    impl = 'TODO';
    if isfield(c, 'implementation'), impl = upper(char(c.implementation)); end
    claims = '';
    if isfield(c, 'claims') && ~isempty(c.claims)
        claims = strjoin(cellstr(string(c.claims)), ',');
    end
    rows(end+1, :) = {c.experiment_id, char(string(getdef(c, 'task', '-'))), st, impl, ...
        char(string(getdef(c, 'entry_function', '-'))), claims}; %#ok<AGROW>
end

fprintf('\nconfigured experiments in %s\n', cfgDir);
fprintf('%-26s %-10s %-9s %-13s %-26s %s\n', 'experiment_id', 'task', 'status', 'implementation', 'entry_function', 'claims');
for k = 1:size(rows, 1)
    fprintf('%-26s %-10s %-9s %-13s %-26s %s\n', rows{k, :});
end
fprintf(['\nstatus         = authorisation/gate  (READY runs; TODO/not enabled need ''require_ready'',false)\n' ...
         'implementation = is the entry function written yet\n' ...
         '=> READY + TODO means "cleared methodologically, code not written yet".\n']);
fprintf('\nrun one with  reproduce(''<experiment_id>'')  or validate with  reproduce(''<experiment_id>'',''dry_run'',true)\n');

out = struct();
out.experiments = rows;
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
