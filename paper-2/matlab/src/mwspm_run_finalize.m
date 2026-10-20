function man = mwspm_run_finalize(ctx, status, payload, errInfo, outputs)
%MWSPM_RUN_FINALIZE  Write manifest.json, close the log, return the manifest.
%
%   man = mwspm_run_finalize(ctx, 'SUCCESS', payload, [], outputs)
%   man = mwspm_run_finalize(ctx, 'FAILED',  payload, errInfo, outputs)
%
%   Protocol requirements implemented here (PROJECT_RULES.md Sec.4.4/4.5):
%     * status is ONLY 'SUCCESS' or 'FAILED'; a failed run NEVER leaves a
%       success flag, and the error message plus stack are recorded;
%     * input hashes (config, model card, data files) and code hashes
%       (every .m under src/matlab and scripts) are recorded;
%     * the environment snapshot and the output list (with hashes) are recorded;
%     * RNG information is recorded for stochastic experiments and is absent
%       for deterministic ones.

if nargin < 3 || isempty(payload)
    payload = struct();
end
if nargin < 4
    errInfo = [];
end
if nargin < 5 || isempty(outputs)
    outputs = {};
end

status = upper(char(status));
assert(any(strcmp(status, {'SUCCESS', 'FAILED'})), 'mwspm:run', ...
    'manifest status must be SUCCESS or FAILED, got "%s"', status);

man = struct();
man.status        = status;
man.run_id        = ctx.runId;
man.experiment_id = ctx.expId;
man.task          = char(string(getdef(ctx.cfg, 'task', '-')));
man.title         = char(string(getdef(ctx.cfg, 'title', '')));
man.claims        = cellstr(string(getdef(ctx.cfg, 'claims', {})));
man.model_id      = ctx.modelCard.model_id;
man.model_version = ctx.modelCard.model_version;
man.spec_version  = ctx.modelCard.spec_version;
man.started_at    = ctx.startedAt;
man.finished_at   = datestr(now, 'yyyy-mm-ddTHH:MM:SS');      %#ok<TNOW1,DATST>
if ~isempty(ctx.tStart)
    man.duration_seconds = toc(ctx.tStart);
else
    man.duration_seconds = NaN;
end

man.environment = mwspm_env();

% ---- hashes --------------------------------------------------------------
inFiles = { ...
    ctx.cfgPath, ...
    ctx.modelCardPath, ...
    fullfile(ctx.root, ctx.modelCard.forcing.T_daily.file), ...
    fullfile(ctx.root, ctx.modelCard.forcing.Lambda_daily.file), ...
    fullfile(ctx.root, ctx.modelCard.forcing.observations_weekly.file), ...
    fullfile(ctx.root, 'data/raw/upstream_cleaned/guangzhou_population.csv')};
man.input_hashes = hashList(ctx.root, inFiles);

codeFiles = [listFiles(fullfile(ctx.root, 'src', 'matlab'), '*.m'); ...
             listFiles(fullfile(ctx.root, 'scripts'), '*.m')];
man.code_hashes = hashList(ctx.root, codeFiles);

% ---- outputs -------------------------------------------------------------
outRec = cell(numel(outputs), 1);
for k = 1:numel(outputs)
    f = outputs{k};
    if ~isAbsolute(f)
        f = fullfile(ctx.dataDir, f);
    end
    e = struct('file', '', 'exists', false, 'bytes', NaN, 'sha256', '');
    if exist(f, 'file') == 2
        d = dir(f);
        e.file = relPath(ctx.root, f);
        e.exists = true;
        e.bytes = d.bytes;
        e.sha256 = mwspm_sha256(f);
    else
        e.file = relPath(ctx.root, f);
        e.exists = false;
    end
    outRec{k} = e;
end
man.outputs = outRec;

% ---- payload (solver settings, termination reasons, metrics, rng, ...) ----
man.result = payload;

if ~isempty(errInfo)
    man.error = struct( ...
        'identifier', errInfo.identifier, ...
        'message', errInfo.message, ...
        'stack', errInfo.stack);
    try
        fprintf(ctx.logFid, '\n[FAILED] %s : %s\n', errInfo.identifier, errInfo.message);
        for k = 1:numel(errInfo.stack)
            fprintf(ctx.logFid, '  at %s (line %d)\n', errInfo.stack(k).name, errInfo.stack(k).line);
        end
    catch
        % logging must never mask the original failure
    end
end

mwspm_writejson(fullfile(ctx.runDir, 'manifest.json'), man);

try
    fprintf(ctx.logFid, '\nstatus: %s\nfinished_at: %s\nduration_seconds: %.3f\n', ...
        status, man.finished_at, man.duration_seconds);
    fclose(ctx.logFid);
catch
end
end

% ---------------------------------------------------------------------------
function L = listFiles(dirPath, pattern)
L = {};
if exist(dirPath, 'dir') ~= 7
    return
end
d = dir(fullfile(dirPath, '**', pattern));
for k = 1:numel(d)
    L{end+1, 1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
end
L = sort(L);
end

% ---------------------------------------------------------------------------
function R = hashList(root, files)
R = cell(numel(files), 1);
for k = 1:numel(files)
    f = files{k};
    e = struct('file', relPath(root, f), 'exists', false, 'sha256', '');
    if exist(f, 'file') == 2
        e.exists = true;
        e.sha256 = mwspm_sha256(f);
    end
    R{k} = e;
end
end

% ---------------------------------------------------------------------------
function p = relPath(root, f)
p = char(f);
rootSlash = root;
if ~endsWith(rootSlash, filesep)
    rootSlash = [rootSlash filesep];
end
if startsWith(p, rootSlash)
    p = p(numel(rootSlash) + 1:end);
end
p = strrep(p, filesep, '/');
end

function tf = isAbsolute(p)
tf = ~isempty(regexp(char(p), '^([A-Za-z]:[\\/]|\\\\|/)', 'once'));
end

function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
