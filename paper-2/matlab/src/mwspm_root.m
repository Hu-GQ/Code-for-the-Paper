function root = mwspm_root()
%MWSPM_ROOT  Locate the project root.
%
%   root = mwspm_root() returns the absolute path of the directory that
%   contains config/project.json.  The search starts at the directory of
%   THIS file (<root>/src/matlab) and walks upwards, so the result does NOT
%   depend on the current working directory, on the MATLAB base workspace or
%   on any GUI state (the project record Sec.4.1/4.2).
%
%   The answer is cached for the duration of the session because the file
%   location cannot change while MATLAB is running.

persistent cachedRoot

if ~isempty(cachedRoot) && exist(fullfile(cachedRoot, 'config', 'project.json'), 'file') == 2
    root = cachedRoot;
    return
end

anchor = fileparts(mfilename('fullpath'));
cand   = anchor;

for k = 1:10
    if exist(fullfile(cand, 'config', 'project.json'), 'file') == 2
        cachedRoot = cand;
        root = cand;
        return
    end
    parent = fileparts(cand);
    if strcmp(parent, cand)          % reached the filesystem root
        break
    end
    cand = parent;
end

error('mwspm:rootNotFound', ...
    'could not locate the project root (a directory containing config/project.json) above %s', anchor);
end
