function out = pp_util(action, varargin)
%PP_UTIL  Shared helpers for the T17 paper-figure scripts.
%
%   pp_util('root')                    project root (char)
%   pp_util('dirs')                    struct: main/supp/preview/source_data paths
%   pp_util('runDir', runId)           absolute path of results/runs/<runId>
%   pp_util('snap', files)             read-only snapshot: {files, sha256} (cellstr)
%   pp_util('assert', snap)            re-hash; error if any file changed
%   pp_util('col', S, name)            numeric column from mwspm_readcsv struct (NaN for '', NA)
%   pp_util('style')                   apply unified figure/axes style (call before plotting)
%   pp_util('palette')                 struct of RGB triplets used across figures
%   pp_util('savefig', fig, figId)     export preview PNG (300 dpi) + vector PDF
%
%   Read-only contract: plotting scripts read ONLY saved files under
%   results/runs/<runId>/data/ (or explicitly annotated historical saves).
%   After reading, pp_util('assert', snap) proves the files were not modified.

switch lower(action)
    case 'root'
        % <root>/scripts/plots/paper/pp_util.m  ->  up 4 levels
        out = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));

    case 'dirs'
        r = pp_util('root');
        out = struct( ...
            'main',   fullfile(r, 'figures', 'main'), ...
            'supp',   fullfile(r, 'figures', 'supp'), ...
            'preview',fullfile(r, 'figures', 'preview'), ...
            'src',    fullfile(r, 'figures', 'source_data'));

    case 'rundir'
        out = fullfile(pp_util('root'), 'results', 'runs', char(varargin{1}));

    case 'snap'
        files = varargin{1};
        if ischar(files) || isstring(files); files = cellstr(files); end
        h = cell(size(files));
        for k = 1:numel(files)
            h{k} = mwspm_sha256(files{k});
        end
        out = struct('files', {files}, 'sha256', {h});

    case 'assert'
        snap = varargin{1};
        for k = 1:numel(snap.files)
            h2 = mwspm_sha256(snap.files{k});
            assert(strcmp(h2, snap.sha256{k}), 'pp:notReadOnly', ...
                'file "%s" changed during plotting - read-only contract violated', snap.files{k});
        end
        fprintf('pp_util: read-only check passed (%d files unchanged)\n', numel(snap.files));

    case 'col'
        S = varargin{1}; name = varargin{2};
        idx = find(strcmp(S.names, name), 1);
        assert(~isempty(idx), 'pp:col', 'column "%s" not found', name);
        v = NaN(S.n_rows, 1);
        for k = 1:S.n_rows
            c = strtrim(S.rows{k}{idx});
            if isempty(c) || strcmpi(c, 'na') || strcmpi(c, 'nan'); continue; end
            v(k) = str2double(c);
        end
        out = v;

    case 'strcol'
        S = varargin{1}; name = varargin{2};
        idx = find(strcmp(S.names, name), 1);
        assert(~isempty(idx), 'pp:col', 'column "%s" not found', name);
        out = cell(S.n_rows, 1);
        for k = 1:S.n_rows
            out{k} = strtrim(S.rows{k}{idx});
        end

    case 'style'
        set(groot, 'defaultFigureColor', 'w');
        set(groot, 'defaultAxesFontName', 'Helvetica');
        set(groot, 'defaultTextFontName', 'Helvetica');
        set(groot, 'defaultAxesFontSize', 8);
        set(groot, 'defaultAxesLineWidth', 0.75);
        set(groot, 'defaultAxesTickDir', 'out');
        set(groot, 'defaultAxesBox', 'on');
        set(groot, 'defaultLegendBox', 'off');
        set(groot, 'defaultAxesColorOrder', pp_util('palette').order);
        out = true;

    case 'palette'
        % restrained family: neutrals + one signal + one accent
        out = struct();
        out.blue   = [0.121 0.310 0.635];   % signal (main series / T0 caliber)
        out.red    = [0.700 0.150 0.150];   % accent (danger / threshold caliber)
        out.gray   = [0.45 0.45 0.45];
        out.lgray  = [0.80 0.80 0.80];
        out.blue_l = [0.62 0.72 0.86];      % light signal
        out.red_l  = [0.88 0.70 0.70];      % light accent
        out.green  = [0.18 0.50 0.30];
        out.orange = [0.90 0.60 0.20];
        out.order  = [out.blue; out.red; out.orange; out.green; out.gray];

    case 'savefig'
        fig = varargin{1}; figId = varargin{2};
        d = pp_util('dirs');
        % strip interactive axes toolbars so exports stay clean
        axs = findall(fig, 'Type', 'axes');
        for k = 1:numel(axs)
            try
                axs(k).Toolbar = [];
            catch
            end
        end
        pngPath = fullfile(d.preview, [figId '.png']);
        pdfPath = fullfile(d.main,   [figId '.pdf']);
        exportgraphics(fig, pngPath, 'Resolution', 300);
        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        fprintf('pp_util: wrote %s\n           %s\n', pngPath, pdfPath);
        out = struct('png', pngPath, 'pdf', pdfPath);

    case 'savefig_supp'
        fig = varargin{1}; figId = varargin{2};
        d = pp_util('dirs');
        axs = findall(fig, 'Type', 'axes');
        for k = 1:numel(axs)
            try
                axs(k).Toolbar = [];
            catch
            end
        end
        pngPath = fullfile(d.preview, [figId '.png']);
        pdfPath = fullfile(d.supp,    [figId '.pdf']);
        exportgraphics(fig, pngPath, 'Resolution', 300);
        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        fprintf('pp_util: wrote %s\n           %s\n', pngPath, pdfPath);
        out = struct('png', pngPath, 'pdf', pdfPath);

    otherwise
        error('pp_util:action', 'unknown action "%s"', action);
end
end
