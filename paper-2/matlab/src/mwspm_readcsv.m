function S = mwspm_readcsv(path)
%MWSPM_READCSV  Minimal quote-aware CSV reader (UTF-8, BOM tolerant).
%
%   S = mwspm_readcsv(path) returns a struct with
%       S.names : 1xM cellstr of the header names (trimmed, quotes removed)
%       S.rows  : (N-1)x1 cell of 1xM cellstr data rows
%       S.n_rows: number of data rows
%
%   Written in-house on purpose: readtable() would have to be trusted to
%   guess delimiter/encoding/missing rules, and the project rule is to
%   define the import contract explicitly.  Fields are returned as text;
%   numeric conversion happens in the caller so that missing values stay
%   distinguishable from zeros.

fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid > 0, 'mwspm:readcsv', 'cannot open %s', path);
raw = fread(fid, Inf, '*uint8');
fclose(fid);

% Orientation: a row of bytes decodes to a row of characters.  native2unicode
% preserves the input shape, and passing a column of characters to strrep()
% later raises "arguments must be character vectors" - see mwspm_readjson.
raw = raw(:).';

% strip an optional UTF-8 byte-order mark
if numel(raw) >= 3 && all(raw(1:3) == uint8([239 187 191]))
    raw = raw(4:end);
end

txt = native2unicode(raw, 'UTF-8');
txt = txt(:).';
txt = strrep(txt, char(13), '');
lines = strsplit(txt, char(10), 'CollapseDelimiters', false);
while ~isempty(lines) && isempty(strtrim(lines{end}))
    lines(end) = [];                                   %#ok<AGROW>
end
assert(~isempty(lines), 'mwspm:readcsv', 'empty file %s', path);

parsed = cell(numel(lines), 1);
for i = 1:numel(lines)
    parsed{i} = mwspm_parse_csv_line(lines{i});
end

width = max(cellfun(@numel, parsed));
for i = 1:numel(parsed)
    if numel(parsed{i}) < width
        parsed{i}(end+1:width) = {''};                 %#ok<AGROW>
    end
end

names = strtrim(regexprep(parsed{1}, '^"(.*)"$', '$1'));
assert(numel(unique(names)) == numel(names), ...
    'mwspm:readcsv', 'duplicate column names in %s', path);

S = struct('names', {names}, 'rows', {parsed(2:end)}, 'n_rows', numel(lines)-1, 'path', path);
end

% ---------------------------------------------------------------------------
function p = mwspm_parse_csv_line(ln)
% One-pass RFC4180-ish parser: "..." quoting, "" escape, commas inside quotes.
p   = {};
cur = '';
inq = false;
k   = 1;
while k <= numel(ln)
    ch = ln(k);
    if inq
        if ch == '"'
            if k < numel(ln) && ln(k+1) == '"'
                cur(end+1) = '"';                      %#ok<AGROW>
                k = k + 2;
                continue
            end
            inq = false;
        else
            cur(end+1) = ch;                           %#ok<AGROW>
        end
    else
        if ch == '"'
            inq = true;
        elseif ch == ','
            p{end+1} = cur;                            %#ok<AGROW>
            cur = '';
        else
            cur(end+1) = ch;                           %#ok<AGROW>
        end
    end
    k = k + 1;
end
p{end+1} = cur;
end
