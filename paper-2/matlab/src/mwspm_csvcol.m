function v = mwspm_csvcol(S, colname, kind)
%MWSPM_CSVCOL  Extract one column of a struct produced by mwspm_readcsv.
%
%   v = mwspm_csvcol(S, colname)            -> cellstr of the column
%   v = mwspm_csvcol(S, colname, 'numeric') -> double column, blanks -> NaN
%
%   The column is looked up by exact (trimmed) header name; a missing column
%   is an error rather than a silent empty result.

if nargin < 3
    kind = 'text';
end

idx = find(strcmp(strtrim(S.names), colname), 1);
assert(~isempty(idx), 'mwspm:csvcol', 'column "%s" not found in %s', colname, S.path);

v = cellfun(@(r) r{idx}, S.rows, 'UniformOutput', false);

if strcmpi(kind, 'numeric')
    v = cellfun(@mwspm_str2double, v);
elseif ~strcmpi(kind, 'text')
    error('mwspm:csvcol', 'unknown conversion kind "%s"', kind);
end
end
