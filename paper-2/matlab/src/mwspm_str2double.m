function x = mwspm_str2double(c)
%MWSPM_STR2DOUBLE  Text -> double with explicit missing-value rules.
%
%   Blank, "NA", "NaN", "null" and "." are all treated as MISSING and become
%   NaN (NOT zero).  Callers decide what a missing value means; the loader
%   asserts that the analysis window contains no missing forcing values.

if iscell(c)
    x = cellfun(@mwspm_str2double, c);
    return
end
if isstring(c)
    c = char(c);
end
if isnumeric(c)
    x = double(c);
    return
end

s = strtrim(char(c));
if isempty(s) || any(strcmpi(s, {'NA', 'NaN', 'null', 'NULL', '.', 'n/a'}))
    x = NaN;
    return
end
x = str2double(s);
end
