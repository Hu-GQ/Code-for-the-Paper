function s = mwspm_readjson(path)
%MWSPM_READJSON  Read a UTF-8 JSON file (BOM tolerant) into a struct.
%
%   Two properties matter here and both are easy to get wrong:
%     1. The bytes are decoded as UTF-8 explicitly.  The project's JSON files
%        contain Chinese annotations (parameter provenance notes), so decoding
%        through the platform code page would silently corrupt them.
%     2. The decoded text must be a ROW character vector.  native2unicode()
%        preserves the orientation of its input, so feeding it a column of
%        bytes yields a column of characters, and jsondecode() rejects that
%        with "JSON text must be a character vector" - a message that points
%        nowhere near the real cause.  Hence the explicit reshape.
%
%   mwspm_writejson() is the exact inverse.

fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid > 0, 'mwspm:readjson', 'cannot open %s', path);
raw = fread(fid, Inf, '*uint8');
fclose(fid);

% Orientation first: a row of bytes decodes to a row of characters.
raw = raw(:).';

% Strip a UTF-8 byte order mark if present (some Windows editors add one).
if numel(raw) >= 3 && all(raw(1:3) == uint8([239 187 191]))
    raw = raw(4:end);
end

txt = native2unicode(raw, 'UTF-8');
txt = txt(:).';
txt(txt == char(0)) = [];

assert(~isempty(txt), 'mwspm:readjson', 'file is empty: %s', path);
s = jsondecode(txt);
end
