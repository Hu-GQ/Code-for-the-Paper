function mwspm_writejson(path, s)
%MWSPM_WRITEJSON  Write a struct to a pretty-printed UTF-8 JSON file.
%
%   NaN / Inf are encoded as JSON null (jsonencode default).  Deliberately no
%   custom injection: the file is exactly what jsondecode(jsonencode(x))
%   reconstructs.

txt = jsonencode(s, 'PrettyPrint', true);
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'mwspm:writejson', 'cannot write %s', path);
fwrite(fid, unicode2native(txt, 'UTF-8'), 'uint8');
fclose(fid);
end
