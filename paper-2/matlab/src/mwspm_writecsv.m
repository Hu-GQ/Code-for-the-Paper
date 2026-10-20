function mwspm_writecsv(path, header, data)
%MWSPM_WRITECSV  Write a numeric matrix to CSV with an explicit header.
%
%   mwspm_writecsv(path, headerCellstr, dataMatrix)
%
%   The header is written verbatim (no BOM) so that downstream tooling sees
%   plain ASCII/UTF-8; the project's own readers tolerate a BOM but never
%   require one.  Column count is validated to prevent silently misaligned
%   tables - a class of error that is very hard to see afterwards.
%
%   PRECISION (P2, fixed 2026-10-07)
%   -------------------------------
%   Values used to be written with '%.10g'.  Measured consequence: the T12.01
%   weekly_fit.csv printed weeks 36 and 37 as BITWISE EQUAL counts, which reads
%   like a real plateau but was a rounding pseudo-tie.  Values are now written
%   with '%.17g', which round-trips any IEEE double exactly, so a printed tie is
%   a real tie at double precision.  This changes the TEXT of every CSV this
%   function produces; runs archived before 2026-10-07 keep their own frozen
%   files and must not be compared byte-wise with newer ones.  The VALUES are
%   unchanged.

assert(iscell(header) || isstring(header), 'mwspm:writecsv', 'header must be cellstr');
header = cellstr(header);
data = double(data);
assert(size(data, 2) == numel(header), 'mwspm:writecsv', ...
    'header has %d columns but data has %d', numel(header), size(data, 2));

fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'mwspm:writecsv', 'cannot write %s', path);

fprintf(fid, '%s', header{1});
for k = 2:numel(header)
    fprintf(fid, ',%s', header{k});
end
fprintf(fid, '\n');

for i = 1:size(data, 1)
    fprintf(fid, '%.17g', data(i, 1));
    for k = 2:size(data, 2)
        fprintf(fid, ',%.17g', data(i, k));
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
