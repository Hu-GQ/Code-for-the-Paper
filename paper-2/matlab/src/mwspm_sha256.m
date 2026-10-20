function h = mwspm_sha256(x)
%MWSPM_SHA256  SHA-256 hex digest of bytes, text or a file.
%
%   h = mwspm_sha256(bytes)   bytes: numeric/char array treated as raw bytes
%   h = mwspm_sha256(path)    existing file -> raw bytes of the file
%   h = mwspm_sha256('text')  string that is not a file path -> UTF-8 bytes
%
%   Used for the provenance hashes required by the project record
%   Sec.4.4 (input hashes, code hashes, run config hash).
%
%   Implementation note: the digest is computed with the JVM
%   (java.security.MessageDigest).  On this installation usejava('jvm') is 1
%   and MessageDigest('SHA-256') is available even though version('-java')
%   reports "Java is not enabled" under -batch.  If the JVM is ever
%   unavailable this function errors out loudly instead of silently
%   degrading to a weaker hash.

persistent md

if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
    if exist(s, 'file') == 2
        fid = fopen(s, 'r', 'n', 'UTF-8');
        assert(fid > 0, 'mwspm:sha256', 'cannot open %s', s);
        b = fread(fid, Inf, '*uint8');
        fclose(fid);
    else
        b = uint8(unicode2native(s, 'UTF-8'));
    end
else
    b = uint8(x(:));
end

assert(usejava('jvm'), 'mwspm:sha256', ...
    'JVM is not available; SHA-256 provenance hashes cannot be computed');

if isempty(md)
    md = java.security.MessageDigest.getInstance('SHA-256');
end
md.reset();
d = md.digest(b);
h = lower(reshape(dec2hex(typecast(d, 'uint8'))', 1, []));
end
