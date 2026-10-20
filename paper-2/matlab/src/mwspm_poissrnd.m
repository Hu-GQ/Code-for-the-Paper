function k = mwspm_poissrnd(lam, stream)
%MWSPM_POISSNRND  Poisson deviates WITHOUT the Statistics Toolbox.
%
%   k = mwspm_poissrnd(lam)            uses the global random stream
%   k = mwspm_poissrnd(lam, stream)    draws from an explicit RandStream
%
%   The project's MATLAB installation has NO Statistics Toolbox (recorded in
%   every manifest as statistics_toolbox_present: false), so poissrnd() is not
%   available.  This is a self-contained replacement:
%
%     lam <= LAM_SWITCH : Knuth's product-of-uniforms method (exact, subject
%                         to floating point);
%     lam >  LAM_SWITCH : normal approximation round(lam + sqrt(lam)*z).
%
%   The switch exists because Knuth's method costs O(lam) iterations and
%   underflows once exp(-lam) reaches 0 (lam >~ 745).  The caller is expected
%   to keep lam small; every use is therefore accompanied by a counter of how
%   often the approximate branch was taken, and that counter is written into
%   the run output.  The approximation is never used silently.
%
%   Shape follows the input.

LAM_SWITCH = 25.0;

if nargin < 2 || isempty(stream)
    randU = @(n) rand(n, 1);
    randN = @(n) randn(n, 1);
else
    randU = @(n) rand(stream, n, 1);
    randN = @(n) randn(stream, n, 1);
end

sz = size(lam);
lam = double(lam(:));
k = zeros(size(lam));

big = lam > LAM_SWITCH;
if any(big)
    % Declared approximation branch.  lam is large here, so the absolute error
    % of a round()-ed normal is O(1) against a mean of O(25+); it is used only
    % as a safety valve for population sizes the design is not meant to reach.
    k(big) = max(0, round(lam(big) + sqrt(lam(big)) .* randN(nnz(big))));
end

idx = find(~big);
if ~isempty(idx)
    L = exp(-lam(idx));
    p = ones(numel(idx), 1);
    done = false(numel(idx), 1);
    kk = zeros(numel(idx), 1);
    guard = 0;
    while ~all(done) && guard < 2000
        cur = find(~done);
        p(cur) = p(cur) .* randU(numel(cur));
        kk(cur) = kk(cur) + 1;
        done(cur) = p(cur) <= L(cur);
        guard = guard + 1;
    end
    k(idx) = max(kk - 1, 0);
end

k = reshape(k, sz);
end
