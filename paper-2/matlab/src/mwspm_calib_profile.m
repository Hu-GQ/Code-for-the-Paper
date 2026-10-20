function P = mwspm_calib_profile(pre, C, name, grid, opts)
%MWSPM_CALIB_PROFILE  Profile likelihood / profile objective of one parameter.
%
%   P = mwspm_calib_profile(pre, C, 'rho', 0.02:0.01:1.00)
%   P = mwspm_calib_profile(pre, C, 'rho', grid, struct('crit', 3.8415))
%
%   For every value g in `grid` the parameter `name` is HELD FIXED at g and all
%   parameters in C.free are re-optimised from several starts (mwspm_calib_fit).
%   The result is a profile objective, not a slice: a slice would understate the
%   curvature and is exactly how a flat likelihood gets mistaken for a peak.
%
%   P fields
%     name, grid, value, J (objective), dnll = J - min(J),
%     nll, nll_full, d_weeks, theta, pearson, dispersion, peak_week,
%     true_total, reported_total, family, n_eval
%     min_at (value minimising J), crit, ci_low, ci_high, ci_is_connected
%
%   Interval rule: a Wilks (profile-likelihood) interval at level p is
%   { g : 2*(J(g) - min J) <= chi2_1(p) }.  mwspm_chisq1_quantile supplies the
%   threshold without the Statistics Toolbox.

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'crit') || isempty(opts.crit)
    opts.crit = mwspm_chisq1_quantile(0.95);   % 3.8415
end

name = char(name);
grid = double(grid(:));
assert(iscellstr(cellstr(string(C.free))), 'mwspm:calib', 'C.free must be a cellstr');
free = cellstr(string(C.free));
assert(~any(strcmp(free, name)), 'mwspm:calib', ...
    '"%s" is in C.free and cannot also be profiled', name);

n = numel(grid);
rows = repmat(struct('value', NaN, 'J', NaN, 'nll', NaN, 'nll_full', NaN, ...
                     'd_weeks', NaN, 'theta', struct(), 'pearson', NaN, ...
                     'dispersion', NaN, 'peak_week', NaN, 'true_total', NaN, ...
                     'reported_total', NaN, 'family', '', 'n_eval', 0, ...
                     'n_starts_used', 0), n, 1);

warm = {};
if isfield(opts, 'warm') && ~isempty(opts.warm)
    warm = opts.warm;
end

% A PERMANENT start set, kept at every grid point.  Continuation alone is not
% enough: it walks the profile along one local branch, and if the branch at the
% previous grid point is a poor local optimum the next point inherits it.  That
% happened in the first full T12 run: the lag-2 profile reported nll = -2066.80
% at rho = 0.80 while -2067.17 was already known to be achievable at the same
% rho (a 0.36 unit discrepancy measured against the anchor).
% Passing the anchor here makes a MEASURABLE contract hold: no grid point can be
% reported worse than what the anchor theta itself achieves at that value.  The
% start cannot bias the result, because every candidate objective is evaluated
% and the best one wins.
anchor = {};
if isfield(opts, 'anchor') && ~isempty(opts.anchor)
    anchor = opts.anchor;
    if ~iscell(anchor), anchor = {anchor}; end
    anchor = anchor(:)';
end

for i = 1:n
    Cg = C;
    Cg.fixed.(name) = grid(i);
    oi = opts;
    oi.warm = warm;
    [b, info] = mwspm_calib_fit(pre, Cg, oi);
    % continuation: the next grid point starts from this one's optimum, which is
    % legitimate because the profile is a sequence of nearby constrained fits
    warm = [{b.theta}, anchor];
    rows(i).value = grid(i);
    rows(i).J = b.J;
    rows(i).nll = b.nll;
    rows(i).nll_full = b.nll_full;
    rows(i).d_weeks = b.d_weeks;
    rows(i).theta = b.theta;
    rows(i).pearson = b.pearson;
    rows(i).dispersion = b.dispersion;
    rows(i).peak_week = b.peak_week;
    rows(i).true_total = b.true_total;
    rows(i).reported_total = b.reported_total;
    rows(i).family = b.family;
    rows(i).n_eval = info.n_eval;
    rows(i).n_starts_used = info.n_starts;
end

J = [rows.J]';
nll = [rows.nll]';
dnll = nll - min(nll);
[jmin, imin] = min(J);   % jmin -> nll_min, imin -> location of the minimum

P = struct();
P.name = name;
P.grid = grid;
P.J = J;
P.nll = nll;
P.dnll = dnll;
P.min_at = grid(imin);
P.nll_min = jmin;
P.crit = opts.crit;
P.rows = rows;
[P.ci_low, P.ci_high, P.ci_is_connected] = wilksBand(grid, dnll, imin, opts.crit);
P.family = C.family;
P.free = free;
P.fixed = C.fixed;
end

% ---------------------------------------------------------------------------
function [lo, hi, connected] = wilksBand(grid, dnll, imin, crit)
% Largest contiguous run around the optimum with nll - min nll <= crit/2.
% A band is reported as NOT connected if any grid point between its ends lies
% outside the band - in that case the point estimate is not even a local peak of
% the profile, which is a finding in itself and must not be smoothed away.
half = crit / 2;
ok = dnll <= half;
n = numel(grid);
k = imin;
while k > 1 && ok(k - 1)
    k = k - 1;
end
lo = grid(k);
k = imin;
while k < n && ok(k + 1)
    k = k + 1;
end
hi = grid(k);
connected = all(ok(grid >= lo & grid <= hi));
end
