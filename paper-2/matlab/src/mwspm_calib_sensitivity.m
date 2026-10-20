function S = mwspm_calib_sensitivity(pre, th, names, opts)
%MWSPM_CALIB_SENSITIVITY  Structural and practical identifiability diagnostics
%   of the calibration mean model.
%
%   S = mwspm_calib_sensitivity(pre, th, names)
%   S = mwspm_calib_sensitivity(pre, th, names, struct('rel_step', 1e-4))
%
%   The mean model is  mu = rho * dC_{w-d}(theta).  Its local behaviour is
%   described by the sensitivity (Jacobian) matrix
%       S_{wj} = d mu_w / d log(theta_j)  =  theta_j * d mu_w / d theta_j
%   evaluated by central differences.  Working in log-units makes every column
%   "counts per 100% change", so columns are comparable and the diagnostics do
%   not simply rank parameters by their physical size.
%
%   Diagnostics returned
%     S_raw            (47 x p) sensitivity in log-scale units
%     colnorm          ||column||
%     Sh               column-normalised sensitivity (unit norm columns)
%     svd              singular values of Sh (descending)
%     cond_number      sigma_max / sigma_min
%     cond_indices     Belsley condition indices sigma_max / sigma_i
%     var_prop         Belsley variance-decomposition proportions (p x p)
%     collinear_sets   groups of parameters with var_prop > 0.5 at some
%                      condition index above 30 (classic collinearity detector)
%     corr_estimates   correlation of the least-squares ESTIMATES from
%                      (Sh'Sh)^{-1}, normalised to a correlation matrix
%     practical_rank   number of singular values above `tol`
%     tol, names, rel_step, th
%
%   Structural vs practical: a parameter is STRUCTURALLY non-identifiable here
%   if its column is exactly (numerically) a linear combination of the others,
%   i.e. if a singular value collapses to machine level.  Everything above that
%   is a PRACTICAL statement about this 47-week data set.

if nargin < 4 || isempty(opts)
    opts = struct();
end
relStep = 1e-4;
if isfield(opts, 'rel_step') && ~isempty(opts.rel_step)
    relStep = double(opts.rel_step);
end
tol = 1e-10;
if isfield(opts, 'tol') && ~isempty(opts.tol)
    tol = double(opts.tol);
end

names = cellstr(string(names));
p = numel(names);

mu0 = mwspm_calib_kernel(pre, th, struct()).Yhat;
n = numel(mu0);

Sraw = zeros(n, p);
usedStep = zeros(1, p);
for j = 1:p
    f = names{j};
    v = double(th.(f));
    h = relStep * max(abs(v), 1e-3);
    % central difference in the NATURAL parameter, then scaled by theta
    tp = th; tp.(f) = v + h;
    tm = th; tm.(f) = v - h;
    d = (mwspm_calib_kernel(pre, tp, struct()).Yhat ...
       - mwspm_calib_kernel(pre, tm, struct()).Yhat) / (2 * h);
    Sraw(:, j) = v * d;
    usedStep(j) = h;
end

colnorm = sqrt(sum(Sraw.^2, 1))';
Sh = Sraw ./ max(colnorm', 1e-300);

[~, Sig, V] = svd(Sh, 'econ');
sig = diag(Sig);
sigMax = max(sig);
condIndices = sigMax ./ max(sig, 1e-300);

% Belsley variance-decomposition proportions
denom = sum((V.^2) ./ (sig.^2)', 2);         % p x 1
varProp = ((V.^2) ./ (sig.^2)') ./ max(denom, 1e-300);

% collinearity sets: for every condition index > 30, the parameters whose
% variance proportion exceeds 0.5
sets = {};
for k = 1:numel(sig)
    if condIndices(k) <= 30
        continue
    end
    idx = find(varProp(:, k) > 0.5);
    if numel(idx) >= 2
        sets{end + 1} = names(idx); %#ok<AGROW>
    end
end

% correlation matrix of the least-squares estimates
XtX = Sh' * Sh;
XtXi = pinv(XtX);
d = sqrt(max(diag(XtXi), 1e-300));
corrEst = XtXi ./ (d * d');

S = struct();
S.names = names;
S.th = th;
S.rel_step = relStep;
S.used_step = usedStep;
S.S_raw = Sraw;
S.colnorm = colnorm;
S.Sh = Sh;
S.svd = sig;
S.cond_number = sigMax / max(min(sig), 1e-300);
S.cond_indices = condIndices;
S.var_prop = varProp;
S.collinear_sets = sets;
S.corr_estimates = corrEst;
S.practical_rank = sum(sig > tol * sigMax);
S.tol = tol;
S.n_obs = n;
end
