function m = mwspm_metrics(obs, pred)
%MWSPM_METRICS  Goodness-of-fit metrics with hand-implemented likelihoods.
%
%   m = mwspm_metrics(obs, pred)
%
%   The Statistics and Machine Learning Toolbox is NOT installed on this
%   installation (recorded in environment/ENVIRONMENT.md and re-verified by
%   mwspm_env).  Everything here is therefore implemented from first
%   principles; no statistical distribution object is used.
%
%   Metrics returned:
%     R2_log, RMSE_log   on the log1p scale (legacy project convention)
%     R2_raw             on the raw count scale
%     poisson_nll        sum(pred - obs.*log(pred)) (parameter-free terms dropped)
%     peak_week_model, peak_week_observed, peak_value_model, peak_value_observed
%
%   WARNING (project memory, methodology pitfall 1): R2_log alone can be
%   inflated by a temporal shift (a solution mis-phased by 11 weeks scored
%   0.71 -> 0.84), so the phase soft penalty must stay in the objective and
%   R2_raw must always be reported next to R2_log.  This function returns BOTH
%   and never a single number.

obs  = double(obs(:));
pred = double(pred(:));
assert(numel(obs) == numel(pred), 'mwspm:metrics', 'length mismatch');

eps0 = 1e-12;
o = log1p(obs);
q = log1p(max(pred, 0));
ssResLog = sum((o - q).^2);
ssTotLog = sum((o - mean(o)).^2);
ssResRaw = sum((obs - pred).^2);
ssTotRaw = sum((obs - mean(obs)).^2);

m = struct();
m.R2_log  = 1 - ssResLog / max(ssTotLog, eps0);
m.RMSE_log = sqrt(mean((o - q).^2));
m.R2_raw  = 1 - ssResRaw / max(ssTotRaw, eps0);
m.RMSE_raw = sqrt(mean((obs - pred).^2));
m.poisson_nll = sum(max(pred, eps0) - obs .* log(max(pred, eps0)));
m.n_weeks = numel(obs);
m.sum_observed = sum(obs);
m.sum_predicted = sum(pred);
m.peak_week_observed = find(obs == max(obs), 1);
m.peak_week_model    = find(pred == max(pred), 1);
m.peak_value_observed = max(obs);
m.peak_value_model    = max(pred);
m.peak_week_shift     = m.peak_week_model - m.peak_week_observed;
end
