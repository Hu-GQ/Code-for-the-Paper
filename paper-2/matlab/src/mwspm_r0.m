function [R0, tr] = mwspm_r0(T, m, mu, gamma_h, opts)
%MWSPM_R0  Pointwise basic reproduction ratio of the active model.
%
%   [R0, tr] = mwspm_r0(T, m, mu, gamma_h)
%   [R0, tr] = mwspm_r0(T, m, mu, gamma_h, opts)
%
%     T        temperature(s), degrees C
%     m        vector-to-human ratio(s) used in the NGM  (m = m0*s(t))
%     mu       adult mosquito mortality rate(s), 1/day   (mu = mu_s*mu_v_raw)
%     gamma_h  human recovery rate, 1/day
%     opts.coef  optional thermal-coefficient override forwarded to
%                mwspm_thermal_traits (the analysis trait sensitivity); default OFF
%
%     R0 = sqrt( m * a^2 * b * c * sigma_v / ( gamma_h * mu * (mu + sigma_v) ) )
%
%   CONVENTIONS THAT MUST BE DECLARED WHENEVER THIS IS REPORTED
%   (theory/model_spec.md Sec.11.2, obligation O6 and the "square root" pitfall):
%     1. the sqrt convention is used ("two full cycles"); the threshold R0=1
%        is unchanged by dropping the sqrt, but the VALUE is not comparable
%        across conventions;
%     2. mu must be the CONTINUOUS mortality rate (variant 'raw'), never the
%        daily-corrected variant;
%     3. the vector-abundance caliber (N_v vs V) must be stated - see
%        mwspm_r0_timeseries and the statement/the analysis.
%
%   sigma_h cancels out of R0 (the exposed human eventually becomes
%   infectious), which is verified symbolically in
%   src/python/sym_model_spec_checks.py.

if nargin < 5 || isempty(opts)
    opts = struct();
end
trOpts = struct();
if isfield(opts, 'coef') && ~isempty(opts.coef)
    trOpts.coef = opts.coef;
end
if isfield(opts, 'alt') && ~isempty(opts.alt)
    trOpts.alt = opts.alt;   % C1: alternative functional forms (hybrid source)
end
tr = mwspm_thermal_traits(T, 'raw', trOpts);
T  = double(T(:));
m  = double(m(:));
if nargin < 3 || isempty(mu)
    mu = tr.mu;
end
mu = double(mu(:));

assert(numel(m) == numel(T) || isscalar(m), 'mwspm:r0', 'm must be scalar or match T');
assert(numel(mu) == numel(T) || isscalar(mu), 'mwspm:r0', 'mu must be scalar or match T');

num = m .* tr.a.^2 .* tr.b .* tr.c .* tr.sigma_v;
den = max(gamma_h .* mu .* (mu + tr.sigma_v), 0);
R0  = sqrt(max(num ./ max(den, realmin), 0));
end
