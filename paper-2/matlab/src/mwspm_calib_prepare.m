function pre = mwspm_calib_prepare(D, modelCard, opts)
%MWSPM_CALIB_PREPARE  Precompute the parameter-INDEPENDENT pieces of the
%   observation-consistent forward map used by the T12 calibration audit.
%
%   pre = mwspm_calib_prepare(D, modelCard)
%   pre = mwspm_calib_prepare(D, modelCard, struct('mu_variant','day'))
%
%   Why this exists
%   ---------------
%   mwspm_make_forcing rebuilds the whole thermal-trait block on every call
%   (~78 ms on this installation).  A profile likelihood needs thousands of
%   forward solves in which only (m0, mu_s, onset, offset, gamma, beta, rho, d)
%   change.  The temperature-dependent traits a(T), b(T), c(T), sigma_v(T),
%   mu_v(T) depend ONLY on the data, so they are computed once here and reused.
%
%   Not a second model definition: the trait curves still come from
%   mwspm_thermal_traits (the single source), and the season factor still comes
%   from mwspm_season inside mwspm_calib_kernel.  Whether the fast kernel is a
%   faithful replacement for the canonical path (mwspm_rhs +
%   mwspm_integrate_fixed + mwspm_observe) is NOT assumed - it is MEASURED by
%   an equivalence test in exp_t12_calibration (recorded in the run data).
%
%   Fields
%   ------
%     mu_variant        'day' (source-faithful) | 'raw'
%     a,b,c,sigma_v     (n_days x 1) daily traits
%     mu_variant_daily  the SELECTED mortality variant (rate per day)
%     mu_raw_daily      the continuous rate (variant independent)
%     Lambda            (n_days x 1) importation forcing
%     doy,T_daily,N_h,sigma_h,gamma_h,kappa_days,s_floor,min_window_days
%     n_days,n_weeks,Y  observations (weekly reported local cases)

if nargin < 3 || isempty(opts)
    opts = struct();
end
muVariant = 'day';
if isfield(opts, 'mu_variant') && ~isempty(opts.mu_variant)
    muVariant = lower(char(opts.mu_variant));
end
assert(any(strcmp(muVariant, {'raw', 'day'})), 'mwspm:calib', ...
    'unknown mu variant "%s"', muVariant);

p0 = mwspm_params(modelCard);
trOpts = struct('mu_cap', p0.mu_v_cap);
if isfield(opts, 'coef') && ~isempty(opts.coef)
    trOpts.coef = opts.coef;   % T14.02 trait-parameterisation sensitivity (opt-in)
end
if isfield(opts, 'alt') && ~isempty(opts.alt)
    trOpts.alt = opts.alt;     % C1: alternative functional forms (hybrid source)
end
tr = mwspm_thermal_traits(D.T_daily, muVariant, trOpts);

pre = struct();
pre.mu_variant        = muVariant;
pre.a                 = tr.a;
pre.b                 = tr.b;
pre.c                 = tr.c;
pre.sigma_v           = tr.sigma_v;
pre.mu_variant_daily  = tr.mu;
pre.mu_raw_daily      = tr.mu_raw;
pre.coef              = tr.coef;             % provenance: which triples were used
pre.coef_override     = tr.coef_override;
pre.forms             = tr.forms;            % C1 provenance: named alternative forms
pre.form_override     = tr.form_override;
pre.mu_source         = tr.mu_source;
pre.Lambda            = D.Lambda_daily(:);
pre.doy               = D.doy(:);
pre.T_daily           = D.T_daily(:);
pre.N_h               = p0.N_h;
pre.sigma_h           = p0.sigma_h;
pre.gamma_h           = p0.gamma_h;
pre.kappa_days        = p0.kappa_days;
pre.s_floor           = p0.s_floor;
pre.min_window_days   = p0.min_window_days;
pre.n_days            = D.n_days;
pre.n_weeks           = D.n_weeks;
pre.Y                 = D.Y_weekly(:);
pre.Y_import_weekly   = D.Yimp_weekly(:);
pre.week_index        = (1:D.n_weeks)';
assert(pre.n_weeks * 7 == pre.n_days, 'mwspm:calib', ...
    'n_weeks*7 must equal n_days for exact weekly edges');
end
