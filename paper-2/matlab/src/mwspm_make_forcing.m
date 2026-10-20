function F = mwspm_make_forcing(D, p, opts)
%MWSPM_MAKE_FORCING  Pre-compute the time-dependent coefficients of the RHS.
%
%   F = mwspm_make_forcing(D, p)              ZOH sampling, mu variant 'raw'
%   F = mwspm_make_forcing(D, p, opts)
%       opts.interp      'zoh' (default) | 'linear'
%       opts.mu_variant  'raw' (default) | 'day'
%
%   Sampling contract
%   -----------------
%   The repository's forcing is defined on a DAILY grid (day 1 .. day 329 of
%   the analysis window, t = 0 .. 328 in days elapsed).  A continuous-time ODE
%   needs values in between, so an explicit sampling rule is declared:
%     'zoh'    -> t in [k, k+1) uses day k+1 (zero-order hold).  With
%                 dt = 1 this reproduces the source project's daily map
%                 exactly, which is what makes the C5 comparison meaningful.
%     'linear' -> piecewise-linear interpolation of the daily values.
%   The rule is recorded in manifest.json; it is NOT a free tuning knob.
%
%   Traits are pre-computed on the daily grid and then sampled with the same
%   rule (faster, and for ZOH identical to computing them from the sampled
%   temperature).
%
%   F has NO function handles, so it can be partially serialised for audit.
%
%   Fields that are VARIANT DEPENDENT : mu_daily, mu_eff_daily (these hold the
%   variant selected by opts.mu_variant).
%   Fields that are VARIANT INDEPENDENT : mu_raw_daily, mu_eff_raw_daily (the
%   continuous rate).  Diagnostics that must be comparable across constructions
%   must read the latter pair - see mwspm_stiffness.

if nargin < 3 || isempty(opts)
    opts = struct();
end
interpRule = 'zoh';
if isfield(opts, 'interp') && ~isempty(opts.interp)
    interpRule = lower(char(opts.interp));
end
muVariant = 'raw';
if isfield(opts, 'mu_variant') && ~isempty(opts.mu_variant)
    muVariant = lower(char(opts.mu_variant));
end
assert(any(strcmp(interpRule, {'zoh', 'linear'})), ...
    'mwspm:forcing', 'unknown interpolation rule "%s"', interpRule);
assert(any(strcmp(muVariant, {'raw', 'day'})), ...
    'mwspm:forcing', 'unknown mu variant "%s"', muVariant);

tr = mwspm_thermal_traits(D.T_daily, muVariant, struct('mu_cap', p.mu_v_cap));

s = ones(D.n_days, 1);
if strcmp(p.season_mode, 'window')
    % normalised on the 329-day window; mwspm_season returns the exact factor
    % so that any other caller recomputing s on the same grid gets the same
    % curve (p is passed by value and is deliberately NOT mutated here).
    [s, sInfo] = mwspm_season(D.doy, p.season);
else
    sInfo = struct('norm_factor', 1, 'normalised', false, 'mean_shat', 1, ...
                   'onset_doy', NaN, 'offset_doy', NaN, 'span_days', NaN, 'kappa_days', NaN);
end

F = struct();
F.interp      = interpRule;
F.mu_variant  = muVariant;
F.tGrid       = D.t_days(:);                      % 0 .. n_days-1
F.n           = D.n_days;
F.T_daily     = D.T_daily(:);
F.a_daily     = tr.a;
F.b_daily     = tr.b;
F.c_daily     = tr.c;
F.sigma_v_daily = tr.sigma_v;
F.mu_daily    = tr.mu;                            % the SELECTED variant
F.mu_eff_daily = p.mu_s * tr.mu;                  % mu_s scaling is applied HERE
% The continuous (raw) rate is kept separately and is INDEPENDENT of
% opts.mu_variant, because tr.mu itself is not: under 'day' the array stored in
% F.mu_daily already holds 1-exp(-mu_raw) (above 0.5).  Any diagnostic that must
% be comparable across constructions (mwspm_stiffness, T11 checks V51/V55) has
% to read these two fields rather than F.mu_daily.  Purely additive: no equation
% reads them, so the dynamics are unchanged.
F.mu_raw_daily     = tr.mu_raw;
F.mu_eff_raw_daily = p.mu_s * tr.mu_raw;
F.s_daily     = s;
F.Nv_daily    = p.m0 * p.N_h * s;                 % external supply N_v(t)
F.Lambda_daily = D.Lambda_daily(:);
F.N_h         = p.N_h;
F.m0          = p.m0;
F.mu_s        = p.mu_s;

F.diagnostics = struct();
F.diagnostics.max_mu_eff            = max(F.mu_eff_daily);
F.diagnostics.min_mu_eff            = min(F.mu_eff_daily);
F.diagnostics.max_a                 = max(F.a_daily);
F.diagnostics.max_sigma_v           = max(F.sigma_v_daily);
F.diagnostics.Nv_peak               = max(F.Nv_daily);
F.diagnostics.Nv_peak_day           = find(F.Nv_daily == max(F.Nv_daily), 1);
F.diagnostics.s_mean                = mean(s);
F.diagnostics.s_max                 = max(s);
F.diagnostics.s_norm_factor         = sInfo.norm_factor;
F.diagnostics.s_onset_doy           = sInfo.onset_doy;
F.diagnostics.s_offset_doy          = sInfo.offset_doy;
F.diagnostics.s_span_days           = sInfo.span_days;
end
