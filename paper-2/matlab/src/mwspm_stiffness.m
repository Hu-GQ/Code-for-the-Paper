function s = mwspm_stiffness(F, p, dtProbe)
%MWSPM_STIFFNESS  Stiffness diagnostics of the active model (solver selection input).
%
%   s = mwspm_stiffness(F, p, dtProbe)
%
%   Two complementary indicators are returned, both computed from the FORCED
%   COEFFICIENTS only (no state, no random numbers), so they are reproducible
%   and auditable:
%
%   1) relaxation-rate spread.  On the disease-free state the linearisation
%      decouples into the human block {-sigma_h, -gamma_h, 0} and the mosquito
%      block {-mu_eff, -(mu_eff+sigma_v), -mu_eff}.  The spread
%          ratio = max(rates) / min(rates)
%      is the classical stiffness ratio proxy for this model class.  The
%      CANONICAL ratio is the one built from the CONTINUOUS (raw) mortality
%      rate, because that is the rate that appears in R0 and in the daily-map
%      stability condition; the day-map variant is reported alongside.
%
%   2) explicit-Euler stability metric for the mosquito mortality channel,
%          m_dt = max_t |1 - mu_eff(t)*dtProbe|
%      which is the quantity that the C5 sub-claim is about (CL-06).  It is
%      reported for BOTH mortality variants because the source project's daily
%      map uses the 1-exp corrected variant and then multiplies by mu_s, so
%      the two variants have different stability behaviour (see the model
%      card field mu_variant_note).  Both counts are returned and must be
%      reported together.
%
%   INVARIANCE TO HOW F WAS CONSTRUCTED (T11 check V55)
%   --------------------------------------------------
%   Both variants are derived from F.mu_raw_daily, i.e. from the CONTINUOUS
%   rate that mwspm_make_forcing stores independently of opts.mu_variant, and
%   from F.mu_s.  They are therefore independent of how F was built.
%
%   This is deliberate and was a defect before, in two stages:
%     (1) the original revision read the "raw" numbers out of F.mu_eff_daily,
%         which IS the day-variant value when F was built with
%         mu_variant='day', yet still labelled them "raw": the same call
%         returned stiffness_ratio 587.984 for a raw-built F and 118.722 for a
%         day-built F;
%     (2) the first repair still used F.mu_daily, which is ALSO variant
%         dependent for the same reason, so the spread persisted (11.2215 vs
%         60.522) - which is why mwspm_make_forcing now carries
%         mu_raw_daily / mu_eff_raw_daily.
%   T11.01 and T07.02 quote these numbers, so a construction-dependent value
%   was not admissible.  The raw-variant numbers are unchanged by the fix (they
%   are what T10's E0 run recorded: ratio 587.984, dt_euler_stable_max
%   0.032376 d, max|1-mu_eff*dt| 60.522, 17/329 days outside the domain).
%
%   This function draws NO conclusion: it only measures.  The C5 verdict
%   belongs to T11.01.

assert(nargin >= 3 && ~isempty(dtProbe), 'mwspm:stiffness', 'dtProbe is required');
assert(isfield(F, 'mu_s'), 'mwspm:stiffness', ...
    'F must carry mu_s (see mwspm_make_forcing)');
assert(isfield(F, 'mu_raw_daily'), 'mwspm:stiffness', ...
    ['F must carry mu_raw_daily, the continuous rate independent of mu_variant; ' ...
     'a forcing built before this field existed cannot support a variant-invariant diagnostic']);
assert(abs(F.mu_s - p.mu_s) <= 1e-12 * max(1, abs(p.mu_s)), 'mwspm:stiffness', ...
    'F.mu_s (%.17g) and p.mu_s (%.17g) disagree; the forcing was built with other parameters', ...
    F.mu_s, p.mu_s);

muRaw = F.mu_raw_daily(:);                          % continuous rate, variant-free
muDay = 1 - exp(-muRaw);
muDay(muRaw <= 0.5) = muRaw(muRaw <= 0.5);          % model-card day mapping
muEffRaw = F.mu_s * muRaw;
muEffDay = F.mu_s * muDay;

sv = F.sigma_v_daily(:);

% ---- 1) relaxation-rate spread ------------------------------------------
ratesRaw = [muEffRaw; muEffRaw + sv; p.sigma_h; p.gamma_h];
ratesDay = [muEffDay; muEffDay + sv; p.sigma_h; p.gamma_h];
ratesRaw = ratesRaw(ratesRaw > 0);
ratesDay = ratesDay(ratesDay > 0);
assert(~isempty(ratesRaw) && ~isempty(ratesDay), 'mwspm:stiffness', ...
    'no positive relaxation rate found; the stiffness ratio would be meaningless');

s = struct();
s.rate_fast_per_day   = max(ratesRaw);
s.rate_slow_per_day   = min(ratesRaw);
s.stiffness_ratio     = s.rate_fast_per_day / s.rate_slow_per_day;
s.dt_euler_stable_max = 2 / s.rate_fast_per_day;
s.dt_probe_days       = dtProbe;

s.rate_fast_day_per_day   = max(ratesDay);
s.rate_slow_day_per_day   = min(ratesDay);
s.stiffness_ratio_day     = s.rate_fast_day_per_day / s.rate_slow_day_per_day;

s.mu_eff_raw_min = min(muEffRaw);
s.mu_eff_raw_max = max(muEffRaw);
s.mu_eff_day_min = min(muEffDay);
s.mu_eff_day_max = max(muEffDay);

% ---- 2) explicit-Euler stability metric ---------------------------------
oneMinusRaw = abs(1 - muEffRaw * dtProbe);
oneMinusDay = abs(1 - muEffDay * dtProbe);

s.max_abs_1_minus_mu_eff_dt_raw = max(oneMinusRaw);
s.n_days_outside_stability_raw  = sum(oneMinusRaw > 1);
s.max_abs_1_minus_mu_eff_dt_day = max(oneMinusDay);
s.n_days_outside_stability_day  = sum(oneMinusDay > 1);

% ---- provenance ---------------------------------------------------------
s.canonical_ratio_variant = 'raw';
s.mu_variant_of_F = '';
if isfield(F, 'mu_variant') && ~isempty(F.mu_variant)
    s.mu_variant_of_F = char(F.mu_variant);
end
s.invariant_to_F_construction = true;
s.invariance_note = ['Both variants are recomputed from F.mu_raw_daily (the continuous rate, ' ...
    'stored independently of mu_variant) and F.mu_s, so the numbers do not depend on the ' ...
    'mu_variant used to build F (verified by T11 check V55).'];
s.legacy_reference = ['theory/derivations/EDGE-01 (T07) reports max|1-mu_eff*dt|=11.22 and ' ...
    '17/329 days outside the domain, computed with the continuous (raw) rate; the raw ' ...
    'entries above are the same quantities and are invariant to the forcing construction.'];
end
