function tr = mwspm_thermal_traits(T, variant, opts)
%MWSPM_THERMAL_TRAITS  Temperature-dependent mosquito traits (Mordecai 2017).
%
%   tr = mwspm_thermal_traits(T)              continuous-rate variant ('raw')
%   tr = mwspm_thermal_traits(T, 'day')       source-faithful daily variant
%   tr = mwspm_thermal_traits(T, 'raw', opts) with opts.mu_cap override
%   tr = mwspm_thermal_traits(T, 'raw', opts) with opts.coef override (T14.02)
%
%   Source of the curves: Mordecai et al. 2017, PLoS NTD 11(4):e0005568,
%   Table A, Ae. albopictus (registered in data/metadata/parameters.csv).
%
%   THE TWO MORTALITY VARIANTS ARE NOT INTERCHANGEABLE:
%     'raw'  mu_v_raw = min(1/lf, mu_cap)              -> continuous-time rate
%            used for R0(t), T0(t) and for any adaptive ODE reference.
%     'day'  mu_v_day = 1-exp(-mu_v_raw) where raw>0.5 -> the exact daily
%            survival fraction; used ONLY to reproduce the source project's
%            daily stepping (source matlab/dengue_rhs.m calls
%            thermal_traits_day).  Note that the gen3 fit multiplies this by
%            mu_s AFTERWARDS, which destroys the |1-mu*dt|<=1 safety; the two
%            variants must therefore always be simulated and reported
%            separately (see config/experiments/model_mwspm_gen3_rhohat.json,
%            field mu_variant_note, and task T11.01).
%
%   COEFFICIENT OVERRIDE (T14.02, opt-in, default OFF)
%   --------------------------------------------------
%   opts.coef.<trait> = [c, T0, Tm] replaces the fitted coefficient triple of
%   that trait ONLY; every other trait keeps its published value.  This exists
%   so that the sensitivity of R0(t) to the thermal-TRAIT PARAMETERISATION can
%   be MEASURED (T14, EXPERIMENT_PLAN Sec.5-4).  It is NOT a second model
%   definition and NOT an alternative trait source: the functional forms and
%   the single source of the curves stay here, in this function.
%
%   ALTERNATIVE FUNCTIONAL FORMS (C1 hybrid source, 2026-10-07, opt-in)
%   -------------------------------------------------------------------
%   A second PUBLICATION may parameterise a trait with a DIFFERENT functional
%   form, which cannot be expressed as a [c,T0,Tm] triple.  opts.alt.<trait>
%   therefore declares a named form explicitly:
%
%     opts.alt.a       = struct('form','linear','par',[slope intercept], ...
%                               'range',[lo hi],'outside','clip')
%     opts.alt.sigma_v = struct('form','eip_exp','par',[floor c1 c2])
%                        => sigma_v = 1/(floor + exp(c1 - c2*T))   (EIP in days)
%     opts.alt.mu      = struct('form','poly_direct_mu','par',[c0 c1 c2 c3 c4])
%                        => mu_raw = polyval(par, T)  (mortality rate DIRECTLY;
%                           there is no lifespan intermediate)
%
%   'outside' (default 'clip') says what happens outside 'range': 'clip' holds
%   the boundary value, 'zero' returns 0, 'extrapolate' uses the formula as is.
%   The rule must be PRE-REGISTERED in the experiment config: choosing it after
%   seeing the result would be a post-hoc degree of freedom.
%
%   tr.forms (the alt spec actually used), tr.form_override and tr.mu_source
%   record which parameterisation produced the numbers.  When opts.alt.mu is
%   used, tr.lf is DERIVED as 1./mu_raw for reporting only (tr.lf_derived=true).
%
%   Recognised fields (all optional):
%     a, b, c, sigma_v   Briere    [c, T0, Tm]  defaults below
%     lf                 quadratic [c, T0, Tm]  (c is negative)
%   The coefficients actually used are returned in tr.coef, and tr.coef_override
%   is true iff at least one triple was overridden, so a run can PROVE which
%   parameterisation produced its numbers.
%
%   Output fields: a, b, c, sigma_v (=PDR=1/EIP), mu, lf, mu_variant,
%                  coef, coef_override, mu_cap, forms, form_override, mu_source.

if nargin < 2 || isempty(variant)
    variant = 'raw';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
muCap = 5.0;
if isfield(opts, 'mu_cap') && ~isempty(opts.mu_cap)
    muCap = opts.mu_cap;
end
alt = struct();
if isfield(opts, 'alt') && ~isempty(opts.alt) && isstruct(opts.alt)
    alt = opts.alt;
end

% ---- published Mordecai 2017 (Ae. albopictus) coefficient triples ----------
coef = struct();
coef.a       = [1.93e-4, 10.25, 38.32];   % biting rate, 1/day
coef.b       = [7.35e-4, 15.84, 36.40];   % mosquito -> human
coef.c       = [4.39e-4,  3.62, 36.82];   % human -> mosquito
coef.sigma_v = [1.09e-4, 10.39, 43.05];   % EIP progression rate (=PDR)
coef.lf      = [-1.43,   13.41, 31.51];   % adult lifespan, days

isQuadratic = struct('a', false, 'b', false, 'c', false, 'sigma_v', false, 'lf', true);

% ---- opt-in override: only whole triples, only known trait names ----------
override = false;
if isfield(opts, 'coef') && ~isempty(opts.coef)
    given = opts.coef;
    if isstruct(given)
        fn = fieldnames(given);
    else
        fn = {};
    end
    for i = 1:numel(fn)
        f = fn{i};
        assert(isfield(coef, f), 'mwspm:thermal', ...
            'unknown trait "%s" in opts.coef (known: %s)', f, strjoin(fieldnames(coef), ', '));
        v = double(given.(f));
        v = v(:)';
        assert(numel(v) == 3, 'mwspm:thermal', ...
            'opts.coef.%s must be a 3-vector [c T0 Tm]', f);
        assert(all(isfinite(v)), 'mwspm:thermal', 'opts.coef.%s must be finite', f);
        coef.(f) = v;
        override = true;
    end
end

T = double(T(:));

if isQuadratic.a, tr.a = mwspm_quadratic(T, coef.a(1), coef.a(2), coef.a(3));
else,             tr.a = mwspm_briere(T,     coef.a(1), coef.a(2), coef.a(3)); end
if isQuadratic.b, tr.b = mwspm_quadratic(T, coef.b(1), coef.b(2), coef.b(3));
else,             tr.b = mwspm_briere(T,     coef.b(1), coef.b(2), coef.b(3)); end
if isQuadratic.c, tr.c = mwspm_quadratic(T, coef.c(1), coef.c(2), coef.c(3));
else,             tr.c = mwspm_briere(T,     coef.c(1), coef.c(2), coef.c(3)); end
if isQuadratic.sigma_v, tr.sigma_v = mwspm_quadratic(T, coef.sigma_v(1), coef.sigma_v(2), coef.sigma_v(3));
else,                   tr.sigma_v = mwspm_briere(T,     coef.sigma_v(1), coef.sigma_v(2), coef.sigma_v(3)); end
if isQuadratic.lf, tr.lf = mwspm_quadratic(T, coef.lf(1), coef.lf(2), coef.lf(3));
else,              tr.lf = mwspm_briere(T,     coef.lf(1), coef.lf(2), coef.lf(3)); end
tr.PDR = tr.sigma_v;                                    % source-project name

% ---- opt-in alternative functional forms (C1) ------------------------------
formOverride = false;
muAlt = [];
if isfield(alt, 'a')
    tr.a = mwspm_alt(T, alt.a, 'a'); formOverride = true;
end
if isfield(alt, 'b')
    tr.b = mwspm_alt(T, alt.b, 'b'); formOverride = true;
end
if isfield(alt, 'c')
    tr.c = mwspm_alt(T, alt.c, 'c'); formOverride = true;
end
if isfield(alt, 'sigma_v')
    tr.sigma_v = mwspm_alt(T, alt.sigma_v, 'sigma_v'); formOverride = true;
    tr.PDR = tr.sigma_v;
end
if isfield(alt, 'lf')
    tr.lf = mwspm_alt(T, alt.lf, 'lf'); formOverride = true;
end
if isfield(alt, 'mu')
    muAlt = mwspm_alt(T, alt.mu, 'mu'); formOverride = true;
end

mu = zeros(size(T));
if ~isempty(muAlt)
    % A source that publishes MORTALITY DIRECTLY: no lifespan intermediate.
    mu = min(max(muAlt, 0), muCap);
    tr.mu_source = 'poly_direct_mu';
    tr.lf = 1 ./ max(mu, 1e-6);          % derived, for reporting only
    tr.lf_derived = true;
else
    pos = tr.lf > 0;
    mu(pos)  = min(1 ./ tr.lf(pos), muCap);
    mu(~pos) = muCap;                                       % out-of-range temperature
    tr.mu_source = 'lifespan_reciprocal';
    tr.lf_derived = false;
end
tr.mu_raw = max(mu, 1e-6);

switch lower(char(variant))
    case 'raw'
        mu = tr.mu_raw;
    case 'day'
        mu = tr.mu_raw;
        idx = mu > 0.5;
        mu(idx) = 1 - exp(-mu(idx));
    otherwise
        error('mwspm:thermal', 'unknown mu variant "%s"', char(variant));
end
tr.mu = mu;
tr.mu_variant = lower(char(variant));

% defensive floor on the transmission-relevant traits
tr.a = max(tr.a, 0);
tr.b = max(tr.b, 0);
tr.c = max(tr.c, 0);
tr.sigma_v = max(tr.sigma_v, 0);
tr.PDR = tr.sigma_v;

tr.coef = coef;
tr.coef_override = override;
tr.mu_cap = muCap;
tr.forms = alt;
tr.form_override = formOverride;
end

% ---------------------------------------------------------------------------
function q = mwspm_alt(T, spec, traitName)
%MWSPM_ALT  Evaluate one named alternative functional form (C1 hybrid source).
assert(isstruct(spec) && isfield(spec, 'form') && isfield(spec, 'par'), ...
    'mwspm:thermal', 'opts.alt.%s must be a struct with fields form and par', traitName);
form = lower(char(string(spec.form)));
par  = double(spec.par(:))';
assert(all(isfinite(par)), 'mwspm:thermal', 'opts.alt.%s.par must be finite', traitName);

switch form
    case 'linear'
        assert(numel(par) == 2, 'mwspm:thermal', ...
            'opts.alt.%s linear needs par = [slope intercept]', traitName);
        q = par(1) * T + par(2);
    case 'eip_exp'
        % EIP in days = par(1) + exp(par(2) - par(3)*T); the trait is its RATE.
        assert(numel(par) == 3, 'mwspm:thermal', ...
            'opts.alt.%s eip_exp needs par = [floor_days c1 c2]', traitName);
        eip = par(1) + exp(par(2) - par(3) * T);
        q = 1 ./ max(eip, 1e-6);
    case 'poly_direct_mu'
        assert(numel(par) >= 1, 'mwspm:thermal', ...
            'opts.alt.%s poly_direct_mu needs par = [c0 c1 ... cn] (descending)', traitName);
        q = polyval(par, T);
    otherwise
        error('mwspm:thermal', 'unknown alternative form "%s" for trait %s', form, traitName);
end

% ---- declared validity range and the PRE-REGISTERED out-of-range rule ------
if isfield(spec, 'range') && ~isempty(spec.range)
    rng_ = double(spec.range(:))';
    assert(numel(rng_) == 2 && rng_(1) < rng_(2), 'mwspm:thermal', ...
        'opts.alt.%s.range must be [lo hi]', traitName);
    outside = 'clip';
    if isfield(spec, 'outside') && ~isempty(spec.outside)
        outside = lower(char(string(spec.outside)));
    end
    below = T < rng_(1);
    above = T > rng_(2);
    switch outside
        case 'clip'
            loV = mwspm_altAt(spec, form, par, rng_(1));
            hiV = mwspm_altAt(spec, form, par, rng_(2));
            q(below) = loV;
            q(above) = hiV;
        case 'zero'
            q(below | above) = 0;
        case 'extrapolate'
            % keep q as computed
        otherwise
            error('mwspm:thermal', 'unknown outside rule "%s" for trait %s', outside, traitName);
    end
end
q = q(:);
end

% ---------------------------------------------------------------------------
function v = mwspm_altAt(spec, form, par, Tat)   %#ok<INUSD>
% Evaluate the raw (un-ranged) form at a single temperature, for 'clip'.
switch form
    case 'linear',          v = par(1) * Tat + par(2);
    case 'eip_exp',         v = 1 / max(par(1) + exp(par(2) - par(3) * Tat), 1e-6);
    case 'poly_direct_mu',  v = polyval(par, Tat);
    otherwise,              error('mwspm:thermal', 'unknown form "%s"', form);
end
end

% ---------------------------------------------------------------------------
function q = mwspm_briere(T, c, T0, Tm)
q = zeros(size(T));
msk = (T > T0) & (T < Tm);
q(msk) = c * T(msk) .* (T(msk) - T0) .* sqrt(Tm - T(msk));
end

function q = mwspm_quadratic(T, c, T0, Tm)
q = zeros(size(T));
msk = (T > T0) & (T < Tm);
q(msk) = c * (T(msk) - T0) .* (T(msk) - Tm);
end
