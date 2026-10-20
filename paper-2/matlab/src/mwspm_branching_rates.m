function R = mwspm_branching_rates(D, pBase, th, muVariant)
%MWSPM_BRANCHING_RATES  Daily rates of the 4-stage vector-host branching process.
%
%   R = mwspm_branching_rates(D, pBase, th, muVariant)
%
%   Returns the event rates of the continuous-time multi-type branching process
%   that is the linearisation of M-WSPM gen3_rhohat about the disease-free
%   state, on the SAME daily grid and with the SAME forcing as the ODE model.
%
%   Particle types (order is fixed and used by every caller):
%       1 = E_h   exposed (not yet infectious) human
%       2 = I_h   infectious human
%       3 = E_v   exposed (latent) mosquito
%       4 = I_v   infectious mosquito
%
%   Events (all rates are per DAY):
%       E_h -> I_h                  sigma_h
%       I_h -> removed              gamma_h
%       I_h -> I_h + E_v            rhv(t) = a(t) c(t) m(t)     [m = m0 s(t)]
%       E_v -> I_v                  sigma_v(t)
%       E_v -> removed              mu(t)                       [competing]
%       I_v -> I_v + E_h            rvh(t) = a(t) b(t)
%       I_v -> removed              mu(t)
%
%   Mean-field identity used as a check: the expectation of the branching
%   process obeys dX/dt = A(t) X with
%        A = [ -sigma_h      0            0          rvh       ;
%               sigma_h     -gamma_h      0          0         ;
%               0            rhv         -(sv+mu)    0         ;
%               0            0            sv        -mu        ]
%   which is EXACTLY the linearisation of mwspm_rhs with S_h = N_h,
%   S_v = N_v(t) = m0 N_h s(t) and Lambda = 0.  That is the branching
%   approximation, declared, not proved (obligation O7/O8 territory).
%
%   IMPORTATION IS ZEROED HERE.  The fitted model uses Lambda(t) > 0 as the
%   seeding mechanism; a branching process of a SINGLE introduction must have
%   Lambda = 0 or extinction is impossible by construction.  This is a
%   modelling choice and is recorded in raw_choices.csv.

if nargin < 4 || isempty(muVariant)
    muVariant = 'day';
end
muVariant = lower(char(muVariant));

assert(isfield(th, 'm0') && isfield(th, 'mu_s') && isfield(th, 'gamma') && isfield(th, 'beta'), ...
    'mwspm:branching', 'theta must carry m0, mu_s, gamma, beta');

sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pBase.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pBase.s_floor, 'min_window_days', pBase.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(D.doy, sp);

tr = mwspm_thermal_traits(D.T_daily, muVariant, struct('mu_cap', pBase.mu_v_cap));
trRaw = mwspm_thermal_traits(D.T_daily, 'raw', struct('mu_cap', pBase.mu_v_cap));

m = th.m0 * s;
mu = th.mu_s * tr.mu;
muRaw = th.mu_s * trRaw.mu;

R = struct();
R.mu_variant   = muVariant;
R.n_days       = numel(D.T_daily);
R.n_weeks      = D.n_weeks;
R.sigma_h      = pBase.sigma_h;
R.gamma_h      = pBase.gamma_h;
R.s_daily      = s(:);
R.m_daily      = m(:);
R.rhv_daily    = tr.a(:) .* tr.c(:) .* m(:);      % human -> mosquito
R.rvh_daily    = tr.a(:) .* tr.b(:);              % mosquito -> human
R.sigma_v_daily = tr.sigma_v(:);
R.mu_daily     = mu(:);
R.mu_raw_daily = muRaw(:);

% R0(t) with the declared convention: continuous (raw) mortality, sqrt form,
% N_v supply calibre.  mwspm_r0 recomputes the traits from T itself.
R0d = mwspm_r0(D.T_daily, m(:), muRaw(:), pBase.gamma_h);
R.R0_daily = R0d(:);

% Malthusian parameter r(t) = dominant eigenvalue of the FROZEN mean matrix.
% This is the quantity the single-type reduction (Carmona & Gandon 2020 Eq.2)
% is matched on, so it is computed here and nowhere else.
nd = R.n_days;
r = zeros(nd, 1);
for k = 1:nd
    A = zeros(4, 4);
    sv = R.sigma_v_daily(k);
    muK = R.mu_daily(k);
    A(1, 1) = -R.sigma_h;  A(1, 4) = R.rvh_daily(k);
    A(2, 1) =  R.sigma_h;  A(2, 2) = -R.gamma_h;
    A(3, 2) =  R.rhv_daily(k); A(3, 3) = -(sv + muK);
    A(4, 3) =  sv;         A(4, 4) = -muK;
    ev = eig(A);
    r(k) = max(real(ev));
end
R.r_daily = r;
end
