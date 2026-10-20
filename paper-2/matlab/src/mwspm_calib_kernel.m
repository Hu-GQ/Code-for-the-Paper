function out = mwspm_calib_kernel(pre, th, opts)
%MWSPM_CALIB_KERNEL  Forward map of the calibration audit: parameters ->
%   weekly reported counts, for a bounded number of likelihood evaluations.
%
%   out = mwspm_calib_kernel(pre, th)
%   out = mwspm_calib_kernel(pre, th, struct('sub_weeks', 1:37))
%
%   th fields (observation-consistent parameterisation of T12):
%     m0, rho, d_weeks, onset_doy, offset_doy, mu_s, gamma, beta
%   optional (initial-value / importation uncertainty):
%     I_h0        initial infectious humans (subtracted from S_h)
%     E_v0, I_v0  initial exposed / infectious mosquitoes
%     lambda_scale  multiplicative factor on the importation forcing
%                   (1.0 = imported cases are taken at face value, which
%                    ASSUMES they are fully reported; <1 models their
%                    under-reporting, which is the same rho-scale critique)
%
%   opts.sub_weeks  subset of week indices used to score (default all).  The
%                   integration is ALWAYS run on the full window: holding out
%                   later weeks is a held-out SCORING block, not a shorter
%                   epidemic, so no information from the held-out block enters
%                   the trajectory (see the leakage note in the T12 report).
%
%   opts.y_obs      optional replacement for the observed weekly series.  Used
%                   by the synthetic-recovery block, which must score against
%                   SIMULATED counts rather than the Guangzhou observations.
%                   The default is pre.Y, i.e. the real data.
%   opts.solver     'euler' (historical default), 'rk4', or independent 'ode15s'.
%   opts.dt         RK4 step in days, default 0.025; must divide one day.
%                   Stable paths use daily ZOH inputs without state clipping.
%   opts.fast       opt-in inlined RK4 stages (default false).  Mathematically
%                   identical to the dailyRhs form, used to make long profile
%                   sweeps affordable; the caller must still cross-check one
%                   point against the default path and ode15s.
%
%   Relation to the canonical path
%   ------------------------------
%   The algebra below mirrors mwspm_rhs.m line by line (same equations, same
%   state order, standard incidence with denominator N_h, Euler dt = 1 day with
%   the source-faithful max(.,0) clamp).  The only differences are mechanical:
%   the forcing coefficients are read from the precomputed daily vectors in
%   `pre` instead of being rebuilt per step, and the zeros(8,1) allocation per
%   step is replaced by scalars.  exp_t12_calibration MEASURES the equivalence
%   against mwspm_rhs + mwspm_integrate_fixed + mwspm_observe and records the
%   maximum discrepancy; the number is evidence, not an assumption.

if nargin < 3 || isempty(opts)
    opts = struct();
end
solver = 'euler';
if isfield(opts, 'solver'), solver = char(opts.solver); end
assert(any(strcmp(solver, {'euler','rk4','ode15s'})), 'mwspm:calib', ...
    'unknown calibration solver "%s"', solver);
useFast = strcmp(solver, 'rk4') && isfield(opts, 'fast') && ~isempty(opts.fast) ...
    && logical(opts.fast);

nD = pre.n_days;
nW = pre.n_weeks;

% ---- season factor: single source is mwspm_season -------------------------
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
            'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
            's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
            'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);

m0   = th.m0;
mu_s = th.mu_s;
Nv   = m0 * pre.N_h * s;                 % external supply N_v(t)
mue  = mu_s * pre.mu_variant_daily;      % mu_eff(t)

a = pre.a; b = pre.b; c = pre.c; sv = pre.sigma_v;
Nh = pre.N_h; sh = pre.sigma_h; gh = pre.gamma_h;
lamScale = 1.0;
if isfield(th, 'lambda_scale') && ~isempty(th.lambda_scale)
    lamScale = double(th.lambda_scale);
end
Lam = lamScale * pre.Lambda;

% ---- initial condition (mwspm_initial_state + optional infection seeds) ---
Ih0 = 0;
if isfield(th, 'I_h0') && ~isempty(th.I_h0), Ih0 = double(th.I_h0); end
Ev0 = 0;
if isfield(th, 'E_v0') && ~isempty(th.E_v0), Ev0 = double(th.E_v0); end
Iv0 = 0;
if isfield(th, 'I_v0') && ~isempty(th.I_v0), Iv0 = double(th.I_v0); end

Sh = Nh - Ih0; Eh = 0; Ih = Ih0; Rh = 0;
Sv = Nv(1);
Ev = Ev0; Iv = Iv0; Cm = 0;

Cum = zeros(nD + 1, 1);
Cum(1) = Cm;

if strcmp(solver, 'euler')
for k = 1:nD
    lam_h = a(k) * b(k) * (Iv / Nh);
    lam_v = a(k) * c(k) * (Ih / Nh);
    mk    = mue(k);

    dSh = -lam_h * Sh;
    dEh =  lam_h * Sh - sh * Eh;
    dIh =  sh * Eh - gh * Ih + Lam(k);
    dRh =  gh * Ih;
    dSv =  mk * (Nv(k) - Sv) - lam_v * Sv;
    dEv =  lam_v * Sv - (mk + sv(k)) * Ev;
    dIv =  sv(k) * Ev - mk * Iv;
    dCm =  sh * Eh;

    Sh = max(Sh + dSh, 0);
    Eh = max(Eh + dEh, 0);
    Ih = max(Ih + dIh, 0);
    Rh = max(Rh + dRh, 0);
    Sv = max(Sv + dSv, 0);
    Ev = max(Ev + dEv, 0);
    Iv = max(Iv + dIv, 0);
    Cm = Cm + dCm;

    Cum(k + 1) = Cm;
end
minState = NaN; % Legacy clamping is preserved for exact historical replay.
else
    h = 0.025;
    if isfield(opts, 'dt'), h = double(opts.dt); end
    assert(isscalar(h) && isfinite(h) && h > 0 && h <= 1, 'mwspm:calib', ...
        'dt must be a finite scalar in (0,1] days');
    nSub = round(1 / h);
    assert(abs(nSub * h - 1) < 1e-12, 'mwspm:calib', 'dt must divide one day');
    x = [Sh; Eh; Ih; Rh; Sv; Ev; Iv; Cm];
    minState = min(x);
    odeOpt = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.125);
    for k = 1:nD
        % Daily ZOH coefficients stay fixed at all RK stages and solver endpoints.
        rates = [a(k)*b(k)/Nh, a(k)*c(k)/Nh, sh, gh, mue(k), sv(k), Nv(k), Lam(k)];
        if strcmp(solver, 'rk4') && useFast
            r1 = rates(1); r2 = rates(2); r3 = rates(3); r4 = rates(4);
            r5 = rates(5); r6 = rates(6); r7 = rates(7); r8 = rates(8);
            for j = 1:nSub
                x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
                x5 = x(5); x6 = x(6); x7 = x(7); x8 = x(8);
                % stage 1
                lh = r1*x7*x1; lv = r2*x3*x5;
                y1 = -lh;
                y2 = lh - r3*x2;
                y3 = r3*x2 - r4*x3 + r8;
                y4 = r4*x3;
                y5 = r5*(r7 - x5) - lv;
                y6 = lv - (r5 + r6)*x6;
                y7 = r6*x6 - r5*x7;
                y8 = r3*x2;
                k1_1 = y1;
                k1_2 = y2;
                k1_3 = y3;
                k1_4 = y4;
                k1_5 = y5;
                k1_6 = y6;
                k1_7 = y7;
                k1_8 = y8;
                % stage 2
                x1 = x(1) + h*k1_1/2;
                x2 = x(2) + h*k1_2/2;
                x3 = x(3) + h*k1_3/2;
                x4 = x(4) + h*k1_4/2;
                x5 = x(5) + h*k1_5/2;
                x6 = x(6) + h*k1_6/2;
                x7 = x(7) + h*k1_7/2;
                x8 = x(8) + h*k1_8/2;
                lh = r1*x7*x1; lv = r2*x3*x5;
                y1 = -lh;
                y2 = lh - r3*x2;
                y3 = r3*x2 - r4*x3 + r8;
                y4 = r4*x3;
                y5 = r5*(r7 - x5) - lv;
                y6 = lv - (r5 + r6)*x6;
                y7 = r6*x6 - r5*x7;
                y8 = r3*x2;
                k2_1 = y1;
                k2_2 = y2;
                k2_3 = y3;
                k2_4 = y4;
                k2_5 = y5;
                k2_6 = y6;
                k2_7 = y7;
                k2_8 = y8;
                % stage 3
                x1 = x(1) + h*k2_1/2;
                x2 = x(2) + h*k2_2/2;
                x3 = x(3) + h*k2_3/2;
                x4 = x(4) + h*k2_4/2;
                x5 = x(5) + h*k2_5/2;
                x6 = x(6) + h*k2_6/2;
                x7 = x(7) + h*k2_7/2;
                x8 = x(8) + h*k2_8/2;
                lh = r1*x7*x1; lv = r2*x3*x5;
                y1 = -lh;
                y2 = lh - r3*x2;
                y3 = r3*x2 - r4*x3 + r8;
                y4 = r4*x3;
                y5 = r5*(r7 - x5) - lv;
                y6 = lv - (r5 + r6)*x6;
                y7 = r6*x6 - r5*x7;
                y8 = r3*x2;
                k3_1 = y1;
                k3_2 = y2;
                k3_3 = y3;
                k3_4 = y4;
                k3_5 = y5;
                k3_6 = y6;
                k3_7 = y7;
                k3_8 = y8;
                % stage 4
                x1 = x(1) + h*k3_1;
                x2 = x(2) + h*k3_2;
                x3 = x(3) + h*k3_3;
                x4 = x(4) + h*k3_4;
                x5 = x(5) + h*k3_5;
                x6 = x(6) + h*k3_6;
                x7 = x(7) + h*k3_7;
                x8 = x(8) + h*k3_8;
                lh = r1*x7*x1; lv = r2*x3*x5;
                y1 = -lh;
                y2 = lh - r3*x2;
                y3 = r3*x2 - r4*x3 + r8;
                y4 = r4*x3;
                y5 = r5*(r7 - x5) - lv;
                y6 = lv - (r5 + r6)*x6;
                y7 = r6*x6 - r5*x7;
                y8 = r3*x2;
                k4_1 = y1;
                k4_2 = y2;
                k4_3 = y3;
                k4_4 = y4;
                k4_5 = y5;
                k4_6 = y6;
                k4_7 = y7;
                k4_8 = y8;
                % stage combination
                x1 = x(1) + h*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
                x2 = x(2) + h*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
                x3 = x(3) + h*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3)/6;
                x4 = x(4) + h*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4)/6;
                x5 = x(5) + h*(k1_5 + 2*k2_5 + 2*k3_5 + k4_5)/6;
                x6 = x(6) + h*(k1_6 + 2*k2_6 + 2*k3_6 + k4_6)/6;
                x7 = x(7) + h*(k1_7 + 2*k2_7 + 2*k3_7 + k4_7)/6;
                x8 = x(8) + h*(k1_8 + 2*k2_8 + 2*k3_8 + k4_8)/6;
                x = [x1; x2; x3; x4; x5; x6; x7; x8];
                mn = x1;
                if x2 < mn, mn = x2; end
                if x3 < mn, mn = x3; end
                if x4 < mn, mn = x4; end
                if x5 < mn, mn = x5; end
                if x6 < mn, mn = x6; end
                if x7 < mn, mn = x7; end
                if x8 < mn, mn = x8; end
                minState = min(minState, mn);
            end
        elseif strcmp(solver, 'rk4')
            for j = 1:nSub
                k1 = dailyRhs(x, rates);
                k2 = dailyRhs(x + h*k1/2, rates);
                k3 = dailyRhs(x + h*k2/2, rates);
                k4 = dailyRhs(x + h*k3, rates);
                x = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
                minState = min(minState, min(x));
            end
        else
            [~, values] = ode15s(@(~, z) dailyRhs(z, rates), [0 1], x, odeOpt);
            minState = min(minState, min(values(:)));
            x = values(end,:)';
        end
        assert(all(isfinite(x)) && minState >= -1e-8, 'mwspm:calib:unstable', ...
            'nonfinite or negative state in unclamped integration');
        Cum(k + 1) = x(8);
    end
    Sh=x(1); Eh=x(2); Ih=x(3); Rh=x(4); Sv=x(5); Ev=x(6); Iv=x(7); Cm=x(8);
end

% ---- weekly aggregation (exact: the accumulator is part of the state) -----
% Cum holds the cumulative true incidence on every integer day 0..n_days, so
% the weekly total is the difference of the accumulator at the weekly EDGES
% t = 0, 7, ..., 7*n_weeks (identical to mwspm_observe, which interpolates the
% same edges).  Cum has n_days+1 entries indexed 1..n_days+1 for t = 0..n_days.
edgeIdx = (0:nW)' * 7 + 1;
dC = diff(Cum(edgeIdx));

d = round(th.d_weeks);
d = min(max(d, 0), nW);
if d > 0
    dCd = [zeros(d, 1); dC(1:nW - d)];
else
    dCd = dC;
end
Yhat = th.rho * dCd;

% ---- scoring subset -------------------------------------------------------
sub = (1:nW)';
if isfield(opts, 'sub_weeks') && ~isempty(opts.sub_weeks)
    sub = double(opts.sub_weeks(:));
end
yAll = pre.Y;
if isfield(opts, 'y_obs') && ~isempty(opts.y_obs)
    yAll = double(opts.y_obs(:));
    assert(numel(yAll) == nW, 'mwspm:calib', ...
        'y_obs has %d entries but the window has %d weeks', numel(yAll), nW);
end

out = struct();
out.Cum_full      = Cum;
out.dC            = dC;
out.dC_delayed    = dCd;
out.Yhat          = Yhat;
out.sub_weeks     = sub;
out.Yhat_sub      = Yhat(sub);
out.Y_obs_sub     = yAll(sub);
out.true_total    = sum(dC);
out.reported_total = sum(Yhat);
out.peak_week     = find(Yhat == max(Yhat), 1);
out.d_weeks_used  = d;
out.s_daily       = s;
out.Nv_daily      = Nv;
out.Iv_final      = Iv;
out.x_final       = [Sh; Eh; Ih; Rh; Sv; Ev; Iv; Cm];
out.solver        = solver;
if strcmp(solver,'rk4') && useFast, out.integrator_path = 'rk4-inline';
else, out.integrator_path = solver; end
out.min_state     = minState;
end

function dx = dailyRhs(x, r)
lh = r(1)*x(7)*x(1);
lv = r(2)*x(3)*x(5);
dx = [-lh; lh-r(3)*x(2); r(3)*x(2)-r(4)*x(3)+r(8); r(4)*x(3); ...
      r(5)*(r(7)-x(5))-lv; lv-(r(5)+r(6))*x(6); ...
      r(6)*x(6)-r(5)*x(7); r(3)*x(2)];
end
