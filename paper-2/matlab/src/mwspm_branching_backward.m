function out = mwspm_branching_backward(R, opts)
%MWSPM_BRANCHING_BACKWARD  Exact extinction / no-detection probabilities.
%
%   out = mwspm_branching_backward(R, opts)
%
%   Backward Kolmogorov equations of the time-inhomogeneous 4-type branching
%   process whose rates come from mwspm_branching_rates.  This is the EXACT
%   answer for the estimands below; the stochastic simulator is validated
%   against it, never the other way round.
%
%   Two independent backward systems are integrated:
%
%     q_i(t) = P(no particle alive at t_end | one type-i particle at t)
%        dq_1/dt = sigma_h (q_1 - q_2)
%        dq_2/dt = gamma_h (q_2 - 1) + rhv(t) (q_2 - q_2 q_3)
%        dq_3/dt = sigma_v(t) (q_3 - q_4) + mu(t) (q_3 - 1)
%        dq_4/dt = mu(t) (q_4 - 1) + rvh(t) (q_4 - q_4 q_1)
%        terminal q(t_end) = 0        (a live particle is not extinct)
%
%     u_i(t) = P(no DETECTED autochthonous human infection is ever produced
%                 | one type-i particle at t), each new human infection being
%                 detected independently with probability rho
%        du_1/dt = sigma_h (u_1 - u_2)
%        du_2/dt = gamma_h (u_2 - 1) + rhv(t) (u_2 - u_2 u_3)
%        du_3/dt = sigma_v(t) (u_3 - u_4) + mu(t) (u_3 - 1)
%        du_4/dt = mu(t) (u_4 - 1) + rvh(t) (u_4 - (1-rho) u_1 u_4)
%        terminal u(t_end) = 1
%
%   Both are ODEs in the probability simplex: the solution is a probability
%   vector, so q, u stay in [0,1] and the backward integration is stable near
%   the relevant fixed point (1/R0 in the supercritical case, 1 in the
%   subcritical case).  A max-violation of [0,1] is still measured and
%   reported; it is never assumed.
%
%   opts fields
%     rho         vector of detection fractions (default: [])
%     t_eval      times (days elapsed) at which to report (default: week starts)
%     static      if true, freeze the rates at each t_eval (quasi-static
%                 variant used to quantify obligation O8); default false
%     max_step    solver MaxStep in days (default 0.1 for the time-varying
%                 case, 5 for the frozen case)
%
%   out fields
%     t_eval, rho, q (nEval x 4), u (nEval x 4 x nRho),
%     max_violation, n_solves, solver

if nargin < 2 || isempty(opts)
    opts = struct();
end
rhoVec = [];
if isfield(opts, 'rho') && ~isempty(opts.rho)
    rhoVec = double(opts.rho(:)');
end
static = false;
if isfield(opts, 'static') && ~isempty(opts.static)
    static = logical(opts.static);
end
nd = R.n_days;
tEnd = double(nd);
if isfield(opts, 't_eval') && ~isempty(opts.t_eval)
    tEval = double(opts.t_eval(:));
else
    tEval = (0 : 7 : 7 * (R.n_weeks - 1))';
end
maxStep = 0.1;
if isfield(opts, 'max_step') && ~isempty(opts.max_step)
    maxStep = double(opts.max_step);
end

nR = numel(rhoVec);
nE = numel(tEval);
Q = nan(nE, 4);
U = nan(nE, 4, nR);
nSolves = 0;
worst = 0;

odeOpt = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', maxStep, ...
    'NormControl', 'off');

if ~static
    % ---------------------------------------------------------- time-varying
    y0 = zeros(4 * (1 + nR), 1);
    for j = 1:nR
        y0(4 + (j - 1) * 4 + (1:4)) = 1;      % u terminal = 1
    end
    solTV = ode15s(@(t, y) rhsTV(t, y, R, rhoVec), [tEnd, 0], y0, odeOpt);
    nSolves = nSolves + 1;
    Y = deval(solTV, tEval)';
    for i = 1:nE
        Q(i, :) = Y(i, 1:4);
        for j = 1:nR
            U(i, :, j) = Y(i, 4 + (j - 1) * 4 + (1:4));
        end
    end
else
    % ------------------------------------------------- frozen (quasi-static)
    for i = 1:nE
        Rf = freezeAt(R, tEval(i));
        y0 = zeros(4 * (1 + nR), 1);
        for j = 1:nR
            y0(4 + (j - 1) * 4 + (1:4)) = 1;
        end
        span = [tEnd, tEval(i)];
        sol = ode15s(@(t, y) rhsStatic(y, Rf, rhoVec), span, y0, odeOpt);
        nSolves = nSolves + 1;
        Y = deval(sol, tEval(i))';
        Q(i, :) = Y(1:4);
        for j = 1:nR
            U(i, :, j) = Y(4 + (j - 1) * 4 + (1:4));
        end
    end
end

% honesty check: the solution is a probability vector, so clamp violations are
% a numerical defect, not a modelling outcome.  They are MEASURED here.
allVals = [Q(:); U(:)];
allVals = allVals(isfinite(allVals));
worst = max([0, max(allVals - 1), max(-allVals)]);

out = struct();
out.t_eval = tEval;
out.rho = rhoVec;
out.q = Q;
out.u = U;
out.max_violation_of_unit_interval = worst;
out.n_solves = nSolves;
out.solver = 'ode15s';
out.static = static;
out.max_step_days = maxStep;
out.seed_type = 'I_h (type 2): one imported infectious human';
out.note = ['q_i and u_i are probabilities; max_violation_of_unit_interval is a ' ...
    'MEASURED numerical residual and must be ~0.'];
end

% ---------------------------------------------------------------------------
function dy = rhsTV(t, y, R, rhoVec)
f = sampleDay(R, t);
dy = coreRHS(y, f, rhoVec);
end

% ---------------------------------------------------------------------------
function dy = rhsStatic(y, Rf, rhoVec)
f = struct('sh', Rf.sigma_h, 'gh', Rf.gamma_h, ...
    'rhv', Rf.rhv_daily(1), 'sv', Rf.sigma_v_daily(1), ...
    'mu', Rf.mu_daily(1), 'rvh', Rf.rvh_daily(1));
dy = coreRHS(y, f, rhoVec);
end

% ---------------------------------------------------------------------------
function dy = coreRHS(y, f, rhoVec)
n = 4 * (1 + numel(rhoVec));
dy = zeros(n, 1);

q1 = y(1); q2 = y(2); q3 = y(3); q4 = y(4);
dy(1) = f.sh * (q1 - q2);
dy(2) = f.gh * (q2 - 1) + f.rhv * (q2 - q2 * q3);
dy(3) = f.sv * (q3 - q4) + f.mu * (q3 - 1);
dy(4) = f.mu * (q4 - 1) + f.rvh * (q4 - q4 * q1);

for j = 1:numel(rhoVec)
    b = 4 + (j - 1) * 4;
    u1 = y(b + 1); u2 = y(b + 2); u3 = y(b + 3); u4 = y(b + 4);
    keep = 1 - rhoVec(j);
    dy(b + 1) = f.sh * (u1 - u2);
    dy(b + 2) = f.gh * (u2 - 1) + f.rhv * (u2 - u2 * u3);
    dy(b + 3) = f.sv * (u3 - u4) + f.mu * (u3 - 1);
    dy(b + 4) = f.mu * (u4 - 1) + f.rvh * (u4 - keep * u1 * u4);
end
end

% ---------------------------------------------------------------------------
function f = sampleDay(R, t)
idx = min(max(floor(double(t)) + 1, 1), R.n_days);
f = struct('sh', R.sigma_h, 'gh', R.gamma_h, ...
    'rhv', R.rhv_daily(idx), 'sv', R.sigma_v_daily(idx), ...
    'mu', R.mu_daily(idx), 'rvh', R.rvh_daily(idx));
end

% ---------------------------------------------------------------------------
function Rf = freezeAt(R, t)
idx = min(max(floor(double(t)) + 1, 1), R.n_days);
Rf = struct('sigma_h', R.sigma_h, 'gamma_h', R.gamma_h, ...
    'rhv_daily', R.rhv_daily(idx), 'rvh_daily', R.rvh_daily(idx), ...
    'sigma_v_daily', R.sigma_v_daily(idx), 'mu_daily', R.mu_daily(idx));
end
