function S = mwspm_integrate_fixed(rhs, t0, tf, x0, spec, reportGrid, events)
%MWSPM_INTEGRATE_FIXED  Fixed-step integrators (explicit Euler, classic RK4).
%
%   S = mwspm_integrate_fixed(rhs, t0, tf, x0, spec, reportGrid, events)
%
%   spec.method             'euler' | 'rk4'
%   spec.dt_days            INTERNAL step.  This IS the integration step - it
%                           is never derived from the reporting grid.
%   spec.clamp_nonnegative  logical; source-faithful max(.,0) clamping
%   reportGrid              times at which the solution is reported
%   events                  struct array with fields name, kind, state_index,
%                           direction, threshold
%
%   Termination reasons: 't_end_reached' | 'event:<name>' | 'nonfinite_state'.
%
%   The internal mesh is kept in full (t_internal/x_internal) so that the
%   number of steps and the step-size distribution can be audited later; the
%   reporting grid is applied afterwards by interpolation and never influences
%   the integration.

if nargin < 7
    events = struct('name', {}, 'kind', {}, 'state_index', {}, 'direction', {}, 'threshold', {});
end

method = lower(char(spec.method));
assert(any(strcmp(method, {'euler', 'rk4'})), 'mwspm:fixed', 'unsupported method "%s"', method);
dt = double(spec.dt_days);
assert(isfinite(dt) && dt > 0, 'mwspm:fixed', 'dt_days must be positive');

clamp = false;
if isfield(spec, 'clamp_nonnegative') && ~isempty(spec.clamp_nonnegative)
    clamp = logical(spec.clamp_nonnegative);
end

nState  = numel(x0);
nSteps  = ceil((tf - t0) / dt);

% Internal times are built by MULTIPLICATION, never by accumulating
% t = t + min(dt, tf - t).  Accumulation drifts: over 13161 steps of dt = 0.025
% the running total ended at 328.999999999975, i.e. BELOW tf, so the last
% reporting grid point (t = 329) lay beyond the reached horizon and interp1
% returned NaN there while termination_reason still said t_end_reached.  The
% endpoint is now pinned to tf exactly.
tNom = t0 + (0:nSteps)' * dt;
tNom = min(tNom, tf);
tNom(end) = tf;

tInt    = zeros(nSteps + 1, 1);
Xint    = zeros(nSteps + 1, nState);
tInt(1) = t0;
Xint(1, :) = x0(:)';

t = t0;
x = x0(:);
nRhs = 0;
termination = 't_end_reached';
eventsHit = struct('name', {}, 'time', [], 'index', []);

for k = 1:nSteps
    tNext = tNom(k + 1);
    h = tNext - t;
    if h <= 0
        % The tail of tNom collapsed onto tf (dt did not divide tf - t0): drop
        % the duplicate node.  tInt(end) is already exactly tf.
        tInt = tInt(1:k);
        Xint = Xint(1:k, :);
        break
    end

    f1 = rhs(t, x); nRhs = nRhs + 1;
    switch method
        case 'euler'
            xn = x + h * f1;
        case 'rk4'
            f2 = rhs(t + h/2, x + (h/2) * f1); nRhs = nRhs + 1;
            f3 = rhs(t + h/2, x + (h/2) * f2); nRhs = nRhs + 1;
            f4 = rhs(t + h,   x + h     * f3); nRhs = nRhs + 1;
            xn = x + (h/6) * (f1 + 2*f2 + 2*f3 + f4);
    end

    if clamp
        xn = max(xn, 0);
    end

    t = tNext;
    x = xn;
    tInt(k + 1) = t;
    Xint(k + 1, :) = x';

    if ~all(isfinite(x))
        termination = 'nonfinite_state';
        tInt = tInt(1:k+1);
        Xint = Xint(1:k+1, :);
        break
    end

    hit = false;
    for e = 1:numel(events)
        ev = events(e);
        if ~strcmp(ev.kind, 'threshold')
            continue
        end
        idx = ev.state_index;
        prev = Xint(k,   idx) - ev.threshold;
        cur  = Xint(k+1, idx) - ev.threshold;
        crossed = false;
        if ev.direction > 0
            crossed = (prev < 0) && (cur >= 0);
        elseif ev.direction < 0
            crossed = (prev > 0) && (cur <= 0);
        else
            crossed = (prev < 0) ~= (cur < 0);
        end
        if crossed
            termination = ['event:' ev.name];
            eventsHit(end+1) = struct('name', ev.name, 'time', t, 'index', idx); %#ok<AGROW>
            hit = true;
            break
        end
    end
    if hit
        tInt = tInt(1:k+1);
        Xint = Xint(1:k+1, :);
        break
    end
end

S = struct();
S.kind        = 'fixed_step';
S.method      = method;
S.spec        = spec;
S.t_internal  = tInt;
S.x_internal  = Xint;
S.n_steps     = numel(tInt) - 1;
steps         = diff(tInt);
if isempty(steps)
    S.min_internal_step = NaN;
    S.max_internal_step = NaN;
else
    S.min_internal_step = min(steps);
    S.max_internal_step = max(steps);
end
S.n_rhs_evaluations = nRhs;
% Contract consistency with the adaptive branch of mwspm_solve: if the loop
% finished without an event but the mesh stops short of tf, reporting
% 't_end_reached' would be a false statement.
if strcmp(termination, 't_end_reached') && tInt(end) < tf - 1e-9 * max(1, abs(tf))
    termination = 'solver_terminated_early';
end
S.termination_reason = termination;
S.events_hit  = eventsHit;
S.t_end_actual = tInt(end);

rg = double(reportGrid(:));
S.t = rg;
S.x = interp1(tInt, Xint, rg, 'linear', NaN);      % NaN beyond the reached horizon
end
