function S = mwspm_solve(opts)
%MWSPM_SOLVE  Solver wrapper: stiff-aware selection, event handling, provenance.
%
%   S = mwspm_solve(opts)
%
%   opts.rhs        function handle @(t,x) -> dx/dt
%   opts.t0,tf      integration interval (days)
%   opts.x0         initial state (column)
%   opts.spec       solver specification, one entry of config solvers[]
%        .kind       'adaptive' | 'fixed_step' | 'auto'
%        .method     'auto' | ode15s | ode23s | ode23t | ode23tb | ode45 |
%                    ode113 | ode78 | ode89 | euler | rk4
%        .RelTol,.AbsTol,.NonNegative     (adaptive)
%        .dt_days,.clamp_nonnegative,.mu_variant (fixed_step)
%   opts.reportGrid report times (interpolated after the fact - the grid does
%                   NOT define the internal step)
%   opts.events     event struct array (see mwspm_integrate_fixed / README)
%   opts.stiff      precomputed mwspm_stiffness output (used by 'auto')
%   opts.pickCfg    pass-through config for mwspm_pick_solver
%   opts.logFid     optional log file id
%
%   Everything the protocol asks to be saved is returned here and written to
%   manifest.json by the caller: method, tolerances, report grid, INTERNAL
%   mesh statistics, event settings and the termination reason.

if nargin < 1 || ~isstruct(opts)
    error('mwspm:solve', 'opts (a struct) is required');
end

t0 = double(getfielddef(opts, 't0', 0));
tf = double(getfielddef(opts, 'tf', NaN));
x0 = opts.x0(:);
rhs = opts.rhs;
reportGrid = getfielddef(opts, 'reportGrid', (t0:tf)');
events = getfielddef(opts, 'events', struct('name', {}, 'kind', {}, 'state_index', {}, 'direction', {}, 'threshold', {}));
spec = opts.spec;
logFid = getfielddef(opts, 'logFid', 1);

assert(isfinite(tf) && tf > t0, 'mwspm:solve', 'invalid interval');

method = lower(char(getfielddef(spec, 'method', 'auto')));
kind   = lower(char(getfielddef(spec, 'kind', 'adaptive')));

% ------------------------------------------------------------------ 'auto'
selInfo = struct('method', method, 'rule', 'explicit method requested (no auto selection)', ...
                 'stiffness_ratio', NaN, 'is_stiff', NaN, 'method_note', '');
if strcmp(method, 'auto')
    assert(isfield(opts, 'stiff') && ~isempty(opts.stiff), 'mwspm:solve', ...
        'method ''auto'' requires opts.stiff (run mwspm_stiffness first)');
    selInfo = mwspm_pick_solver(opts.stiff, getfielddef(opts, 'pickCfg', struct()));
    method = selInfo.method;
end

fixedMethods = {'euler', 'rk4'};
if any(strcmp(method, fixedMethods))
    assert(strcmp(kind, 'fixed_step'), 'mwspm:solve', ...
        'method "%s" requires spec.kind = ''fixed_step''', method);
    spec.dt_days = getfielddef(spec, 'dt_days', NaN);
    S = mwspm_integrate_fixed(rhs, t0, tf, x0, spec, reportGrid, events);
else
    % ------------------------------------------------------- adaptive solver
    assert(any(strcmp(method, {'ode15s','ode23s','ode23t','ode23tb','ode45','ode113','ode78','ode89'})), ...
        'mwspm:solve', 'unsupported method "%s"', method);

    odeOpts = odeset('RelTol', getfielddef(spec, 'RelTol', 1e-8), ...
                     'AbsTol', getfielddef(spec, 'AbsTol', 1e-10));
    % MaxStep is opt-in and NO caller before 2026-10-06 set it, so adding this
    % line cannot change any previously recorded result.  It became necessary
    % because the temperature / importation forcing is piecewise constant on a
    % DAILY grid: without a cap the adaptive solver steps straight over the
    % jumps.  Measured failure (T11.01, baseline model): ode15s at
    % RelTol 1e-10 / AbsTol 1e-12 returned only 18 internal steps with steps up
    % to 32.9 days and a cumulative incidence of exactly zero, while the same
    % model integrated with dt = 0.01 gave ~380 infections.  A reference that
    % resolves nothing is not a reference.
    mx = getfielddef(spec, 'MaxStep', []);
    if ~isempty(mx) && isfinite(mx) && mx > 0
        odeOpts = odeset(odeOpts, 'MaxStep', double(mx));
    end
    nn = getfielddef(spec, 'NonNegative', []);
    if ~isempty(nn)
        odeOpts = odeset(odeOpts, 'NonNegative', nn);
    end

    evList = events;
    useEvents = ~isempty(evList);
    if useEvents
        for e = 1:numel(evList)
            assert(strcmp(evList(e).kind, 'zero_crossing'), 'mwspm:solve', ...
                ['adaptive solvers only support zero_crossing events; "' evList(e).name ...
                 '" is a "%s" event - use a fixed-step integrator for threshold events'], evList(e).kind);
        end
        odeOpts = odeset(odeOpts, 'Events', @(t, y) mwspm_event_fcn(t, y, evList));
    end

    % NOTE: the tspan passed to the solver is the TWO-element interval.  The
    % report grid is applied afterwards with deval, so the adaptive internal
    % step is never tied to the reporting density.
    %
    % The single-output form is used deliberately: the solution struct carries
    % the adaptive mesh (sol.x/sol.y) AND, when an Events option was supplied,
    % the event log (sol.xe/sol.ye/sol.ie).  The multi-output form
    % [t,y,te,ye,ie] must not be mixed with sol, because the return signature
    % differs between solvers (ode15s rejects a 6-output request).
    lastwarn('');
    sol = feval(method, rhs, [t0 tf], x0, odeOpts);
    [wmsg, wid] = lastwarn;

    % sol.ie / sol.xe exist only when 'Events' was set; keep them empty
    % otherwise so the downstream logic stays branch-free.
    if useEvents && isfield(sol, 'ie')
        ie = sol.ie(:);
        xe = sol.xe(:);
    else
        ie = [];
        xe = [];
    end

    S = struct();
    S.kind            = 'adaptive';
    S.method          = method;
    S.spec            = spec;
    S.t_internal      = sol.x(:);
    S.x_internal      = sol.y';
    S.n_steps         = numel(sol.x) - 1;
    steps             = diff(sol.x);
    if isempty(steps)
        S.min_internal_step = NaN;
        S.max_internal_step = NaN;
    else
        S.min_internal_step = min(steps);
        S.max_internal_step = max(steps);
    end
    S.t_end_actual    = sol.x(end);
    if S.t_end_actual < tf - 1e-12 && isempty(ie)
        S.termination_reason = 'solver_terminated_early';
    else
        S.termination_reason = 't_end_reached';
    end
    S.events_hit = struct('name', {}, 'time', [], 'index', []);
    if ~isempty(ie)
        for k = 1:numel(ie)
            S.events_hit(end+1) = struct('name', evList(ie(k)).name, 'time', xe(k), 'index', ie(k)); %#ok<AGROW>
        end
        S.termination_reason = ['event:' evList(ie(1)).name];
    end

    rg = double(reportGrid(:));
    Xr = nan(numel(rg), numel(x0));
    inRange = rg <= S.t_end_actual + 1e-12;
    if any(inRange)
        Xr(inRange, :) = deval(sol, rg(inRange))';
    end
    S.t = rg;
    S.x = Xr;

    if ~isempty(wid)
        S.solver_warning_id = wid;
        S.solver_warning_message = wmsg;
        fprintf(logFid, '[warn] solver %s: %s (%s)\n', method, wmsg, wid);
    else
        S.solver_warning_id = '';
        S.solver_warning_message = '';
    end
end

S.selection = selInfo;
if isfield(spec, 'mu_variant')
    S.mu_variant = spec.mu_variant;
else
    S.mu_variant = 'raw';
end
if strcmp(method, 'ode15s')
    S.solver_realisation = ['EXPERIMENT_PLAN names "LSODA"; MATLAB does not ship LSODA, so ' ...
        'ode15s (variable-order BDF, stiff) is used and this substitution is recorded here.'];
else
    S.solver_realisation = '';
end
end

% ---------------------------------------------------------------------------
function [value, isterminal, direction] = mwspm_event_fcn(~, y, evList)
value = zeros(numel(evList), 1);
isterminal = ones(numel(evList), 1);
direction = zeros(numel(evList), 1);
for k = 1:numel(evList)
    value(k) = y(evList(k).state_index) - evList(k).threshold;
    direction(k) = evList(k).direction;
end
end

% ---------------------------------------------------------------------------
function v = getfielddef(s, f, dflt)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
