function sel = mwspm_pick_solver(stiff, cfg)
%MWSPM_PICK_SOLVER  Choose a solver from the measured stiffness (rule, not taste).
%
%   sel = mwspm_pick_solver(stiff, cfg)
%
%   The rule is explicit and logged, so that "solver was chosen according to
%   stiffness and model type" is an auditable statement rather than a statement:
%
%     stiffness_ratio > cfg.stiffness_ratio_threshold   ->  stiff method
%     otherwise                                         ->  non-stiff method
%
%   For this project the forced coefficients already span a wide rate range
%   (mu_eff up to ~12/day versus sigma_h ~0.17/day), i.e. ratio ~ 10^2-10^3,
%   so the default best performs as a stiff run.  A non-stiff method is still
%   available and is used automatically when a model variant is genuinely
%   non-stiff.
%
%   Model-type note: the forcing is zero-order held on a daily grid, i.e. it
%   has jump discontinuities.  Stiff BDF solvers (ode15s) handle such jumps
%   without reducing the order adversarially, which is a second, independent
%   reason to prefer them for the production runs.

if nargin < 2 || isempty(cfg)
    cfg = struct();
end
thr = 1e3;
if isfield(cfg, 'stiffness_ratio_threshold') && ~isempty(cfg.stiffness_ratio_threshold)
    thr = double(cfg.stiffness_ratio_threshold);
end
stiffMethod = 'ode15s';
if isfield(cfg, 'stiff_method') && ~isempty(cfg.stiff_method)
    stiffMethod = char(cfg.stiff_method);
end
nonStiffMethod = 'ode45';
if isfield(cfg, 'nonstiff_method') && ~isempty(cfg.nonstiff_method)
    nonStiffMethod = char(cfg.nonstiff_method);
end

sel = struct();
sel.stiffness_ratio = stiff.stiffness_ratio;
sel.threshold       = thr;
sel.rule            = sprintf('stiffness_ratio %.3g %s %.3g -> %s method', ...
    stiff.stiffness_ratio, ternary(stiff.stiffness_ratio > thr, '>', '<='), thr, ...
    ternary(stiff.stiffness_ratio > thr, 'stiff', 'non-stiff'));
if stiff.stiffness_ratio > thr
    sel.method = stiffMethod;
    sel.is_stiff = true;
else
    sel.method = nonStiffMethod;
    sel.is_stiff = false;
end
sel.method_note = ['MATLAB does not ship LSODA; the plan''s "LSODA" is realised as ' sel.method ...
    ' and this substitution is recorded in manifest.json (solver_realisation).'];
end

% ---------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
