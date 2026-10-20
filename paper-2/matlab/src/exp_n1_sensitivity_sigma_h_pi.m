function out = exp_n1_sensitivity_sigma_h_pi(ctx)
%EXP_N1_SENSITIVITY_SIGMA_H_PI  N1 -- sensitivity of the warning week tau*
%   to sigma_h (human incubation rate) and to pi (phase-soft-penalty weight).
%
%   Preregistered: protocols/decisions/2026-10-10_N1_sensitivity_preregistration.md
%   (D2026-10-10-01).  Two one-dimensional refit scans, both at integer lag d=1:
%
%     Axis A: sigma_h in {1/4, 1/5, 1/5.9, 1/8} per day  (IIP 4-8 d, pi=300 fixed)
%     Axis B: pi      in {0, 30, 100, 300, 1000, 3000}   (sigma_h = 1/5.9 fixed)
%
%   Each variant refits theta = (log m0, rho, gamma, beta, t_on, t_off, mu_s)
%   with the SAME objective as T12.01 (Poisson nll + phase soft penalty), then
%   recomputes tau* in the gate-G1 caliber (identical to exp_c4_uncertainty's
%   unc_tau_star):  tau* = min{ w : cum_w >= T0(w) },  T0(w) = rho / log R0(w)
%   for R0(w) > 1, cum_w = cumulative DETECTED model cases, weekly R0 = MEAN of
%   daily R0 (continuous 'raw' mortality).  tau* is reported at BOTH lag 1 and
%   lag 2 (the second lag is the same refit theta re-evaluated with d_weeks=2).
%
%   Honest constraints (config "honesty"):
%     * sigma_h changes phase, but the pi penalty pins the peak at week 37, so
%       the phase effect of sigma_h can be partially absorbed; Axis A is thus a
%       CONDITIONAL sensitivity at pi=300, not an unconditional sigma_h claim.
%     * the pi=0 branch is reported SEPARATELY (optimiser unconstrained), never
%       merged into a single "robust" narrative.
%     * no new data evidence; the T12.01 main report point is unchanged.
%     * a sensitive result (tau* spread >= 2 weeks) is reported as-is.

logFid = ctx.logFid;
root   = ctx.root;
modelCard = ctx.modelCard;
cfg    = ctx.cfg;

fprintf(logFid, '=== exp_n1_sensitivity_sigma_h_pi (N1: sigma_h & pi sensitivity of tau*) ===\n');
fprintf(logFid, 'preregistered: %s\n', char(string(getdef(cfg, 'preregistered_by', '-'))));

D = mwspm_load_data(root, modelCard);
pre0 = mwspm_calib_prepare(D, modelCard, struct('mu_variant', 'day'));
box  = mwspm_calib_theta('box');

% ------------------------------------------------------------- parent fit ---
fitRun = n1_find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 'n1:noParent', ...
    ['no SUCCESS run of exp_c4_rho_fit with data/fit_best.json found under results/runs/. ' ...
     'N1 must refit from the real T12.01 main report point; it may not invent one.']);
fit = mwspm_readjson(fitRun.path);
th0 = n1_theta_from_fit(fit.main_report_point.theta);
fprintf(logFid, 'parent fit run : %s\n', fitRun.runId);
fprintf(logFid, 'theta0: m0=%.6f rho=%.6f mu_s=%.6f onset=%.3f offset=%.3f gamma=%.4f beta=%.4f d=%d\n', ...
    th0.m0, th0.rho, th0.mu_s, th0.onset_doy, th0.offset_doy, th0.gamma, th0.beta, th0.d_weeks);

PD = cfg.parameter_domain;
sigmaHvals = double(PD.axis_A_sigma_h.values(:)');
piVals     = double(PD.axis_B_pi.values(:)');
sigmaH0    = double(PD.axis_A_sigma_h.baseline);
pi0        = double(PD.axis_B_pi.baseline);
targetPeak = double(PD.target_peak_week);

o = cfg.optimiser;
fitOpt = struct('n_screen', double(o.n_screen), 'seeds', double(o.seeds(:)'), ...
    'n_seed', double(o.n_seed), 'restarts', double(o.restarts), ...
    'maxFunEvals', double(o.maxFunEvals), 'maxIter', double(o.maxIter), ...
    'TolX', 1e-8, 'TolFun', 1e-9, 'verbose', false);

freeAll = {'m0', 'rho', 'onset_doy', 'offset_doy', 'mu_s', 'gamma', 'beta'};
warm0   = {th0};
nWeeks  = pre0.n_weeks;
kParam  = 8;

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_sigma_h', numel(sigmaHvals), 'n_pi', numel(piVals), ...
        'sigma_h_baseline', sigmaH0, 'pi_baseline', pi0, 'parent_run_id', fitRun.runId);
    fprintf(logFid, '[dry-run] %d sigma_h variants + %d pi variants; parent readable; no refit\n', ...
        numel(sigmaHvals), numel(piVals));
    return
end

% ================================================================ Axis A ----
fprintf(logFid, '\n--- Axis A: sigma_h sensitivity (pi = %g fixed, d = 1) ---\n', pi0);
rowsA = {};
tauA  = nan(numel(sigmaHvals), 2);
for i = 1:numel(sigmaHvals)
    sh  = sigmaHvals(i);
    pre = pre0; pre.sigma_h = sh;
    C = n1_baseC(freeAll, box, pi0, targetPeak, 1);
    b  = mwspm_calib_fit(pre, C, n1_merge(fitOpt, struct('warm', {warm0})));
    th = b.theta;
    bic = 2 * b.nll_full + kParam * log(nWeeks);
    [t1, ~, ~, ~] = n1_tau_star(pre, D, th, th.rho, pre.gamma_h);
    th2 = th; th2.d_weeks = 2;
    [t2, ~, ~, ~] = n1_tau_star(pre, D, th2, th2.rho, pre.gamma_h);
    tauA(i, :) = [t1, t2];
    atB = (th.rho >= box.rho(2) - 1e-6) || (th.rho <= box.rho(1) + 1e-6);
    mm  = mwspm_metrics(pre.Y, b.Yhat);
    rowsA{end+1, 1} = {sh, 1/sh, b.J, b.nll, b.nll_full, bic, th.m0, th.rho, th.mu_s, ...
        th.onset_doy, th.offset_doy, th.gamma, th.beta, b.peak_week, b.true_total, ...
        b.reported_total, mm.R2_raw, mm.R2_log, t1, t2, double(atB)}; %#ok<AGROW>
    fprintf(logFid, '  sigma_h=1/%.4g (IIP %.4g d): J=%.4f rho=%.6f%s m0=%.6f mu_s=%.4f BIC=%.3f R2log=%.4f | tau*=[%s,%s]\n', ...
        sh, 1/sh, b.J, th.rho, n1_bnd(atB), th.m0, th.mu_s, bic, mm.R2_log, n1_s(t1), n1_s(t2));
end

% ================================================================ Axis B ----
fprintf(logFid, '\n--- Axis B: pi sensitivity (sigma_h = 1/%.4g fixed, d = 1) ---\n', sigmaH0);
rowsB = {};
tauB  = nan(numel(piVals), 2);
for i = 1:numel(piVals)
    piV = piVals(i);
    pre = pre0; pre.sigma_h = sigmaH0;
    C = n1_baseC(freeAll, box, piV, targetPeak, 1);
    b  = mwspm_calib_fit(pre, C, n1_merge(fitOpt, struct('warm', {warm0})));
    th = b.theta;
    bic = 2 * b.nll_full + kParam * log(nWeeks);
    [t1, ~, ~, ~] = n1_tau_star(pre, D, th, th.rho, pre.gamma_h);
    th2 = th; th2.d_weeks = 2;
    [t2, ~, ~, ~] = n1_tau_star(pre, D, th2, th2.rho, pre.gamma_h);
    tauB(i, :) = [t1, t2];
    atB = (th.rho >= box.rho(2) - 1e-6) || (th.rho <= box.rho(1) + 1e-6);
    mm  = mwspm_metrics(pre.Y, b.Yhat);
    rowsB{end+1, 1} = {piV, b.J, b.nll, b.nll_full, bic, th.m0, th.rho, th.mu_s, ...
        th.onset_doy, th.offset_doy, th.gamma, th.beta, b.peak_week, b.true_total, ...
        b.reported_total, mm.R2_raw, mm.R2_log, t1, t2, double(atB)}; %#ok<AGROW>
    fprintf(logFid, '  pi=%g: J=%.4f rho=%.6f%s m0=%.6f mu_s=%.4f peak=%d BIC=%.3f | tau*=[%s,%s]\n', ...
        piV, b.J, th.rho, n1_bnd(atB), th.m0, th.mu_s, b.peak_week, bic, n1_s(t1), n1_s(t2));
end

% ============================================================ verdicts ------
rangeA1 = n1_range(tauA(:, 1));  rangeA2 = n1_range(tauA(:, 2));
rangeB1 = n1_range(tauB(:, 1));  rangeB2 = n1_range(tauB(:, 2));
verdA1  = n1_verdict(rangeA1);   verdA2  = n1_verdict(rangeA2);
verdB1  = n1_verdict(rangeB1);   verdB2  = n1_verdict(rangeB2);

fprintf(logFid, '\n--- tau* spread (weeks) ---\n');
fprintf(logFid, '  Axis A (sigma_h): lag1 [%s] spread=%s -> %s ; lag2 [%s] spread=%s -> %s\n', ...
    n1_s(min(tauA(:,1))), n1_s(rangeA1), verdA1, n1_s(min(tauA(:,2))), n1_s(rangeA2), verdA2);
fprintf(logFid, '  Axis B (pi):      lag1 [%s] spread=%s -> %s ; lag2 [%s] spread=%s -> %s\n', ...
    n1_s(min(tauB(:,1))), n1_s(rangeB1), verdB1, n1_s(min(tauB(:,2))), n1_s(rangeB2), verdB2);

% baseline self-consistency: the sigma_h=1/5.9 variant (Axis A) and the pi=300
% variant (Axis B) are the SAME refit (baseline) and must reproduce the T12.01
% main report point tau* = 20 at lag 1.
[~, iA0] = min(abs(sigmaHvals - sigmaH0));
[~, iB0] = min(abs(piVals - pi0));
baselineLag1 = unique([tauA(iA0, 1), tauB(iB0, 1)]);
selfConsistent = (numel(baselineLag1) == 1) && isfinite(baselineLag1);
fprintf(logFid, '  baseline tau* (lag1) from Axis A and Axis B: %s ; self-consistent=%s\n', ...
    mat2str(baselineLag1), n1_yes(selfConsistent));

% =============================================================== outputs ----
hdrA = {'sigma_h_per_day', 'iip_days', 'J', 'nll', 'nll_full', 'BIC', 'm0', 'rho', ...
    'mu_s', 'onset_doy', 'offset_doy', 'gamma', 'beta', 'peak_week', 'true_total', ...
    'reported_total', 'R2_raw', 'R2_log', 'tau_star_lag1', 'tau_star_lag2', 'rho_at_bound'};
n1_write_rows(fullfile(ctx.dataDir, 'axisA_sigma_h.csv'), hdrA, rowsA, '%.17g,%d');

hdrB = {'pi', 'J', 'nll', 'nll_full', 'BIC', 'm0', 'rho', 'mu_s', 'onset_doy', ...
    'offset_doy', 'gamma', 'beta', 'peak_week', 'true_total', 'reported_total', ...
    'R2_raw', 'R2_log', 'tau_star_lag1', 'tau_star_lag2', 'rho_at_bound'};
n1_write_rows(fullfile(ctx.dataDir, 'axisB_pi.csv'), hdrB, rowsB, '%.17g,%d');

% combined tau* table (variant, axis, value, tau_star_lag1, tau_star_lag2)
fid = fopen(fullfile(ctx.dataDir, 'tau_star_table.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'axis,variant_value,tau_star_lag1,tau_star_lag2\n');
for i = 1:numel(sigmaHvals)
    fprintf(fid, 'sigma_h,%.17g,%s,%s\n', sigmaHvals(i), n1_s(tauA(i,1)), n1_s(tauA(i,2)));
end
for i = 1:numel(piVals)
    fprintf(fid, 'pi,%.17g,%s,%s\n', piVals(i), n1_s(tauB(i,1)), n1_s(tauB(i,2)));
end
fclose(fid);

verdict = struct();
verdict.axis_A_sigma_h = struct('tau_star_range_lag1_weeks', rangeA1, ...
    'tau_star_range_lag2_weeks', rangeA2, 'verdict_lag1', verdA1, 'verdict_lag2', verdA2, ...
    'condition', sprintf('pi = %g fixed, d = 1', pi0));
verdict.axis_B_pi = struct('tau_star_range_lag1_weeks', rangeB1, ...
    'tau_star_range_lag2_weeks', rangeB2, 'verdict_lag1', verdB1, 'verdict_lag2', verdB2, ...
    'condition', sprintf('sigma_h = 1/%.4g fixed, d = 1', sigmaH0));
verdict.baseline_selfconsistency = struct('baseline_tau_star_lag1', baselineLag1, ...
    'self_consistent', selfConsistent, ...
    'expectation', 'must equal the T12.01 main report point tau* = 20 (lag 1)');
verdict.rule = 'tau* spread <= 1 week => insensitive; >= 2 weeks => sensitive (preregistered D2026-10-10-01)';
mwspm_writejson(fullfile(ctx.dataDir, 'verdict.json'), verdict);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.preregistered_by = char(string(getdef(cfg, 'preregistered_by', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.main_report_point = fit.main_report_point;
summary.deterministic = true;
summary.rng = 'none (deterministic refit; fixed optimiser seeds for reproducible starts)';
summary.axis_A_sigma_h = struct('values', sigmaHvals, 'iip_days', 1 ./ sigmaHvals, ...
    'pi_fixed', pi0, 'tau_star_lag1', tauA(:, 1)', 'tau_star_lag2', tauA(:, 2)', ...
    'range_lag1_weeks', rangeA1, 'range_lag2_weeks', rangeA2);
summary.axis_B_pi = struct('values', piVals, 'sigma_h_fixed', sigmaH0, ...
    'tau_star_lag1', tauB(:, 1)', 'tau_star_lag2', tauB(:, 2)', ...
    'range_lag1_weeks', rangeB1, 'range_lag2_weeks', rangeB2);
summary.verdict = verdict;
summary.honesty = cfg.honesty;
summary.tau_star_caliber = ['tau* = min{ w : cum_w >= T0(w) }, T0(w) = rho/log R0(w) for R0(w) > 1; ' ...
    'cum_w = cumulative DETECTED model cases (mwspm_calib_kernel, Euler dt=1 d, mu_variant=day); ' ...
    'weekly R0 = mean of daily R0 (continuous raw mortality); lag2 = same refit theta re-evaluated at d_weeks=2.'];
summary.value_type_note = ['Every number is a deterministic point: the refit objective is a deterministic ODE ' ...
    'trajectory scored against the weekly observations. No interval is produced: these are OAT/one-axis ' ...
    'perturbations, NOT a joint uncertainty region.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.deterministic = true;
out.rng = 'none';
out.tau_star_range_axisA_lag1_weeks = rangeA1;
out.tau_star_range_axisA_lag2_weeks = rangeA2;
out.tau_star_range_axisB_lag1_weeks = rangeB1;
out.tau_star_range_axisB_lag2_weeks = rangeB2;
out.verdict_axisA_lag1 = verdA1;
out.verdict_axisB_lag1 = verdB1;
out.baseline_selfconsistency_pass = selfConsistent;
out.outputs = {'axisA_sigma_h.csv', 'axisB_pi.csv', 'tau_star_table.csv', 'verdict.json', 'summary.json'};
out.note = ['N1 (preregistered D2026-10-10-01). Axis A is a CONDITIONAL sensitivity at pi=300; ' ...
    'the pi=0 branch is reported separately, never merged into a robust narrative.'];
end

% ---------------------------------------------------------------------------
function C = n1_baseC(freeAll, box, piV, targetPeak, dList)
C = struct('free', {freeAll}, 'fixed', struct(), 'box', box, ...
    'family', 'poisson', 'family_opts', struct(), 'sub_weeks', [], 'y_obs', [], ...
    'phase_penalty', true, 'target_peak', targetPeak, 'pen', piV, 'd_list', dList);
end

% ---------------------------------------------------------------------------
function [tauStar, r0StartWeek, T0atTau, cumAtTau] = n1_tau_star(pre, D, th, rhoUsed, gammaH)
% Identical caliber to exp_c4_uncertainty's unc_tau_star (gate-G1 / T14.01).
nW = pre.n_weeks;
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
trRaw = mwspm_thermal_traits(D.T_daily, 'raw');
muEff = th.mu_s * trRaw.mu;
R0daily = mwspm_r0(D.T_daily, th.m0 * s, muEff, gammaH);
wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);
gt1 = find(R0daily > 1, 1);
if isempty(gt1)
    r0StartWeek = NaN;
else
    r0StartWeek = wk(gt1);
end
R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);
preK = pre;
preK.gamma_h = gammaH;
kern = mwspm_calib_kernel(preK, th, struct());
cumDet = cumsum(kern.Yhat(:));
T0w = nan(nW, 1);
ok = R0w > 1;
T0w(ok) = rhoUsed ./ log(R0w(ok));
cross = find(isfinite(T0w) & cumDet >= T0w, 1);
if isempty(cross)
    tauStar = NaN; T0atTau = NaN; cumAtTau = NaN;
else
    tauStar = cross; T0atTau = T0w(cross); cumAtTau = cumDet(cross);
end
end

% ---------------------------------------------------------------------------
function th = n1_theta_from_fit(tf)
th = struct('m0', double(tf.m0), 'rho', double(tf.rho), 'onset_doy', double(tf.onset_doy), ...
    'offset_doy', double(tf.offset_doy), 'mu_s', double(tf.mu_s), 'gamma', double(tf.gamma), ...
    'beta', double(tf.beta), 'd_weeks', double(tf.d_weeks));
end

% ---------------------------------------------------------------------------
function P = n1_find_parent(root, wantId, expId, relPath)
runsDir = fullfile(root, 'results', 'runs');
d = dir(runsDir);
P = [];
for k = 1:numel(d)
    if ~d(k).isdir || strcmp(d(k).name, '.') || strcmp(d(k).name, '..')
        continue
    end
    if ~isempty(wantId) && ischar(wantId) && ~strcmp(char(wantId), d(k).name)
        continue
    end
    f = fullfile(runsDir, d(k).name, relPath);
    if exist(f, 'file') ~= 2
        continue
    end
    manF = fullfile(runsDir, d(k).name, 'manifest.json');
    if exist(manF, 'file') == 2
        man = mwspm_readjson(manF);
        if isfield(man, 'status') && ~strcmp(char(string(man.status)), 'SUCCESS')
            continue
        end
        if isfield(man, 'experiment_id') && ~strcmp(char(string(man.experiment_id)), expId)
            continue
        end
    end
    try
        j = mwspm_readjson(f);
    catch
        continue
    end
    if ~isfield(j, 'provenance') || ~isfield(j.provenance, 'experiment_id') ...
            || ~strcmp(char(string(j.provenance.experiment_id)), expId)
        continue
    end
    if isempty(P) || strcmp(d(k).name, P.runId) > 0
        P = struct('runId', d(k).name, 'path', f);
    end
end
end

% ---------------------------------------------------------------------------
function o = n1_merge(base, over)
o = base;
fn = fieldnames(over);
for k = 1:numel(fn)
    o.(fn{k}) = over.(fn{k});
end
end

% ---------------------------------------------------------------------------
function r = n1_range(v)
f = isfinite(v);
if ~any(f)
    r = NaN;
else
    r = max(v(f)) - min(v(f));
end
end

function v = n1_verdict(r)
if isnan(r)
    v = 'NA';
elseif r <= 1
    v = 'insensitive';
else
    v = 'sensitive';
end
end

function s = n1_s(v)
if isfinite(v)
    s = sprintf('%.6g', v);
else
    s = 'NA';
end
end

function s = n1_bnd(tf)
if tf, s = '*'; else s = ' '; end
end

function s = n1_yes(c)
if c, s = 'YES'; else s = 'NO'; end
end

function n1_write_rows(path, hdr, rows, numFmt)
fid = fopen(path, 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:numel(rows)
    rw = rows{i};
    nNum = numel(rw);
    for k = 1:nNum
        if k > 1, fprintf(fid, ','); end
        if k == nNum   % last column is the at-bound flag (0/1 integer)
            fprintf(fid, '%d', rw{k});
        else
            fprintf(fid, '%.17g', rw{k});
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
