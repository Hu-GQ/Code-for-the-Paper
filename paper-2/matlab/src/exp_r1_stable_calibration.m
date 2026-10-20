function payload = exp_r1_stable_calibration(ctx)
%EXP_R1_STABLE_CALIBRATION  Unpenalized conditional fits; never a rho estimate.
cfg = ctx.cfg;
if isfield(ctx,'dryRun') && ctx.dryRun
    payload = struct('outputs',{{}},'dry_run',true);
    return
end
diary(fullfile(ctx.runDir,'logs','stdout.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
snap = fullfile(ctx.dataDir,'code_snapshot');
mkdir(snap);
copyfile(fullfile(ctx.root,'src','matlab','*.m'), snap);
copyfile(fullfile(ctx.root,'scripts','reproduce.m'), snap);
copyfile(fullfile(ctx.root,'tests','check_revision_stable.m'), snap);
addpath(fullfile(ctx.root,'tests'));
verification = check_revision_stable(ctx.root);
mwspm_writejson(fullfile(ctx.dataDir,'verification.json'),verification);
assert(~cfg.phase_penalty, 'new revision must not use a peak penalty');
D = mwspm_load_data(ctx.root,ctx.modelCard);
pre = mwspm_calib_prepare(D,ctx.modelCard,struct('mu_variant',cfg.mu_variant));
old = mwspm_readjson(fullfile(ctx.root,'results','runs',cfg.legacy_run,'data','fit_best.json'));
warm = mwspm_readjson(fullfile(ctx.root,'results','runs',cfg.stable_warm_run,'data','refine.json'));
C = struct('free',{cellstr(string(cfg.free))},'fixed',struct('rho',cfg.fixed_rho), ...
    'box',mwspm_calib_theta('box'),'family','poisson','phase_penalty',false, ...
    'kernel_opts',cfg.solver);
fits = struct();
allTrace = table();
for d = cfg.lag_grid(:)'
    key = sprintf('lag%d',d);
    C.d_list = d;
    opt = cfg.optimiser;
    opt.warm = {old.free_fits.(key).theta, warm.(key).theta};
    [best,info] = mwspm_calib_fit(pre,C,opt);
    a = mwspm_calib_kernel(pre,best.theta,cfg.solver);
    b = mwspm_calib_kernel(pre,best.theta,struct('solver','rk4','dt',cfg.solver.dt/2));
    c = mwspm_calib_kernel(pre,best.theta,struct('solver','ode15s'));
    Lc = mwspm_calib_nll(pre.Y,c.Yhat,'poisson');
    solverCheck = struct('step_halving_max_weekly',max(abs(a.Yhat-b.Yhat)), ...
        'ode15s_max_weekly',max(abs(a.Yhat-c.Yhat)), ...
        'ode15s_nll_difference',abs(best.nll-Lc.nll));
    assert(solverCheck.step_halving_max_weekly<cfg.verification.max_weekly_solver_difference && ...
        solverCheck.ode15s_max_weekly<cfg.verification.max_weekly_solver_difference && ...
        solverCheck.ode15s_nll_difference<cfg.verification.max_nll_solver_difference, ...
        'candidate solver verification failed');
    assert(abs(best.J-best.nll)<1e-10, 'phase penalty leaked into new fit');
    fits.(key) = struct('best',best,'info',info,'solver_check',solverCheck);
    allTrace = [allTrace; struct2table(info.trace)]; %#ok<AGROW>
    mwspm_writejson(fullfile(ctx.dataDir,'stable_fit.json'),fits);
    writetable(allTrace,fullfile(ctx.dataDir,'fit_trace.csv'));
    fprintf('R1 lag%d nll=%.9f peak=%d flag=%d solver error=%.3g\n', ...
        d,best.nll,best.peak_week,best.exitflag,solverCheck.ode15s_max_weekly);
end
summary = struct('interpretation','conditional pure-likelihood candidates; rho fixed, no inference acceptance', ...
    'fixed_rho',cfg.fixed_rho,'phase_penalty',false,'solver',cfg.solver, ...
    'all_start_flags',allTrace.exitflag','all_starts_converged',all(allTrace.exitflag==1));
mwspm_writejson(fullfile(ctx.dataDir,'summary.json'),summary);
payload = summary;
payload.outputs = cfg.outputs;
end
