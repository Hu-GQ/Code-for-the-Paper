function out = mwspm_branching_sim(R, tauWeek, rho, nSim, opts)
%MWSPM_BRANCHING_SIM  Stochastic tau-leap replicates of the 4-type process.
%
%   out = mwspm_branching_sim(R, tauWeek, rho, nSim, opts)
%
%   N_sim independent replicates of the time-inhomogeneous branching process
%   whose rates are R (mwspm_branching_rates), seeded by ONE imported
%   infectious human (type 2, I_h) at the start of week tauWeek.
%
%   Estimand that is actually simulated
%   -----------------------------------
%   established = "the lineage produces at least one DETECTED autochthonous
%   human infection inside the analysis window".  Each newly created human
%   infection is detected independently with probability rho; the imported seed
%   itself is NOT counted, because counting it would make the estimand trivial
%   (P >= rho for every tau at every tau).  Replicates stop at the first
%   detection, so the run also yields the empirical distribution of the
%   detection week, which is the genuine REPLICATE DISTRIBUTION the summary
%   must keep separate from parameter intervals.
%
%   The survival-to-t_end estimand is NOT simulated here: surviving lineages
%   grow without bound and the run cost explodes.  It is available EXACTLY from
%   mwspm_branching_backward, so simulating it would add cost without adding
%   information.  This choice is recorded in raw_choices.csv.
%
%   Discretisation
%   --------------
%   Fixed sub-steps inside each day (the forcing is zero-order-hold on the
%   daily grid, so rates are constant within a day and the sub-stepping is
%   pure time discretisation).  The number of sub-steps per day is chosen so
%   that max_rate * dt <= opts.p_target; the resulting discretisation error is
%   MEASURED by re-running with p_target/4, never assumed away.
%
%   opts fields
%     p_target    max per-capita event probability per sub-step (default 0.05)
%     seed        RandStream seed (REQUIRED; no time-based seeds)
%     cap_live    live-population cap (default 400)
%     cap_cum     cumulative human infections cap (default 5000)
%
%   out fields
%     p_est, mc_se, n_established, n_extinct, n_alive_at_end, n_cap
%     det_week_q05/q25/q50/q75/q95  (among established replicates)
%     det_week_all_nan_fraction
%     n_substeps_total, n_overshoot, n_normal_fallback, elapsed_seconds
%     seed, generator, p_target

if nargin < 5 || isempty(opts)
    opts = struct();
end
pTarget = 0.05;
if isfield(opts, 'p_target') && ~isempty(opts.p_target)
    pTarget = double(opts.p_target);
end
capLive = 400;
if isfield(opts, 'cap_live') && ~isempty(opts.cap_live)
    capLive = double(opts.cap_live);
end
capCum = 5000;
if isfield(opts, 'cap_cum') && ~isempty(opts.cap_cum)
    capCum = double(opts.cap_cum);
end
assert(isfield(opts, 'seed') && ~isempty(opts.seed), ...
    'mwspm:branching', 'opts.seed is required; time-based seeds are not allowed');
seed = double(opts.seed);
stream = RandStream('mt19937ar', 'Seed', seed);

nDays = R.n_days;
dStart = 7 * (double(tauWeek) - 1) + 1;      % 1-based day index of introduction
dStart = min(max(dStart, 1), nDays);

sh = R.sigma_h;
gh = R.gamma_h;

Eh = zeros(nSim, 1);
Ih = ones(nSim, 1);          % seed: one imported infectious human
Ev = zeros(nSim, 1);
Iv = zeros(nSim, 1);
cumH = zeros(nSim, 1);

detected    = false(nSim, 1);
extinctFlag = false(nSim, 1);
aliveEnd    = false(nSim, 1);
capFlag     = false(nSim, 1);
firstDetDay = nan(nSim, 1);

act = (1:nSim)';
nSubTotal = 0;
nOver = 0;
t0 = tic;

for k = dStart:nDays
    if isempty(act)
        break
    end
    rhv = R.rhv_daily(k);
    rvh = R.rvh_daily(k);
    sv  = R.sigma_v_daily(k);
    mu  = R.mu_daily(k);

    rmax = max([sh, gh, rhv, sv + mu, mu, rvh]);
    nsub = max(1, ceil(rmax / pTarget));
    dt = 1 / nsub;
    nSubTotal = nSubTotal + nsub;

    for s = 1:nsub
        K = numel(act);
        eh = Eh(act); ih = Ih(act); ev = Ev(act); iv = Iv(act);

        bV = zeros(K, 1);
        if rhv > 0 && any(ih)
            bV = mwspm_poissrnd(ih * rhv * dt, stream);
        end
        bH = zeros(K, 1);
        if rvh > 0 && any(iv)
            bH = mwspm_poissrnd(iv * rvh * dt, stream);
        end
        pH = zeros(K, 1);
        if any(eh)
            pH = mwspm_poissrnd(eh * sh * dt, stream);
        end
        dH = zeros(K, 1);
        if any(ih)
            dH = mwspm_poissrnd(ih * gh * dt, stream);
        end
        pV = zeros(K, 1);
        if sv > 0 && any(ev)
            pV = mwspm_poissrnd(ev * sv * dt, stream);
        end
        dV = zeros(K, 1);
        if mu > 0 && any(ev)
            dV = mwspm_poissrnd(ev * mu * dt, stream);
        end
        dI = zeros(K, 1);
        if mu > 0 && any(iv)
            dI = mwspm_poissrnd(iv * mu * dt, stream);
        end

        % The tau-leap draws the competing exits independently, so their sum can
        % exceed the compartment.  Clamping is a defect of the scheme; it is
        % COUNTED, not hidden, and the count shrinks as p_target shrinks.
        pH = min(pH, eh);
        dH = min(dH, ih);
        nOver = nOver + nnz((pV + dV) > ev);
        pV = min(pV, ev);
        dV = min(dV, ev - pV);
        dI = min(dI, iv);

        eh2 = eh - pH + bH;
        ih2 = ih + pH - dH;
        ev2 = ev + bV - pV - dV;
        iv2 = iv + pV - dI;

        if any(bH)
            % P(at least one of bH new infections is detected) = 1-(1-rho)^bH.
            % One uniform per replicate is exact for the indicator we need.
            pdet = 1 - (1 - rho) .^ bH;
            u = rand(stream, K, 1);
            newDet = (bH > 0) & (u < pdet);
        else
            newDet = false(K, 1);
        end

        cumH(act) = cumH(act) + bH;
        justDet = act(newDet & ~detected(act));
        if ~isempty(justDet)
            detected(justDet) = true;
            firstDetDay(justDet) = k;
        end

        Eh(act) = eh2; Ih(act) = ih2; Ev(act) = ev2; Iv(act) = iv2;

        tot = eh2 + ih2 + ev2 + iv2;
        dead = tot <= 0;
        if any(dead)
            idx = act(dead);
            extinctFlag(idx) = true;
        end

        hitCap = (tot > capLive) | (cumH(act) > capCum);
        if any(hitCap)
            idx = act(hitCap);
            capFlag(idx) = true;
            % Reaching the cap without a detection has probability (1-rho)^cap,
            % which is <= 0.87^400 ~ 1e-24 on this grid.  Declaring the lineage
            % established here is therefore numerically safe, and the number of
            % times it happens is reported.
            detected(idx) = true;
            needDay = idx(isnan(firstDetDay(idx)));
            if ~isempty(needDay)
                firstDetDay(needDay) = k;
            end
        end

        keep = ~dead & ~hitCap & ~detected(act);
        act = act(keep);
        if isempty(act)
            break
        end
    end
end

if ~isempty(act)
    aliveEnd(act) = true;
end

el = toc(t0);
pEst = mean(detected);
se = sqrt(pEst * (1 - pEst) / nSim);

dw = firstDetDay(detected & ~isnan(firstDetDay));
q = @(p) tern(isempty(dw), NaN, quantile(dw, p));

out = struct();
out.tau_week = double(tauWeek);
out.rho = double(rho);
out.n_sim = nSim;
out.p_est = pEst;
out.mc_se = se;
out.n_established = nnz(detected);
out.n_extinct = nnz(extinctFlag);
out.n_alive_at_end = nnz(aliveEnd);
out.n_cap_forced = nnz(capFlag);
out.det_week_q05 = q(0.05);
out.det_week_q25 = q(0.25);
out.det_week_q50 = q(0.50);
out.det_week_q75 = q(0.75);
out.det_week_q95 = q(0.95);
out.n_det_week_observed = numel(dw);
out.n_substeps_total = nSubTotal;
out.n_overshoot = nOver;
out.cap_live = capLive;
out.cap_cum = capCum;
out.p_target = pTarget;
out.seed = seed;
out.generator = 'mt19937ar';
out.elapsed_seconds = el;
out.first_detection_day = firstDetDay;
end

% ---------------------------------------------------------------------------
function v = tern(c, a, b)
if c, v = a; else v = b; end
end
