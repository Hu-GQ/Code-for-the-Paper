% t04_pilot_wangzhao_example2.m
% -------------------------------------------------------------------------
% T04 pilot (SMALL-SCALE, method validation only).
%
% Purpose: validate the periodic basic reproduction ratio implementation
%          R0 = rho(L) (Wang & Zhao 2008, Theorem 2.1/2.2) against the
%          NUMBERS ALREADY PUBLISHED in Wang & Zhao (2008) Example 2
%          (J Dyn Diff Equat 20:699-717, PDF p.10-11, Fig.1 and Fig.2):
%            (a) delta = 1 : [R0] = 1 at k = 0.70, but true R0 = 1 at k = 0.614
%            (b) k = 0.65 : [R0] is always 0.9636 as delta varies in [0,1],
%                           while R0 > 1 for 0.75 < delta < 1
%
% This is a NUMERICAL SIGN ONLY. It is not a proof, and it says nothing
% about the Guangzhou 2025 application. Its only job is to confirm that the
% monodromy / root-finding machinery we plan to use in T08 reproduces a
% published benchmark before we trust it on our own model.
%
% Protocol: standalone script; reads its run config from
%           results/runs/<run_id>/config.json; writes outputs into the same
%           run directory; never touches the source project.
% -------------------------------------------------------------------------
clear; clc; close all;

fprintf('=== T04 pilot: Wang-Zhao periodic R0 validation ===\n');

% ---- locate project root (the dir containing config/project.json) --------
here = fileparts(mfilename('fullpath'));
root = here;
for kk = 1:6
    if exist(fullfile(root, 'config', 'project.json'), 'file'), break; end
    root = fileparts(root);
end
assert(exist(fullfile(root, 'config', 'project.json'), 'file') == 2, ...
    'project root not found from %s', here);
fprintf('project root: %s\n', root);

% ---- load run config -----------------------------------------------------
runId  = '20261003T024149_ddbae9e8';
runDir = fullfile(root, 'results', 'runs', runId);
cfg    = jsondecode(fileread(fullfile(runDir, 'config.json')));
fprintf('run_id: %s\n', runId);

P     = cfg.params;
beta0 = P.beta0; b = P.b; gam = P.gamma; c = P.c; omega = P.omega;
RelTol = cfg.RelTol; AbsTol = cfg.AbsTol; btol = cfg.bisection_tol;

V = diag([b + gam, c]);   % constant transition matrix; Vbar = V

% ================= TEST 1: Fig.1 -- delta = 1, find k where R0 = 1 ========
fprintf('\n--- Test 1 (Fig.1): delta = 1, solve k such that R0 = 1 ---\n');
delta = 1.0;

% analytic: [R0]=1 <=> k = (b+gamma)*c/beta0   (since int_0^1 cos(2 pi t) dt = 0)
kAvg = (b + gam) * c / beta0;
fprintf('[R0] = 1 at k = %.4f   (published: 0.70)\n', kAvg);

% true periodic threshold: R0 = 1 <=> rho(Phi_{F-V}(omega)) = 1 (Theorem 2.2)
lo = 0.30; hi = 1.20;
rhoLo = rhoMono(1.0, lo, delta, beta0, V, omega, RelTol, AbsTol);
rhoHi = rhoMono(1.0, hi, delta, beta0, V, omega, RelTol, AbsTol);
fprintf('rho(Phi_{F-V}) at k=%.2f : %.6f ; at k=%.2f : %.6f\n', lo, rhoLo, hi, rhoHi);
assert((rhoLo - 1) < 0 && (rhoHi - 1) > 0, 'Test 1: no sign change on [0.30,1.20]');
while (hi - lo) > btol
    mid = 0.5*(lo + hi);
    if (rhoMono(1.0, mid, delta, beta0, V, omega, RelTol, AbsTol) - 1) < 0
        lo = mid;
    else
        hi = mid;
    end
end
kTrue = 0.5*(lo + hi);
fprintf('R0  = 1 at k = %.4f   (published: 0.614)\n', kTrue);
if kAvg > kTrue
    fprintf('=> time-average UNDERESTIMATES risk (kAvg %.3f > kTrue %.3f)\n', kAvg, kTrue);
else
    fprintf('=> time-average OVERESTIMATES risk (kAvg %.3f < kTrue %.3f)\n', kAvg, kTrue);
end

% ================= TEST 2: Fig.2 -- k = 0.65, sweep delta =================
fprintf('\n--- Test 2 (Fig.2): k = 0.65, sweep delta in [0,1] ---\n');
kFix      = 0.65;
deltaList = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0];
R0avgList  = zeros(size(deltaList));
R0trueList = zeros(size(deltaList));
for i = 1:numel(deltaList)
    d            = deltaList(i);
    R0avgList(i)  = avgR0(kFix, beta0, V);
    R0trueList(i) = solveR0(kFix, d, beta0, V, omega, RelTol, AbsTol, btol);
    if R0trueList(i) > 1, flag = 'R0>1'; else, flag = 'R0<1'; end
    fprintf('  delta=%.2f : [R0]=%.4f   R0=%.4f   %s\n', ...
        d, R0avgList(i), R0trueList(i), flag);
end

% ---- write results -------------------------------------------------------
outFile = fullfile(runDir, 'data', 'pilot_results.csv');
fid = fopen(outFile, 'w');
fprintf(fid, 'test,param,param_value,metric,value,published\n');
fprintf(fid, 'fig1,delta,1.0,k_for_R0avg_eq_1,%.6f,0.70\n', kAvg);
fprintf(fid, 'fig1,delta,1.0,k_for_R0_eq_1,%.6f,0.614\n', kTrue);
for i = 1:numel(deltaList)
    fprintf(fid, 'fig2,delta,%.2f,R0avg,%.6f,0.9636\n', deltaList(i), R0avgList(i));
    fprintf(fid, 'fig2,delta,%.2f,R0true,%.6f,NA\n', deltaList(i), R0trueList(i));
end
fclose(fid);
fprintf('\nwrote %s\n', outFile);

% ---- manifest ------------------------------------------------------------
man = struct();
man.run_id           = runId;
man.status           = 'SUCCESS';
man.purpose          = 'validate Wang-Zhao periodic R0 implementation against published Example 2';
man.k_for_R0avg_eq_1 = kAvg;
man.k_for_R0_eq_1    = kTrue;
man.delta_list       = deltaList;
man.R0avg_list       = R0avgList;
man.R0true_list      = R0trueList;
man.note             = 'numerical sign only; NOT a proof; benchmark reproduction, not a new result';
man.matlab_version   = version();
fid = fopen(fullfile(runDir, 'manifest.json'), 'w');
fwrite(fid, jsonencode(man));
fclose(fid);

fprintf('\n=== pilot finished ===\n');

% ============================ LOCAL FUNCTIONS ============================
function F = Fmat(t, k, delta, beta0)
    bs = k     * (1 + delta * cos(2*pi*t));
    bm = beta0 * (1 + delta * cos(2*pi*t));
    F  = [0, bs; bm, 0];
end

function rho = rhoMono(lam, k, delta, beta0, V, omega, RelTol, AbsTol)
    % Spectral radius of the monodromy matrix Phi(omega,0,lam) of
    %   dW/dt = (-V + F(t)/lam) W ,  W(0) = I
    A   = @(t) (-V + Fmat(t, k, delta, beta0) / lam);
    rhs = @(t, y) [A(t)*y(1:2); A(t)*y(3:4)];   % two columns stacked
    opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
    [~, Y] = ode45(rhs, [0, omega], [1; 0; 0; 1], opts);
    W   = [Y(end,1:2).', Y(end,3:4).'];
    rho = max(abs(eig(W)));
end

function lam = solveR0(k, delta, beta0, V, omega, RelTol, AbsTol, btol)
    % Theorem 2.1: R0 is the unique lambda with rho(W(omega,0,lambda)) = 1;
    % rho is non-increasing in lambda.
    lo = 1e-4; hi = 50;
    if (rhoMono(lo, k, delta, beta0, V, omega, RelTol, AbsTol) - 1) <= 0
        lam = lo; return;
    end
    if (rhoMono(hi, k, delta, beta0, V, omega, RelTol, AbsTol) - 1) >= 0
        lam = hi; return;
    end
    while (hi - lo) > btol
        mid = 0.5*(lo + hi);
        if (rhoMono(mid, k, delta, beta0, V, omega, RelTol, AbsTol) - 1) > 0
            lo = mid;
        else
            hi = mid;
        end
    end
    lam = 0.5*(lo + hi);
end

function R0avg = avgR0(k, beta0, Vbar)
    % [R0] = rho( Fbar * Vbar^{-1} );  int_0^1 cos(2 pi t) dt = 0
    %   =>  Fbar = [0, k; beta0, 0]
    Fbar  = [0, k; beta0, 0];
    R0avg = max(abs(eig(Fbar / Vbar)));
end
