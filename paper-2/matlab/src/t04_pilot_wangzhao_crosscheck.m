% t04_pilot_wangzhao_crosscheck.m
% -------------------------------------------------------------------------
% T04 pilot, part 2. The first pilot reproduced Wang & Zhao (2008) Example 2's
% TIME-AVERAGED number exactly ([R0] = 0.9636 at k=0.65; [R0]=1 at k=0.70)
% but NOT their periodic number (we get R0=1 at k=0.6697; they publish 0.614).
%
% This script checks whether that gap is an integration artifact of ode45 by
% recomputing the same monodromy with an INDEPENDENT method: a fine
% product-of-exponentials approximation,
%   Phi(omega,0,lam) ~= prod_j expm(A(t_j) * h),   h = omega/N
% which converges as N -> infinity but is completely different in structure.
%
% If both methods agree, the discrepancy against the published 0.614 is NOT
% our integration error and must be recorded as an open question (do not use
% their Example 2 magnitude as a calibration target until resolved).
% -------------------------------------------------------------------------
clear; clc; close all;

beta0 = 0.3; b = 2.0; gam = 0.1; c = 0.1; omega = 1.0;
V = diag([b + gam, c]);

fprintf('=== cross-check: ode45 vs product-of-exponentials ===\n');

function F = Fmat2(t, k, delta, beta0)
    bs = k     * (1 + delta * cos(2*pi*t));
    bm = beta0 * (1 + delta * cos(2*pi*t));
    F  = [0, bs; bm, 0];
end

function rho = rho_expmprod(lam, k, delta, beta0, V, omega, N)
    % Phi ~= prod_j expm(A(t_j) h)  (forward product, left-to-right ordering)
    h  = omega / N;
    Wm = eye(2);
    for j = 0:(N-1)
        t = j * h;
        A = -V + Fmat2(t, k, delta, beta0) / lam;
        Wm = expm(A * h) * Wm;
    end
    rho = max(abs(eig(Wm)));
end

% ---- case A: delta = 1, k = 0.65 (the Fig.2 endpoint) --------------------
k = 0.65; delta = 1.0;
fprintf('\ncase A: delta=1, k=0.65\n');
for N = [200, 1000, 5000, 20000]
    r = rho_expmprod(1.0, k, delta, beta0, V, omega, N);
    fprintf('  N=%6d : rho(Phi_{F-V}) = %.8f\n', N, r);
end
% ode45 reference (high accuracy) for the same case
Afun = @(t) (-V + Fmat2(t, k, delta, beta0));
rhs  = @(t, y) [Afun(t)*y(1:2); Afun(t)*y(3:4)];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, Y] = ode45(rhs, [0, omega], [1; 0; 0; 1], opts);
W  = [Y(end,1:2).', Y(end,3:4).'];
r45 = max(abs(eig(W)));
fprintf('  ode45 (RelTol 1e-12): rho = %.8f\n', r45);

% ---- case B: solve R0 for delta=1, k=0.65 via lambda root-find -----------
% (lambda root-find gives R0 directly)
fL = @(lam) rho_expmprod(lam, k, delta, beta0, V, omega, 5000) - 1;
lo = 0.5; hi = 2.0;
fprintf('\ncase B: R0 at delta=1, k=0.65 by lambda root-find (N=5000)\n');
fprintf('  f(%.2f)=%+.6f  f(%.2f)=%+.6f\n', lo, fL(lo), hi, fL(hi));
if fL(lo) > 0 && fL(hi) < 0
    while (hi - lo) > 1e-6
        mid = 0.5*(lo+hi);
        if fL(mid) > 0, lo = mid; else, hi = mid; end
    end
    fprintf('  R0 = %.6f\n', 0.5*(lo+hi));
else
    fprintf('  no sign change in [0.5,2.0]\n');
end

% ---- case C: what k makes R0 = 1 (delta=1), expm-product route -----------
fprintf('\ncase C: k such that rho(Phi_{F-V})=1, delta=1 (N=5000)\n');
g = @(kk) rho_expmprod(1.0, kk, delta, beta0, V, omega, 5000) - 1;
lo = 0.40; hi = 0.90;
fprintf('  g(%.2f)=%+.6f  g(%.2f)=%+.6f\n', lo, g(lo), hi, g(hi));
if g(lo) < 0 && g(hi) > 0
    while (hi - lo) > 1e-6
        mid = 0.5*(lo+hi);
        if g(mid) < 0, lo = mid; else, hi = mid; end
    end
    fprintf('  k(R0=1) = %.6f   (published 0.614; ode45 gave 0.6697)\n', 0.5*(lo+hi));
else
    fprintf('  no sign change in [0.40,0.90]\n');
end

fprintf('\n=== cross-check finished ===\n');
