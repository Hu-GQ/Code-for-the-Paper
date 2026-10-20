function out = mwspm_floquet_metrics(S, opts)
%MWSPM_FLOQUET_METRICS  Three periodic criteria for a 329-day daily series.
%
%   out = mwspm_floquet_metrics(S)
%   out = mwspm_floquet_metrics(S, struct('doLinear', true, 'h', 1/40, 'hCheck', 1/80))
%
%   S fields (all 329x1 daily vectors, one window):
%     m    mosquito-to-host ratio m0*s(t)
%     a,b,c,sigma_v   temperature traits (continuous-rate)
%     mu   effective adult mortality  mu_s*mu_v_raw(T)
%
%   Criteria
%     (i)   pointwise   sup_t R0(t)          -- insufficient
%     (ii)  envelope    R_env on the elementwise Metzler envelope -- sufficient
%     (iii) Floquet     rho(Phi_T) for a 329-day zero-order-held (ZOH)
%                       artificial periodic repetition -- necessary and sufficient
%                       within the periodic setting
%
%   kappa scales the vector-to-human term a*c*m; kappa_F solves
%   log rho(Phi_T(kappa)) = 0 and kappa* = 1/R_env^2.
%
%   The implementation is the one used by
%   reviews/internal/20261016_followup_audit/gzhou_floquet_audit.m, moved here
%   so that production runs and the audit cannot drift apart.  The optional
%   linear-interpolation branch reproduces the older 328-day calculation and is
%   kept only as a regression check.

if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'doLinear'), opts.doLinear = false; end
if ~isfield(opts, 'h'),       opts.h = 1/40;     end
if ~isfield(opts, 'hCheck'),  opts.hCheck = 1/80; end

sigmaH = 1/5.9;
gammaH = 1/5;

n = numel(S.m);
assert(n == 329 && all(isfinite(S.m)) && all(S.mu > 0), 'floquet:input', ...
    'expected 329 finite days with positive mortality');

instant = sqrt(max(S.m .* S.a.^2 .* S.b .* S.c .* S.sv ./ ...
    (gammaH .* S.mu .* (S.mu + S.sv)), 0));
[peak, day] = max(instant);

P = S.a .* S.b;
Q = S.a .* S.c * max(S.m);
bar = [-sigmaH, 0, 0, max(P);
        sigmaH, -gammaH, 0, 0;
        0, max(Q), -min(S.mu + S.sv), 0;
        0, 0, max(S.sv), -min(S.mu)];
den = prod(-diag(bar));
loop = bar(1,4) * bar(2,1) * bar(3,2) * bar(4,3);
renv = sqrt(loop / den);

zoh = logrhoZoh(S, 1, sigmaH, gammaH);
kapZ = fzero(@(k) logrhoZoh(S, k, sigmaH, gammaH), [1e-3, 1e3]);

out = struct();
out.n_days = n;
out.sup_R0 = peak;
out.peak_day = day;
out.peak_week = floor((day - 1) / 7) + 1;
out.R0_envelope = renv;
out.kappa_star = 1 / renv^2;
out.envelope_spectral_abscissa = max(real(eig(bar)));
out.zoh_period_days = n;
out.zoh_logrho = zoh;
out.zoh_rho = exp(zoh);
out.zoh_kappa_F = kapZ;
out.zoh_kappa_root_logrho = logrhoZoh(S, kapZ, sigmaH, gammaH);
out.safety_factor = renv / peak;

if opts.doLinear
    lin = logrhoLinear(S, 1, opts.h, sigmaH, gammaH);
    fine = logrhoLinear(S, 1, opts.hCheck, sigmaH, gammaH);
    kapL = fzero(@(k) logrhoLinear(S, k, opts.h, sigmaH, gammaH), [1e-3, 1e3]);
    out.linear_period_days = n - 1;
    out.linear_logrho = lin;
    out.linear_rho = exp(lin);
    out.linear_logrho_step_difference = lin - fine;
    out.linear_kappa_F = kapL;
    out.linear_kappa_root_logrho = logrhoLinear(S, kapL, opts.hCheck, sigmaH, gammaH);
end
end

% ---------------------------------------------------------------------------
function lr = logrhoZoh(S, k, sigmaH, gammaH)
Y = eye(4); acc = 0;
for i = 1:numel(S.m)
    M = matrixAt(S, i, k, sigmaH, gammaH);
    Y = expm(M) * Y;
    c = max(abs(Y(:)));
    assert(c > 0 && isfinite(c), 'floquet:overflow', 'non-finite fundamental matrix');
    Y = Y / c; acc = acc + log(c);
end
lr = acc + log(max(abs(eig(Y))));
end

% ---------------------------------------------------------------------------
function lr = logrhoLinear(S, k, h, sigmaH, gammaH)
n = round(1/h);
assert(abs(n*h - 1) < 1e-12, 'floquet:step', 'h must divide one day');
Y = eye(4); acc = 0;
for day = 1:numel(S.m)-1
    for j = 1:n
        u = (j-1) / n;
        A = matrixBetween(S, day, u, k, sigmaH, gammaH);
        B = matrixBetween(S, day, u + h/2, k, sigmaH, gammaH);
        C = matrixBetween(S, day, u + h, k, sigmaH, gammaH);
        y1 = A * Y; y2 = B * (Y + h*y1/2); y3 = B * (Y + h*y2/2); y4 = C * (Y + h*y3);
        Y = Y + h*(y1 + 2*y2 + 2*y3 + y4)/6;
    end
    c = max(abs(Y(:)));
    assert(c > 0 && isfinite(c), 'floquet:overflow', 'non-finite fundamental matrix');
    Y = Y / c; acc = acc + log(c);
end
lr = acc + log(max(abs(eig(Y))));
end

% ---------------------------------------------------------------------------
function M = matrixAt(S, i, k, sigmaH, gammaH)
M = [-sigmaH, 0, 0, S.a(i)*S.b(i);
      sigmaH, -gammaH, 0, 0;
      0, k*S.a(i)*S.c(i)*S.m(i), -(S.mu(i) + S.sv(i)), 0;
      0, 0, S.sv(i), -S.mu(i)];
end

% ---------------------------------------------------------------------------
function M = matrixBetween(S, i, u, k, sigmaH, gammaH)
a  = S.a(i)  + (S.a(i+1)  - S.a(i))  * u;
b  = S.b(i)  + (S.b(i+1)  - S.b(i))  * u;
c  = S.c(i)  + (S.c(i+1)  - S.c(i))  * u;
m  = S.m(i)  + (S.m(i+1)  - S.m(i))  * u;
mu = S.mu(i) + (S.mu(i+1) - S.mu(i)) * u;
sv = S.sv(i) + (S.sv(i+1) - S.sv(i)) * u;
M = [-sigmaH, 0, 0, a*b;
      sigmaH, -gammaH, 0, 0;
      0, k*a*c*m, -(mu + sv), 0;
      0, 0, sv, -mu];
end
