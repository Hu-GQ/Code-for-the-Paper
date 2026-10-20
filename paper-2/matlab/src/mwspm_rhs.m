function dX = mwspm_rhs(t, X, p, F)
%MWSPM_RHS  Right-hand side of the active model M-WSPM gen3_rhohat.
%
%   dX = mwspm_rhs(t, X, p, F)
%
%   X = [S_h; E_h; I_h; R_h; S_v; E_v; I_v; Cum]        (8x1)
%   t   elapsed days since 2025-02-03 (continuous for adaptive solvers,
%       grid values for fixed-step integrators)
%   p   parameter struct from mwspm_params
%   F   forcing struct from mwspm_make_forcing
%
%   This function is the MODEL layer only: it is deliberately free of any
%   experiment logic (no tolerances, no grids, no output paths, no stopping
%   rules).  Experiment entry points in src/matlab/exp_*.m own all of that.
%
%   Equations (theory/model_spec.md Sec.6, standard/frequency-dependent
%   incidence with denominator N_h):
%     dS_h = -lambda_h S_h
%     dE_h =  lambda_h S_h - sigma_h E_h
%     dI_h =  sigma_h E_h - gamma_h I_h + Lambda(t)
%     dR_h =  gamma_h I_h
%     dS_v =  mu_eff(N_v(t) - S_v) - lambda_v S_v
%     dE_v =  lambda_v S_v - (mu_eff + sigma_v) E_v
%     dI_v =  sigma_v E_v - mu_eff I_v
%     dCum =  sigma_h E_h                      (observational accumulator)
%   with
%     lambda_h = a(T) b(T) I_v / N_h ,  lambda_v = a(T) c(T) I_h / N_h
%     N_v(t)   = m0 N_h s(t) ,          mu_eff   = mu_s mu_v_variant(T)
%
%   Cum feeds back into NO equation, therefore the 7x7 Jacobian, the
%   next-generation matrix and every T08/T09 stability statement are
%   unaffected by its presence (recorded in the model card, state_note).

X = X(:);
assert(numel(X) == 8, 'mwspm:rhs', 'state must have 8 entries, got %d', numel(X));

f = mwspm_forcing_eval(F, t);

Sh = X(1); Eh = X(2); Ih = X(3); Rh = X(4);      %#ok<NASGU>  (Rh is recorded, not differentiated)
Sv = X(5); Ev = X(6); Iv = X(7);

lam_h = f.a * f.b * (Iv / p.N_h);
lam_v = f.a * f.c * (Ih / p.N_h);

dX = zeros(8, 1);
dX(1) = -lam_h * Sh;
dX(2) =  lam_h * Sh - p.sigma_h * Eh;
dX(3) =  p.sigma_h * Eh - p.gamma_h * Ih + f.Lambda;
dX(4) =  p.gamma_h * Ih;
dX(5) =  f.mu_eff * (f.Nv - Sv) - lam_v * Sv;
dX(6) =  lam_v * Sv - (f.mu_eff + f.sigma_v) * Ev;
dX(7) =  f.sigma_v * Ev - f.mu_eff * Iv;
dX(8) =  p.sigma_h * Eh;
end
