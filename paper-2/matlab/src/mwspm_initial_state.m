function X0 = mwspm_initial_state(D, p, F)
%MWSPM_INITIAL_STATE  Initial condition of the active model (theory/model_spec.md Sec.5.2).
%
%   X0 = mwspm_initial_state(D, p, F)
%
%   S_h(0)=N_h, E_h=I_h=R_h=0, S_v(0)=m0 N_h s(t0), E_v=I_v=0, Cum(0)=0.
%
%   Meaning: the entire population is susceptible (2025 Guangzhou is not an
%   endemic setting and there is no large 2024 wave to leave immunity), the
%   vector pool is initialised at the t0 abundance and is fully susceptible,
%   and every local infection must be seeded by the importation forcing
%   Lambda(t).  When the season is constant (baseline model) s(t0) = 1.

s0 = 1.0;
if strcmp(p.season_mode, 'window')
    f0 = mwspm_forcing_eval(F, 0);
    s0 = f0.s;
end

X0 = zeros(8, 1);
X0(1) = p.N_h;
X0(2) = 0;
X0(3) = 0;
X0(4) = 0;
X0(5) = p.m0 * p.N_h * s0;
X0(6) = 0;
X0(7) = 0;
X0(8) = 0;

% sanity: the susceptible vector pool must be non-negative and finite
assert(all(isfinite(X0)) && all(X0 >= 0), 'mwspm:initial', 'invalid initial state');
end
