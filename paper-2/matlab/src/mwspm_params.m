function p = mwspm_params(modelCard, overrides)
%MWSPM_PARAMS  Build the parameter struct of the active model from its card.
%
%   p = mwspm_params(modelCard)             values straight from the card
%   p = mwspm_params(modelCard, overrides)  scalar fields overridden
%
%   The model card is config/experiments/model_mwspm_gen3_rhohat.json; it is
%   the machine-readable mirror of theory/model_spec.md (spec v1.0) and
%   data/metadata/parameters.csv.  Overrides are the ONLY way an experiment
%   changes a parameter, and they are always recorded in the run config.
%
%   The returned struct is plain (no objects, no handles) so that it can be
%   written verbatim into manifest.json.

assert(isstruct(modelCard), 'mwspm:params', 'modelCard must be a struct');
assert(strcmp(modelCard.model_version, 'gen3_rhohat'), ...
    'mwspm:params', 'unexpected model version "%s"', modelCard.model_version);

f  = modelCard.parameters_fitted;
fx = modelCard.parameters_fixed;

p = struct();
p.model_id       = modelCard.model_id;
p.model_version  = modelCard.model_version;
p.spec_version   = modelCard.spec_version;

% fitted (main report point)
p.m0             = double(f.m0);
p.rho            = double(f.rho);
p.d_weeks        = double(f.d_weeks);
p.onset_doy      = double(f.onset_doy);
p.offset_doy     = double(f.offset_doy);
p.mu_s           = double(f.mu_s);
p.gamma          = double(f.gamma);
p.beta           = double(f.beta);

% fixed
p.N_h            = double(fx.N_h);
p.sigma_h        = double(fx.sigma_h);
p.gamma_h        = double(fx.gamma_h);
p.kappa_days     = double(fx.kappa_days);
p.s_floor        = double(fx.s_floor);
p.mu_v_cap       = double(fx.mu_v_cap_per_day);
p.min_window_days = 30.0;

% window
p.t0_date        = modelCard.window.t0_date;
p.t_end_date     = modelCard.window.t_end_date;
p.n_days         = double(modelCard.window.n_days);
p.n_weeks        = double(modelCard.window.n_weeks);

% season switch: 'window' (default) or 'constant' (s(t) == 1, source baseline)
p.season_mode    = 'window';

if nargin >= 2 && ~isempty(overrides)
    assert(isstruct(overrides), 'mwspm:params', 'overrides must be a struct');
    fn = fieldnames(overrides);
    for k = 1:numel(fn)
        p.(fn{k}) = overrides.(fn{k});
    end
end

% derived convenience values (recomputed after overrides)
p.days_elapsed = p.n_days;
p.span_days    = double(days(datetime(p.t_end_date) - datetime(p.t0_date))) + 1;

% season spec passed to mwspm_season
sp = struct();
sp.onset_doy       = p.onset_doy;
sp.offset_doy      = p.offset_doy;
sp.kappa_days      = p.kappa_days;
sp.gamma           = p.gamma;
sp.beta            = p.beta;
sp.s_floor         = p.s_floor;
sp.min_window_days = p.min_window_days;
sp.normalize       = strcmp(p.season_mode, 'window');
sp.norm_factor     = NaN;
p.season = sp;
end
