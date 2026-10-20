function R = mwspm_r0_timeseries(D, p, opts)
%MWSPM_R0_TIMESERIES  Daily R0(t) of the active model plus window diagnostics.
%
%   R = mwspm_r0_timeseries(D, p)
%   R = mwspm_r0_timeseries(D, p, opts)
%       opts.mu_variant  'raw' (default, REQUIRED for R0) | 'day'
%       opts.caliber     'Nv' (default) | 'V'
%       opts.V           in-model vector abundance V(t) per day; REQUIRED for
%                        caliber 'V' (it must come from a solve, not be assumed)
%       opts.coef        optional thermal-coefficient override (the analysis trait
%                        sensitivity); forwarded to mwspm_thermal_traits, OFF by
%                        default.  R0(t) IS affected by the trait curves, so the
%                        override changes a genuinely T-dependent quantity.
%
%   CALIBER (the analysis, resolved by T08 as: default N_v):
%     'Nv'  m(t) = m0*s(t)          -> "vector-supply caliber" (consistent with
%                                      the initial version of the paper)
%     'V'   m(t) = V(t)/N_h         -> "in-model abundance caliber"; needs the
%                                      solved V(t), and using it is a MODEL
%                                      DEFINITION CHANGE (model_spec.md Sec.16).
%
%   R fields: R0 (n_days x 1), max_R0, peak_week, gt1_weeks, gt1_first_week,
%             gt1_last_week, caliber, mu_variant, warning (cellstr)

if nargin < 3 || isempty(opts)
    opts = struct();
end
muVariant = 'raw';
if isfield(opts, 'mu_variant') && ~isempty(opts.mu_variant)
    muVariant = lower(char(opts.mu_variant));
end
caliber = 'Nv';
if isfield(opts, 'caliber') && ~isempty(opts.caliber)
    caliber = char(opts.caliber);
end

warningText = {};
if ~strcmp(muVariant, 'raw')
    warningText{end+1} = ['R0 was computed with the daily-corrected mortality variant (' ...
        muVariant '); this is NOT the continuous-rate quantity and must not be reported as R0.'];
end

s = ones(D.n_days, 1);
if strcmp(p.season_mode, 'window')
    s = mwspm_season(D.doy, p.season);
end

switch caliber
    case 'Nv'
        mVec = p.m0 * s;
    case 'V'
        assert(isfield(opts, 'V') && ~isempty(opts.V), 'mwspm:r0', ...
            ['caliber ''V'' requires opts.V (the solved in-model vector abundance). ' ...
             'V(t) is a first-order lag filter of N_v(t) and must not be assumed equal to it (the statement).']);
        V = double(opts.V(:));
        assert(numel(V) == D.n_days, 'mwspm:r0', 'opts.V must have n_days entries');
        mVec = V / p.N_h;
        warningText{end+1} = ['caliber ''V'' changes the model definition (model_spec.md Sec.16): ' ...
            'downstream tasks depending on R0 must be reopened.'];
    otherwise
        error('mwspm:r0', 'unknown caliber "%s"', caliber);
end

trRawOpts = struct();
if isfield(opts, 'coef') && ~isempty(opts.coef)
    trRawOpts.coef = opts.coef;
end
if isfield(opts, 'alt') && ~isempty(opts.alt)
    trRawOpts.alt = opts.alt;   % C1: alternative functional forms
end
trRaw  = mwspm_thermal_traits(D.T_daily, 'raw', trRawOpts);
muEff  = p.mu_s * trRaw.mu;
R0     = mwspm_r0(D.T_daily, mVec, muEff, p.gamma_h, trRawOpts);

R = struct();
R.R0            = R0;
R.max_R0        = max(R0);
R.peak_week     = mwspm_day2week(find(R0 == max(R0), 1), p);
R.above1        = R0 > 1;
R.gt1_weeks     = unique(mwspm_day2week(find(R.above1)', p));
if any(R.above1)
    R.gt1_first_week = min(R.gt1_weeks);
    R.gt1_last_week  = max(R.gt1_weeks);
else
    R.gt1_first_week = NaN;
    R.gt1_last_week  = NaN;
end
R.n_days_above1 = sum(R.above1);
R.caliber       = caliber;
R.mu_variant    = muVariant;
R.s             = s;
R.m             = mVec;
R.coef          = trRaw.coef;
R.coef_override = trRaw.coef_override;
R.warning       = warningText;
end

% ---------------------------------------------------------------------------
function w = mwspm_day2week(dayIdx, p)
% day 1 -> week 1, day 7 -> week 1, day 8 -> week 2 ...
if isempty(dayIdx) || ~isfinite(dayIdx(1))
    w = NaN;
    return
end
w = floor((double(dayIdx) - 1) / 7) + 1;
w = min(max(w, 1), p.n_weeks);
end
