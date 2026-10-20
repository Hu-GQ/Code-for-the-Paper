function L = mwspm_calib_nll(y, mu, family, opts)
%MWSPM_CALIB_NLL  Count-data likelihoods for the T12 calibration audit.
%
%   L = mwspm_calib_nll(y, mu, family)
%   L = mwspm_calib_nll(y, mu, family, struct('df', 41, 'theta', 5))
%
%   y      (n x 1) observed weekly counts
%   mu     (n x 1) model mean (rho * weekly true incidence at lag d)
%   family 'poisson' | 'quasipoisson' | 'negbin' | 'gaussian_sqrt'
%
%   The Statistics and Machine Learning Toolbox is NOT installed on this
%   installation (SKILL_MAP.md Sec.6, re-verified by mwspm_env), so no
%   distribution object is used: every log-likelihood is written out and the
%   normalising constants are kept explicit.
%
%   Why the family matters here (the point of the check)
%   ---------------------------------------------------
%   Whether weekly dengue counts are OVERDISPERSED relative to Poisson is the
%   QUESTION this function exists to answer, so it is deliberately NOT assumed
%   in either direction.  Poisson (or, worse, "errors are independent Gaussian")
%   understates uncertainty whenever the counts really are overdispersed; but
%   the opposite failure - reporting an inflated variance the data do not
%   support - is just as wrong.  L therefore returns the Pearson dispersion
%   X2/df alongside the log-likelihood, so the adequacy of the Poisson
%   assumption is a MEASURED quantity rather than an assumption.
%
%   MEASURED (T12, results/summaries/calibration.md): at the reference estimate
%   X2/df = 1.046 with df = 46, and the NB2 size MLE is theta = 71.6 with only a
%   0.86 nll improvement.  On THIS data set the counts are therefore very close
%   to Poisson, and choosing a heavier-tailed family would not change the
%   conclusions.  The contrast families are kept because that had to be
%   established, not presumed.
%
%   Fields
%   ------
%     family, n, mu (clamped), y
%     nll        negative log-likelihood in the parameter-dependent part only
%                (the legacy project convention, so it is comparable with the
%                 objective in output/results/refit_v3.json)
%     nll_full   negative log-likelihood INCLUDING the parameter-free
%                gammaln(y+1) term, so that AIC/BIC are comparable ACROSS
%                families (Poisson is the theta->inf limit of NB2)
%     pearson, df, dispersion = pearson/df
%     theta      NB2 size (only for 'negbin')
%     scale      quasi-Poisson scale (only for 'quasipoisson')
%     notes

if nargin < 4 || isempty(opts)
    opts = struct();
end

y  = double(y(:));
mu = double(mu(:));
assert(numel(y) == numel(mu), 'mwspm:calib', 'y and mu length mismatch');

tiny = 1e-12;
muC  = max(mu, tiny);
n    = numel(y);

L = struct();
L.family = lower(char(family));
L.n      = n;
L.y      = y;
L.mu     = mu;
L.theta  = NaN;
L.scale  = NaN;
L.n_extra = 0;                        % parameters beyond the model fit (see gaussian_sqrt)
L.notes  = '';

const = sum(gammaln(y + 1));          % parameter-free for every family below

switch L.family
    case 'poisson'
        L.nll      = sum(muC - y .* log(muC));
        L.nll_full = L.nll + const;
        L.notes    = ['Poisson NLL = sum(mu - y*log mu); var = mu assumed. ' ...
                      'Overdispersion is NOT modelled.'];

    case 'quasipoisson'
        % scale is supplied (or estimated by the caller from the Pearson X2);
        % the point estimate of mu is unchanged, only the curvature is.
        scale = 1.0;
        if isfield(opts, 'scale') && ~isempty(opts.scale) && isfinite(opts.scale)
            scale = double(opts.scale);
        end
        L.scale    = scale;
        L.nll      = sum(muC - y .* log(muC)) / scale;
        L.nll_full = L.nll + const / scale;
        L.notes    = ['quasi-Poisson: same mean model, variance inflated by the ' ...
                      'estimated scale; NLL divided by the scale.'];

    case 'negbin'
        % NB2: var = mu + mu^2/theta.  theta -> inf recovers Poisson.
        th = NaN;
        if isfield(opts, 'theta') && ~isempty(opts.theta)
            th = double(opts.theta);
        end
        assert(isfinite(th) && th > 0, 'mwspm:calib', ...
            'negbin requires a positive size theta (got %g)', th);
        L.theta = th;
        L.nll = -sum( gammaln(y + th) - gammaln(th) ...
                      + th .* (log(th) - log(th + muC)) ...
                      + y .* (log(muC) - log(th + muC)) );
        L.nll_full = L.nll + const;
        L.notes    = 'NB2 via gammaln; var = mu + mu^2/theta.';

    case 'gaussian_sqrt'
        % Anscombe-variance-stabilised Gaussian contrast.  Kept ONLY as a
        % contrast so that "we did not silently default to an independent
        % Gaussian error model" is a verifiable statement.
        %
        % TWO corrections are needed before this row may be compared with the
        % count families, and BOTH were missing in the first full run, where the
        % raw version produced AIC 112 against the Poisson 182 - i.e. it made an
        % independent-Gaussian error model look BEST, which is exactly the
        % conclusion the check exists to test:
        %   (a) a sqrt-transformed Gaussian is a density for z = sqrt(y+3/8),
        %       NOT for y, so the change-of-variables Jacobian
        %       sum(log(2*sqrt(y+3/8))) is added back;
        %   (b) the residual variance is profiled by maximum likelihood
        %       (sigma2 = RSS/n) instead of being fixed at 1, which costs one
        %       extra parameter and is reported in L.n_extra so the caller's
        %       AIC/BIC parameter count stays honest.
        z   = sqrt(y + 3/8);
        m   = sqrt(muC + 3/8);
        rss = sum((z - m).^2);
        s2  = max(rss / n, tiny);   % guard: a perfect fit would give log(0)
        jac = sum(log(2 * sqrt(y + 3/8)));
        L.nll      = 0.5 * rss;               % raw sum of squares, NOT comparable
        L.nll_full = 0.5 * n * log(2 * pi * s2) + 0.5 * n + jac;
        L.n_extra  = 1;
        L.scale    = s2;
        L.notes    = ['Anscombe sqrt + i.i.d. Gaussian error with the residual ' ...
                      'variance profiled by MLE and the sqrt Jacobian restored, so ' ...
                      'nll_full is a density for y and its AIC/BIC ARE comparable ' ...
                      'with the count families (plus one parameter for sigma2).'];

    otherwise
        error('mwspm:calib', 'unknown likelihood family "%s"', L.family);
end

% Pearson dispersion of the mean model (family independent)
L.pearson    = sum((y - mu).^2 ./ muC);
L.df         = n - 1;
if isfield(opts, 'df') && ~isempty(opts.df)
    L.df = double(opts.df);
end
L.dispersion = L.pearson / max(L.df, 1);
end
