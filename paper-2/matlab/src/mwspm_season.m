function [s, info] = mwspm_season(doy, sp)
%MWSPM_SEASON  Phenomenological relative vector abundance s(t).
%
%   [s, info] = mwspm_season(doy, sp)
%
%   sp fields: onset_doy, offset_doy, kappa_days, gamma, beta, s_floor,
%              min_window_days, normalize (logical), norm_factor (optional)
%
%   s_hat(t) = W(t) * exp(gamma*(clip(u,-0.5,1.5)-0.5)) * exp(beta*4*uhat*(1-uhat))
%   W(t)     = sigmoid((doy-onset)/kappa) * (1 - sigmoid((doy-offset)/kappa))
%   u        = (doy-onset)/(offset-onset)
%   s        = max(s_hat/mean(s_hat), s_floor) with the mean taken over the
%              analysis window (the caller passes the full daily grid once).
%
%   IMPORTANT (theory/model_spec.md Sec.7.3, the statement): s(t) is the EXTERNALLY
%   SET vector supply/driver.  The in-model abundance V(t) is a first-order
%   lag filter of N_v(t)=m0*N_h*s(t) and is NOT equal to it.  Never describe
%   s(t) as "the number of mosquitoes in the model".
%
%   Normalisation contract: when sp.normalize is true the mean is computed on
%   the vector `doy` passed in.  Continuous-time use must therefore pass the
%   full analysis-window grid (or a precomputed sp.norm_factor) so that the
%   normalisation does not depend on the sampling of t.

doy = double(doy(:));

on  = double(sp.onset_doy);
off = max(double(sp.offset_doy), on + double(sp.min_window_days));
k   = double(sp.kappa_days);

z1 = (doy - on)  / k;
z2 = (doy - off) / k;
W  = mwspm_sigmoid(z1) .* (1 - mwspm_sigmoid(z2));

span = max(off - on, 1);
u    = (doy - on) / span;
uc   = min(max(u, -0.5), 1.5);
uh   = min(max(u,  0.0), 1.0);

shat = W .* exp(sp.gamma * (uc - 0.5)) .* exp(sp.beta * 4 * uh .* (1 - uh));

info = struct();
info.onset_doy  = on;
info.offset_doy = off;
info.span_days  = off - on;
info.kappa_days = k;

normFactor = 1.0;
normalised = false;
if isfield(sp, 'normalize') && sp.normalize
    if isfield(sp, 'norm_factor') && ~isempty(sp.norm_factor) && isfinite(sp.norm_factor)
        normFactor = double(sp.norm_factor);
    else
        m = mean(shat);
        if m > 0
            normFactor = 1 / m;
        end
    end
    normalised = true;
end

info.norm_factor = normFactor;
info.normalised  = normalised;
info.mean_shat   = mean(shat);

s = shat * normFactor;
floorVal = 1e-8;
if isfield(sp, 's_floor') && ~isempty(sp.s_floor)
    floorVal = double(sp.s_floor);
end
s = max(s, floorVal);
end

% ---------------------------------------------------------------------------
function y = mwspm_sigmoid(z)
y = 1 ./ (1 + exp(-min(max(z, -500), 500)));
end
