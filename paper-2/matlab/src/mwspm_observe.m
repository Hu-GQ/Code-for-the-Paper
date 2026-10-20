function O = mwspm_observe(tDaily, CumVec, p, D)
%MWSPM_OBSERVE  Map true incidence to weekly reported counts and score the fit.
%
%   O = mwspm_observe(tDaily, CumVec, p, D)
%
%   tDaily  report grid including every integer day 0..n_days
%   CumVec  cumulative TRUE local incidence (state 8) on that grid
%
%   Weekly aggregation is done by differencing the cumulative accumulator at
%   the exact weekly edges t = 0,7,...,329.  This is quadrature-EXACT for the
%   accumulator and avoids the trapezoidal error that would come from
%   integrating sigma_h*E_h numerically after the fact.
%
%   Observation model (theory/model_spec.md Sec.12):
%       Y_w | dC_{w-d} ~ Binomial(dC_{w-d}, rho),  E[Y_w | .] = rho*dC_{w-d}
%   The scalar d-week lag is applied by shifting the weekly true-incidence
%   series; it is an OBSERVATION-layer lag and is NOT a delay inside the ODE
%   (obligation O9).
%
%   Output O fields: dC, dC_delayed, Yhat, Y, residuals, and the metrics
%   R2_log, R2_raw, RMSE_log, poisson_nll, totals, peak weeks.

tDaily = double(tDaily(:));
CumVec = double(CumVec(:));
assert(numel(tDaily) == numel(CumVec), 'mwspm:observe', 'time/state length mismatch');

assert(p.n_weeks * 7 == p.n_days, 'mwspm:observe', ...
    'n_weeks*7 must equal n_days for exact weekly edges (got %d*7 vs %d)', p.n_weeks, p.n_days);

edges = (0:p.n_weeks) * 7;                        % 48 edges for 47 weeks
CumEdge = interp1(tDaily, CumVec, edges, 'linear');
assert(~any(isnan(CumEdge)), 'mwspm:observe', 'report grid does not cover all weekly edges');

% Force a column: interp1 echoes the orientation of the query points, and the
% lag shift below vertically concatenates, which requires a column here.
dC = diff(CumEdge(:));                            % 47x1 true weekly incidence

d = round(p.d_weeks);
d = min(max(d, 0), p.n_weeks);
dCdelayed = [zeros(d, 1); dC(1:p.n_weeks - d)];

Yhat = p.rho * dCdelayed;
Y    = D.Y_weekly(:);

O = struct();
O.week_index  = (1:p.n_weeks)';
O.dC          = dC;
O.dC_delayed  = dCdelayed;
O.Yhat        = Yhat;
O.Y           = Y;
O.residual    = Y - Yhat;
O.d_weeks     = d;

O.metrics = mwspm_metrics(Y, Yhat);
O.metrics.true_total   = sum(dC);
O.metrics.reported_total_model = sum(Yhat);
O.metrics.reported_total_observed = sum(Y);
end
