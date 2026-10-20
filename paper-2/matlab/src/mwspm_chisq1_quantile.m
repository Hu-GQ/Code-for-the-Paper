function q = mwspm_chisq1_quantile(p)
%MWSPM_CHISQ1_QUANTILE  Quantile of the chi-square distribution with 1 d.o.f.
%
%   q = mwspm_chisq1_quantile(0.95)   -> 3.8415
%
%   Implemented from the identity chi2_1 = Z^2 and
%       P(chi2_1 <= q) = 2*Phi(sqrt(q)) - 1 = p
%   so that
%       sqrt(q) = Phi^{-1}((1+p)/2) = -sqrt(2) * erfcinv(1+p)   =>   q = 2*erfcinv(1+p)^2.
%
%   erfcinv is part of base MATLAB (specfun); the Statistics and Machine
%   Learning Toolbox, which owns chi2inv, is NOT installed on this
%   installation.  Using a self-contained formula keeps the profile-likelihood
%   threshold reproducible without a hidden toolbox dependency.

assert(isscalar(p) && p > 0 && p < 1, 'mwspm:stats', 'p must be in (0,1)');
q = 2.0 * erfcinv(1.0 + p)^2;
end
