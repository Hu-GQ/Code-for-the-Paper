function varargout = mwspm_calib_theta(action, varargin)
%MWSPM_CALIB_THETA  Parameterisation of the T12 calibration audit.
%
%   th   = mwspm_calib_theta('reference')            existing estimate (model card)
%   box  = mwspm_calib_theta('box')                  search boxes + provenance
%   nm   = mwspm_calib_theta('names')                canonical parameter names
%   z    = mwspm_calib_theta('pack',   th, names, box)
%   th   = mwspm_calib_theta('unpack', z,  names, box, base)
%
%   Single source for "which parameters exist, over what range, and how a
%   bounded problem is mapped onto an unbounded one".  Every optimiser in T12
%   goes through pack/unpack, so the boxes cannot drift between the profile,
%   the holdout and the synthetic-recovery runs.
%
%   Bounded -> unbounded map: for a field with box [lo, hi],
%       u = lo + (hi - lo) * sigmoid(z),    z = logit((u - lo)/(hi - lo)).
%   The map is monotone and never leaves the box, so box constraints are
%   enforced by construction rather than by a penalty.
%
%   Provenance of the boxes: code/dengue_fit.py Spec.bounds() as expanded by
%   code/step12_rho_free.py _expansion(); rho is replaced by the source model's
%   hard constraint [0.02, 0.80] (run_01_parameter_estimation.m:
%   rho = 0.02 + sigmoid(x2)*0.78).  The profile scan relaxes the upper bound to
%   1.00 for the SHAPE of the likelihood only (never as the main caliber).

switch lower(char(action))
    case 'names'
        varargout{1} = {'m0', 'rho', 'd_weeks', 'onset_doy', 'offset_doy', ...
                        'mu_s', 'gamma', 'beta'};

    case 'box'
        % Names are the kernel's own field names, so a packed/unpacked theta can
        % be handed straight to mwspm_calib_kernel.  The legacy optimiser used
        % an unbounded log m0 in [-7, 3]; the equivalent box is m0 in
        % [exp(-7), exp(3)] = [9.12e-4, 20.09].
        b = struct();
        b.m0        = [9.118819655545162e-04, 20.085536923187668];
        b.rho       = [0.02, 0.80];    % source hard constraint
        b.d_weeks   = [0.0, 3.0];      % discretised by round() in the kernel
        b.onset_doy = [40.0, 280.0];
        b.offset_doy = [80.0, 365.0];
        b.mu_s      = [0.25, 15.0];
        b.gamma     = [-3.0, 2.0];
        b.beta      = [-3.0, 3.0];
        b.provenance = ['code/dengue_fit.py Spec.bounds() with log_m0 in [-7,3] ' ...
            'mapped to m0 = exp(log_m0); rho from the source parameterisation ' ...
            'rho = 0.02 + sigmoid(x2)*0.78 (code/step12_rho_free.py)'];
        varargout{1} = b;

    case 'reference'
        % The existing parameter estimate under audit: the main report point of
        % output/results/refit_v3.json (headline), mirrored in the model card
        % config/experiments/model_mwspm_gen3_rhohat.json.
        th = struct();
        th.m0         = 1.6637;
        th.rho        = 0.8000;
        th.d_weeks    = 2.0;
        th.onset_doy  = 180.668;
        th.offset_doy = 253.456;
        th.mu_s       = 12.3044;
        th.gamma      = -1.5338;
        th.beta       = -1.8797;
        varargout{1} = th;

    case 'pack'
        [th, names, box] = deal(varargin{1:3});
        names = cellstr(string(names));
        n = numel(names);
        z = zeros(n, 1);
        for k = 1:n
            f = names{k};
            assert(isfield(th, f), 'mwspm:calib', 'field "%s" missing from theta', f);
            z(k) = zfromu(double(th.(f)), box.(f));
        end
        varargout{1} = z;

    case 'unpack'
        [z, names, box, base] = deal(varargin{1:4});
        names = cellstr(string(names));
        z = double(z(:));
        assert(numel(z) == numel(names), 'mwspm:calib', ...
            'z has %d entries but %d names were given', numel(z), numel(names));
        th = base;
        for k = 1:numel(names)
            th.(names{k}) = ufromz(z(k), box.(names{k}));
        end
        varargout{1} = th;

    otherwise
        error('mwspm:calib', 'unknown action "%s"', char(action));
end
end

% ---------------------------------------------------------------------------
function z = zfromu(u, box)
lo = box(1); hi = box(2);
w = (u - lo) / (hi - lo);
w = min(max(w, 1e-9), 1 - 1e-9);
z = log(w / (1 - w));
end

function u = ufromz(z, box)
lo = box(1); hi = box(2);
s = 1.0 / (1.0 + exp(-min(max(double(z), -60), 60)));
u = lo + (hi - lo) * s;
end
