function e = mwspm_env()
%MWSPM_ENV  Bounded environment snapshot for manifest.json.
%
%   Records the MATLAB release, platform, JVM availability, the products the
%   project actually relies on and whether the solver / optimiser entry points
%   exist.  It does NOT dump the full environment, PATH or any credential
%   (the project record hard constraint 5).

e.matlab_version = version();
e.release        = version('-release');
e.computer       = computer();
e.matlabroot     = matlabroot;
e.jvm            = usejava('jvm');
e.captured_at    = datestr(now, 'yyyy-mm-ddTHH:MM:SS');   %#ok<TNOW1,DATST>
e.timezone       = 'local (UTC+08:00 assumed for run_id)';

v = ver;
e.toolbox_count = numel(v);
names = cell(numel(v), 1);
for k = 1:numel(v)
    names{k} = sprintf('%s|%s|%s', v(k).Name, v(k).Version, v(k).Release);
end
e.toolboxes = names;

% entry points the project depends on (existence only, never executed here)
probe = { ...
    'ode45', 'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'ode113', 'ode78', 'ode89', ...
    'odeset', 'deval', 'interp1', 'cumtrapz', 'rng', 'RandStream'};
e.entry_points = struct();
for k = 1:numel(probe)
    e.entry_points.(probe{k}) = (exist(probe{k}, 'file') > 0) || (exist(probe{k}, 'class') == 8);
end

% products the project must NOT assume (recorded as a fact, per SKILL_MAP)
st = false;
for k = 1:numel(v)
    if contains(lower(v(k).Name), 'statistics and machine learning')
        st = true;
    end
end
e.statistics_toolbox_present = st;
e.statistics_toolbox_note = ['Recorded for honesty: statistical distributions, bootstrap and ' ...
    'negative-binomial likelihoods must be hand-implemented in src/matlab/ when this is false.'];

lic = license('inuse');
if isstruct(lic)
    e.licenses_in_use = {lic.feature};
else
    e.licenses_in_use = {};
end
end
