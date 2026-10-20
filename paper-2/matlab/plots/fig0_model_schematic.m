function fig0_model_schematic()
%FIG0_MODEL_SCHEMATIC  Model schematic for the M-WSPM manuscript (T19 prereq O-03).
%
%   Deterministic drawing (no run data, no solver, no randomness; AI-generated
%   artwork is not used). Shows: human SEIR column + mosquito SEI column, the
%   two cross-infection fluxes (standard incidence), external forcing
%   (imported cases Lambda(t); seasonal supply N_v = m0*N_h*s(t); daily mean
%   temperature T(t) driving the trait curves), and the observation layer
%   (reporting ratio rho, lag d).
%
%   Exports: figures/main/fig0_model_schematic.pdf + figures/preview/fig0_model_schematic.png
%   (via pp_util('savefig')). Reads no saved results => no read-only snapshot.

pp_util('style');
pal = pp_util('palette');

W = 16; H = 9;                                  % canvas data units
fig = figure('Units', 'centimeters', 'Position', [1 1 17.0 9.6]); %#ok<NASGU>
ax = axes('Units', 'normalized', 'Position', [0 0 1 1]); hold(ax, 'on');
axis(ax, 'off'); xlim(ax, [0 W]); ylim(ax, [0 H]);

blue = pal.blue; red = pal.red; gr = pal.gray;
blueL = [0.90 0.94 0.99]; grL = [0.93 0.93 0.93]; redL = [0.985 0.92 0.92];

xb = 1.15;  w = 3.9;  hb = 1.02;               % human boxes
xm = 10.95; wm = 3.9;                           % mosquito boxes

% ---------- human column (top to bottom) ----------
rbox(xb, 6.55, w, hb, blueL, blue); txt(xb+w/2, 7.06, 'S_h  (susceptible)', blue, 9);
rbox(xb, 5.05, w, hb, blueL, blue); txt(xb+w/2, 5.56, 'E_h  (exposed)',     blue, 9);
rbox(xb, 3.55, w, hb, blueL, blue); txt(xb+w/2, 4.06, 'I_h  (infectious)',  blue, 9);
rbox(xb, 2.05, w, hb, blueL, blue); txt(xb+w/2, 2.56, 'R_h  (recovered)',   blue, 9);
ar(xb+w/2, 6.55, xb+w/2, 6.09, blue, 0.9); txt(xb+w/2+1.60, 6.32, '\lambda_h S_h', blue, 7.5);
ar(xb+w/2, 5.05, xb+w/2, 4.59, blue, 0.9); txt(xb+w/2+1.15, 4.82, '\sigma_h E_h',  blue, 7.5);
ar(xb+w/2, 3.55, xb+w/2, 3.09, blue, 0.9); txt(xb+w/2+1.10, 3.32, '\gamma_h I_h',  blue, 7.5);
txt(xb+w/2, 8.05, 'Humans  (SEIR)', blue, 10);

% imported cases into I_h
ar(0.30, 4.06, 1.12, 4.06, red, 1.1);
txt(0.85, 4.80, {'Imported', 'cases  \Lambda(t)'}, red, 7);

% ---------- mosquito column ----------
rbox(xm, 6.55, wm, hb, grL, gr); txt(xm+wm/2, 7.06, 'S_v  (susceptible)', gr, 9);
rbox(xm, 5.05, wm, hb, grL, gr); txt(xm+wm/2, 5.56, 'E_v  (exposed)',     gr, 9);
rbox(xm, 3.55, wm, hb, grL, gr); txt(xm+wm/2, 4.06, 'I_v  (infectious)',  gr, 9);
ar(xm+wm/2, 6.55, xm+wm/2, 6.09, gr, 0.9); txt(xm+wm/2+1.60, 6.32, '\lambda_v S_v', gr, 7.5);
ar(xm+wm/2, 5.05, xm+wm/2, 4.59, gr, 0.9); txt(xm+wm/2+1.15, 4.82, '\sigma_v E_v',  gr, 7.5);
txt(12.00, 8.05, 'Mosquitoes  (SEI; no recovery)', gr, 10);

% seasonal supply (birth) and deaths
ar(13.90, 8.36, 13.90, 7.62, gr, 1.1);
txt(13.90, 8.72, 'N_v(t) = m_0 N_h s(t)', gr, 7.5);
txt(13.90, 8.56, '(seasonal supply)', gr, 7);
for yy = [7.06 5.56 4.06]
    ar(xm+wm+0.06, yy, xm+wm+0.50, yy, gr, 0.9);
end
text(ax, 15.68, 5.56, 'deaths  \mu_v^{eff}', 'Rotation', 90, 'Color', gr, ...
    'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% ---------- cross infections (standard incidence) ----------
ar(xm+0.02, 4.06, xb+w+0.02, 6.87, blue, 1.2);   % I_v -> S_h
txt(8.30, 6.45, '\lambda_h = a(T) b(T) I_v / N_h', blue, 7.5);
ar(xb+w+0.02, 4.06, xm+0.02, 6.87, red, 1.2);    % I_h -> S_v
txt(8.35, 4.30, '\lambda_v = a(T) c(T) I_h / N_h', red, 7.5);

% ---------- temperature note (bottom left) ----------
txt(0.85, 1.06, {'Daily mean temperature  T(t)  (ERA5) drives', ...
    'a(T), b(T), c(T), \sigma_v(T), \mu_v(T):', ...
    '\mu_v^{eff} = \mu_s \mu_v(T)   (\mu_s: fitted mortality multiplier)'}, gr, 7, 'left');

% ---------- observation layer ----------
rbox(9.55, 2.30, 5.95, 1.12, redL, red);
txt(12.53, 2.86, {'Observed weekly cases  Y_w', ...
    'Y_w \sim Binomial(\Delta C_{w-d}, \rho)  \approx  Poisson(\rho \Delta C_{w-d})', ...
    '\rho: reporting ratio (fitted);   d: reporting lag (weeks)'}, red, 7);
% dashed path: true new infections -> observation box
plot(ax, [1.13 0.55], [3.80 3.80], '--', 'Color', red, 'LineWidth', 0.9);
plot(ax, [0.55 0.55], [3.80 1.90], '--', 'Color', red, 'LineWidth', 0.9);
plot(ax, [0.55 9.30], [1.90 1.90], '--', 'Color', red, 'LineWidth', 0.9);
plot(ax, [9.30 9.30], [1.90 2.86], '--', 'Color', red, 'LineWidth', 0.9);
ar(9.30, 2.86, 9.53, 2.86, red, 0.9);
txt(6.85, 2.14, 'true new infections  C(t) = \int \sigma_h E_h dt', red, 7);

pp_util('savefig', fig, 'fig0_model_schematic');
fprintf('fig0: done (deterministic schematic; no saved results read).\n');

% =====================================================================
% nested helpers (axes `ax` and canvas constants assigned before first use)
% =====================================================================
    function rbox(x1, y1, ww, hh, fc, ec)
        rectangle('Position', [x1 y1 ww hh], 'Curvature', 0.14, ...
            'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', 0.9, 'Parent', ax);
    end
    function ar(x1, y1, x2, y2, col, lw)
        annotation(fig, 'arrow', [x1 x2] / W, [y1 y2] / H, ...
            'Color', col, 'LineWidth', lw, 'HeadWidth', 7, 'HeadLength', 8);
    end
    function txt(x, y, s, col, fs, align)
        if nargin < 6; align = 'center'; end
        if strcmp(align, 'left'); ha = 'left'; else; ha = 'center'; end
        text(ax, x, y, s, 'Color', col, 'FontSize', fs, 'FontWeight', 'normal', ...
            'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
    end
end
