clear; clc;

lambda = 1.55e-6;
k0 = 2*pi/lambda;
n_co = 1.537;
n_sub = 1.4442;
n_clad = 1.511;
a_vec = linspace(0.2e-6, 5.0e-6, 200);
maxLat = 5;
maxVer = 5;

modeCurves = {};
for m_lat = 0:maxLat
    for m_ver = 0:maxVer
         Neff_mode = nan(size(a_vec));
         for idx = 1:length(a_vec)
             a = a_vec(idx);
             b = a;
             fun_lat = @(N_x) lateral_TE(N_x, n_co, n_clad, a, k0, m_lat);
             try
                 N_x = fzero(fun_lat, [n_clad, n_co]);
             catch
                 N_x = NaN;
             end
             if isnan(N_x) || N_x < n_clad || N_x > n_co
                 Neff_mode(idx) = NaN;
                 continue;
             end
             N_lower = max(n_sub, n_clad);
             fun_ver = @(N) vertical_TM(N, N_x, n_sub, n_clad, b, k0, m_ver);
             try
                 N_final = fzero(fun_ver, [N_lower, N_x]);
             catch
                 N_final = NaN;
             end
             Neff_mode(idx) = N_final;
         end
         modeCurves{end+1} = struct('m_lat', m_lat, 'm_ver', m_ver, 'a', a_vec, 'Neff', Neff_mode);
     end
end

modeCount = zeros(size(a_vec));
for j = 1:numel(modeCurves)
    Neff_j = modeCurves{j}.Neff;
    modeCount(~isnan(Neff_j)) = modeCount(~isnan(Neff_j)) + 1;
end

single_mode_indices = find(modeCount == 1);
dual_mode_indices = find(modeCount == 2);

if ~isempty(single_mode_indices)
    a_single_min = a_vec(min(single_mode_indices)) * 1e6;
    a_single_max = a_vec(max(single_mode_indices)) * 1e6;
    fprintf('Single-mode range: a ∈ [%.3f, %.3f] µm\n', a_single_min, a_single_max);
else
    fprintf('No single-mode region found.\n');
end

if ~isempty(dual_mode_indices)
    a_dual_min = a_vec(min(dual_mode_indices)) * 1e6;
    a_dual_max = a_vec(max(dual_mode_indices)) * 1e6;
    fprintf('Dual-mode range: a ∈ [%.3f, %.3f] µm\n', a_dual_min, a_dual_max);
else
    fprintf('No dual-mode region found.\n');
end

if ~isempty(single_mode_indices)
    a_single_mid_idx = round(mean(single_mode_indices));
    a_single_mid = a_vec(a_single_mid_idx) * 1e6;
    Neff_single_mid = modeCurves{1}.Neff(a_single_mid_idx);
    fprintf('Waveguide size at middle of single-mode range: a = %.3f µm\n', a_single_mid);
    fprintf('Effective index at this point: N_eff = %.4f\n', Neff_single_mid);
else
    fprintf('Cannot calculate center effective index due to absence of single-mode region.\n');
end

figure;
hold on;
colors = lines(numel(modeCurves));
legend_entries = {};
for j = 1:numel(modeCurves)
    curve = modeCurves{j};
    if any(~isnan(curve.Neff))
       plot(curve.a*1e6, curve.Neff, 'Color', colors(j,:), 'LineWidth', 2);
       legend_entries{end+1} = sprintf('m_{ver}=%d, m_{lat}=%d', curve.m_ver, curve.m_lat);
    end
end

num_modes_per_a = zeros(size(a_vec));
for j = 1:numel(modeCurves)
    curve = modeCurves{j};
    valid_indices = ~isnan(curve.Neff);
    num_modes_per_a(valid_indices) = num_modes_per_a(valid_indices) + 1;
end

yl = ylim;
single_mode_indices = find(num_modes_per_a == 1);
if ~isempty(single_mode_indices)
    x_patch = [a_vec(single_mode_indices)*1e6, fliplr(a_vec(single_mode_indices)*1e6)];
    y_patch = [ones(size(single_mode_indices))*yl(1), fliplr(ones(size(single_mode_indices))*yl(2))];
    fill(x_patch, y_patch, [0.9 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    x_mid = mean(a_vec(single_mode_indices)) * 1e6;
    y_text = yl(1) + 0.015;
    text(x_mid, y_text, 'Single-mode region', 'Color', 'b', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

dual_mode_indices = find(num_modes_per_a == 2);
if ~isempty(dual_mode_indices)
    x_patch = [a_vec(dual_mode_indices)*1e6, fliplr(a_vec(dual_mode_indices)*1e6)];
    y_patch = [ones(size(dual_mode_indices))*yl(1), fliplr(ones(size(dual_mode_indices))*yl(2))];
    fill(x_patch, y_patch, [1.0 0.85 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    x_mid = mean(a_vec(dual_mode_indices)) * 1e6;
    y_text = yl(1) + 0.015;
    text(x_mid, y_text, 'Dual-mode region', 'Color', [0.85 0.33 0.1], 'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Square waveguide side length a (\mum)', 'FontSize', 12);
ylabel('Effective refractive index N_{eff}', 'FontSize', 12);
title('Effective refractive index dispersion curves for square waveguide modes', 'FontSize', 10, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'eastoutside');
grid on;

% --- Functions ---
function F = lateral_TE(N_x, n_co, n_clad, a, k0, m_lat)
    if N_x < n_clad || N_x > n_co
       F = NaN;
       return;
    end
    beta = k0 * sqrt(n_co^2 - N_x^2);
    gamma = k0 * sqrt(N_x^2 - n_clad^2);
    F = k0*a*sqrt(n_co^2 - N_x^2) - ( m_lat*pi + 2*atan( sqrt((N_x^2 - n_clad^2)/(n_co^2 - N_x^2)) ) );
end

function F = vertical_TM(N, N_x, n_sub, n_clad, b, k0, m_ver)
    if N < max(n_sub, n_clad) || N > N_x
       F = NaN;
       return;
    end
    beta = k0 * sqrt(N_x^2 - N^2);
    gamma_sub = k0 * sqrt(N^2 - n_sub^2);
    gamma_clad = k0 * sqrt(N^2 - n_clad^2);
    F = k0*b*sqrt(N_x^2 - N^2) - ( m_ver*pi + ...
        atan( (N_x^2/n_sub^2)*sqrt((N^2 - n_sub^2)/(N_x^2 - N^2) ) ) + ...
        atan( (N_x^2/n_clad^2)*sqrt((N^2 - n_clad^2)/(N_x^2 - N^2) ) ) );
end
