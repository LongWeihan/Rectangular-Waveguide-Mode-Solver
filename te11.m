clear; clc;
lambda = 1.55e-6;
k0 = 2*pi/lambda;
n_co   = 1.5375;
n_sub  = 1.4442;
n_clad = 1.512;
a_vec = linspace(0.2e-6, 5.0e-6, 200);
maxVer = 5;
maxLat = 5;
modeCurves = {};

for m_ver = 0:maxVer
    for m_lat = 0:maxLat
         Neff_mode = nan(size(a_vec));
         for idx = 1:length(a_vec)
             a = a_vec(idx);  
             b = a;
             fun_vert = @(Ny) vertical_TE_TE(Ny, n_co, n_sub, n_clad, b, k0, m_ver);
             lower_vert = max(n_sub, n_clad);
             try
                 Ny = fzero(fun_vert, [lower_vert, n_co]);
             catch
                 Ny = NaN;
             end
             if isnan(Ny) || Ny < lower_vert || Ny > n_co
                 Neff_mode(idx) = NaN;
                 continue;
             end
             fun_horiz = @(N) horizontal_TM_TM(N, Ny, n_clad, a, k0, m_lat);
             try
                 N_final = fzero(fun_horiz, [n_clad, Ny]);
             catch
                 N_final = NaN;
             end
             Neff_mode(idx) = N_final;
         end
         modeCurves{end+1} = struct('m_ver', m_ver, 'm_lat', m_lat, 'a', a_vec, 'Neff', Neff_mode);
     end
end

modeCount = zeros(size(a_vec));
for j = 1:numel(modeCurves)
    Neff_j = modeCurves{j}.Neff;
    modeCount(~isnan(Neff_j)) = modeCount(~isnan(Neff_j)) + 1;
end

single_mode_indices = find(modeCount == 1);
dual_mode_indices   = find(modeCount == 2);

if ~isempty(single_mode_indices)
    a_single_min = a_vec(min(single_mode_indices)) * 1e6;
    a_single_max = a_vec(max(single_mode_indices)) * 1e6;
    fprintf('Single-mode transmission region: a ∈ [%.3f, %.3f] µm\n', a_single_min, a_single_max);
else
    fprintf('No single-mode transmission region found.\n');
end

if ~isempty(dual_mode_indices)
    a_dual_min = a_vec(min(dual_mode_indices)) * 1e6;
    a_dual_max = a_vec(max(dual_mode_indices)) * 1e6;
    fprintf('Dual-mode transmission region: a ∈ [%.3f, %.3f] µm\n', a_dual_min, a_dual_max);
else
    fprintf('No dual-mode transmission region found.\n');
end

if ~isempty(single_mode_indices)
    a_single_mid_idx = round(mean(single_mode_indices));
    a_single_mid = a_vec(a_single_mid_idx) * 1e6;
    Neff_single_mid = modeCurves{1}.Neff(a_single_mid_idx);
    fprintf('Mid-point of single-mode region: a = %.3f µm\n', a_single_mid);
    fprintf('Effective refractive index at this point: N_eff = %.4f\n', Neff_single_mid);
else
    fprintf('No single-mode region found. Cannot calculate mid-point effective index.\n');
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
title('Dispersion curves of effective refractive indices for different modes in square waveguide', 'FontSize', 10, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'eastoutside');
grid on;

function F = vertical_TE_TE(Ny, n_co, n_sub, n_clad, b, k0, m_ver)
    if Ny < max(n_sub, n_clad) || Ny > n_co
       F = NaN;
       return;
    end
    beta_vert = k0 * sqrt(n_co^2 - Ny^2);
    gamma_sub = k0 * sqrt(Ny^2 - n_sub^2);
    gamma_clad = k0 * sqrt(Ny^2 - n_clad^2);
    F = k0 * b * sqrt(n_co^2 - Ny^2) - ( m_ver*pi + atan( sqrt((Ny^2 - n_sub^2)/(n_co^2 - Ny^2)) ) + atan( sqrt((Ny^2 - n_clad^2)/(n_co^2 - Ny^2)) ) );
end

function F = horizontal_TM_TM(N, Ny, n_clad, a, k0, m_lat)
    if N < n_clad || N > Ny
       F = NaN;
       return;
    end
    beta_horiz = k0 * sqrt(Ny^2 - N^2);
    F = k0 * a * sqrt(Ny^2 - N^2) - ( m_lat*pi + 2*atan( (Ny^2/n_clad^2)*sqrt((N^2 - n_clad^2)/(Ny^2 - N^2)) ) );
end
