clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'third party codes', filesep, 'Fornberg'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
cd(current_directory)
warning off  % we turn off warning for numerical Vandermonde inverse
%              calculation in Savitzky-Golay method

%%
frequency_range = [.05:.1:1];
possible_l = [1: 10];

%%
min_rng = inf;
max_rng = -inf;
iteration_delta_x = 0;
for delta_x = 0;[-0.04 -0.02 0 0.02 0.04];
    iteration_delta_x = iteration_delta_x + 1;
    iteration_frequency = 0;
    for frq = frequency_range
        iteration_frequency = iteration_frequency + 1;
        [f, grad, hess, third_order_tensor, fourth_order_tensor, params] = sinusoidal_gratting_2D(frq, delta_x);
        N_x = params.N_x;
        N_y = params.N_y;
        h_x = params.h_x;
        h_y = params.h_y;
        
        %%
        partial_f{1}.staggered = grad.staggered{1};
        partial_f{1}.centralized = grad.centralized{1};
        partial_f{2}.staggered = hess.staggered{1};
        partial_f{2}.centralized = hess.centralized{1};
        partial_f{3}.staggered = third_order_tensor.staggered{1};
        partial_f{3}.centralized = third_order_tensor.centralized{1};
        partial_f{4}.staggered = fourth_order_tensor.staggered{1};
        partial_f{4}.centralized = fourth_order_tensor.centralized{1};
        
        %%  recall derivative matrices
        iteration = 0;
        for l = possible_l
            iteration = iteration + 1;
            [iteration_delta_x, iteration_frequency, iteration]
            %
            boundary_mask = logical(ones(128, 128));
            boundary_mask(l+1: end-l, l+1: end-l) = false;
            
            for n_ord = 1: 4;
                flag_symbolic = false;                
                if 2*l >= n_ord+1
                    %%  MaxPol
                    nod = 'staggered';
                    P = 2*l-1;
                    D = derivmtx(l, P, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.staggered{n_ord, l} = f*D'/h_x^n_ord;
                    nod = 'centralized';
                    P = 2*l;
                    D = derivmtx(l, P, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.centralized{n_ord, l} = f*D'/h_x^n_ord;
                    %
                    error_frame_staggered = partial_f{n_ord}.staggered - partial_f_approximation.staggered{n_ord, l};
                    error_weight = partial_f{n_ord}.staggered;
                    err.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(boundary_mask))/norm(error_weight(boundary_mask));
                    err.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(~boundary_mask))/norm(error_weight(~boundary_mask));
                    %
                    error_frame_centralized = partial_f{n_ord}.centralized - partial_f_approximation.centralized{n_ord, l};
                    error_weight = partial_f{n_ord}.centralized;
                    err.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(boundary_mask))/norm(error_weight(boundary_mask));
                    err.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(~boundary_mask))/norm(error_weight(~boundary_mask));
                    
                    %%  Fornberg
                    nod = 'staggered';
                    D = derivmtx_Fornberg(l, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.staggered{n_ord, l} = f*D'/h_x^n_ord;
                    nod = 'centralized';
                    D = derivmtx_Fornberg(l, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.centralized{n_ord, l} = f*D'/h_x^n_ord;
                    %
                    error_frame_staggered = partial_f{n_ord}.staggered - partial_f_approximation.staggered{n_ord, l};
                    error_weight = partial_f{n_ord}.staggered;
                    err_Fornberg.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(boundary_mask))/norm(error_weight(boundary_mask));
                    err_Fornberg.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(~boundary_mask))/norm(error_weight(~boundary_mask));
                    %
                    error_frame_centralized = partial_f{n_ord}.centralized - partial_f_approximation.centralized{n_ord, l};
                    error_weight = partial_f{n_ord}.centralized;
                    err_Fornberg.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(boundary_mask))/norm(error_weight(boundary_mask));
                    err_Fornberg.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(~boundary_mask))/norm(error_weight(~boundary_mask));
                    
                    %%  Savtizky-Golay
                    nod = 'staggered';
                    P = 2*l-1;
                    D = derivmtx_SavitzkyGolay(l, P, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.staggered{n_ord, l} = f*D'/h_x^n_ord;
                    nod = 'centralized';
                    P = 2*l;
                    D = derivmtx_SavitzkyGolay(l, P, n_ord, N_x, nod, true, flag_symbolic);
                    partial_f_approximation.centralized{n_ord, l} = f*D'/h_x^n_ord;
                    %
                    error_frame_staggered = partial_f{n_ord}.staggered - partial_f_approximation.staggered{n_ord, l};
                    error_weight = partial_f{n_ord}.staggered;
                    err_SavitzyGolay.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(boundary_mask))/norm(error_weight(boundary_mask));
                    err_SavitzyGolay.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_staggered(~boundary_mask))/norm(error_weight(~boundary_mask));
                    %
                    error_frame_centralized = partial_f{n_ord}.centralized - partial_f_approximation.centralized{n_ord, l};
                    error_weight = partial_f{n_ord}.centralized;
                    err_SavitzyGolay.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(boundary_mask))/norm(error_weight(boundary_mask));
                    err_SavitzyGolay.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = norm(error_frame_centralized(~boundary_mask))/norm(error_weight(~boundary_mask));
                    
                else
                    err.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    %
                    err_Fornberg.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_Fornberg.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_Fornberg.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_Fornberg.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    %
                    err_SavitzyGolay.staggered.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_SavitzyGolay.staggered.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_SavitzyGolay.centralized.boundary(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                    err_SavitzyGolay.centralized.interior(iteration_frequency, l, n_ord, iteration_delta_x) = nan;
                end
                min_val = min([...
                    err.centralized.interior(:);...
                    err.centralized.boundary(:);...
                    err.staggered.interior(:);...
                    err.staggered.boundary(:); ...
                    err_Fornberg.centralized.interior(:);...
                    err_Fornberg.centralized.boundary(:);...
                    err_Fornberg.staggered.interior(:);...
                    err_Fornberg.staggered.boundary(:);...
                    err_SavitzyGolay.centralized.interior(:);...
                    err_SavitzyGolay.centralized.boundary(:);...
                    err_SavitzyGolay.staggered.interior(:);...
                    err_SavitzyGolay.staggered.boundary(:)]);
                max_val = max([...
                    err.centralized.interior(:);...
                    err.centralized.boundary(:);...
                    err.staggered.interior(:);...
                    err.staggered.boundary(:); ...
                    err_Fornberg.centralized.interior(:);...
                    err_Fornberg.centralized.boundary(:);...
                    err_Fornberg.staggered.interior(:);...
                    err_Fornberg.staggered.boundary(:);...
                    err_SavitzyGolay.centralized.interior(:);...
                    err_SavitzyGolay.centralized.boundary(:);...
                    err_SavitzyGolay.staggered.interior(:);...
                    err_SavitzyGolay.staggered.boundary(:)]);
                if min_val>0
                    min_rng = min(min_rng, min_val);
                    max_rng = max(max_rng, max_val);
                end
            end
        end
    end
end

%%
err.staggered.boundary = mean(err.staggered.boundary, 4);
err.staggered.interior = mean(err.staggered.interior, 4);
err.centralized.boundary = mean(err.centralized.boundary, 4);
err.centralized.interior = mean(err.centralized.interior, 4);
%
err_Fornberg.staggered.boundary = mean(err_Fornberg.staggered.boundary, 4);
err_Fornberg.staggered.interior = mean(err_Fornberg.staggered.interior, 4);
err_Fornberg.centralized.boundary = mean(err_Fornberg.centralized.boundary, 4);
err_Fornberg.centralized.interior = mean(err_Fornberg.centralized.interior, 4);
%
err_SavitzyGolay.staggered.boundary = mean(err_SavitzyGolay.staggered.boundary, 4);
err_SavitzyGolay.staggered.interior = mean(err_SavitzyGolay.staggered.interior, 4);
err_SavitzyGolay.centralized.boundary = mean(err_SavitzyGolay.centralized.boundary, 4);
err_SavitzyGolay.centralized.interior = mean(err_SavitzyGolay.centralized.interior, 4);

%%  Color specification
col = colormap(hot);
col = col(10:end-10, :);
stp = length(col)/(3+1);
col = col(round(1: stp: length(col)), :);
close all

font_size = 8;
font_size_contour = 10;

figure('rend','painters','pos',[50, 50, 350, 575]);
contour_levels = 10.^[-1, -2, -3, -4, -5, -7, -9, -11 -13];
n_ord = 1;
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err.staggered.interior(:, 2:end,n_ord), contour_levels, '-', ...
    'Color', [.6 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
hold on
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_Fornberg.staggered.interior(:, 2:end,n_ord), contour_levels, '--', ...
    'Color', col(3, :));
clabel(c,h,'FontSize',font_size_contour)
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_SavitzyGolay.staggered.interior(:, 2:end,n_ord), contour_levels, '-.', ...
    'Color', [0 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
legend('MaxPol', 'Fornberg', 'Savitzky-Golay')
xlabel('Tap-length polynomial (l)')
ylabel('Zone-Plate harmonic frequency (\omega)')
set(gca, 'FontSize', font_size_contour)
title('Staggered first order derivative interior approximation')

figure('rend','painters','pos',[450, 50, 350, 575]);
contour_levels = 10.^[-1, -2, -3, -4, -6, -8, 10];
n_ord = 1;
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err.staggered.boundary(:, 2:end,n_ord), contour_levels, '-', ...
    'Color', [.6 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
hold on
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_Fornberg.staggered.boundary(:, 2:end,n_ord), contour_levels, '--', ...
    'Color', col(3, :));
clabel(c,h,'FontSize',font_size_contour)
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_SavitzyGolay.staggered.boundary(:, 2:end,n_ord), contour_levels, '-.', ...
    'Color', [0 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
legend('MaxPol', 'Fornberg', 'Savitzky-Golay')
xlabel('Tap-length polynomial (l)')
ylabel('Zone-Plate harmonic frequency (\omega)')
set(gca, 'FontSize', font_size)
title('Staggered first order derivative boundary approximation')

%
figure('rend','painters','pos',[850, 50, 350, 575]);
contour_levels = 10.^[-1, -2, -3, -4, -6, -8, -10, -12];
n_ord = 2;
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err.centralized.interior(:, 2:end,n_ord), contour_levels, '-', ...
    'Color', [.6 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
hold on
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_Fornberg.centralized.interior(:, 2:end,n_ord), contour_levels, '--', ...
    'Color', col(3, :));
clabel(c,h,'FontSize',font_size_contour)
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_SavitzyGolay.centralized.interior(:, 2:end,n_ord), contour_levels, '-.', ...
    'Color', [0 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
legend('MaxPol', 'Fornberg', 'Savitzky-Golay')
xlabel('Tap-length polynomial (l)')
ylabel('Zone-Plate harmonic frequency (\omega)')
set(gca, 'FontSize', font_size)
title('Centralized first order derivative interior approximation')


figure('rend','painters','pos',[1250, 50, 350, 575]);
contour_levels = 10.^[-1, -2, -3, -4, -6, -8, -10];
n_ord = 2;
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err.centralized.boundary(:, 2:end,n_ord), contour_levels, '-', ...
    'Color', [.6 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
hold on
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_Fornberg.centralized.boundary(:, 2:end,n_ord), contour_levels, '--', ...
    'Color', col(3, :));
clabel(c,h,'FontSize',font_size_contour)
[c,h] = contour(possible_l(2:end), frequency_range*2*pi, ...
    err_SavitzyGolay.centralized.boundary(:, 2:end,n_ord), contour_levels, '-.', ...
    'Color', [0 .6 0]);
clabel(c,h,'FontSize',font_size_contour)
legend('MaxPol', 'Fornberg', 'Savitzky-Golay')
xlabel('Tap-length polynomial (l)')
ylabel('Zone-Plate harmonic frequency (\omega)')
set(gca, 'FontSize', font_size)
title('Centralized first order derivative boundary approximation')

warning on