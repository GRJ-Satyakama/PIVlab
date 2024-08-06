clear all
close all
clc

%%  load image library
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'data'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
addpath([cd, filesep, 'third party codes'])
cd(current_directory)
load Farid_Simoncelli_TIP2004.mat

%%  parameter for frequency response (fft) calculation
n_fft = 512;
stp = 2*pi/(n_fft);
omega = [-pi: stp: pi-stp]';
omega_half = omega(n_fft/2+2:end);

%%  setup parameters
n_ord = input('Please enter the order of differentiation (0, 1, 2, ...): ');
l = 15;
sym_flag = true;

%%  filter measurement setup
prc_w_a = .9;
prc_w_c = 1/2;
prc_w_b = 0.1;

%%  Gaussian filter analaysis
% sigma_values = exp([1:-.05:-.5]);
sigma_values = exp([.7:-.1:-.5]);
iteration_Gaussian = 0;
for sigma_value = sigma_values
    iteration_Gaussian = iteration_Gaussian + 1;
    [c_Gaussian{iteration_Gaussian}] = gaussian_kernel(l, sigma_value, n_ord);
    c_Gaussian{iteration_Gaussian} = c_Gaussian{iteration_Gaussian}(:);
    h_response_Gaussian{iteration_Gaussian} = fftshift(abs(fft(c_Gaussian{iteration_Gaussian}, n_fft)));
    h_quality_Gaussian{iteration_Gaussian} = h_response_Gaussian{iteration_Gaussian}./abs(omega).^n_ord;
    %
    h_response_Gaussian{iteration_Gaussian} = h_response_Gaussian{iteration_Gaussian}(n_fft/2+2:end);
    h_quality_Gaussian{iteration_Gaussian} = h_quality_Gaussian{iteration_Gaussian}(n_fft/2+2:end);
    %
    indx_a = sum(h_quality_Gaussian{iteration_Gaussian} > prc_w_a);
    indx_c = sum(h_quality_Gaussian{iteration_Gaussian} > prc_w_c);
    indx_b = sum(h_quality_Gaussian{iteration_Gaussian} > prc_w_b);
    %
    w_a_Gaussian(iteration_Gaussian) = omega_half(indx_a);
    w_b_Gaussian(iteration_Gaussian) = omega_half(indx_b);
    w_c_Gaussian(iteration_Gaussian) = omega_half(indx_c);
    %
    %%  passband period
    S_passband_Gaussian(iteration_Gaussian) = h_response_Gaussian{iteration_Gaussian}(1: indx_c)' * omega_half(1: indx_c) / indx_c / 2;
    
    %%  Rolloff transition period
    S_transition_Gaussian(iteration_Gaussian) = h_response_Gaussian{iteration_Gaussian}(indx_a: indx_b)' * omega_half(indx_a: indx_b) / (indx_b - indx_a + 1) / 2;
    
    %%  Stopband period
    S_stopband_Gaussian(iteration_Gaussian) = h_response_Gaussian{iteration_Gaussian}(indx_b: end)' * omega_half(indx_b: end) / (numel(omega_half) - indx_b + 1) / 2;
    
end

%%  Savtizky-Golay fitler analysis
P_values = [2: 2: 2*l-2];
iteration_SavitzkyGolay = 0;
for P = P_values
    iteration_SavitzkyGolay = iteration_SavitzkyGolay + 1
    c_SavitzkyGolay{iteration_SavitzkyGolay} = derivstag_SavitzkyGolay(l, P, 0, n_ord, sym_flag);
    h_response_SavitzkyGolay{iteration_SavitzkyGolay} = fftshift(abs(fft(c_SavitzkyGolay{iteration_SavitzkyGolay}', n_fft)));
    h_quality_SavitzkyGolay{iteration_SavitzkyGolay} = h_response_SavitzkyGolay{iteration_SavitzkyGolay}./abs(omega).^n_ord;
    %
    h_response_SavitzkyGolay{iteration_SavitzkyGolay} = h_response_SavitzkyGolay{iteration_SavitzkyGolay}(n_fft/2+2:end);
    h_quality_SavitzkyGolay{iteration_SavitzkyGolay} = h_quality_SavitzkyGolay{iteration_SavitzkyGolay}(n_fft/2+2:end);
    %
    indx_a = sum(h_quality_SavitzkyGolay{iteration_SavitzkyGolay} > prc_w_a);
    indx_c = sum(h_quality_SavitzkyGolay{iteration_SavitzkyGolay} > prc_w_c);
    indx_b = sum(h_quality_SavitzkyGolay{iteration_SavitzkyGolay} > prc_w_b);
    %
    w_a_SavitzkyGolay(iteration_SavitzkyGolay) = omega_half(indx_a);
    w_b_SavitzkyGolay(iteration_SavitzkyGolay) = omega_half(indx_b);
    w_c_SavitzkyGolay(iteration_SavitzkyGolay) = omega_half(indx_c);
    %
    %%  passband period
    S_passband_SavitzkyGolay(iteration_SavitzkyGolay) = h_response_SavitzkyGolay{iteration_SavitzkyGolay}(1: indx_c)' * omega_half(1: indx_c) / indx_c / 2;
    
    %%  Rolloff transition period
    S_transition_SavitzkyGolay(iteration_SavitzkyGolay) = h_response_SavitzkyGolay{iteration_SavitzkyGolay}(indx_a: indx_b)' * omega_half(indx_a: indx_b) / (indx_b - indx_a + 1) / 2;
    
    %%  Stopband period
    S_stopband_SavitzkyGolay(iteration_SavitzkyGolay) = h_response_SavitzkyGolay{iteration_SavitzkyGolay}(indx_b: end)' * omega_half(indx_b: end) / (numel(omega_half) - indx_b + 1) / 2;
end

%%  MaxPol fitler analysis
iteration_MaxPol = 0;
for P = P_values
    iteration_MaxPol = iteration_MaxPol + 1
    c_MaxPol{iteration_MaxPol} = derivstag(l, P, 0, n_ord, sym_flag);
    h_response_MaxPol{iteration_MaxPol} = fftshift(abs(fft(c_MaxPol{iteration_MaxPol}', n_fft)));
    h_quality_MaxPol{iteration_MaxPol} = h_response_MaxPol{iteration_MaxPol}./abs(omega).^n_ord;
    %
    h_response_MaxPol{iteration_MaxPol} = h_response_MaxPol{iteration_MaxPol}(n_fft/2+2:end);
    h_quality_MaxPol{iteration_MaxPol} = h_quality_MaxPol{iteration_MaxPol}(n_fft/2+2:end);
    %
    indx_a = sum(h_quality_MaxPol{iteration_MaxPol} > prc_w_a);
    indx_c = sum(h_quality_MaxPol{iteration_MaxPol} > prc_w_c);
    indx_b = sum(h_quality_MaxPol{iteration_MaxPol} > prc_w_b);
    %
    w_a_MaxPol(iteration_MaxPol) = omega_half(indx_a);
    w_b_MaxPol(iteration_MaxPol) = omega_half(indx_b);
    w_c_MaxPol(iteration_MaxPol) = omega_half(indx_c);
    %
    %%  passband period
    S_passband_MaxPol(iteration_MaxPol) = h_response_MaxPol{iteration_MaxPol}(1: indx_c)' * omega_half(1: indx_c) / indx_c / 2;
    
    %%  Rolloff transition period
    S_transition_MaxPol(iteration_MaxPol) = h_response_MaxPol{iteration_MaxPol}(indx_a: indx_b)' * omega_half(indx_a: indx_b) / (indx_b - indx_a + 1) / 2;
    
    %%  Stopband period
    S_stopband_MaxPol(iteration_MaxPol) = h_response_MaxPol{iteration_MaxPol}(indx_b: end)' * omega_half(indx_b: end) / (numel(omega_half) - indx_b + 1) / 2;
end

%%  Farid-Simoncelli fitler analysis
l_values = [2:9];
iteration_Simoncelli = 0;
for l_value = l_values
    iteration_Simoncelli = iteration_Simoncelli + 1
    c_Simoncelli{iteration_Simoncelli} = dlowpass{l_value-1, n_ord+1};
    h_response_Simoncelli{iteration_Simoncelli} = fftshift(abs(fft(c_Simoncelli{iteration_Simoncelli}, n_fft)));
    h_quality_Simoncelli{iteration_Simoncelli} = h_response_Simoncelli{iteration_Simoncelli}./abs(omega).^n_ord;
    %
    h_response_Simoncelli{iteration_Simoncelli} = h_response_Simoncelli{iteration_Simoncelli}(n_fft/2+2:end);
    h_quality_Simoncelli{iteration_Simoncelli} = h_quality_Simoncelli{iteration_Simoncelli}(n_fft/2+2:end);
    %
    indx_a = sum(h_quality_Simoncelli{iteration_Simoncelli} > prc_w_a);
    indx_c = sum(h_quality_Simoncelli{iteration_Simoncelli} > prc_w_c);
    indx_b = sum(h_quality_Simoncelli{iteration_Simoncelli} > prc_w_b);
    %
    w_a_Simoncelli(iteration_Simoncelli) = omega_half(indx_a);
    w_b_Simoncelli(iteration_Simoncelli) = omega_half(indx_b);
    w_c_Simoncelli(iteration_Simoncelli) = omega_half(indx_c);
    %
    %%  passband period
    S_passband_Simoncelli(iteration_Simoncelli) = h_response_Simoncelli{iteration_Simoncelli}(1: indx_c)' * omega_half(1: indx_c) / indx_c / 2;
    
    %%  Rolloff transition period
    S_transition_Simoncelli(iteration_Simoncelli) = h_response_Simoncelli{iteration_Simoncelli}(indx_a: indx_b)' * omega_half(indx_a: indx_b) / (indx_b - indx_a + 1) / 2;
    
    %%  Stopband period
    S_stopband_Simoncelli(iteration_Simoncelli) = h_response_Simoncelli{iteration_Simoncelli}(indx_b: end)' * omega_half(indx_b: end) / (numel(omega_half) - indx_b + 1) / 2;
end

%%
col = [.6 0 0;
    0 0 .6;
    0 .6 0;
    .6 .6 0];

%%
eng = [S_passband_Gaussian, S_transition_Gaussian, S_stopband_Gaussian, ...
    S_passband_Simoncelli, S_transition_Simoncelli, S_stopband_Simoncelli, ...
    S_passband_SavitzkyGolay, S_transition_SavitzkyGolay, S_stopband_SavitzkyGolay, ...
    S_passband_MaxPol, S_transition_MaxPol, S_stopband_MaxPol];
w_c = [w_c_Gaussian, w_c_Simoncelli, w_c_SavitzkyGolay, w_c_MaxPol];

%%  excute Gaussian results
figure('rend','painters','pos', [50, 550, 400, 400]);
plot(w_c_Gaussian, S_passband_Gaussian, 'LineWidth', 1, 'Color', col(1, :))
hold on
plot(w_c_Gaussian, S_stopband_Gaussian, '--', 'LineWidth', 1, 'Color', col(1, :))
% axis([0 pi min(eng) max(eng)])
axis([0 pi 6e-4 5e-0])
set(gca, 'FontSize', 14, 'YScale', 'log')
grid
xlabel('Cutoff Frequency - (\omega)')
ylabel('Spectral Power')
legend('Pass-band', 'Stop-band', 'Location', 'SouthEast')
%
figure('rend','painters','pos', [50, 50, 400, 400]);
plot(sigma_values, w_c_Gaussian, 'LineWidth', 1, 'Color', col(1, :))
hold on
plot(sigma_values, w_b_Gaussian - w_a_Gaussian, '--', 'LineWidth', 1, 'Color', col(1, :))
axis([min(sigma_values) max(sigma_values) 0 pi])
set(gca, 'FontSize', 14, 'YScale', 'linear')
grid
xlabel('Variance (\sigma)')
ylabel('Frequency - (\omega)')
legend('Cutoff frequency', 'Transition band', 'Location', 'NorthEast')

%%  excute Farid-Simoncelli results
figure('rend','painters','pos', [500, 550, 400, 400]);
plot(w_c_Simoncelli, S_passband_Simoncelli, 'LineWidth', 1, 'Color', col(2, :))
hold on
plot(w_c_Simoncelli, S_stopband_Simoncelli, '--', 'LineWidth', 1, 'Color', col(2, :))
% axis([0 pi min(eng) max(eng)])
axis([0 pi 6e-4 5e-0])
set(gca, 'FontSize', 14, 'YScale', 'log')
grid
xlabel('Cutoff Frequency - (\omega)')
ylabel('Spectral Power')
legend('Pass-band', 'Stop-band', 'Location', 'SouthEast')
%
figure('rend','painters','pos', [500, 50, 400, 400]);
plot(l_values, w_c_Simoncelli, 'LineWidth', 1, 'Color', col(2, :))
hold on
plot(l_values, w_b_Simoncelli - w_a_Simoncelli, '--', 'LineWidth', 1, 'Color', col(2, :))
axis([min(l_values) max(l_values) 0 pi])
set(gca, 'FontSize', 14, 'YScale', 'linear')
grid
xlabel('Tap-length polynomial (l)')
ylabel('Frequency - (\omega)')
legend('Cutoff frequency', 'Transition band', 'Location', 'NorthEast')

%% excute Savitzky-Golay results
figure('rend','painters','pos', [950, 550, 400, 400]);
plot(w_c_SavitzkyGolay, S_passband_SavitzkyGolay, 'LineWidth', 1, 'Color', col(3, :))
hold on
plot(w_c_SavitzkyGolay, S_stopband_SavitzkyGolay, '--', 'LineWidth', 1, 'Color', col(3, :))
% axis([0 pi min(eng) max(eng)])
axis([0 pi 6e-4 5e-0])
set(gca, 'FontSize', 14, 'YScale', 'log')
grid
xlabel('Cutoff Frequency - (\omega)')
ylabel('Spectral Power')
legend('Pass-band', 'Stop-band', 'Location', 'SouthEast')
%
figure('rend','painters','pos', [950, 50, 400, 400]);
plot(P_values, w_c_SavitzkyGolay, 'LineWidth', 1, 'Color', col(3, :))
hold on
plot(P_values, w_b_SavitzkyGolay - w_a_SavitzkyGolay, '--', 'LineWidth', 1, 'Color', col(3, :))
axis([min(P_values) max(P_values) 0 pi])
set(gca, 'FontSize', 14, 'YScale', 'linear')
grid
xlabel('Least-Square polynomial degree (P)')
ylabel('Frequency - (\omega)')
legend('Cutoff frequency', 'Transition band', 'Location', 'NorthEast')

%%  excute MaxPol results
figure('rend','painters','pos', [1400, 550, 400, 400]);
plot(w_c_MaxPol, S_passband_MaxPol, 'LineWidth', 1, 'Color', col(4, :))
hold on
plot(w_c_MaxPol, S_stopband_MaxPol, '--', 'LineWidth', 1, 'Color', col(4, :))
% axis([0 pi min(eng) max(eng)])
axis([0 pi 6e-4 5e-0])
set(gca, 'FontSize', 14, 'YScale', 'log')
grid
xlabel('Cutoff Frequency - (\omega)')
ylabel('Spectral Power')
legend('Pass-band', 'Stop-band', 'Location', 'SouthEast')
%
figure('rend','painters','pos', [1400, 50, 400, 400]);
plot(P_values, w_c_MaxPol, 'LineWidth', 1, 'Color', col(4, :))
hold on
plot(P_values, w_b_MaxPol - w_a_MaxPol, '--', 'LineWidth', 1, 'Color', col(4, :))
axis([min(P_values) max(P_values) 0 pi])
set(gca, 'FontSize', 14, 'YScale', 'linear')
grid
xlabel('Maxflat polynomial degree (P)')
ylabel('Frequency - (\omega)')
legend('Cutoff frequency', 'Transition band', 'Location', 'NorthEast')