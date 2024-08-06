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

%%  setup parameters
n_ord = input('Please enter the order of differentiation (0, 1, 2, ...): ');
l = 25;
sym_flag = true;
if mod(n_ord, 2)
    x_polynomial = [-l+1/2: l-1/2];
else
    x_polynomial = [-l: l];
end

%% Farid-Simoncelli
[FS_Kernel] = dlowpass{2-1, n_ord+1};
if mod(n_ord, 2)
    x_Simoncelli = [-2+1/2: 2-1/2];    
else
    x_Simoncelli = [-2: 2];
end

%% Gaussian kernel
sgm = 1.22;
[G_Kernel, x_Gaussian] = gaussian_kernel(l, sgm, n_ord);

%%  Savtizky-Golay
P_SG = 22;
if mod(n_ord, 2)
    [SG_Kernel] = -derivstag_SavitzkyGolay(l, P_SG, 0, n_ord, sym_flag);
else
    [SG_Kernel] = derivcent_SavitzkyGolay(l, P_SG, 0, n_ord, sym_flag);
end


%%  MaxPol (Hosseini-Plataniotis)
P_HP = 12;
if mod(n_ord, 2)
    [MaxPol_Kernel] = -derivstag(l, P_HP, 0, n_ord, sym_flag);
else
    [MaxPol_Kernel] = derivcent(l, P_HP, 0, n_ord, sym_flag);
end


%%  parameter for frequency response (fft) calculation
n_fft = 512;
stp = 2*pi/(n_fft);
omega = [-pi: stp: pi-stp];

%%
h_G = fftshift(abs(fft(G_Kernel, n_fft)));
h_SG = fftshift(abs(fft(SG_Kernel, n_fft)));
h_FS = fftshift(abs(fft(FS_Kernel, n_fft)));
h_HP = fftshift(abs(fft(MaxPol_Kernel, n_fft)));
legend_string = {'Gaussian', 'Savitzky-Golay', 'Farid-Simoncelli', 'MaxFlat (Proposed)', ['|j\omega|^', num2str(n_ord)]};

%%
figure('rend','painters','pos', [50, 200, 500, 500]);
plot(omega, h_G, 'LineWidth', 1, 'Color', [.6 0 0])
hold on
plot(omega, h_SG, 'LineWidth', 1, 'Color', [0 .6 0])
plot(omega, h_FS, 'LineWidth', 1, 'Color', [0 0 .6])
plot(omega, h_HP, 'LineWidth', 1, 'Color', [.6 .6 0])
plot(omega, abs(omega).^n_ord, '--', 'LineWidth', 1, 'Color', [0 0 0])
range_plot = pi;
axis([-range_plot, range_plot, 0, (range_plot)^n_ord*1.05])
axis square
ylabel(['|H^' num2str(n_ord), '(\omega)|'])
xlabel('\omega')
set(gca, 'FontSize', 12)
legend(legend_string, 'Location', 'North')
%
figure('rend','painters','pos', [600, 200, 500, 500]);
plot(omega, h_G, 'LineWidth', 1, 'Color', [.6 0 0])
hold on
plot(omega, h_SG, 'LineWidth', 1, 'Color', [0 .6 0])
plot(omega, h_FS, 'LineWidth', 1, 'Color', [0 0 .6])
plot(omega, h_HP, 'LineWidth', 1, 'Color', [.6 .6 0])
plot(omega, abs(omega).^n_ord, '--', 'LineWidth', 1, 'Color', [0 0 0])
range_plot = pi;
axis([-range_plot, range_plot, 0, (range_plot)^n_ord*1.05])
axis square
ylabel(['|H^' num2str(n_ord), '(\omega)|'])
xlabel('\omega')
set(gca,'YScale','log', 'FontSize', 12)
legend(legend_string, 'Location', 'Best')

figure('rend','painters','pos', [1150, 200, 1000, 320]);
stem(x_Gaussian, G_Kernel, 'Color', [.6 0 0])
hold on
stem(x_polynomial, SG_Kernel, 'Color', [0 .6 0])
stem(x_Simoncelli, FS_Kernel, 'Color', [0 0 .6])
stem(x_polynomial, MaxPol_Kernel, 'Color', [.6 .6 0])
for iteration = 1: (2*l+1)
        label_kernel{iteration} = [num2str(iteration-l-1)];
end
set(gca,'Xtick',[-l:l],'XTickLabel',label_kernel)
axis tight
grid
xlabel('discrete node - (x)')
legend(legend_string(1:end-1), 'Location', 'NorthEast')

