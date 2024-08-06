clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
addpath([cd, filesep, 'third party codes'])
cd(current_directory)

%%  lowpass filter with cutoff level P_x and P_y
n_ord = input('Please enter the order of differentiation (0, 1, 2, ...): ');  % order of differentiations to be evaluated
l = 15; % polynomial degree of the filter
theta = pi/4;   % rotation degree
N_fft = 513;
sym_flag = true;

%%  parameter for frequency response (fft) calculation
stp = 2*pi/(N_fft);
omega = [-pi: stp: pi-stp];

%%
p_x = 6;   % maxpol degree (cutoff) in x-axis
p_y = 6;   % maxpol degree (cutoff) in y-axis
[G_maxpol] = derivdirec(l, p_x, p_y, n_ord, theta, sym_flag);
%
p_x = 12;   % maxpol degree (cutoff) in x-axis
p_y = 12;   % maxpol degree (cutoff) in y-axis
[G_savtizky_golay] = derivdirec_SavitzkyGolay(l, p_x, p_y, n_ord, theta, sym_flag);
%
sigma_x = 1.35;
sigma_y = 1.35;
[G_Gaussian] = derivdirec_Gaussian(l, sigma_x, sigma_y, n_ord, theta);
%
a = max([G_maxpol(:); G_savtizky_golay(:); G_Gaussian(:)]);
b = min([G_maxpol(:); G_savtizky_golay(:); G_Gaussian(:)]);
if a>abs(b)
    max_range = [-a, a];
else
    max_range = [b, -b];
end

%%  excute 2D fitler plot
figure('rend','painters','pos',[50 50 300 300]);
imagesc(G_maxpol)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Impulse Response (MaxPol)')
%
figure('rend','painters','pos',[50 400 300 300]);
imagesc(G_savtizky_golay)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Impulse Response (Savtizky-Golay)')
%
figure('rend','painters','pos',[50 750 300 300]);
imagesc(G_Gaussian)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Impulse Response (Gaussian)')

%%  excute 2D fitler plot (Fourier Response)
F_maxpol = fftshift(abs(fft2(G_maxpol, N_fft, N_fft)));
F_savtizky_golay = fftshift(abs(fft2(G_savtizky_golay, N_fft, N_fft)));
F_Gaussian = fftshift(abs(fft2(G_Gaussian, N_fft, N_fft)));
a = max([F_maxpol(:); F_savtizky_golay(:); F_Gaussian(:)]);
b = min([F_maxpol(:); F_savtizky_golay(:); F_Gaussian(:)]);
if a>abs(b)
    max_range = [0, a];
else
    max_range = [0, -b];
end

%%
figure('rend','painters','pos',[400 50 300 300]);
imagesc(F_maxpol)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Filter Response (MaxPol)')
%
figure('rend','painters','pos',[400 400 300 300]);
imagesc(F_savtizky_golay)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Filter Response (Savtizky-Golay)')
%
figure('rend','painters','pos',[400 750 300 300]);
imagesc(F_Gaussian)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('Filter Response (Gaussian)')
%

%%  contour level
figure('rend','painters','pos',[750 50 300 300]);
contour_levels = [1, .5, 10.^[-1 -3 -8 -15]];
[c,h] = imcontour(omega, omega, F_maxpol, contour_levels, 'LineColor', [.6 .6 0]);
clabel(c,h)
axis image
xlabel('\omega_x')
ylabel('\omega_y')
title('Contour Level (MaxPol)')
%
figure('rend','painters','pos',[750 400 300 300]);
contour_levels = [1, .5, 10.^[-1 -2 -8 -15]];
[c,h] = imcontour(omega, omega, F_savtizky_golay, contour_levels, 'LineColor', [0 .6 0]);
clabel(c,h)
axis image
xlabel('\omega_x')
ylabel('\omega_y')
title('Contour Level (Savitzky-Golay)')
%
figure('rend','painters','pos',[750 750 300 300]);
contour_levels = [1, .5, 10.^[-1 -3 -8 -15]];
[c,h] = imcontour(omega, omega, F_Gaussian, contour_levels, 'LineColor', [.6 0 0]);
clabel(c,h)
axis image
xlabel('\omega_x')
ylabel('\omega_y')
title('Contour Level (Gaussian)')
