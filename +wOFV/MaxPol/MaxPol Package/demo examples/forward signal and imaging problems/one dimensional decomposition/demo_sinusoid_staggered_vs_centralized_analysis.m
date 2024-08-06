clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%
x_start = -1;
x_end = 1;
N = 128;    % number of discrete sampling
h_T = (x_end - x_start)/(N - 1);
x = [x_start: h_T: x_end]';
L = [1: 10];    % polynomial degree
F = [.1: .1: 10];   % sinusoid frequency

iteration_frequency = 0;
for frequency = F;    % Harmonic frequency
    iteration_frequency = iteration_frequency + 1;
    omega = 2*pi*frequency;
    
    %%  creat analytical vector-valued functions
    y = sin(omega*x);
    dy_centralized = omega * cos(omega*x);
    dy_staggered   = omega * cos(omega*(x - h_T/2));
    
    %%  approximate numerical derivative
    iteration_l = 0
    for l = L;        
        iteration_l = iteration_l + 1;
        [iteration_frequency, iteration_l]
        %   staggered approximation
        P = 2*l - 1;
        n_ord = 1;
        nod = 'staggered';
        sparse_flag = true;
        sym_flag = false;
        [D] = derivmtx(l, P, n_ord, N, nod, sparse_flag, sym_flag);
        dy_approximate = D*y/h_T^n_ord;
        err_staggered = abs(dy_staggered - dy_approximate);
        err.staggered.boundary(iteration_frequency, iteration_l) = mean(err_staggered([1:l, end-l+1: end]));
        err.staggered.interior(iteration_frequency, iteration_l) = mean(err_staggered(l+1: end-l));
        %   centralized approximation
        P = 2*l;
        n_ord = 1;
        nod = 'centralized';
        sparse_flag = true;
        sym_flag = false;
        [D] = derivmtx(l, P, n_ord, N, nod, sparse_flag, sym_flag);
        dy_approximate = D*y/h_T^n_ord;
        err_centralized = abs(dy_centralized - dy_approximate);
        err.centralized.boundary(iteration_frequency, iteration_l) = mean(err_centralized([1:l, end-l+1: end]));
        err.centralized.interior(iteration_frequency, iteration_l) = mean(err_centralized(l+1: end-l));
    end
end

%%
figure('rend','painters','pos',[100, 200, 300, 300]);
[c,h] = imcontour(L, F, log10(err.staggered.interior), [-1:-1:-12]);
clabel(c,h)
title('log_{10}|approx. Error|, interior, staggered')
xlabel('polynomial degree (l)')
ylabel('frequency (\omega)')
%
figure('rend','painters','pos',[425, 200, 300, 300]);
[c,h] = imcontour(L, F, log10(err.centralized.interior), [-1:-1:-12]);
clabel(c,h)
title('log_{10}|Error|, interior, centralized')
xlabel('polynomial degree (l)')
ylabel('frequency (\omega)')
%
figure('rend','painters','pos',[750, 200, 300, 300]);
[c,h] = imcontour(L, F, log10(err.staggered.boundary), [-1:-1:-9]);
clabel(c,h)
title('log_{10}|Error|, boundary, staggered')
xlabel('polynomial degree (l)')
ylabel('frequency (\omega)')
%
figure('rend','painters','pos',[1075, 200, 300, 300]);
[c,h] = imcontour(L, F, log10(err.centralized.boundary), [-1:-1:-9]);
clabel(c,h)
title('log_{10}|Error|, boundary, centralized')
xlabel('polynomial degree (l)')
ylabel('frequency (\omega)')