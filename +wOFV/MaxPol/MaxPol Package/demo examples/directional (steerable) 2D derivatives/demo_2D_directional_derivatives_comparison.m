clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%  lowpass filter with cutoff level (p_h_lowpass)
n = [0:4];
l = 13;
p_x = 6;
p_y = 6;
deg_stp = pi/4;
theta = [0: deg_stp: pi-deg_stp];   % rotation degree
N_fft = 512;

for ord = n
    deg_iteration = 0;
    for deg = theta
        fprintf(['derivative order = ', num2str(ord), ...
            ', cutoff degree p_x = p_y = 6, steering angle = ', ...
            num2str(deg/pi) 'pi\n'])
        
        %%  call MaxPol directional derivative kernel
        deg_iteration = deg_iteration + 1;
        [G] = derivdirec(l, p_x, p_y, ord, deg, true);
        a = max(G(:));
        b = min(G(:));
        if a>abs(b)
            max_range = [-a, a];
        else
            max_range = [b, -b];
        end
        
        %%  excute 2D impulse response plots
        figure('rend','painters','pos',[(deg_iteration-1)*175+50 ...
            (ord)*175+50 150 150]);
        imagesc(G, max_range)
        axis image
        axis off
        colormap gray
        set(gca, 'Xtick', [], 'Ytick', [])
                
        %%  excute 2D fitler response plots
        F_G = fftshift(abs(fft2(G, N_fft, N_fft)));
        a = max(F_G(:));
        b = min(F_G(:));
        if a>abs(b)
            max_range = [0, a];
        else
            max_range = [0, -b];
        end
        figure('rend','painters','pos',[(deg_iteration-1)*175+900 ...
            (ord)*175+50 150 150]);
        imagesc(F_G, max_range)
        axis image
        axis off
        colormap gray
        set(gca, 'Xtick', [], 'Ytick', [])
                
        %%
        if ord == 0
            break
        end
    end
end
