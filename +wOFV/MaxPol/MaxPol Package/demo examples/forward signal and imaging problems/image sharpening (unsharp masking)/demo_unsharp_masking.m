close all
clear all
clc

%%
current_directory = pwd;
cd ..
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'data'])
addpath([cd, filesep, 'third party codes'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
addpath([cd, filesep, 'third party codes', filesep, 'SI'])
cd(current_directory)

%%
params.lambda = 1;
params.cut = 13;
do_histogram_correction = false;

%%  load image library
image_name = 'flower';
reference_RGB = imread([image_name, '.png']);
if true
    %     imshow(reference_RGB)
    %     rect = getrect;
    rect = [164    32   479   501];
    reference_RGB = reference_RGB(rect(2): rect(2) + rect(4)-1, rect(1): rect(1) + rect(3)-1, :);
    close all
end

blurred_image_RGB = reference_RGB;
blurred_image_RGB_post = blurred_image_RGB(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
reference_RGB_post = reference_RGB(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);

%%  convert RGB to Lab and apply unsharp masking only on Lightness channel
[m, n, q] = size(reference_RGB_post);
reference_Lab = rgb2lab(reference_RGB);
blurred_image_Lab = rgb2lab(blurred_image_RGB);
reference_image = reference_Lab(:, :, 1);
blurred_image = blurred_image_Lab(:, :, 1);
params.ref_image = reference_image;

%%  parameter setup
%
params.l = 6;
params.sig = 0.67;
params.degree_step = pi/32;
params.method = 'gaussian';
[out_gaussian] = directional_sharpening(blurred_image, params);
params.l = 2;
params.degree_step = pi/32;
params.method = 'Farid_Simoncelli';
[out_Farid] = directional_sharpening(blurred_image, params);
params.l = 2;
params.px = 10;
params.py = 10;
params.degree_step = pi/24;
params.method = 'Savitzky_Golay';
[out_savtizky_golay] = directional_sharpening(blurred_image, params);
params.l = 6;
params.px = 10;
params.py = 10;
params.degree_step = pi/22;
params.method = 'MaxPol';
[out_maxpol] = directional_sharpening(blurred_image, params);

params.l = 6;
params.px = 12;
params.py = 12;
params.degree_step = pi/12;
params.method = 'MaxPol';
[out_maxpol_fullband] = directional_sharpening(blurred_image, params);

%%  crop the boundary edges to elimiante boundary condition effect
MaxPol_fullband_sharp_image = out_maxpol_fullband.sharpened_image(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
MaxPol_sharp_image = out_maxpol.sharpened_image(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Savtizky_Golay_sharp_image = out_savtizky_golay.sharpened_image(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Gaussian_sharp_image = out_gaussian.sharpened_image(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Farid_Simoncelli_sharp_image = out_Farid.sharpened_image(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);

%% conver lab to rgb for sharpened image
Gaussian_sharpened_image_Lab = blurred_image_Lab(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Gaussian_sharpened_image_Lab(:,:,1) = Gaussian_sharp_image;
Gaussian_sharpened_image_RGB = lab2rgb(Gaussian_sharpened_image_Lab, 'OutputType', 'uint8');
%
Farid_Simoncelli_sharpened_image_Lab = blurred_image_Lab(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Farid_Simoncelli_sharpened_image_Lab(:,:,1) = Farid_Simoncelli_sharp_image;
Farid_Simoncelli_sharpened_image_RGB = lab2rgb(Farid_Simoncelli_sharpened_image_Lab, 'OutputType', 'uint8');
%
Savtizky_Golay_sharpened_image_Lab = blurred_image_Lab(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
Savtizky_Golay_sharpened_image_Lab(:,:,1) = Savtizky_Golay_sharp_image;
Savtizky_Golay_sharpened_image_RGB = lab2rgb(Savtizky_Golay_sharpened_image_Lab, 'OutputType', 'uint8');
%
MaxPol_sharpened_image_Lab = blurred_image_Lab(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
MaxPol_sharpened_image_Lab(:,:,1) = MaxPol_sharp_image;
Maxflat_sharpened_image_RGB = lab2rgb(MaxPol_sharpened_image_Lab, 'OutputType', 'uint8');
%
MaxPol_fullband_sharpened_image_Lab = blurred_image_Lab(params.cut+1: end-params.cut, params.cut+1: end-params.cut, :);
MaxPol_fullband_sharpened_image_Lab(:,:,1) = MaxPol_fullband_sharp_image;
Maxflat_fullband_sharpened_image_RGB = lab2rgb(MaxPol_fullband_sharpened_image_Lab, 'OutputType', 'uint8');


%%  Sharpnes quality index registration
%
si_index_Bluured_image = si_rgb(blurred_image_RGB_post);
si_index_Gaussian = si_rgb(Gaussian_sharpened_image_RGB);
si_index_Farid_Simoncelli = si_rgb(Farid_Simoncelli_sharpened_image_RGB);
si_index_Savtizky_Golay = si_rgb(Savtizky_Golay_sharpened_image_RGB);
si_index_MaxPol = si_rgb(Maxflat_sharpened_image_RGB);
si_index_MaxPol_fullband = si_rgb(Maxflat_fullband_sharpened_image_RGB);

%%
[m, n, q] = size(reference_RGB_post);
scale = 1;
%
figure('rend','painters','pos',[50, m*scale+50, [n, m]*scale]);
img(double(blurred_image_RGB_post)/255)
title(['Original Image (Blurry Observation), SI = ', num2str(round(si_index_Bluured_image))])

figure('rend','painters','pos',[50+1*(n*scale+25), m*scale+50, [n, m]*scale]);
img(double(Gaussian_sharpened_image_RGB)/255)
title(['Gaussian Unsharp Masking, SI = ', num2str(round(si_index_Gaussian))])

figure('rend','painters','pos',[50+2*(n*scale+25), m*scale+50, [n, m]*scale]);
img(double(Farid_Simoncelli_sharpened_image_RGB)/255)
title(['Farid-Simoncelli Unsharp Masking, SI = ', num2str(round(si_index_Farid_Simoncelli))])

figure('rend','painters','pos',[50, 50, [n, m]*scale]);
img(double(Savtizky_Golay_sharpened_image_RGB)/255)
title(['Savtizky-Golay Unsharp Masking, SI = ', num2str(round(si_index_Savtizky_Golay))])

figure('rend','painters','pos',[50+1*(n*scale+25), 50, [n, m]*scale]);
img(double(Maxflat_sharpened_image_RGB)/255)
title(['MaxPol Unsharp Masking, SI = ', num2str(round(si_index_MaxPol))])

figure('rend','painters','pos',[50+2*(n*scale+25), 50, [n, m]*scale]);
img(double(Maxflat_fullband_sharpened_image_RGB)/255)
title(['MaxPol (fullband) Unsharp Masking, SI = ', num2str(round(si_index_MaxPol_fullband))])