clear all
close all
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
cd(current_directory)

%%
load Farid_Simoncelli_TIP2004.mat

%%  load image data
image_scan = imread('pirate.tif');
[m, n, q] = size(image_scan);
if q > 1
    image_scan = rgb2lab(image_scan);
    image_scan = image_scan(:,:,1);
end
image_scan = double(image_scan);

%% initialized parameters
l = 8;
%
sgm = sqrt(2);  % Gaussian
%
P_smoothing_Savitzky_Goaly = 2;
P_derivative_Savitzky_Goaly = 4;
%
P_smoothing_MaxPol = 0;
P_derivative_MaxPol = 2;

%%  Canny via Gaussian kernels (MATLAB recommendation)
[gaussKernel] = gaussian_derivatives(l, sgm);
smoothing_kernel = gaussKernel{1};
derivative_kernel = gaussKernel{2};
[segmented_Gaussian, dx_Gaussian, dy_Gaussian] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%  Canny via Simocelli
smoothing_kernel = dlowpass{l-1, 1};
derivative_kernel = -dlowpass{l-1, 2};
[segmented_Simoncelli, dx_Simoncelli, dy_Simoncelli] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%  Canny via Savitzky-Golay
[smoothing_kernel] = derivcent_SavitzkyGolay(l, P_smoothing_Savitzky_Goaly, 0, 0, true);
[derivative_kernel] = -derivstag_SavitzkyGolay(l, P_derivative_Savitzky_Goaly, 0, 1, true);
[segmented_Savitzky_Golay, dx_Savitzky_Golay, dy_Savitzky_Golay] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%  Canny via MaxPol
[smoothing_kernel] = derivcent(l, P_smoothing_MaxPol, 0, 0, true);
[derivative_kernel] = -derivstag(l, P_derivative_MaxPol, 0, 1, true);
[segmented_MaxPol, dx_MaxPol, dy_MaxPol] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%
dx_images = [dx_Gaussian(:); dx_Simoncelli(:); dx_Savitzky_Golay(:); dx_MaxPol(:)];
dy_images = [dy_Gaussian(:); dy_Simoncelli(:); dy_Savitzky_Golay(:); dy_MaxPol(:)];

%%  convert to 3D
image_scan_org = image_scan;
image_scan = repmat(image_scan, [1,1,3]);
segmented_Gaussian = repmat(segmented_Gaussian, [1,1,3]);
segmented_Simoncelli = repmat(segmented_Simoncelli, [1,1,3]);
segmented_Savitzky_Golay = repmat(segmented_Savitzky_Golay, [1,1,3]);
segmented_MaxPol = repmat(segmented_MaxPol, [1,1,3]);

%%
col = [.6 0 0;
    0 0 .6;
    0 .6 0;
    .6 .6 0];
vector_size = 1.3;

%%  color box patch
[m, n, q] = size(image_scan);
pos = [231, 292];
siz = [32, 32];
i_start = pos(1)-siz(1)/2+1;
i_end = pos(1)+siz(1)/2;
j_start = pos(2)-siz(2)/2+1;
j_end = pos(2)+siz(2)/2;
strip_color = [0, 1, 0];
for channel = 1: 3
    image_scan(i_start: i_end, j_start, channel) = 255*strip_color(channel);
    image_scan(i_start: i_end, j_end, channel) = 255*strip_color(channel);
    image_scan(i_start, j_start: j_end, channel) = 255*strip_color(channel);
    image_scan(i_end, j_start: j_end, channel) = 255*strip_color(channel);
    %
    segmented_Gaussian(i_start: i_end, j_start, channel) = 255*strip_color(channel);
    segmented_Gaussian(i_start: i_end, j_end, channel) = 255*strip_color(channel);
    segmented_Gaussian(i_start, j_start: j_end, channel) = 255*strip_color(channel);
    segmented_Gaussian(i_end, j_start: j_end, channel) = 255*strip_color(channel);
    %
    segmented_Simoncelli(i_start: i_end, j_start, channel) = 255*strip_color(channel);
    segmented_Simoncelli(i_start: i_end, j_end, channel) = 255*strip_color(channel);
    segmented_Simoncelli(i_start, j_start: j_end, channel) = 255*strip_color(channel);
    segmented_Simoncelli(i_end, j_start: j_end, channel) = 255*strip_color(channel);
    %
    segmented_Savitzky_Golay(i_start: i_end, j_start, channel) = 255*strip_color(channel);
    segmented_Savitzky_Golay(i_start: i_end, j_end, channel) = 255*strip_color(channel);
    segmented_Savitzky_Golay(i_start, j_start: j_end, channel) = 255*strip_color(channel);
    segmented_Savitzky_Golay(i_end, j_start: j_end, channel) = 255*strip_color(channel);
    %
    segmented_MaxPol(i_start: i_end, j_start, channel) = 255*strip_color(channel);
    segmented_MaxPol(i_start: i_end, j_end, channel) = 255*strip_color(channel);
    segmented_MaxPol(i_start, j_start: j_end, channel) = 255*strip_color(channel);
    segmented_MaxPol(i_end, j_start: j_end, channel) = 255*strip_color(channel);
    %
end

%%
close all
figure('rend','painters','pos',[50, 50, [n, m]]);
img(image_scan/255)
title('Origianl image')
%
figure('rend','painters','pos',[60, m+150, [n, m]*.45]);
imagesc(image_scan(i_start: i_end, j_start: j_end, :)/255, [0, 1])
set(gca, 'Xtick', [], 'Ytick', [])
axis image
axis off

%%  excute Gaussian results
% close all
figure('rend','painters','pos',[50+1*n, 50, [n, m]]);
img(segmented_Gaussian)
title('Gaussian kernel method')
%
figure('rend','painters','pos',[50+1.025*n, m+150, [n, m]*.45]);
imagesc(segmented_Gaussian(i_start: i_end, j_start: j_end, :), [0, 1])
set(gca, 'Xtick', [], 'Ytick', [])
axis image
axis off
%
figure('rend','painters','pos',[1.625*n, m+150, [n, m]*.45]);
img(image_scan_org(i_start: i_end,j_start: j_end)/255)
hold on
[X, Y] = meshgrid([1:n], [1:m]);
i_threshold = 0;
j_threshold = 0;
quiver(X(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-j_start+1,...
    Y(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-i_start+1,...
    dx_Gaussian(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    dy_Gaussian(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    vector_size,'Color', col(1,:))
axis off

%%  excute Simoncelli results
% close all
figure('rend','painters','pos',[50+2*n, 50, [n, m]]);
img(segmented_Simoncelli)
title('Farid-Simoncelli kernel method')

figure('rend','painters','pos',[50+2.025*n, m+150, [n, m]*.45]);
imagesc(segmented_Simoncelli(i_start: i_end, j_start: j_end, :), [0, 1])
set(gca, 'Xtick', [], 'Ytick', [])
axis image
axis off
%
figure('rend','painters','pos',[2.625*n, m+150, [n, m]*.45]);
img(image_scan_org(i_start: i_end,j_start: j_end)/255)
hold on
[X, Y] = meshgrid([1:n], [1:m]);
i_threshold = 0;
j_threshold = 0;
quiver(X(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-j_start+1,...
    Y(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-i_start+1,...
    dx_Simoncelli(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    dy_Simoncelli(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    vector_size,'Color', col(2,:))
axis off

%%  excute Savitzky-Golay results
% close all
figure('rend','painters','pos',[50+3*n, 50, [n, m]]);
img(segmented_Savitzky_Golay)
title('Savitzky-Golay kernel method')
%
figure('rend','painters','pos',[50+3.025*n, m+150, [n, m]*.45]);
imagesc(segmented_Savitzky_Golay(i_start: i_end, j_start: j_end, :), [0, 1])
set(gca, 'Xtick', [], 'Ytick', [])
axis image
axis off
%
figure('rend','painters','pos',[3.625*n, m+150, [n, m]*.45]);
img(image_scan_org(i_start: i_end,j_start: j_end)/255)
hold on
[X, Y] = meshgrid([1:n], [1:m]);
i_threshold = 0;
j_threshold = 0;
quiver(X(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-j_start+1,...
    Y(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-i_start+1,...
    dx_Savitzky_Golay(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    dy_Savitzky_Golay(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    vector_size,'Color', col(3,:))
axis off

%%  excute MaxPol results
% close all
figure('rend','painters');%,'pos',[50+4*n, 50, [n, m]]);
img(segmented_MaxPol)
title('MaxPol kernel method (proposed)')

figure('rend','painters');%,'pos',[50+4.025*n, m+150, [n, m]*.45]);
imagesc(segmented_MaxPol(i_start: i_end, j_start: j_end, :), [0, 1])
set(gca, 'Xtick', [], 'Ytick', [])
axis image
axis off
%
figure('rend','painters');%,'pos',[4.625*n, m+150, [n, m]*.45]);
img(image_scan_org(i_start: i_end,j_start: j_end)/255)
hold on
[X, Y] = meshgrid([1:n], [1:m]);
i_threshold = 0;
j_threshold = 0;
quiver(X(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-j_start+1,...
    Y(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1)-i_start+1,...
    dx_MaxPol(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    dy_MaxPol(i_start+i_threshold: i_end-i_threshold+1, j_start+j_threshold: j_end-j_threshold+1),...
    vector_size,'Color', col(4,:))
axis off
