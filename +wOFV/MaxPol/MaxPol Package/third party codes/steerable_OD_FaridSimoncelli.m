function [S] = steerable_OD_FaridSimoncelli(l, theta)

%%  Farid-Simoncelli OD
% [dlowpass] = FaridSimoncelli_OD_LP(l);
ll = l;
load Farid_Simoncelli_TIP2004
dlowpass_library = dlowpass;
clear dlowpass l
l = ll;
dlowpass = dlowpass_library(l-1, :);
dlowpass{2} = [0; dlowpass{2}];

dx  = dlowpass{1} * dlowpass{2}';
dy  = dlowpass{2} * dlowpass{1}';
dxx = dlowpass{1} * dlowpass{3}';
dyy = dlowpass{3} * dlowpass{1}';
dxy = dlowpass{2} * dlowpass{2}';

S{1} = dlowpass{1} * dlowpass{1}';
%   rotate
S{2} = cos(theta)*dx - sin(theta)*dy;
S{3} = cos(theta)^2*dxx + sin(theta)^2*dyy - 2*cos(theta)*sin(theta)*dxy;
