function Iw=wOFV_warpP(I0,V)

%Periodic warping operator that distorts an image to Iw according to the
%vector field V using periodic boundary conditions. I0 is a 
%griddedInterpolant object of the image. The original image and V must be 
%the same size.

n=size(V,1);

%Initial grid
[Y,X]=ndgrid(1:n,1:n);

%Define motion matrix with periodic boundary conditions
motmat=cat(3,mod(X-V(:,:,2)-1,n)+1,mod(Y-V(:,:,1)-1,n)+1);

Iw=I0(motmat(:,:,2),motmat(:,:,1));