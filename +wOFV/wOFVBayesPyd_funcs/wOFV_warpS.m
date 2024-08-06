function Iw=wOFV_warpS(I0,V)

%Symmetric warping operator that distorts an image to Iw according to the
%vector field V using symmetric boundary conditions. I0 is a 
%griddedInterpolant object of the image. The original image and V must be 
%the same size.

n=size(V,1);
m=size(V,2);

%Initial grid
[X,Y]=meshgrid(1:n,1:m);

X=X';
Y=Y';

n=n+1;
m=m+1;

%Define motion matrix with symmetric boundary conditions
motmat=cat(3,n-max(X-V(:,:,1),n)+max(min(X-V(:,:,1),n),2-(X-V(:,:,1))),...
    m-max(Y-V(:,:,2),m)+max(min(Y-V(:,:,2),m),2-(Y-V(:,:,2))));

Iw=I0(motmat(:,:,1),motmat(:,:,2));