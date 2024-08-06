function grd=OFgrad(A)
%Computes the gradient using MATLAB's functionality from a scalar matrix A,
%sorting the two components into their proper places for optical flow
%calculations

[gdy,gdx]=gradient(A);
grd=cat(3,gdx,gdy);