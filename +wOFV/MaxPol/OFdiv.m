function dv=OFdiv(V)
%Divergence of a vector field using MATLAB's computation of the divergence
%but sorted into a way that can be used in optical flow computations
dv=divergence(V(:,:,2),V(:,:,1));