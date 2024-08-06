function grad=OFgradMP(A,lD,PD,lS,PS)
%Compute the gradient of a scalar field A using MaxPol 
%differentiation and smoothing (Hosseini and Plataniotis, 2017(a,b)). The 
%strategy is to differentiate in one direction and smooth in the other 
%using independent kernels, which mitigates noise. A is a scalar field,
%and grad will be in optical flow format (ddx and then ddy concatenated in 
%dimension 3, ddx is positive down, ddy is positive to the right, both in 
%units of pixels). lD, PD, lS, PS are optional parameters passed to MaxPol 
%describing the kernels for the differentiation (lD and PD) and smoothing 
%(lS and PS). If left empty, they default to a more-or-less universally 
%effective kernel used by McManus in his PhD thesis (McManus, 2019). The 
%output is grad, a 3D matrix of the same size as A but with 2 entries along
%dimension 3 containing the gradient.

if nargin==1
    %Set default kernel parameters if not passed to the function
    lD=3;
    PD=2;
    lS=2;
    PS=2;
end

%Use a centralized scheme by default, the other option is 'staggered.'
nod='centralized';
%Compute S and D as sparse matrices for speed
sparseflag=1;
%Compute using the symbolic toolbox. If the symbolic toolbox is not 
%installed, this step will throw an error. If so, set to 'false'. Be aware
%that stability issues may occur as a result.
symflag='true';

%Extend A by 4 points on all sides using symmetric BC's
for t=1:2
    A=[A(1:4,:);A;A(end:-1:end-3,:)];
    A=A';
end

Dx=derivmtx(lD,PD,1,size(A,1),nod,sparseflag,symflag);
Sx=derivmtx(lS,PS,0,size(A,2),nod,sparseflag,symflag);

%If A is square, Dy=Dx' and Sy=Sx'
if size(A,1)==size(A,2)
    Dy=Dx';
    Sy=Sx';
else
    Dy=derivmtx(lD,PD,1,size(A,2),nod,sparseflag,symflag)';
    Sy=derivmtx(lS,PS,0,size(A,1),nod,sparseflag,symflag)';
end

grad=cat(3,Dx*A*Sx,Sy*A*Dy);

%Truncate added boundary points
grad=grad(5:end-4,5:end-4,:);