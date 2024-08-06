function strrt=OFstrainrate(V,lD,PD,lS,PS)
%Compute the strain rate of a velocity field V using MaxPol differentiation
%and smoothing (Hosseini and Plataniotis, 2017(a,b)). The strategy is to
%differentiate in one direction and smooth in the other using independent 
%kernels, which mitigates noise. V is a 2D velocity field in optical flow
%format (u and then v concatenated in dimension 3, u is positive down, v is
%positive to the right, both in units of pixels). lD, PD, lS, PS are 
%optional parameters passed to MaxPol describing the kernels for the 
%differentiation (lD and PD) and smoothing (lS and PS). If left empty, they
%default to a more-or-less universally effective kernel used by McManus in 
%his PhD thesis (McManus, 2019). The output is strrt, a matrix of the same 
%size as u and v containing the strain rate.

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

%Extend V by 4 points on all sides using symmetric BC's
for t=1:2
    V=[V(1:4,:,:);V;V(end:-1:end-3,:,:)];
    V=permute(V,[2,1,3]);
end

Dx=derivmtx(lD,PD,1,size(V,1),nod,sparseflag,symflag);
Sx=derivmtx(lS,PS,0,size(V,2),nod,sparseflag,symflag);

%If V is square, Dy=Dx' and Sy=Sx'
if size(V,1)==size(V,2)
    Dy=Dx';
    Sy=Sx';
else
    Dy=derivmtx(lD,PD,1,size(V,2),nod,sparseflag,symflag)';
    Sy=derivmtx(lS,PS,0,size(V,1),nod,sparseflag,symflag)';
end

strrt=.5*(Sy*V(:,:,1)*Dy+Dx*V(:,:,2)*Sx);

%Truncate added boundary points
strrt=strrt(5:end-4,5:end-4);