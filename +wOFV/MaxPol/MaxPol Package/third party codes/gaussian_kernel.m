function [gaussKernel, x] = gaussian_kernel(l, sgm, n_OD)

if mod(n_OD, 2)
    x = (-l+1/2): (l-1/2);
else
    x = -l: l;
end
c = 1/(sqrt(2*pi)*sgm);
gaussKernel = c * exp(-(x.^2)/(2*sgm^2));
gaussKernel = gaussKernel/sum(gaussKernel);
switch n_OD
    case 0
        d = 1;
    case 1
        d = -x/sgm^2;
    case 2
        d = (x.^2-sgm^2)/sgm^4;
    case 3
        d = -(x.^3 - 3*x*sgm^2)/sgm^6;
    case 4
        d = (x.^4 - 6*x.^2*sgm^2 + 3*sgm^4)/sgm^8;
end
gaussKernel = d .* gaussKernel;
gaussKernel = gaussKernel(:)';