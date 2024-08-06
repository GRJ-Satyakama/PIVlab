function [K, bases] = derivdirec_Gaussian(l, sigma_x, sigma_y, n, theta)

%DERIVDIREC_GAUSSIAN 2D direction (steerable) Gaussian derivative kernel
%
K = 0;
for k = 0: n
    %%  lowpass derivative filter with cutoff level (P_x, P_y)
    if mod(n - k, 2)        
        [d_x] = gaussian_kernel(l, sigma_x, n-k);
        d_x = [d_x, 0];
    else
        [d_x] = gaussian_kernel(l, sigma_x, n-k);  
        d_x = d_x;
    end
    %
    if mod(k, 2)
        [d_y] = gaussian_kernel(l, sigma_y, k);
        d_y = d_y';
    else
        [d_y] = gaussian_kernel(l, sigma_y, k);
        d_y = d_y';
    end
    %
    if mod(k, 2)
        if mod(n, 2)
            if cos(theta)>0
                d_y = [d_y; 0];
            else
                d_y = [0; d_y];
            end
        else
            if cos(theta)<0
                d_y = [d_y; 0];
            else
                d_y = [0; d_y];
            end
        end
    end
    %
    %%  2D basis filter construction
    bases{k+1} = [d_y * d_x];
    %
    %%  bionomial coefficients
    c(k+1) = nchoosek(n, k) * cos(theta)^(n-k) * (-sin(theta))^k;
    %
    %%  accumulate 2D basis
    K = K + c(k+1) * bases{k+1};
end

for k = 0: n
    bases{k+1} = bases{k+1};
end