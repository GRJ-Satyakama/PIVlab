function [out] = directional_sharpening(blurred_image, params)

%%
 warning off
blurred_image = double(blurred_image);
%%
Theta = [0: params.degree_step: pi-params.degree_step];
image_laplacian = 0;
K = 0;
for orienting_iteration = 1: numel(Theta)
    theta = Theta(orienting_iteration);
    switch params.method
        case 'MaxPol'           
            [S{1}] = derivdirec(params.l, params.px, params.py, 2, theta, true);            
            K = K + S{1};
        case 'Savitzky_Golay'
            [S{1}] = derivdirec_SavitzkyGolay(params.l, params.px, params.py, 2, theta, true);
            K = K + S{1};
        case 'gaussian'
            [S{1}] = derivdirec_Gaussian(params.l, params.sig, params.sig, 2, theta);
            K = K + S{1};
        case 'Farid_Simoncelli'
            [S] = steerable_OD_FaridSimoncelli(params.l, theta);
            K = K + S{3};
    end
end

image_laplacian = imfilter(blurred_image, K, 'symmetric', 'conv');
sharpened_image = blurred_image - params.lambda * image_laplacian / sqrt(numel(Theta)) / 2;
warning on
%%  excute output
out.sharpened_image = sharpened_image;
out.blurred_image = blurred_image;
out.image_laplacian = image_laplacian;