function [SI] = si_rgb(rgb_image)

for channel = 1: size(rgb_image, 3)
    SI(channel) = sharpness_index(rgb_image(:, :, channel));
%     SI(channel) = s_index(rgb_image(:, :, channel));
end
SI = mean(SI);