function mask = createmask(im, type)
% create initial for chanvese image segmentation

h = [0 1 0; 1 -4 1; 0 1 0];
temp = imfilter(im, h);
thre = max(max(abs(temp))) * .5;
idx = find(abs(temp) >= thre);
[m, n] = size(im);
[x, y] = meshgrid(1:m, 1:n);
[axis_x, axis_y] = ind2sub([m, n], idx);
cx = round(mean(axis_x));
cy = round(mean(axis_y));
mask = zeros(m, n);

switch lower(type)
    case 'small'
        r = 10;
        mask(((x-cx).^2 + (y-cy).^2) < r.^2) = 1;
    case 'medium'
        r1 = 50;
        r2 = 0.2 * min([m, n]);
        r = min([r1, r2]);
        mask(((x-cx).^2 + (y-cy).^2) < r.^2) = 1;
    case 'large'
        r1 = 150;
        r2 = (1/3) * min([m, n]);
        r = min([r1, r2]);
        mask(((x-cx).^2 + (y-cy).^2) < r.^2) = 1;
    otherwise
        error('type must be small, medium, large.')
end

end

