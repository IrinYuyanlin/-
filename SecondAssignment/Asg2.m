im = imread('imForTest', 'tif');
[im revertclass] = tofloat(im);
im_fft = fftshift(fft2(im));
[mm, nn] = size(im_fft);

%% gaussian noise
g1 = imnoise(im, 'gaussian'); % u = 0; sigma = 0.01;
g2 = imnoise(im, 'gaussian', 0, 0.1); % u = 0; sigma = 0.1;

g1_fft = fftshift(fft2(g1));
g2_fft = fftshift(fft2(g2));


%% salt & pepper noise
s1 = imnoise(im, 'salt & pepper'); % noise density = 0.05
s2 = imnoise(im, 'salt & pepper', 0.2); % density = 0.2

s1_fft = fftshift(fft2(s1));
s2_fft = fftshift(fft2(s2));

figure(1);
subplot(221); imshow(g1); title('gaussian_default');
subplot(222); imshow(g2); title('gaussian, u=0, sigma=0.1');
subplot(223); imshow(s1); title('s&p, default');
subplot(224); imshow(s2); title('s&p, density=0.2');
%% generate spatial filter

spa_aver = fspecial('average'); % default [3 3]
spa_gau = fspecial('gaussian'); % sigma = 0.5


%% results

% firstly, trying different types of filters to gaussian and salt & pepper
% noise respectively.

g1_spa_aver = imfilter(g1, spa_aver);
g1_spa_gau = imfilter(g1, spa_gau);

s1_spa_aver = imfilter(s1, spa_aver);
s1_spa_gau = imfilter(s1, spa_gau);

figure(2); 
subplot(221); imshow(g1_spa_aver); title('average_gau');
subplot(222); imshow(g1_spa_gau); title('gaussian_gau');
subplot(223); imshow(s1_spa_aver); title('average_s&p');
subplot(224); imshow(s1_spa_gau); title('gau_s&p')


%% frequency domain filter

[X, Y] = meshgrid(1:nn, 1:mm);
x_cen = (nn - 1) / 2; y_cen = (mm - 1) / 2;
idealp = sqrt((X-x_cen).^2 + (Y-y_cen).^2);
idealp(idealp <= 150) = 1;
idealp(idealp > 150) = 0;

D0 = 150;
D = sqrt((X-x_cen).^2 + (Y-y_cen).^2);
gaulp = exp(-(D.^2) / (2*(D0^2)));

g1_fft_ideal = g1_fft .* idealp;
g1_fft_gaulp = g1_fft .* gaulp;
s1_fft_ideal = s1_fft .* idealp;
s1_fft_gaulp = s1_fft .* gaulp;

g1_ideal = ifft2(fftshift(g1_fft_ideal));
g1_gau = ifft2(fftshift(g1_fft_gaulp));
s1_ideal = ifft2(fftshift(s1_fft_ideal));
s1_gau = ifft2(fftshift(s1_fft_gaulp));

figure(3);
subplot(221); imshow(g1_ideal); title('gau_ideal');
subplot(222); imshow(g1_gau); title('gau_gau');
subplot(223); imshow(s1_ideal); title('s&p_ideal');
subplot(224); imshow(s1_gau); title('s&p_gau');

%% appendix function

function [out, revertclass] = tofloat(in)
%TOFLOAT Convert image to floating point
%   [OUT, REVERTCLASS] = TOFLOAT(IN) converts the inout image IN to
%   floating-point. If IN is a double or single image, then OUT equals IN.
%   Otherwise, OUT equals IM2SINGLE(IN). REVERCLASS is a function handle
%   that can be used to convert back to the class of IN.

identity = @(x) x;
tosingle = @im2single;

table = {'uint8', tosingle, @im2uint8
    'uint16', tosingle, @im2uint16
    'int16', tosingle, @im2int16
    'logical', tosingle, @logical
    'double', identity, identity
    'single', identity, identity};

classIndex = find(strcmp(class(in), table(:, 1)));

if isempty(classIndex)
    error('Unsupported input image class');
end

out  = table{classIndex, 2}(in);

revertclass = table{classIndex, 3};
end
