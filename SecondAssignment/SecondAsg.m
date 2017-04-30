%% task 1: add different noise model and use imfilter function
im = imread('imForTest', 'tif');
% gaussian noise
g1 = imnoise(im, 'gaussian');
% salt and pepper noise
g2 = imnoise(im, 'salt & pepper');
% multiplicative noise
g3 = imnoise(im, 'speckle');
% poission noise
g4 = imnoise(im, 'poisson');

% take gaussian noise as a example to show on the figure
figure(1); subplot(211); imshow(im);
subplot(212); imshow(g1); title('with ganussian noise')
%% task 2: generate gaussian filter
w = fspecial('gaussian');  % default sigma is 0.5
w1 = fspecial('gaussian', 3);
g1_fil = imfilter(g1, w);
g1_fil_w1 = imfilter(g1, w1);

figure(2); subplot(211); imshow(g1_fil); title('sigma = 0.5')
subplot(212); imshow(g1_fil_w1); title('sigma = 1.5')

%% task 3: implements convolution/correlation
% we just simply assume the convolution factor is symmetric
% the way we treat edge is just set them to zero for the purpose
% of simplicity, and for the convenience, i just use a specific to
% illustrate the process of convolution, not in a general form, but the
% idea is the same.

ww = ones(3); % for demonstration
[h, w] = size(im);
im_edge = zeros(h +2 , w + 2);
im_edge(2:h + 1, 2:w + 1) = double(im(:, :)) / 255;
im_conv = zeros(h, w);
for ii = 2:h + 1
    for jj = 2:w + 1
        margin1 = ww(1, 1) * im_edge(ii - 1, jj - 1) + ww(1, 2) * ...
            im_edge(ii - 1, jj) + ww(1, 3) * im_edge(ii - 1, jj + 1);
        margin2 = ww(2, 1) * im_edge(ii, jj - 1) + ww(2, 2) * ...
            im_edge(ii, jj) + ww(2, 3) * im_edge(ii, jj + 1);
        margin3 = ww(3, 1) * im_edge(ii + 1, jj - 1) + ww(3, 2) * ...
            im_edge(ii + 1, jj) + ww(3, 3) * im_edge(ii + 1, jj + 1);
        im_conv(ii - 1, jj - 1) = (margin1 + margin2 + margin3) / 9;
    end
end

figure(3); imshow(im_conv);

%% task 4: filtering in frequency domain

% 2D ideal low pass filter
[xx, yy] = meshgrid(1:w * 2, 1:h * 2);
x_cen = (2*w - 1) / 2; y_cen = (2*h - 1) / 2;
idealp = sqrt((xx - x_cen).^2 + (yy - y_cen).^2);
idealp(idealp <= 100) = 1;
idealp(idealp > 100) = 0;

% 2D ideal high pass filter
ideahp = sqrt((xx - x_cen).^2 + (yy - y_cen).^2);
ideahp(ideahp <= 100) = 0;
ideahp(ideahp > 100) = 1;

% 2D gaussian low pass filter
sigma = 200^2; % square
cof = 1 / (sqrt(2 * pi * sigma));
margin_G = -((xx - x_cen).^2 + (yy - y_cen).^2) / (2 * sigma);
gauLow = cof * exp(margin_G);

% modify the scalar
max_gau = max(max(gauLow));
gauLow = gauLow / max_gau;

% 2D gaussian high pass filter
gauHigh = ones(2 * h, 2 * w);
gauHigh = gauHigh - gauLow;

%% compare different types of filters

% please check the Asg2.m file

%% task 5: add periodic noise and filter it in frequency domain
C = [2 256; 128 128; 2 64]; % randomly choose some frequences
[r_peri, R_peri, S_peri] = imnoise3(h, w, C);
[im_f, revertclass] = tofloat(im);
im_fft = fftshift(fft2(im_f));
im_peri = ifft2(ifftshift(im_fft + R_peri));
im_peri = revertclass(im_peri);
figure(4); imshow(im_peri);
[MM, NN] = size(im_fft);

% formula H = Hk * H-k (k = 1:Q) express the mathmatical description about
% the filter and principle. More details can see the associated
% documentation.

% use butterworth filter as the base of periodic filter.
im_fil_peri = zeros(MM, NN);
im_fil_peri = butt_peri(C, im_fft);
im_fil_spa = ifft2(ifftshift(im_fft + R_peri));
im_fil_spa = revertclass(im_fil_spa);
figure(5); imshow(im_fil_spa);

% its my implement of periodic filter, but it seems not work...
function im_fil_peri = butt_peri(c, im_fft)
%BUTT_PERI generates a periodic filter specified by the m-by-2 matrix c;

mm = length(c);
[h, w] = size(im_fft);
H = ones(h, w);
for ii = 1:mm
    u0 = round(c(ii, 1));
    v0 = round(c(ii, 2));
    
    if u0 < 3
        u0 = 3;
    end
    if v0 < 3
        v0 = 3;
    end
    
    H1 = ones(h, w);
    H1(u0 - 2:u0 + 2, v0 - 2:v0 + 2) = 0;
    H = H1 .* H;
end

im_fil_peri = H .* im_fft;

end

%% bonus: adaptive filter

function gg = adaptive_me(im, varargin)
%ADAPTIVE_ME computes the filtered results using adaptive filter.
% "varagin" default: 1; represent the estimation of global noisy variance.

error(nargchk(1, 2, nargin))

if isempty(varargin)
    sigma = 1;
else
    sigma = varargin{1};
end

[mm, nn] = size(im);
gg = zeros(mm, nn);
im_full = zeros(mm + 2, nn + 2);
im_full(2:mm + 1, 2:nn + 1) = im(:, :);
for ii = 1:mm
    for jj = 1:nn
        ml = sum(sum(im_full(ii:ii + 2, jj:jj + 2))) / 9;
        sigma2 = sqrt((im_full(ii:ii + 2, jj:jj + 2) - ml).^2);
        gg(ii, jj) = im(ii, jj) - (sigma / sigma2) * (im(ii, jj) - ml);
    end
end

end

%% appendix function

function [r, R, S] = imnoise3(M, N, C, A, B)
%IMNOISE3 generates periodic noise
% r: A*sin(2*pi*u0*(x+Bx)/M + 2*pi*v0*(y+By)/N)
% R: the frequency form of r
% S: spectrum of r

K = size(C, 1);
if nargin < 4
    A = ones(1, K); % default value
end
if nargin < 5
    B = zeros(K, 2); 
end

R = zeros(M, N);
for j = 1:K
    u1 = floor(M/2) + 1 - C(j, 1);
    v1 = floor(N/2) + 1 - C(j, 2);
    cof = i * M * N * (A(j) / 2);
    margin1 = exp(-i * 2 * pi * (C(j, 1) * B(j, 1) / M) + ...
        (C(j, 2) * B(j, 2) / N));
    R(u1, v1) = cof * margin1;
    
    % due to the symetric feature of fft
    u2 = floor(M/2) + 1 - C(j, 1);
    v2 = floor(N/2) + 1 - C(j, 2);
    R(u2, v2) = -cof * margin1;
end

S = abs(R);
r = real(ifft2(ifftshift(R)));

end

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

