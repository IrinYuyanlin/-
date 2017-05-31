function kapa = kappa(I, method)
% this function computes kappa of phi0
% method specified use which method to compute
% default is non-linear method
p = size(I, 3);
if p == 3
    I = double(rgb2gray(I));
else
    I =double(I);
end

if (nargin ~= 2)
    error('Input argument must be 2.')
end

switch method
    case 'formula'
        [m n] = size(I);
        temp = padarray(I, [1, 1], 1, 'pre');
        temp = padarray(temp, [1, 1], 1, 'post');
        
        fy = temp(3:end, 2:n+1) - temp(1:m, 2:n+1);
        fx = temp(2:m+1, 3:end) - temp(2:m+1, 1:n);
        fyy = temp(3:end, 2:n+1) - 2*temp(2:m+1, 2:n+1) + temp(1:m, 2:n+1);
        fxx = temp(2:m+1, 3:end) - 2*temp(2:m+1, 2:n+1) + temp(2:m+1, 1:n);
        fxy = 0.25 * (temp(3:end, 3:end) - temp(3:end, 2:n+1) - temp(2:m+1, 3:end)...
             + temp(1:m, 1:n));
%         fxy = 0.25.*(temp(3:end, 3:end) - temp(1:m, 3:end) + temp(3:end, 1:n)...
%             - temp(1:m, 1:n));
        G = (fx.^2 + fy.^2).^0.5;
        K = (fxx.*fx.^2 + fyy.*fy.^2 - 2*fxy.*fx.*fy) ./ (eps + fx.^2 + fy.^2).^1.5;
        kapa = K.*G;
        kapa(1, :) = eps;
        kapa(end, :) = eps;
        kapa(:, end) = eps;
        kapa(:, 1) = eps;
        kapa = kapa / max(max(abs(G)));
    case 'matrix'
        h = [-1/16, 5/16, -1/16; 5/16, -1, 5/16; -1/16, 5/16, -1/16];
        kapa = imfilter(I, h);
    otherwise
        error('method must be "formula" or "matrix".')
end


end

