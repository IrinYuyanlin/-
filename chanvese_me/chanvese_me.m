function seg = chanvese_me(im, num_iters, mu, delta, varargin)
% im: input image
% num_iters: number of iterations
% varargin{1}: how to set initial mask, auto or hand
% varargin{2}: must specified if first one is hand
segMethod = varargin{1};
disp('You choose: '); disp(segMethod);
if length(varargin) == 2
    segSize = varargin{2};
    disp('The size you choose is: '); disp(segSize);
end


switch length(varargin)
    case 2
        if strcmp(segMethod, 'auto')
            type = varargin{2};
            [m, n, p] = size(im);
            if p == 3
                im = rgb2gray(im);
                im = double(im);
            else
                im = double(im);
            end

            % create initial mask
            if ischar(type)
                switch lower(type)
                    case 'small'
                        mask = createmask(im, 'small');
                    case 'medium'
                        mask = createmask(im, 'medium');
                    case 'large'
                        mask = createmask(im, 'large');
                    otherwise
                        error('type must be small, medium, large.')
                end
            end
        end
        
    case 1
        p = size(im, 3);
        if p == 3
            im = double(rgb2gray(im));
        else
            im = double(im);
        end
        figure(1);
        imshow(im, 'displayrange', [0 255]);
        hold on;
        mask = roipoly;
    otherwise
        error('Input argument must be 4 or 5.')

end


phi0 = initialcontour(mask);
figure(2); title('original contour');
contour(phi0, [0 0], 'b');

for iter = 1:num_iters
    idx_outside = find(phi0 > 0);
    idx_inside = find(phi0 <= 0);
    
    c_outside = sum(sum(im(idx_outside))) / (length(idx_outside)+eps);
    c_inside = sum(sum(im(idx_inside))) / (length(idx_inside)+eps);
    
    lambda1 = 1;
    lambda2 = 1;
    
    var_outside = lambda1*(im-c_outside).^2;
    var_inside = lambda2*(im-c_inside).^2;
    
    kapa = kappa(phi0, 'formula');
    force = mu*kapa/max(max(abs(kapa))) + var_inside - var_outside;
    force = force / max(max(abs(force)));
    
    phi0 = phi0 + delta*force;
    showcontour(im, phi0, iter);
    seg = phi0 <= 0;
end
figure(3);
imshow(seg); title('Final segmetation');


end

