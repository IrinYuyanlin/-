%%  Task 1:Import images

im = imread('strawberry', 'tif');
[y, x, z] = size(im);
figure(1); imshow(im); title('My original picture');

%%  Task 2: counting pixels of different intensities and colors

inten = 256;
n = 1:inten;
m = zeros(inten, z);
for ii = 1:3
    for jj = 1:256
        m(jj, ii) = length(find(im(:, :, ii) == (jj - 1)));
    end
end

figure(2); 
subplot(3, 1, 1); stem(n - 1, m(:, 1)); title('The red channel');
subplot(3, 1, 2); stem(n - 1, m(:, 2)); title('The green channel');
subplot(3, 1, 3); stem(n - 1, m(:, 3)); title('The blue channel');

figure(3);
im_red = im(:, :, 1); subplot(3, 1, 1); imhist(im_red); title('Red');
im_green = im(:, :, 2); subplot(3, 1, 2); imhist(im_green); title('Green');
im_blue = im(:, :, 3); subplot(3, 1, 3);; imhist(im_blue); title('Blue');

%%  Task 3: Implementing the interpolation function for gray-scale using binary interpolation method

im2 = imread('chest_x_ray', 'bmp');
figure(4); imshow(im2); title('My original gray-scale picture');

[y1 x1 z1] = size(im2); im2 = im2double(im2);

scale = 2;  % You can change the scaling factor here

ms = zeros(y1 * scale, x1 * scale);
m1 = zeros(y1 + 2, x1 + 2);
m1(2:y1 + 1, 2:x1 + 1) = im2(1:y1, 1:x1);
m1(1, 1) = im2(1, 1); m1(y1 + 2, x1 + 2) = im2(y1, x1); m1(1, x1 + 2) = im2(1, x1); m1(y1 + 2, 1) = im2(y1, 1);
m1(2:y1 + 1, 1) = im2(1:y1, 1); m1(2:y1 + 1, x1 + 2) = im2(1:y1, x1);
m1(1, 2:x1 + 1) = im2(1, 1:x1); m1(y1 + 2, 2:x1 + 1) = im2(y1, 1:x1);



[h w z] = size(ms);
for mi = 1:h
    for mj = 1:w
        i = mi / scale; j = mj / scale;
        ii = floor(i); jj = floor(j);
        u = i - ii; v = j - jj;
        ii = ii + 1; jj = jj + 1;
        f_margin0 = (1 - v) * m1(ii, jj) + v * m1(ii, jj + 1);
        f_margin1 = (1 - v) * m1(ii + 1, jj) + v * m1(ii + 1, jj + 1);
        ms(mi, mj) = (1 - u) * f_margin0 + u * f_margin1;
    end
end

figure(5); imshow(ms); title('This is the scaling picture');

%%  Task 4: Imaging rotation for a gray-scale image

im3 = imread('chest_x_ray', 'bmp');
figure(6); imshow(im3); title('Original picture');

degree = 45; % You can change the degree that you want to rotate
[h1 w1 z1] = size(im3);
rad = degree * pi / 180;
r = [cos(rad) -sin(rad) 0; sin(rad) cos(rad) 0; 0 0 1];
p_l_u = [1 1 1] * r;
p_r_u = [1 w1 1] * r;
p_l_d = [h1 1 1] * r;
p_r_d = [h1 w1 1] * r;

height = round(max([(abs(p_l_u(1) - p_r_d(1))) + .5 abs(p_r_u(1) - p_l_d(1)) + .5])); % compute height
width = round(max([(abs(p_l_u(2) - p_r_d(2))) + .5 abs(p_r_u(2) - p_l_d(2)) + .5]));
im3_r = zeros(height, width);

delta_y = abs(min([p_l_u(1) p_l_d(1) p_r_u(1) p_r_d(1)]));
delta_x = abs(min([p_l_u(2) p_l_d(2) p_r_u(2) p_r_d(2)]));

for i = (1 - delta_y):(height - delta_y)
    for j = (1 - delta_x):(width - delta_x)
        p = [i j 1] / r;
        u = p(1) - floor(p(1));
        v = p(2) - floor(p(2));
        float_Y = u;
        float_X = v;
        
        if ((p(1) >= 1) && (p(2) >= 1) && (p(1) <= h1) && (p(2) <= w1)) % check out if the matching pixels is out of origial picture
            
%             ind1 = floor(p(1)); % down y
%             ind2 = floor(p(2)); % left x
%             ind3 = ceil(p(1)); % up y
%             ind4 = ceil(p(2)); % right x
%             
%             % now lets apply our binary interpolation implementation
%             f_margin0 = im3(ind1, ind2) * (v - 1) / (0 - 1) + im3(ind1, ind2 + 1) * (v - 0) / (1 - 0);
%             f_margin1 = im3(ind1 + 1, ind2) * (v - 1) / (0 - 1) + im3(ind1 + 1, ind2 + 1) * (v - 0) / (1 - 0);
%             im3_r(i + delta_y, j + delta_x) = f_margin0 * (u - 1) / (0 - 1) + f_margin1 * (u - 0) / (1 - 0);
            
% but my result is really disappointing, so i tried another method
% i still reserve my fail code, because i want to know where the mistake happens.
            
% this code uses weithed neareat interpolation.

            pix_up_left = [floor(p(1)) floor(p(2))];
            pix_up_right = [floor(p(1)) ceil(p(2))];
            pix_down_left = [ceil(p(1)) floor(p(2))];
            pix_down_right = [ceil(p(1)) ceil(p(2))]; 
        
            value_up_left =( 1 - float_X) * (1 - float_Y);
            value_up_right =float_X * (1 - float_Y);
            value_down_left = (1 - float_X) * float_Y;
            value_down_right = float_X * float_Y;
                                                            
            im3_r(i + delta_y, j + delta_x) = value_up_left * im3(pix_up_left(1), pix_up_left(2))+ ...
                                        value_up_right * im3(pix_up_right(1), pix_up_right(2))+ ...
                                        value_down_left * im3(pix_down_left(1), pix_down_left(2))+ ...
                                        value_down_right * im3(pix_down_right(1), pix_down_right(2));

            
        end
    end
end

figure(7); imshow(uint8(im3_r)); title('This is the rotated picture');
