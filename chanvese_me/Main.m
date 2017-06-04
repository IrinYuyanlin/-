im = imread('brain.jpg');
im2 = imread('flowers.jpg');
im3 = imread('r06aug97.erdak.10--1---2.png');
% imshow(im3);
im3 = double(im3)/255;
im3 = imadjust(im3, [0 0.15], []);
im3 = im3 * 255;
imshow(im3, 'displayrange', [0 255]);

% seg1 = chanvese_me(im, 400, 0.5, 0.5, 'hand');
% seg2 = chanvese_me(im2, 400, 0.5, 0.5, 'auto', 'large');
% seg3 = chanvese_me(im3, 400, 0.5, 0.5, 'auto', 'large');