function showcontour(I, phi0, iter)
% this function shows the evolution of contour
imshow(I, 'initialmagnification', 'fit', 'displayrange', [0 255]);
hold on;

contour(phi0, [0 0], 'r', 'LineWidth', 4);
contour(phi0, [0 0], 'g', 'LineWidth', 1.3);
hold off;
title([num2str(iter) ' Iterations']);
drawnow;

end

