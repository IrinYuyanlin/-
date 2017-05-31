function contour_ = initialcontour(bw)
temp = 2 * bw - 1;
contour_ = bwdist(temp > 0) - bwdist(temp <= 0);
end

