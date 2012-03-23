function [x] = softThresh(x,s)
% softThresh performs soft thresholding on the vector x using positive
% threshold s
x(abs(x)<s) = 0;
x(x>=s) = x(x>=s)-s;
x(x<=-s) = x(x<=-s)+s;
end