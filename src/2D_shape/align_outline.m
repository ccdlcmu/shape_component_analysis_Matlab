function [z2] = align_outline(y0,y1)
% output z2 the alignment of y1 with respect to y0
% y0 and y1 are complex column vectors representing landmark shapes

dim = size(y0,1);
H = -gallery('orthog',dim,4);
H(1,:) = [];
z0 = H'*H*y0;
z1 = H'*H*y1;
theta0 = -angle(z0'*z1);
z2 = exp(1i*theta0)*z1;




