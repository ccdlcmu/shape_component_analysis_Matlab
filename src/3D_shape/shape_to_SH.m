function SHs = shape_to_SH(shapes)
% Transhform an array of shape representatives to shperical harmonic
% coefficients
% shapes: 3*Lc*N array of shapes, Lc/2+1 is the total number of SH coefficients, N number of shapes
% SH: 3*(Lc/2+1)*N array of spherical harmonic coefficients

D = size(shapes,1);
Lc = size(shapes,2);
SHs = cat(2, zeros(D,1), shapes(:,1:Lc/2) + 1i * shapes(:,Lc/2 + (1:Lc/2)));