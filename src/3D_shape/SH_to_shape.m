function  shapes = SH_to_shape(SHs)
% Transhform an array of shperical harmonic
% coefficients to shape representatives
% SH: 3 * (Lc/2+1) *  N array of spherical harmonic coefficients
% shapes: 3*Lc*N array of shapes, Lc/2+1 is the total number of SH coefficients, N number of shapes

shapes = zeros([size(SHs,1) 2*size(SHs,2)-2 size(SHs,3)]); 
for ii = 1:size(SHs,3)    
        shapes(:,:,ii) = cat(2, real(SHs(:,2:end,ii)), imag(SHs(:,2:end,ii)));
        shapes(:,:,ii) = shapes(:,:,ii)/norm(shapes(:,:,ii), 'fro');
end