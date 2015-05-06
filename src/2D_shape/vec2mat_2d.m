function mat = vec2mat_2d(vec)
% convert complex vector representation of 2D shapes into matrix representation
% vec: N by S complex matrix, N: number of landmarks, S: number of samples
% mat: 2 by N by S tensor
N = size(vec,1);
S = size(vec,2);

mat = zeros(2, N, S);
for ii = 1:S
    mat(:,:,ii) = [real(vec(:,ii))'; imag(vec(:,ii))'];
end