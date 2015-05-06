function vec = mat2vec_2d(mat)
% convert matrix representation of 2D shapes into complex vector representation
% mat: 2 by N by S tensor, N: number of landmarks, S: number of samples
% vec: N by S complex matrix

vec = reshape(mat(1,:,:) + 1i * mat(2,:,:), size(mat, 2), size(mat,3));

