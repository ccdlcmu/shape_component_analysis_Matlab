function distmap = pdist2_NDriemann(xt, yt)
% compute pairwise distance matrix using Riemannian metric. 
% xt, yt, D by N by S array, xt(:,:,ii) is the i-th example
% if receive one input, output pairwise distance between examples.

if nargin == 1  
    s_numx = size(xt,3);

    distmap = zeros(s_numx,s_numx);
    for ii = 1:s_numx
        for jj = (ii+1):s_numx
            [U S V] = svd(xt(:,:,ii)*xt(:,:,jj)');
            distmap(ii,jj) = abs(acos(sum(abs(diag(S)))));            
        end
    end

    distmap = distmap + distmap';
else
    s_numx = size(xt,3);
    s_numy = size(yt,3);
    
    distmap = zeros(s_numx,s_numy);
    for ii = 1:s_numx
        for jj = 1:s_numy
            [U S V] = svd(xt(:,:,ii)*yt(:,:,jj)');
            distmap(ii,jj) = abs(acos(sum(abs(diag(S)))));            
        end
    end

end
