function [yt U S] = shape_factor_analysis(y0, dim_reduced)
% y0 s-by-d-by-n tensor, 
% where s = dim(landmark), d = # landmarks, n = #examples

s_num = size(y0,1);
dim = size(y0,2);
n_sample = size(y0,3);

for ii = 1:n_sample
    y0(:,:,ii) = y0(:,:,ii)/norm(y0(:,:,ii),'fro');
end

% tic
cov = zeros(dim,dim);
for ii = 1:n_sample 
    ym = y0(:,:,ii);
    cov = cov + ym'*ym;
end

[U S] = eigs(cov, dim_reduced);
S = diag(S);
% toc

yt = zeros(s_num, dim_reduced, n_sample);
for ii = 1:n_sample
    ym = y0(:,:,ii);
    yt(:,:,ii) = ym*U;
    yt(:,:,ii) = yt(:,:,ii)/norm(yt(:,:,ii),'fro');
end
