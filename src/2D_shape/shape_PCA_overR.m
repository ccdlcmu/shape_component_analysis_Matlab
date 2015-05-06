function [mu V score sp sp_cord ev dim_proj] = shape_PCA_overR(y0, cutoff)

mu = mean_shape_riemann(y0,3);
s_num = size(y0,2);

dim = size(y0,1);
dim_sub = dim-1;

H = -gallery('orthog',dim,4);
H(1,:) = [];

y1 = H*y0;
mu = H*mu;
y1 = y1./ sqrt(repmat(diag((y1')*y1)',dim_sub,1));
mu = mu./ sqrt(repmat(diag((mu')*mu)',dim_sub,1));

mulog = Log_map(mu,y1);
mulogR = [real(mulog); imag(mulog)];
[U S V] = svd(mulogR');
ev = diag(S).^2;

dim_proj = find(cumsum(ev)/sum(ev)>cutoff,1,'first');

score = V(:,1:dim_proj)'*mulogR;

[sp_cord v_proj sp] = score_to_cord(score, V(:,1:dim_proj), mu);


