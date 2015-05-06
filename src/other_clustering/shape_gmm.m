function [yy_gmm idx best_model] = shape_gmm(score, Maxiter, max_gnum, reg, mu, V, dim_proj)
options = statset('Display','off','Maxiter',Maxiter);
obj = cell(1,max_gnum);
BIC = zeros(max_gnum,1);
for ii = 1:max_gnum
    obj{ii} = gmdistribution.fit(score',ii,'Options',options,'Regularize',reg,'Replicates',10);
    BIC(ii) = obj{ii}.BIC;
end

[tmp tmp2] = min(BIC);
best_model = obj{tmp2};

v_proj = V(:,1:dim_proj)*best_model.mu';

s0 = sqrt(sum(v_proj.*conj(v_proj)));
sp_cord = [cos(s0);repmat(sin(s0)./(s0+eps),dim_proj,1).*best_model.mu'];

yy_gmm = pull_back([real(mu);imag(mu)], V(:,1:dim_proj), sp_cord);
idx = cluster(best_model, score');