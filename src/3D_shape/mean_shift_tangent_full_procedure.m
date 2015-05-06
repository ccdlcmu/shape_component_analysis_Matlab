function [yc gplabel] = mean_shift_tangent_full_procedure(score, sigma, iter)

yt = mean_shift_Rn(score', sigma, iter);
distmap = pdist2(yt, yt);

% yc_tg = return_example_tanget(yt_tg', mu, V(:,1:dim_proj), modes_tg);
[modes] = return_cluster(distmap);
ind = modes2ind(modes);
yc = yt(ind,:);

distmap = pdist2(score', yc);
[temp gplabel] = min(distmap,[],2);
yc = yc';