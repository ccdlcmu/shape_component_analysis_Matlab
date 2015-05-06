function  CM = contingency_table(cluster1,cluster2)
% cluster1: a cell list of member indes for each cluster.
% cluster2: a cell list of member indes for each cluster.
% CM(i,j) = # of cluster1(i) ^ cluster2(j)
CM = zeros(length(cluster1),length(cluster2));
for ii = 1:length(cluster1)
    for jj = 1:length(cluster2)
        CM(ii,jj) = numel(intersect(cluster1{ii},cluster2{jj}));
    end
end