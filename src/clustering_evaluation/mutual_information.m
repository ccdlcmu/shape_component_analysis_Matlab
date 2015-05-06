function MI = mutual_information(CM)
% calculate the mutual information form a contingency_table.
% MI = I(cluster2; cluster1)
% CM(i,j) = # of cluster1(i) ^ cluster2(j)
N = sum(CM(:));
Pi = sum(CM,1)/N + eps;
Pj = sum(CM,2)/N + eps;
Pij = CM/N + eps;
MI = sum(sum(Pij.*log(Pij./(Pj*Pi))));