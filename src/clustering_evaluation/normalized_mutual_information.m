function NMI = normalized_mutual_information(CM)
% calculate the normalized mutual information form a contingency_table.
% NMI = 2*I(cluster2; cluster1)/ (H(cluster1) + H(cluster2))
% CM(i,j) = # of cluster1(i) ^ cluster2(j)
N = sum(CM(:));
Pi = sum(CM,1)/N + eps;
Pj = sum(CM,2)/N + eps;
Pij = CM/N + eps;
NMI = -2 * sum(sum(Pij.*log(Pij./(Pj*Pi)))) / ( sum(sum(Pi.*log(Pi)))  +  sum(sum(Pj.*log(Pj))) );