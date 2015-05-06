function [V D1] = LapEigMap(DD,knn_num,sigma)

num = size(DD,1);
[kd IDX] = knn(DD,knn_num);

G = zeros(num,num);
for ii = 1:num
    G(ii,IDX(:,ii)) = 1;
end

W = exp(-DD.^2/sigma);
D = diag(sum(W));
L = D-W;
[V D1] = eig(L,D);

D1 = diag(D1);
[temp ind] = sort(D1,'ascend');
V = V(:,ind);