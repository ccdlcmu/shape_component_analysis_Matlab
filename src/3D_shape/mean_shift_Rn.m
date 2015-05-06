function yt = mean_shift_Rn(y0, sigma, iter_max)
% clear 
% close all
% data = randn(2,1000);
% data(1,1:200) = data(1,1:200) + 3;
% data(2,1:200) = data(2,1:200) + 3;
% data(1,201:250) = data(1,201:250) + -2;
% data(2,201:250) = data(2,201:250) + -4;
% yt = data';
yt = y0;
dim = size(y0, 2);
for ii = 1:iter_max
    distmap = pdist2(yt,y0);
    W = exp(-distmap.^2/2/sigma^2);
    Ws = sum(W,2);
    yt = W* yt./repmat(Ws,1,dim);

end