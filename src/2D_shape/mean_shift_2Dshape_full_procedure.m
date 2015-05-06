function [yc] = mean_shift_2Dshape_full_procedure(y0,sigma,iter, V, mu)

% perform nonlinear mean-shift clustering
[yt distmap] = mean_shift_2Dshape(y0, sigma, iter);

% discard duplicated modes
[modes] = return_cluster(distmap);

% discard duplicated modes
if nargin==5 
    yc = return_example_PGA(yt, mu, V, modes);
else
    ind = zeros(1,length(modes)); 
    for ii = 1:length(modes)
        ind(ii) = modes{ii}(1);
    end
    yc = yt(:,ind);
end
