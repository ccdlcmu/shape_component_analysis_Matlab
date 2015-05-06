function yc = return_example_PGA(yt, mu, V, modes)


ind = zeros(1,length(modes)); 
for ii = 1:length(modes)
    ind(ii) = modes{ii}(1);
end

yc = pull_back([real(mu);imag(mu)], V, yt(:,ind));
dim = size(yc,1)/2;
yc = yc(1:dim,:) + yc((dim+1):end,:)*1i;