function yc = return_example_tanget(yt, mu, V, modes)

ind = zeros(1,length(modes)); 
for ii = 1:length(modes)
    ind(ii) = modes{ii}(1);
end

[sp_cord] = score_to_cord(yt(:,ind), V, mu);

yc = pull_back([real(mu);imag(mu)], V, sp_cord);
dim = size(yc,1)/2;
yc = yc(1:dim,:) + yc((dim+1):end,:)*1i;

