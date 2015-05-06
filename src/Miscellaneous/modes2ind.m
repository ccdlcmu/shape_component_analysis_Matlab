function ind = modes2ind(modes)

ind = zeros(1,length(modes)); 
for ii = 1:length(modes)
    ind(ii) = modes{ii}(1);
end