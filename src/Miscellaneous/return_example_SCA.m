function yc = return_example_SCA(yt, U, modes)


% ind = zeros(1,length(modes)); 
yc = zeros([size(yt,1) size(U,1) length(modes)]);
for ii = 1:length(modes)
    yc(:,:,ii) = yt(:,:, modes{ii}(1))*U';
end



