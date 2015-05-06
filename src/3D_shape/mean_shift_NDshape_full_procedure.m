function yc = mean_shift_NDshape_full_procedure(yr, y0, sigma, iter, U)

[yt distmap] = mean_shift_NDshape(yr, y0, sigma, iter);
[modes] = return_cluster(abs(distmap));

if nargin==5 
    yc = return_example_2DPCA(yt, U, modes);
else
    yc = zeros([size(yr,1) size(yr,2) length(modes)]); 
    for ii = 1:length(modes)
        yc(:,:,ii) = yt(:,:,modes{ii}(1));
    end
    
end
