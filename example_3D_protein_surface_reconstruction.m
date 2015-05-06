clear 
%% reproduce table 6

addpath(genpath('src\'))
addpath(genpath('data\'))
addpath(genpath('third_party_packages\'))

load fvec_3proteins
% load fvec_smooth

load sphere_standard
% reshape fvec output from SHARM
deg = 10;

fvec_all(1:(deg+1)^2,:,:);
fvec_tmp = zeros([size(fvec_all,2), size(fvec_all,1), size(fvec_all,3)]);

for ii = 1:size(fvec_all,3)    
    fvec_tmp(:,:,ii) = fvec_all(:,:,ii)';
end

% 
y0 = SH_to_shape(fvec_tmp);


s_num = size(y0,1);
dim = size(y0,2);
n_sample = size(y0,3);


tic
[yt U S] = shape_factor_analysis(y0, 160);

clear yt
toc



dimr = [40 80 160];
for ii = 1:n_sample
    
    ym = y0(:,:,ii);
    y0rec = shape_to_SH(ym);
    Zvert1 = real(Zs(:,1:(deg+1)^2)*y0rec(:,1:(deg+1)^2)');
    
    for rr = 1:length(dimr)
        Ut = U(:,1:dimr(rr));
        yt = ym*Ut;
        yt = yt/norm(yt,'fro');
        ytt = yt*Ut';

        ytrec = shape_to_SH(ytt);        
        Zvert = real(Zs(:,1:(deg+1)^2)*ytrec(:,1:(deg+1)^2)');
        
        error_fvec(ii,rr) = norm(ym - ytt,'fro')/norm(ym,'fro');
        error_surf(ii,rr) = norm(Zvert - Zvert1,'fro')/norm(Zvert1,'fro');
        error_surfarea(ii,rr) = abs(trimeshSurfaceArea(Zvert, fs) - trimeshSurfaceArea(Zvert1, fs))/trimeshSurfaceArea(Zvert1, fs);
    end
    
end

results_mean = [mean(error_surfarea,1) ;mean(error_fvec,1); mean(error_surf,1)]
results_std = [std(error_surfarea,1) ;std(error_fvec,1); std(error_surf,1)]



% print out six ,metrics for simplicity
printmat(results_mean, 'mean value', 'Surface_area_distortion nRMSE_in_SH_coef?cients nRMSE_in_coordinates', '40 80 160' )