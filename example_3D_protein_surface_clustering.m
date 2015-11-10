clear 
%% produce table 5 in the 3D case
addpath(genpath('src\'))
addpath(genpath('data\'))

load fvec_3proteins
% load fvec_smooth
load sphere_standard

% reshape fvec output from SHARM
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
[yr_2d U S] = shape_component_analysis(y0, 40);


display('processing nonlinear meanshift')
tic
yc_2d_full = mean_shift_NDshape_full_procedure(y0, y0, 0.1, 100);
distmap = pdist2_NDriemann(y0, yc_2d_full);
[tmp ind_full] = min(distmap,[],2);
% yc_2d = mean_shift_NDshape_adaptive_full_procedure(yr_2d, yr_2d, 2, 100, U);
toc

display('processing nonlinear meanshift with 2DPCA')
tic
yc_2d = mean_shift_NDshape_full_procedure(yr_2d, yr_2d, 0.1, 100, U);
% yc_2d = mean_shift_NDshape_adaptive_full_procedure(yr_2d, yr_2d, 2, 100, U);
distmap = pdist2_NDriemann(y0, yc_2d);
[tmp ind] = min(distmap,[],2);

toc




% rotation invariant SH features

for ii  = 1:31
    tmp = fvec_all(((ii-1)^2+1):(ii)^2,:,:);
    SHf(ii,:) = sqrt(sum(abs(reshape(tmp,[],60)).^2,1));
end

SHfn = SHf./repmat(sqrt(sum(SHf.^2)),31,1);

display('processing nonlinear meanshift with rotation invariant SH features')
tic
[yc, label_sh] = mean_shift_tangent_full_procedure(SHfn, .0145, 100);
toc

% rotation invariant SH features + PCA

display('processing nonlinear meanshift with rotation invariant SH features + PCA')
[coefs,scores,variances,t2] = princomp(SHfn');
tic
[yc_shpca, label_shpca] = mean_shift_tangent_full_procedure(scores(:,1:6)', .0145, 100);
toc
% Laplacian Eigen map

display('processing nonlinear meanshift + Laplacian Eigen map')
distmap = pdist2_NDriemann(y0);
tic
[V D1] = LapEigMap(distmap,5,1);

Lapf = V(:,2:11)';
[yc_lap, label_lap] = mean_shift_tangent_full_procedure(Lapf, .015, 100);
toc

modes_gt = ind2modes(ceil([1:60]/20));

modes_RMS = ind2modes(ind_full);
CM0 = contingency_table(modes_RMS, modes_gt);
vec0 = clustering_performance(CM0)

modes_RMSSFA = ind2modes(ind);
CM1 = contingency_table(modes_RMSSFA, modes_gt);
vec1 = clustering_performance(CM1);


modes_MSSH = ind2modes(label_sh);
CM2 = contingency_table(modes_MSSH, modes_gt);
vec2 = clustering_performance(CM2);

modes_MSSHPCA = ind2modes(label_shpca);
CM3 = contingency_table(modes_MSSHPCA, modes_gt);
vec3 = clustering_performance(CM3);

modes_MSLAP = ind2modes(label_lap);
CM4 = contingency_table(modes_MSLAP, modes_gt);
vec4 = clustering_performance(CM4);




results = [vec0; vec1; vec2; vec3; vec4;];

% print out six ,metrics for simplicity
printmat(results(:,1:6), 'Clustering metric', 'RMS RMSSFA MSSH MSSHPCA MSLAP', 'purity MI NMI AR MI2 v-measures' )