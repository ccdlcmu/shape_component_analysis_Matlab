clear
%% reproduce table 5 in the 2D case

close all

addpath(genpath('src\'))
addpath(genpath('data\'))
load MPEG7_05
y0 = MPEG7_05;

f_dim = 30;

% set up ground truth labels
group_number = 5;
member_number = 20;

for ii = 1:group_number
    ground_truth{ii} = (ii-1)*member_number + [1:member_number];
end

s_num = size(y0,2);
dim = size(y0,1);


% normalize scaling and translation
H = -gallery('orthog',dim,4);
H(1,:) = [];


y0 = H'*H*y0;
for ii = 1:s_num
    y0(:,ii) = y0(:,ii)/norm(y0(:,ii));
end

% principal geodesic analysis and corresponding projections
[mu V score sp sp_cord ev dim_proj] = shape_PCA_overR(y0, .99);


% nonlinear meanshift without dimension reduction
display('processing nonlinear meanshift without dimension reduction')
tic
[yc] = mean_shift_2Dshape_full_procedure(y0, 0.12, 200);
toc
[modes_RMS ind] = oneNN_2dshape(y0,yc);

% nonlinear meanshift with dimension reduction
display('processing nonlinear meanshift with dimension reduction')
[yc2] = mean_shift_2Dshape_full_procedure(sp_cord, 0.12, 200, V(:,1:dim_proj), mu);
[modes_RMSSFA ind2] = oneNN_2dshape(y0,H'*yc2);



% linear mean-shift with elliptic fourier descriptors

FSD_oriented = zeros(4, f_dim, s_num);
for ii = 1:s_num
     rFSDs = fEfourier([real(y0(:,ii)) imag(y0(:,ii))], f_dim, 1, 1);
     FSD_oriented(:,:,ii) = rFSDs;
end
FSD_oriented(:,1,:) = 0;

FSD_o_vec = reshape(FSD_oriented(:,2:end,:),[],s_num);

display('processing linear mean-shift with elliptic fourier descriptors')
tic
[yc_fs ind_fs] = mean_shift_tangent_full_procedure(FSD_o_vec, 0.2, 100);
toc
modes_MSFD = ind2modes(ind_fs);

% linear mean-shift on the tangent plane
display('linear mean-shift on the tangent plane')
tic
[yc_tg ind_tg] = mean_shift_tangent_full_procedure(score, 0.085, 100);
toc
modes_tMS = ind2modes(ind_tg);


% laplacian Eigen map
display('linear mean-shift + laplacian eigenmap')
distmap = abs(acos(abs(y0'*y0)));
[V D1] = LapEigMap(distmap,5,1);

Lapf = V(:,2:7)';
tic
[yc_lap, label_lap] = mean_shift_tangent_full_procedure(Lapf, .0075, 100);
toc
modes_MSLAP = ind2modes(label_lap);


% output performance 
CM1 = contingency_table(modes_RMS, ground_truth);
CM2 = contingency_table(modes_RMSSFA, ground_truth);
CM3 = contingency_table(modes_tMS, ground_truth);
CM4 = contingency_table(modes_MSFD, ground_truth);
CM5 = contingency_table(modes_MSLAP, ground_truth);

vec1 = clustering_performance(CM1);
vec2 = clustering_performance(CM2);
vec3 = clustering_performance(CM3);
vec4 = clustering_performance(CM4);
vec5 = clustering_performance(CM5);

results = [vec1; vec2; vec3; vec4; vec5;];

% print out six ,metrics for simplicity
printmat(results(:,1:6), 'Clustering metric', 'RMS RMSSFA tMS MSFD MSLAP', 'purity MI NMI AR MI2 v-measures' )
% plot clustering results


