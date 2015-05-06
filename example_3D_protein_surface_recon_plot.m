clear 
%% reproduce figure 5, left

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

[yr_2d05 U S] = shape_factor_analysis(y0, 80);

[vs fs] = sphereMesh([0 0 0 1]);
Zs = calculate_SPHARM_basis(vs, 31);

dimr =  [40 20 10];

close all
ha = tight_subplot(3, 4, [.01 .01],[.1 .01],[.01 .01]);
for jj = 1:3
    ym = y0(:,:,(20*jj - 15) + 1);
    y0rec = shape_to_SH(ym);
    Zvert = real(Zs*y0rec');
    axes(ha((jj-1)*4 + 1))
    patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    axis off
    for ii = 2:4
        
        Ut = U(:,1:dimr(ii-1));
        yt = ym*Ut;
        yt = yt/norm(yt,'fro');
        ytt = yt*Ut';
        y0rec = shape_to_SH(ytt);
        Zvert = real(Zs*y0rec');
        axes(ha((jj-1)*4 + ii))
        patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
        axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
        axis off 
    end
end
