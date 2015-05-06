clear 
%% reproduce figure 5, right

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

[yr_2d05 U S] = shape_factor_analysis(y0, 10);
yc_2d05 = mean_shift_NDshape_full_procedure(yr_2d05, yr_2d05, 0.1, 100, U);

[yr_2d10 U S] = shape_factor_analysis(y0, 20);
yc_2d10 = mean_shift_NDshape_full_procedure(yr_2d10, yr_2d10, 0.1, 100, U);


[yr_2d20 U S] = shape_factor_analysis(y0, 40);
yc_2d20 = mean_shift_NDshape_full_procedure(yr_2d20, yr_2d20, 0.1, 100, U);


[vs fs] = sphereMesh([0 0 0 1]);
Zs = calculate_SPHARM_basis(vs, 31);

close all
ha = tight_subplot(3, 5, [.01 .01],[.1 .01],[.01 .01]);


for jj = 1:15
    axes(ha((jj)))
    axis off

end

for ii = 1:3
    ym = yc_2d05(:,:,ii);
    y0rec = shape_to_SH(ym);
    Zvert = real(Zs*y0rec');
    axes(ha(ii))
    patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    axis off
end

for ii = 1:4
    ym = yc_2d10(:,:,ii);
    y0rec = shape_to_SH(ym);
    Zvert = real(Zs*y0rec');
    axes(ha(5 + ii))
    patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    axis off
end

for ii = 1:5
    ym = yc_2d20(:,:,ii);
    y0rec = shape_to_SH(ym);
    Zvert = real(Zs*y0rec');
    axes(ha(10 + ii))
    patch('vertices', Zvert, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    axis off
end

