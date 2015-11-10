clear
close all
load MPEG7_05
y0 = MPEG7_05;

dim_projs = [10, 20, 40];

for dd = 1:length(dim_projs)

    dim_proj = dim_projs(dd);
    s_num = size(y0,2);
    dim = size(y0,1);

    % dim_sub = 6;

    H = -gallery('orthog',dim,4);
    H(1,:) = [];


    y0 = H'*H*y0;
    for ii = 1:s_num 
        y0(:,ii) = y0(:,ii)/norm(y0(:,ii));
    end

    % principal geodesic analysis
    mu = mean_shape_riemann(y0,10);

    mulog = Log_map(mu,y0);
    mulogR = [real(mulog); imag(mulog)];
    [U S V] = svds(mulogR', dim_proj);
    tgvector_rec = (U*S*V')';



    tic


    ym= zeros(2, dim, s_num);
    for ii = 1:s_num 
        ym(:,:,ii) = [real(y0(:,ii)) imag(y0(:,ii))]';
    end

    % shape component analysis
    [yt U S] = shape_component_analysis(ym, ceil(dim_proj/2));

    toc

    yc1 = zeros(dim,s_num);
    for ii = 1:s_num
        ytt = yt(:,:,ii)*U';
    %     yc1(:,:,ii) = yc1(:,:,ii)/norm(yc1(:,:,ii),'fro');
        yc1(:,ii) = ytt(1,:)' + 1i* ytt(2,:)';
        yc1(:,ii) = yc1(:,ii)/norm(yc1(:,ii));
    end


    dist_SCA(dd) = sqrt(mean(diag(abs(acos(abs(y0'*yc1))).^2)))
    dist_PGA(dd) = sqrt(mean(sum((tgvector_rec - mulogR).^2)))

    weight = norm(abs(acos(abs(y0'*y0))),'fro');
    distortion_SCA(dd) =  norm(abs(acos(abs(yc1'*yc1))) - abs(acos(abs(y0'*y0))),'fro')/weight;
    distortion_PGA(dd) = norm(pdist2(tgvector_rec',tgvector_rec') - abs(acos(abs(y0'*y0))),'fro')/weight;
end

