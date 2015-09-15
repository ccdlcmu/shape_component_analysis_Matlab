clear
close all
load MPEG7_05
y0 = MPEG7_05;

dim_projs = [40, 20, 10];
examples = [1, 21, 41, 61, 81];



s_num = size(y0,2);
dim = size(y0,1);

H = -gallery('orthog',dim,4);
H(1,:) = [];


y0 = H'*H*y0;
for ii = 1:s_num 
    y0(:,ii) = y0(:,ii)/norm(y0(:,ii));
end

for ii = 1:length(examples)
    subplot(length(examples), 4, 1+4*(ii-1))
    plot(y0(:,examples(ii)),'r','linewidth',2)
    axis tight
    axis square
end

for dd = 1:length(dim_projs)

    dim_proj = dim_projs(dd);

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
    [yt U S] = shape_factor_analysis(ym, dim_proj);

    toc

    
    for ii = 1:length(examples)
        ytt = yt(:,:,examples(ii))*U';
        yc1 = ytt(1,:)' + 1i* ytt(2,:)';
        yc1 = yc1/norm(yc1);
        subplot(length(examples), 4, 1+dd+4*(ii-1))
        plot(yc1)
        axis tight
        axis square
    end
    
   
end

