clear
close all

addpath(genpath('src\'))
addpath(genpath('data\'))

load meanshape_moving_average
y0 = shape;


s_num = size(y0,2);
dim = size(y0,1);

H = -gallery('orthog',dim,4);
H(1,:) = [];

y0 = H'*H*y0;
for ii = 1:s_num
    y0(:,ii) = y0(:,ii)/norm(y0(:,ii));
end
y02 = H*y0;


[mu V score sp sp_cord ev dim_proj] = shape_PCA_overR(y0, .99);


tic
[yt distmap] =  mean_shift_2Dshape(y02,0.035,200);
toc
tic
[yt2 distmap2] = mean_shift_2Dshape(sp_cord,0.035,200);
toc

[modes] = return_cluster(distmap);
[modes2] = return_cluster(distmap2);

ind = zeros(1,length(modes));
for ii = 1:length(modes)
    ind(ii) = modes{ii}(1);
end

ind2 = zeros(1,length(modes2));
for ii = 1:length(modes2)
    ind2(ii) = modes2{ii}(1);
end

yy2 = pull_back([real(mu);imag(mu)], V(:,1:dim_proj), yt2);
yy2 = yy2(1:(dim-1),:) + yy2(dim:end,:)*1i;

group_number = 10;
member_number = 20;

y0 = H'*H*y0;
for ii = 1:s_num
    y0(:,ii) = y0(:,ii)/norm(y0(:,ii));
end

y02 = H*y0;



for ii = 1:16
    for jj = 1:8
        plot(4*ell_centrize(H'*y02(:,(ii-1)*80 + jj*10)) + 0.2 + ii + 1i*jj, 'r', 'linewidth',2)
        hold on
    end
end
axis equal
axis([0 17 0 9.])
set(gca,'XTick',[])
set(gca,'YTick',[])

figure
for ii = 1:length(ind)
    subplot(2,22,(ii-1)*2 + [1 2])
    plot(ell_centrize(H'*yy(:,ind(ii))),'linewidth',2)
    axis equal
    axis([-0.12 0.12 -0.12 0.12])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

for ii = 1:length(ind2)
    subplot(2,22,22 + (ii-1)*2 + [1 2])
    plot(ell_centrize(H'*yy2(:,ind2(ii))),'k','linewidth',2)
    axis equal
    axis([-0.12 0.12 -0.12 0.12])
    
    set(gca,'XTick',[])
    set(gca,'YTick',[]) 
end

