function mu = mean_shape_riemann(y0,iter)
% calculate mean shape mu from input data y0.
% y0: mxn complex matrix. n is the number of repetition. m is the number of landmarks.  
if nargin<2
    iter = 200;
end

mu = mean_shape_procrustes(y0');
mu = conj(mu);

s_num = size(y0,2);
dim = size(y0,1);
dim_sub = dim-1;

H = -gallery('orthog',dim,4);
H(1,:) = [];

y0 = H*y0;
mu = H*mu;

y0 = y0./ sqrt(repmat(diag((y0')*y0)',dim_sub,1));

OF = zeros(1,iter);
tic
for tt = 1:iter 
    mu = mu./ sqrt(repmat(diag((mu')*mu)',dim_sub,1));
    for ii = 1:s_num
        theta0 = -angle(mu'*y0);
        z1 = repmat(exp(1i*theta0),dim_sub,1).*y0;
        rho = abs(acos(abs(mu'*z1)));        
        
        mulog = repmat(rho./(sin(rho)+eps),dim_sub,1).*(z1-mu*cos(rho));
        mhy = mean(mulog,2);
        s0 = norm(mhy);
        mu = sin(s0)*mhy/(s0+eps) + mu*cos(s0);
        
    end
    OF(tt) = sum(rho);
%     yt = yt./ sqrt(repmat(diag((yt')*yt)',dim_sub,1));
%     corrmat = yt'*yt;
%     distmap2 =  triu(abs(sqrt(1-corrmat.*conj(corrmat))),1);
%     [nn xx] = hist(distmap2(distmap2>0),0:0.01:1);
%     pp = nn / sum(nn);
%     entro(tt) = -sum(pp.*(log(pp+eps)));
end
plot(OF)
mu = H'*mu;
toc