function [yt distmap] = mean_shift_2Dshape(y0,sigma,iter)


s_num = size(y0,2);
dim = size(y0,1);

yt = y0;


for ii = 1:s_num
    s0 = 1000;
    mu = yt(:,ii);
    tt = 0;
    while tt  < iter && s0 > 1e-6 
        
        theta0 = -angle(mu'*y0);
        z1 = repmat(exp(1i*theta0),dim,1).*y0;
        rho = abs(acos(abs(mu'*z1)));

        W0 = exp(-rho.^2/2/sigma^2);
        W0 = W0 ./sum(W0);       
        mulog = repmat(rho./(sin(rho)+eps),dim,1).*(z1-mu*cos(rho));
        
        mhy = mulog*W0';
        s0 = max(abs(mhy));
        mu = sin(s0)*mhy/(s0+eps) + mu*cos(s0);
        mu = mu/norm(mu);
        
        tt = tt + 1;
    end

    yt(:,ii) = mu;

end
distmap = abs(acos(abs(yt'*yt)));