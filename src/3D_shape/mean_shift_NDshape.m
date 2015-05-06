function [xt distmap] = mean_shift_NDshape(zi, x0, sigma, iter_max)

dim = size(x0,1);
l_num = size(x0,2);
s_numx = size(x0,3);
s_numz = size(zi,3);

xt = zeros(size(x0));

 
for ii = 1:s_numx

    mu = x0(:,:,ii);
    s0 = 1000;
    tt = 0;
    while (tt < iter_max  && s0>1e-6) 
        mhx = zeros(dim, l_num);
        total_weight = 0;
        for jj = 1:s_numz

            [U S V] = svd(zi(:,:,jj)*mu');
            Rj = det(V)*det(U)*V*U';
            rhoj = acos(sum(abs(diag(S))));
            Tj = rhoj/sin(rhoj+eps)*(Rj*zi(:,:,jj) - mu * cos(rhoj));

            W0 = exp(-rhoj.^2/2/sigma^2);
            total_weight = total_weight + W0;
            mhx = mhx + Tj*W0;
        end
        mhx = mhx/(total_weight+eps);       
        s0 = norm(mhx,'fro');
        mu = sin(s0)*mhx/(s0+eps) + mu*cos(s0);
        mu = mu/norm(mu,'fro');
        tt = tt + 1;
    end
    xt(:,:,ii) = mu;

end


 distmap = pdist2_NDriemann(xt);