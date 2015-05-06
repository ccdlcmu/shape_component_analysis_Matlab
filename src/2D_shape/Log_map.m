function mulog = Log_map(mu,y0)
dim_sub = size(y0,1);
y0 = y0./ sqrt(repmat(diag((y0')*y0)',dim_sub,1));
theta0 = -angle(mu'*y0);
z1 = repmat(exp(1i*theta0),dim_sub,1).*y0;
rho = abs(acos(abs(mu'*z1)));
mulog = repmat(rho./(sin(rho)+eps),dim_sub,1).*(z1-mu*cos(rho));
