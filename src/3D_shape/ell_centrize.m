function xxc = ell_centrize(xx)
% calibrate rotation by two dimensional ellipse fitting
% aa = [a1 a2 a3 a4 a5 a6]', coefficient of ellipse,
% xx(k,1) = x(k) + i * y(k), 
% ellipse = a1 x^2 + a2 xy +a3 y^2 +a4 x + a5 y + a6  
D = [real(xx).^2 real(xx).*imag(xx) imag(xx).^2 real(xx) imag(xx) ones(size(xx,1),1)];
A = zeros(6,6);
A(3,1) = 2;
A(1,3) = 2;
A(2,2) = -1;
[U D1] = eig(D'*D,A);
lambda = diag(U'*A*U);

ind = (lambda>1e-4);
aa_temp = U(:,ind).*diag(1/sqrt(lambda(ind)));
[temp ind2] = min(diag(aa_temp'*D'*D*aa_temp));
aa = aa_temp(:,ind2);


ELL = sign(aa(1))*[aa(1) aa(2)/2 aa(4)/2;aa(2)/2 aa(3) aa(5)/2;aa(4)/2 aa(5)/2 aa(6)];
% [V D2] = eig(ELL);
% lambda = diag(D2);
% ind2 = find(sign(lambda)~=prod(sign(lambda)));

% V = -sign(lambda(ind2(1)))*V;
% [temp ind4] = sort(abs(lambda(ind2)),'descend');
% ori = det(V(1:2,[ind2(ind4)]));
% V(:,ind2(ind4(1))) = ori *V(:,ind2(ind4(1)));
% 
% xxc_temp = V(:,[ind2(ind4)])'*[real(xx) imag(xx) ones(length(xx),1)].';
% xxc = xxc_temp(1,:)' + i*xxc_temp(2,:)';

R1 = ELL(1:2,1:2);
t1 = ELL(1:2,3);
xc = inv(R1)*t1;

[U1 S1 V1] = svd(R1);
V1 = sign(V1(:,1)'*[1;0])*V1;
V1(:,2) = sign(det(V1)) * V1(:,2);

xxc_temp = V1'*[real(xx)+xc(1) imag(xx)+xc(2)].';
xxc = xxc_temp(1,:)' + i*xxc_temp(2,:)';