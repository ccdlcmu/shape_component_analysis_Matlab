function [yn] = shape_normalize(y0)

dim = size(y0,1);
dim_sub = dim-1;

H = -gallery('orthog',dim,4);
H(1,:) = [];

yn = H*y0;
yn = yn./ sqrt(repmat(diag((yn')*yn)',dim_sub,1));