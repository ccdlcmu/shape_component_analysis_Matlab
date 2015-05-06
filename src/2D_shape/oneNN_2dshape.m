function [modes ind] = oneNN_2dshape(y0,yc)

distmap = abs(acos(abs(y0'*yc)));
[temp ind] = min(distmap,[],2);
modes = ind2modes(ind);