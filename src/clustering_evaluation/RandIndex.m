function [AR,RI,MI,HI]=RandIndex(CM)
% CM contingency matrix

n=sum(sum(CM));
nis=sum(sum(CM,2).^2);		%sum of squares of sums of rows
njs=sum(sum(CM,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(CM.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)
