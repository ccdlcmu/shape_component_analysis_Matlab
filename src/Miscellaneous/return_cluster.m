function [modes] = return_cluster(distmap,TOL)

if nargin<2
    TOL = 1e-3;
end

s_num = size(distmap,1);
modes = cell(1,50); 
set1 = find(distmap(1,:)<TOL);
modes{1} = set1;

set2 = setdiff(1:s_num,set1); 
ii = 2;
while (~isempty(set2))  
    modes{ii} = find(distmap(set2(1),:)<TOL);
    set1 = union(set1,modes{ii});
    ii = ii + 1;
    set2 = setdiff(1:s_num,set1);
end
modes(ii:end) = [];
 


% for ii = 1:min(35, length(modes))
% subplot(7,5,ii)
% plot(H'*yt(:,modes{ii}(1)))
% axis equal
% set(gca,'XTick',[],'YTick',[])
% end
% figure
% for ii = 1:10
% subplot(4,5,ii)
% plot(H'*y0(:,1 + 20*(ii-1)))
% axis equal
% set(gca,'XTick',[],'YTick',[])
% end
% [distmap2(:); diag(distmap2)]  