function purity = calculate_purity(CM)
% computing purity from a contingency_table.
% assuming CM(i,j) represents contingency of an estimated cluster i and
% groundtruth cluster j.
purity = sum(max(CM,[],2))/sum(sum(CM));
