function [vec metric_name] = clustering_performance(CM)

purity = calculate_purity(CM);
MI = mutual_information(CM);
NMI = normalized_mutual_information(CM);
[AR,RI,MI2,HI2]=RandIndex(CM);


[v,h, c, hc,hk,h_ck,h_kc] = calculate_v_measure (CM);
vec = [purity MI NMI AR,MI2, v, h, c, hc,hk,h_ck,h_kc];
metric_name = {'purity', 'MI', 'NMI', 'AR', 'MI2','v-measures'};