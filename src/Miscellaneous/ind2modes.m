function modes = ind2modes(cluster_label)

for ii = 1:max(cluster_label)
    modes{ii} = find(cluster_label == ii);
end