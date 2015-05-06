function y0 = pull_back(mu,V, sp_cord)

y0 = mu*sp_cord(1,:) + V*sp_cord(2:end,:);