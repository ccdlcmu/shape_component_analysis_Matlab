function [sp_cord v_proj sp] = score_to_cord(score, V, mu)

v_proj = V*score;
% v_proj = v_proj(1:dim_sub/2,:) + v_proj(dim_sub/2+[1:dim_sub/2],:)*1i;


s0 = sqrt(sum(v_proj.*conj(v_proj)));
sp = repmat(sin(s0)./(s0+eps),size(V,1),1).*v_proj +  [real(mu); imag(mu)]*cos(s0);
sp_cord = [cos(s0);repmat(sin(s0)./(s0+eps),size(V,2),1).*score];