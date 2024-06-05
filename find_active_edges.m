function [edg_ind] = find_active_edges(c_aft,tol)
% Input: c_aft matrix to analyze which edges are active (i.e. > tol)
% Output: Active edge indices

c_aft = c_aft';
[N1,nt] = size(c_aft);

c_aft_norm = zeros(N1,1);

for i = 1:N1
   c_aft_norm(i) = norm(c_aft(i,:))/nt;
end

edg_ind = find(c_aft_norm >= tol);

end

