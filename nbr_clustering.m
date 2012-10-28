function [C_cell, n_clst] = nbr_clustering(distM, s, I_want)

% Compute the neighborhood clustering using a scale parameter
% This is also the Rips complex

if nargin < 3
    I_want = 1:length(distM);
end

[n_clst, C] = graphconncomp(sparse(distM(I_want,I_want)<=s),'Weak',true);

C_cell = cell(n_clst,1);
for i = 1:n_clst
    C_cell{i} = I_want((C==i));
end