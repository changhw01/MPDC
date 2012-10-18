
%> @file Idxs4LargeClst.m
%> @brief Find the indices for the larger clusters
%=========================================================================
%> @brief Find the indices for the larger clusters
%>
%> @param c_size The sizes of the clusters
%> @param propThres The proportion threshold
%>
%> @retval idxs_large the idxs for the larger clst
%=========================================================================

function [idxs_large] = Idxs4LargeClst(c_size, propThres)

n_clst = length(c_size);

[c_size_sort idxs_sort] = sort(c_size,'descend');
sizeProp = c_size_sort / sum(c_size_sort);

acc_prop = 0;
up2idx = 0;
for i = 1:n_clst
    if acc_prop > propThres
        up2idx = i;
        break
    else
        acc_prop = acc_prop + sizeProp(i);
    end
end

if up2idx == 0
    up2idx = n_clst;
end

idxs_large = idxs_sort(1:up2idx);
