function asg = assign_target_clst(I_tgtbasis, distM, idx_last)

%% Note

% Assuming the data are sorted according to the density (1 as the most dense), we will assign
% the data to the I_tgtbasis one by one


%% Main routine

n_data = size(distM,1);
asg = zeros(n_data,1);

if nargin < 3
	idx_last = n_data;
end

for i = 1:idx_last
    if ismember(i,I_tgtbasis)
        asg(i) = i;
    else
        % assign according to the nearest assigned neighbor
        [~,idx] = min(distM(i,1:i-1));
        asg(i) = asg(idx);
    end
end
