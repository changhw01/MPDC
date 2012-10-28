function asg_all = construct_ptt_whole_3(I_tgtbasis, distM)

%% Note

% Assuming the data are sorted according to the density (1 as the most dense), we will assign
% the data to the I_tgtbasis one by one


%% Main routine

n_data = size(distM,1);
asg_all = zeros(n_data,1);

for i = 1:n_data
    if ismember(i,I_tgtbasis)
        asg_all(i) = i;
    else
        % assign according to the nearest assigned neighbor
        [~,idx] = min(distM(i,1:i-1));
        asg_all(i) = asg_all(idx);
    end
end
