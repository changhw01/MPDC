function rev = pclustering_fullbasis(distM, f, thres_f, thres_dist, order_f)

% Description: Persistence Clustering
%   Filter the graph nodes using the filter functions and then partition the data using the
%   distance thresholds.
%
%   * The main difference from 'pclustering' is that this one will take each data point as a
%   seperate basis
%
%
% Note: 1. We assume the data are ordered according to the filter function (in the descending order)
%       2. The main trick in this implementation is we need to adjust the bases and representations
%          to the correct dimension using size_levelset
%
%
% Input:
%   distM: The distance matrix of the whole data
%   f: The filter function
%   thres_f: The filter thresholds (in increasing order)
%   thres_d: The distance thresholds (in increasing order)
%   order_f: Indicate whether low of high f values are important
%
% Output: (See the end of this file)
%   subptts: (cell) sub partitions of the data
%   rank_level: (array)
%   spM: (matrix)
%   sM: (matrix)
%   size_levelset: (array)

%% parameters and data structures

n_data      = length(f);
n_h         = length(thres_f);
n_r         = length(thres_dist);

subptts     = {};
rank_level  = zeros(n_h,1);
spM         = zeros(n_h,n_r);
sM          = zeros(n_h,n_r);

iptt_next   = 1;

% default is using density, so higher value is better
if nargin < 5
    order_f = 'high';
end

%% Count how many data within threshold
% make sure the threshold are in right order
f = reshape(f,n_data,1);

if strcmp(order_f,'low')
    thres_f = sort(thres_f,'descend');
    thres_f = reshape(thres_f,n_h,1);
    size_levelset = sum(repmat(f',n_h,1)<=repmat(thres_f,1,n_data),2);
else
    thres_f = sort(thres_f);
    thres_f = reshape(thres_f,n_h,1);
    size_levelset = sum(repmat(f',n_h,1)>=repmat(thres_f,1,n_data),2);
end

%% Main Loop

size_ls_pre = 0;

for ih = n_h:-1:1
    
    %% ---- If no new data appear ---------------------------------------------------------------- %
    if size_ls_pre == size_levelset(ih)
        spM(ih,:) = spM(ih+1,:);
        sM(ih,:) = sM(ih+1,:);
        rank_level(ih) = rank_level(ih+1);
        continue
    end
    
    %% ---- Process ------------------------------------------------------------------------------ %
    size_ls_cur = size_levelset(ih);
    
    for ir = 1:n_r
        fprintf('ih=%d, ir=%d\n',ih,ir);
        
        if ir == 1  
            %-- Compute new bases
            I = 1:size_ls_cur;
            [n_clst_cur, C_cur] = graphconncomp(sparse(distM(I,I)<=thres_dist(1)),'Weak',true);
            sM(ih,ir) = n_clst_cur;
            
            rank_cur = size_ls_cur;
            rank_level(ih) = rank_cur;
            
            %-- compute the representation
            subptt_cur = false(n_clst_cur, rank_cur);
            for i = 1:n_clst_cur
                subptt_cur(i,:) = (C_cur==i);
            end
            
        else
            % Compute merging graph
            mG = false(n_clst_pre);
            for i = 1:n_clst_pre
                for j = i+1:n_clst_pre
                    if any(any(distM(subptt_pre(i,:), subptt_pre(j,:)) <= thres_dist(ir)))
                        mG(i,j) = true;
                    end
                end
            end            
            
            %-- Clusters in the merging graph
            if ~any(mG(:)) % no merge happens
                spM(ih,ir) = spM(ih,ir-1);
                sM(ih,ir) = n_clst_pre;
                continue  % continute to next ir
            else
                [n_clst_cur, C_cur] = graphconncomp(sparse(mG),'Weak',true);
                sM(ih,ir) = n_clst_cur;
                subptt_cur = false(n_clst_cur,rank_cur);
                for i = 1:n_clst_cur
                    subptt_cur(i,:) = any(subptt_pre((C_cur==i),:),1);
                end
            end
                
        end
        
        %-- Check with the subptt below
        % need to check all sb of the same size unless we have found the equivlent one
        ih_to_test = ih+1;
        while ih_to_test <= n_h
            ir_to_test = find(sM(ih_to_test,:)==n_clst_cur,1,'first');
            if ~isempty(ir_to_test)
                subptt_cur_adj = adjust_rank(subptt_cur,rank_level(ih_to_test));
                subptt_pre_adj = adjust_rank(subptts{spM(ih_to_test,ir_to_test)},rank_level(ih_to_test));
                if ~any(any(xor(subptt_cur_adj,subptt_pre_adj)))
                    spM(ih,ir) = spM(ih_to_test,ir_to_test);
                    % need to update the previous subptt
                    subptts{spM(ih_to_test,ir_to_test)} = subptt_cur;
                    break;
                end
            end
            ih_to_test = ih_to_test + 1;
        end
        
        % if still not assigned, must be a new one
        if spM(ih,ir) == 0
            subptts     = cat(1,subptts,subptt_cur);
            spM(ih,ir)   = iptt_next;
            iptt_next   = iptt_next+1;
        end        
        
        %-- Update
        n_clst_pre = n_clst_cur;
        subptt_pre = subptt_cur;        
    end  % End of the ir loop
    
    %% ---- Update -----------------------------------------------------------%
    size_ls_pre = size_ls_cur;
    
end % End of ih loop

%% Run though the subptts to get the life time of each basis
blife = zeros(n_data,1);
for ih = n_h:-1:1
    for ir = 1:n_r
        sptt = subptts{spM(ih,ir)};
        for i = 1:sM(ih,ir)
            idx = find(sptt(i,:),1,'first');
            blife(idx) = blife(idx) + 1;
        end
    end
end


%% Data to return
rev.n_data      = n_data;
rev.n_h         = n_h;
rev.n_r         = n_r;
rev.subptts     = subptts;
rev.rank_level  = rank_level;
rev.spM         = spM;
rev.sM          = sM;
rev.blife       = blife;
rev.size_levelset = size_levelset;
