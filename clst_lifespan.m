function [rev] = clst_lifespan(basis_tgt, pdiag)

% Description: Cluster Life Span
%   Given a persistence diagram and a data, find the life span of that clst.
%
%
% Input:
%   basis_tgt: the cluster we are interested in. Should be a scalar or column vector.
%   pdiag: the persistence diagram
%
% Output:
%   cM: (matrix)
%   sM: (matrix)

%% Find clst_tgt and the first level (densest) that the clst and appear

% Find the level that the data appear
ih_start = find(basis_tgt<=pdiag.size_levelset,1,'last');

cM = zeros(pdiag.n_h,pdiag.n_r);
sM = zeros(pdiag.n_h,pdiag.n_r);

if ih_start > pdiag.n_h
    disp('Invalid input');
    return
end

%% Check it from the lowest level

i_next = 1;
rep_clst = {};
ir_lbound = 1;
ir_rbound = pdiag.n_r;

for ih = ih_start:-1:1
    
    %-- Get the info of this level -------------------------------------------------%
    rank_cur = pdiag.rank_level(ih);
    size_ls  = pdiag.size_levelset(ih);
    flag_found = false;
    
    %-- Move toward right -------------------------------------------------------------------------%
    for ir = ir_lbound:ir_rbound
        % find the clst that contain the tgt basis
        sptt = pdiag.subptts{pdiag.spM(ih,ir)};
        idx_clst = find(sptt(:,basis_tgt),1,'first');
        
        fbasis_cur_clst = find(sptt(idx_clst,1:rank_cur),1,'first');
        if fbasis_cur_clst < basis_tgt  % get merged
            break
        end
        
        % We know it is not merged by other clusters
        rep_tgt = sptt(idx_clst,1:rank_cur);
        sM(ih,ir) = sum(rep_tgt);
        
        if ir>1 && sM(ih,ir)==sM(ih,ir-1)
            cM(ih,ir)=cM(ih,ir-1);
        else
            flag_new = true;
            
            % check below to see whether the clst has been found before
            for ih_l = ih+1:pdiag.n_h
                idx_r_test = find(sM(ih_l,:)==sM(ih,ir),1,'first');
                if isempty(idx_r_test)
                    continue
                end
                rep_v_tgt = adjust_rank(rep_tgt, pdiag.rank_level(ih_l));
                rep_v_test = adjust_rank(rep_clst{cM(ih_l,idx_r_test)}, pdiag.rank_level(ih_l));
                if ~any(xor(rep_v_tgt,rep_v_test)) % same
                    flag_new = false;
                    cM(ih,ir) = cM(ih_l,idx_r_test);
                    rep_clst{cM(ih,ir)} = rep_tgt;  % need to update to get the right rank
                    break;
                end
            end
            
            if flag_new
                rep_clst = cat(1,rep_clst, rep_tgt);
                cM(ih,ir) = i_next;
                i_next = i_next + 1;
            end
        end                   
        
    end % end of ir loop
    
    ir_lbound = find(sM(ih,:),1,'first');
    ir_rbound = find(sM(ih,:),1,'last');
    
end % end of ih loop

%% data to return
rev.cM           = cM;
rev.sM           = sM;
rev.rep_clst = rep_clst';
