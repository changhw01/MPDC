function rep_new = adjust_rank(rep, r)

r_pre = size(rep,2);
if r_pre == r;
    rep_new = rep;
elseif r_pre > r
    rep_new = rep(:,1:r);
else
    rep_new = padarray(rep,[0,r-r_pre],false,'post');    
end
