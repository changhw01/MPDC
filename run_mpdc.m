clear

%% import data
data = load('D_cos_doc');
%pD_raw = data.D_cos_doc;
pD = data.D_cos_doc;

clear data

% pD = exp(-pD_raw.*pD_raw/0.7);
% pD = 1-pD./max(pD(:));

n_data = length(pD);
%% density estimate

G = pD;

n_NN = 10;
h = 0.5;

dsty_est = zeros(n_data,1);
h_sqr = h*h;

for i = 1:n_data
    [v_sorted, idxs_sorting] = sort(nonzeros(G(i,:)));

    if length(v_sorted) < n_NN+1
        dsty_est(i) = sum(exp(-v_sorted(2:end)/h_sqr))/n_NN;
    else
        dsty_est(i) = sum(exp(-v_sorted(2:n_NN+1)/h_sqr))/n_NN;
    end
end

% NOTE: need to sort the data according to density first 
I_large = Idxs4LargeClst(dsty_est, 1);
dsty_est_large = dsty_est(I_large);
distM = pD(I_large,I_large);
%figure;showMatAsImg(distM);
figure; plot(dsty_est_large,'.');

[~,I_large_back] = sort(I_large);

%% mpca

n_thres_dsty = 20;
n_thres_dist = 40; 

n_points = 500;

% thresholds: 2~100
min_dsty = dsty_est_large(n_points); max_dsty = dsty_est_large(5);
min_dist = 0.2; max_dist = 0.7;

thres_dsty = linspace(max_dsty, min_dsty,n_thres_dsty);
thres_dist = linspace(min_dist,max_dist,n_thres_dist);

% persistence clustering
pdiag = pclustering_fullbasis(distM, dsty_est_large, thres_dsty, thres_dist);
%draw_clMat(pdiag.spM,pdiag.sM,1,'jet');


[~ ,I_blife] = sort(pdiag.blife,'descend');

%% Check some diagrams

cls_top = cell(16,1);
figure;
for i = 1:16    
    cls = clst_lifespan(I_blife(i), pdiag);
    
    subplot(4,4,i);
    draw_clMat(cls.cM,cls.sM,false);
    title(sprintf('basis %d',I_blife(i)));
    cls_top{i} = cls;
end
clear cls

% 
% idx = 1;
% idx_plot = 1;
% figure;
% while idx <= length(I)
%     cls = clst_lifespan(I(idx), pdiag);
%     
%     if diag_filter(cls, 10, 10)
%         subplot(4,4,idx_plot);
%         draw_clMat(cls.cM,cls.sM,false);
%         idx_plot = idx_plot+1;
%         title(sprintf('basis %d',I(idx)));
%         
%         if idx_plot > 16
%             break
%         end
%     end
%     idx = idx + 1;
% end
% clear cls

%% Check all

dg_property = zeros(n_points,3);
for i = 1:n_points
    fprintf('process data: %d\n',i);
    
    cls = clst_lifespan(I_blife(i), pdiag);
    [dg_property(i,1),dg_property(i,2),dg_property(i,3)] = diagram_property_2(cls, ones(1,n_points));
end
clear cls
I_dg_large = dg_property(:,3)>0;
figure; hold on
scatter(dg_property(I_dg_large,1),dg_property(I_dg_large,2),2*sqrt(dg_property(I_dg_large,3))+10,'b','filled');
