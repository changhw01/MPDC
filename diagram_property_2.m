function [diag_max_dist, diag_max_dsty, max_population] = diagram_property_2(cls,data_population)

[n_h,n_r] = size(cls.sM);

%% must size must appear at the right most 

% diag_size = max(cls.sM(:));

max_population = 0;
for i = 1:n_h
    idx = find(cls.cM(i,:)>0,1,'last');    
    if isempty(idx)
        continue
    end
    cidx = cls.cM(i,idx);
    size_local = sum(data_population(find(cls.rep_clst{cidx})));
    if size_local > max_population
        max_population = size_local;
    end
end

%% max scale persistence must be in the bottom
diag_max_dist = 0;
for i = n_h:-1:1
    I = find(cls.sM(i,:));
    if ~isempty(I)
        diag_max_dist = max(I)-min(I)+1;
        break;
    end
end

%% max level set persistence must be the right must column
diag_max_dsty = 0;
I = find(cls.sM(:,1));
if ~isempty(I)
    diag_max_dsty = max(I)-min(I)+1;
end
