function draw_clMat(cvM,ncM, flag_new_fig, ctype)

% use 'summer' for clst diagram
% use 'jet' for partitions

[n_h, n_r] = size(cvM);

if nargin < 3
    flag_new_fig = true;
    ctype = 'summer';
end
if nargin < 4
    ctype = 'summer';
end

%% Prepare for color
ncM(ncM>0) = ncM(ncM>0) - min(ncM(ncM>0))+1;
max_ncM = max(ncM(:));

I_color = [1:max_ncM]';
if length(I_color) > 64
    % suppress it
    I_color_log = log2(I_color);
    I_color = floor(63*I_color_log/max(I_color_log(:))) + 1;
% else
%     % scale up
%     I_color = floor(63*I_color/max(I_color(:))) + 1;
end


    
%% base coordinates

centers_x = 1:n_r; 
centers_y = n_h:-1:1;

gap_conn = 0.2;
gap_disc = 0.4;
base_rec_width = gap_disc*2;
base_rec_height = gap_disc*2;

base_coords_x_left = centers_x - gap_disc;
base_coords_x_right = centers_x + gap_disc;
base_coords_y_bottom = centers_y - gap_disc;
base_coords_y_up = centers_y + gap_disc;

%% the main loop
if flag_new_fig   
    figure; 
end
cmap = colormap(ctype);
hold on;

%-- one rectangle for each clustering
for i = 1:n_h
    for j = 1:n_r       
        if cvM(i,j) > 0
            x_left = base_coords_x_left(j);
            y_bottom = base_coords_y_bottom(i);

            c = cmap(I_color(ncM(i,j)),:);
            rectangle('Position',[x_left,y_bottom,base_rec_width,base_rec_height], ...
                      'FaceColor', c,'LineStyle','none');        
        end
    end
end

%-- fill up the vertical gaps if necessary
for i = 1:n_h
    for j = 2:n_r
        
        if cvM(i,j)==0
            continue;
        end
        
        if cvM(i,j-1) == cvM(i,j)
            x_left = base_coords_x_right(j-1);
            y_bottom = base_coords_y_bottom(i);

            c = cmap(I_color(ncM(i,j)),:);
            rectangle('Position',[x_left,y_bottom,gap_conn,base_rec_height], ...
                      'FaceColor', c,'LineStyle','none');   
        end
    end
end

%-- fill up the horizonal gaps if necessary
% care needed for a 2x2 block
for i = 2:n_h
    for j = 1:n_r
                
        if cvM(i,j)==0
            continue;
        end
        
        if cvM(i-1,j) == cvM(i,j)
            x_left = base_coords_x_left(j);
            y_bottom = base_coords_y_up(i);
                        
            if (j<n_r) && (cvM(i,j)==cvM(i,j+1)) && (cvM(i,j)==cvM(i-1,j+1))
                rec_hlen = base_rec_width + gap_conn;
            else
                rec_hlen = base_rec_width;
            end
            
            c = cmap(I_color(ncM(i,j)),:);
            rectangle('Position',[x_left,y_bottom,rec_hlen,gap_conn], ...
                      'FaceColor', c,'LineStyle','none');   
        end
    end
end

%% -- fine ture the figure
title('Cluster Evolution','FontSize',12);
%xlabel('Distance Threshold Index')
% set(gca, 'XTick', 1:n_r);
% labels = cell(n_r,1);
% for i = 1:n_r
%     labels{i} = int2str(i);
% end
% set(gca, 'XTickLabel', labels, 'FontSize',12);
% 
% %ylabel('Density Threshold Index')
% set(gca, 'YTick', 1:n_h);
% labels = cell(n_h,1);
% for i = 1:n_h
%     labels{i} = int2str(n_h-i+1);
% end
% set(gca, 'YTickLabel', labels,'FontSize',12);

% colorbar('location','EastOutside')
% sorted_val = sort(I_color);
% labels = cell(6,1);
% inc = floor(length(sorted_val)/6);
% for i = 1:6
%     labels{i} = int2str(sorted_val(i*inc));
% end
% colorbar('YTickLabel', labels);

axis([0,n_r+1,0,n_h+1])
% grid on
% box on
