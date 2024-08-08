
start_row = 230;%212;
start_col = 214;%228;

start_year = 1980;
end_year = 2021;
num_years = end_year - start_year;
start_day = date_index([start_year 9 15], ymd_time); % sept 15th, 1980 

last_locs = get_end_locs(start_col, start_row, v, u, ymd_time, x_ref, y_ref);
yrs_filtered = last_locs;%vecs(:, logicalMask);

idx = idxtemp;

figure(13);
%geographic_basemap({start_col, start_row}, true, latitude, longitude);
plot_clusters(idx, yrs_filtered, {start_col, start_row}, latitude, longitude);
sgtitle("Beau CLUSTERS")

figure(16)
plot_clusters_hist(idx, yrs_filtered, {start_col, start_row}, latitude, longitude);

%% helpers 👋

function plot_clusters(idx, vecs, start_pixel, latitude, longitude)
    [start_col, start_row] = start_pixel{:};
    num_groups = length(unique(idx));
    lets = 'ABCDEFGHIJKLMNOPQ';
    
    for i = 1:num_groups
        subplot(1, num_groups, i)
        hold on
        geographic_basemap({start_col, start_row}, true, latitude, longitude);

        % if there are outliers, grp id is the val of that outlier
        grpid = i;
        if i == num_groups && sum(idx(find(idx == i))) == 0
            grpid = -1;
        end
        
        
        for j = 1:size(idx,1)

           if idx(j) == grpid
              scatter_vec(vecs(:,j),'b',latitude, longitude);
           end
        end
        hold off
        title(lets(i)+", N="+sum(idx==grpid), 'FontSize',20)
    
    end

end


function plot_clusters_hist(idx, vecs, start_pixel, latitude, longitude)
    [start_col, start_row] = start_pixel{:};
    num_groups = length(unique(idx));
    lets = "ABCDE";
    
    for i = 1:num_groups

        countmap = histomapgram_vecs(vecs(:, idx==i), 361, 361);
        countmap(countmap < 1) =NaN;
        subplot(1, num_groups, i)
        hold on
        geographic_basemap({start_col, start_row}, false, latitude, longitude);        
        pcolorm(latitude, longitude, countmap, 'FaceAlpha', 1); 
        shading flat
        colormap(cbrewer('seq', 'Oranges', 10))
       
        land = readgeotable("landareas.shp");
        geoshow(land, "FaceColor", [0.8 0.8 0.8])

        box = [10, 10]; % size of boxed pixels
        box_lats = [latitude(start_row, start_col), ...
            latitude(start_row, start_col + box(1)), ...
            latitude(start_row + box(1), start_col+box(1)), ...
            latitude(start_row + box(1), start_col), ...
            latitude(start_row, start_col)...
            ];
        
        box_lons = [longitude(start_row, start_col), ...
            longitude(start_row, start_col + box(2)), ...
            longitude(start_row + box(2), start_col + box(2)), ...
            longitude(start_row + box(2), start_col), ...
            longitude(start_row, start_col)...
            ];

        linem(box_lats, box_lons, 'k', 'LineWidth', 1);
        hold off
        caxis([0 15])
        colorbar('southoutside')
    
    end

end

function plot_clusters_same_graph(idx, vecs, start_pixel, latitude, longitude)

    [start_col, start_row] = start_pixel{:};
    num_groups = length(unique(idx));
    colors = 'rbgkyo';
    
    for i = 1:num_groups

        hold on
        geographic_basemap({start_col, start_row}, true, latitude, longitude);

        % if there are outliers, grp id is the val of that outlier
        grpid = i;
        if i == num_groups && sum(idx(find(idx == i))) == 0
            grpid = -1;
        end
        
        
        for j = 1:size(idx,1)

           if idx(j) == grpid
              scatter_vec(vecs(:,j),colors(i),latitude, longitude);
           end
        end
        hold off
    
    end

    t = "";
    for i = 1:num_groups
        t = t + "n_"+colors(i) + "=" + sum(idx==i)+"   ";
    end
    title(t);

end


function plot_components(idx, vecs, id, start_pixel, latitude, longitude)
    grp_vecs = zeros(size(vecs, 1), length(find(idx == id))); % 200 x num in group
    
    % gather points in vector    
    hold on
    i_grp_vecs = 1;
    for j = 1:size(idx,1)

       if idx(j) == id
            grp_vecs(:, i_grp_vecs) = vecs(:,j);
            i_grp_vecs = i_grp_vecs + 1;
       end

    end
    hold off

    num_members = size(grp_vecs, 2);

    % plot all
    figure(20);
    geographic_basemap(start_pixel, true, latitude, longitude);
    hold on
    for i = 1:num_members
        scatter_vec(grp_vecs(:,i), 'b', latitude, longitude);
    end
    hold off

    num_cols = 4;
    num_rows = ceil(num_members/num_cols);

    % plot each group individually
    figure(21);
    for i = 1:num_members
        subplot(num_rows, num_cols, i);
        geographic_basemap(start_pixel, true, latitude, longitude);
        scatter_vec(grp_vecs(:,i), 'b', latitude, longitude);
    end

    
end

function paired_vec = pair(vec)
    paired_vec = zeros(length(vec)/2, 2);
    
    i = 1;
    i_paired = 1;
    while i < 200
        paired_vec(i_paired, 1) = vec(i);
        paired_vec(i_paired, 2) = vec(i + 1);
        i = i + 2;
        i_paired = i_paired + 1;
    end

end

function unpaired_vec = unpair(vec)
    unpaired_vec = zeros(length(vec)*2, 1);
    
    i = 1;
    i_paired = 1;
    while i<length(vec)*2
        unpaired_vec(i) = vec(i_paired,1);
        unpaired_vec(i+1) = vec(i_paired,2);
    
        i = i + 2;
        i_paired = i_paired + 1;
    end

end

% use dbscan with a wide radius to  find outlier pixels, and replace them
% with the nearest non-outlier pixel
function [clean_vec, foundOutliers] = clean(vec, start_col, start_row, latitude, longitude)
    foundOutliers = false;
    
    % convert vector into x y pairs
    paired_vec = pair(vec);
    
    % cluster (-1 is outlier)
    % 5 points minimum to be a core point
    % within a radius of 10 pixels
    dbclust_vec = dbscan(paired_vec, 9, 12);
    
    % plot outliers and clusters
    %figure();
    %if sum(dbclust_vec == -1) > 0
     %   geographic_basemap({start_col, start_row}, true, latitude, longitude);
      %  scatter_vec(vec, dbclust_vec, latitude, longitude);
    %title("IDENTIFIED OUTLIERS")
    %end
    % remove outliers by... replacing with nearest non-outlier pixel
    for i = 1:length(dbclust_vec)
        
        % outlier
        if dbclust_vec(i) == -1
            foundOutliers = true;
    
            % outlier becomes value of nearest nonoutlier
            if (i > 1)
                paired_vec(i,:) = paired_vec(i - 1,:);
            else % edge case where no previous val is outlier
                j = 2; % look ahead for nonoutlier
                while j<length(dbclust_vec)
                    if dbclust_vec(j) ~= -1
                        paired_vec(i,:) = paired_vec(j, :);
                        break;
                    end
                end
            end
        end
    end
    
    % plot where outlier values got mapped to
    % figure();
    % geographic_basemap({start_col, start_row}, true, latitude, longitude);
    % scatter_vec(unpair(paired_vec), dbclust_vec, latitude, longitude);
    
    clean_vec = unpair(paired_vec);

end
