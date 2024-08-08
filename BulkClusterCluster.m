
Ymat = zeros(43);

for MM = 240:244       %columns
    for NN = 142:146   %rows
    
    %if (MM == 256 && NN == 166) || (MM == 254 && NN == 168)...
    %        || (MM == 258 && NN == 168) || (MM == 259 && NN == 168)

    %else
          start_col = MM%182;%203;
          start_row = NN%161;%160;%200;
        try
          kmeans_dbscrub_gen;
          AB_clusterclustering;
        catch
        end
        
        close all;
        end

    %end
end
 

    figure(51)
    imagesc(years, years, Ymat)

    figure(52)
    plot(graph(Ymat, 'omitselfloops'))