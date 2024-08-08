
load BeaufortCluster2.mat
Location(1).name = "Beaufort";
Location(1).startrow = 74;
Location(1).startcol = 213;
Location(1).A = years(idxtemp == 1);%[1980, 1981, 1982, 1984, 1985, 1986,  1987, 1990,  1993, 1998, 1999,  2001,  2011];
Location(1).B = years(idxtemp == 2);%[1997, 2004, 2010, 2015,  2017, 2020];
Location(1).C = years(idxtemp == 3);%[1983, 1989,  1991, 1992,  1994,  1995,  1996,  2000,  2003, 2005, 2019,  2021];
Location(1).D = years(idxtemp == 4);%2000:2016;
Location(1).block = zeros(304, 448);
Location(1).block(Location(1).startrow:Location(1).startrow+10, Location(1).startcol:Location(1).startcol+10) = 1;

load LaptevCluster.mat
Location(2).name = "Laptev";
Location(2).startrow = 160;
Location(2).startcol = 165;
Location(2).A = years(idxtemp == 1);%[1980, 1981, 1982, 1983, 1984, 1985, 1986, 1988, 1989, 1990, 1991, 1993, 1995, 1998, 1999,  2000,  2002,  2005];
Location(2).B = years(idxtemp == 2);%[2006,  2007, 2008, 2013,  2018, 2021];
Location(2).C = years(idxtemp == 3);%[1987, 1992, 1994,  1996,  1997,  2001,  2003, 2004,  2009,  2010, 2011,  2012,  2014, 2015,  2017,  2019, 2020];
Location(2).D = years(idxtemp == 4);%2000:2016;
Location(2).block = zeros(304, 448);
Location(2).block(Location(2).startrow:Location(2).startrow+10, Location(2).startcol:Location(2).startcol+10) = 1;

load ChukchiCluster.mat
Location(3).name = "Chukchi";
Location(3).startrow = 85;
Location(3).startcol = 184;
Location(3).A = years(idxtemp == 1);%[1996, 2000, 2003, 2004, 2006, 2007, 2011, 2014, 2018, 2021];
Location(3).B = years(idxtemp == 2);%[1982, 1983, 1984, 1998, 2017];
Location(3).C = years(idxtemp == 3);%[1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993,  1995,  1997, 1999,  2002,  2005,  2008, 2010,  2012,  2013, 2015,  2016, 2019, 2020];
Location(3).D = years(idxtemp == 4);%2000:2016;
Location(3).block = zeros(304, 448);
Location(3).block(Location(3).startrow:Location(3).startrow+10, Location(3).startcol:Location(3).startcol+10) = 1;

Location(4).name = "Elesmere";
Location(4).startrow = 120;
Location(4).startcol = 240;
Location(4).A = [1982, 1984, 1985, 1986, 1998, 2003, 2004, 2007,   2012,   2017,  2019];
Location(4).B = [1988, 1991, 1992, 1994, 2008, 2010,  2011,2013, 2018];
Location(4).C = [1980,  1981,  1983, 1987, 1989,  1990, 1993,  1995,1996,   1997, 1999, 2000,  2001,    2002,   2005,  2006,   2009, 2014,  2015, 2020,  2021];
Location(4).D = 2000:2016;
Location(4).block = zeros(304, 448);
Location(4).block(Location(4).startrow:Location(4).startrow+10, Location(4).startcol:Location(4).startcol+10) = 1;


%% Find FUD for the box

% calculate averages for each subsection and add that to structure
for i = 2%1:length(Location)
   [Location(i).Ax.doaAvg Location(i).Ax.doano]  = getAvg(doa, Location(i).startrow, Location(i).startcol, Location(i).A);
    [Location(i).Ax.docAvg Location(i).Ax.docno]= getAvg(doc, Location(i).startrow, Location(i).startcol, Location(i).A);

    [Location(i).Bx.doaAvg Location(i).Bx.doano]= getAvg(doa, Location(i).startrow, Location(i).startcol, Location(i).B);
    [Location(i).Bx.docAvg Location(i).Bx.docno]= getAvg(doc, Location(i).startrow, Location(i).startcol, Location(i).B);

    [Location(i).Cx.doaAvg Location(i).Cx.doano] = getAvg(doa, Location(i).startrow, Location(i).startcol, Location(i).C);
    [Location(i).Cx.docAvg Location(i).Cx.docno] = getAvg(doc, Location(i).startrow, Location(i).startcol, Location(i).C);

    [Location(i).Dx.doaAvg Location(i).Dx.doano] = getAvg(doa, Location(i).startrow, Location(i).startcol, Location(i).D);
    [Location(i).Dx.docAvg Location(i).Dx.docno] = getAvg(doc, Location(i).startrow, Location(i).startcol, Location(i).D);

end



%% make figures

M = 4*2;
for L = 2%length(Location)
figure(L)
set(gcf, 'Position', [25 375 1312 782]);

Tx = tight_subplot(4, M, 0.02, 0.05, 0.05);

    %% Cluster A/DOC
        % Row 1: DOO
           
            X = Location(L).Ax.docAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year 
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+1));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Ax.docAvg;
            Z = Location(L).A-1978; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+1));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Ax.docAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+1));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Ax.docAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+1));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

        % Row 1: DOO
           
            X = Location(L).Ax.doaAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+2));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Ax.doaAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+2));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Ax.doaAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+2));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Ax.doaAvg;
            Z = Location(L).A-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+2));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

            %% Cluster B
        % Row 1: DOO
           
            X = Location(L).Bx.docAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+3));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Bx.docAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+3));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Bx.docAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+3));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Bx.docAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+3));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

        % Row 1: DOO
           
            X = Location(L).Bx.doaAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+4));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Bx.doaAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+4));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Bx.doaAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+4));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Bx.doaAvg;
            Z = Location(L).B-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+4));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));


    %% Location C

            % Row 1: DOO
           
            X = Location(L).Cx.docAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+5));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Cx.docAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+5));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Cx.docAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+5));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Cx.docAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+5));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

        % Row 1: DOO
           
            X = Location(L).Cx.doaAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+6));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Cx.doaAvg;
            Z = Location(L).C-1978; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+6));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Cx.doaAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+6));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Cx.doaAvg;
            Z = Location(L).C-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+6));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));


              %% Location D

            % Row 1: DOO
           
            X = Location(L).Dx.docAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+7));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Dx.docAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+7));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Dx.docAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+7));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Dx.docAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+7));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

        % Row 1: DOO
           
            X = Location(L).Dx.doaAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = doo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(0*M+8));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
        % Row 2: DOR
  
            X = Location(L).Dx.doaAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(1*M+8));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 3: CMO
            X = Location(L).Dx.doaAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(2*M+8));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));
    
    
        % Row 4: DOR - CMO
            X = Location(L).Dx.doaAvg;
            Z = Location(L).D-1978 ; % make this the next calendar year
            Y = dor(:,:,Z(Z<=37)+1) - cmo(:,:,Z(Z<=37)+1);
            
            corrmap = corrmapfreezemelt(X, Y, yesdatalist);
            axes(Tx(3*M+8));
            subimageplot_cc(corrmap, Location(L).block, land, "PRGn", [-1.04 1], sum(~isnan(X)));

end

%% FUNCTIONS

function [avgs, no] = getAvg(data, startrow, startcol, years)

    years = years(years <= 2016);
    subset = data(startrow:startrow+10, startcol:startcol+10, years-1978);
    subset(subset == 0) = NaN;
    avgs = squeeze(mean(subset, [1 2], 'omitnan'));
    no = sum(~isnan(subset))/121;
end


function [corrmap] = corrmapfreezemelt(X, Y, yesdatalist)

Y(Y<= 0) = NaN;

corrmap = nan(304, 448);

if sum(~isnan(X)) > 5

    for m = 1:length(yesdatalist)
        [r, c] = ind2sub([304 448], yesdatalist(m));

        a = squeeze(Y(r,c,:));

        [coef, p] = corrcoef(X(~isnan(a)), a(~isnan(a)) );
        count = length(X(~isnan(a)));
if length(p) > 1
        if p(2) < .05 && count >= 4
            corrmap(r,c) = coef(2);
        end
end


    end

end


end
