addpath /Users/acb8/Documents/Research/MatlabUseful/
% uncomment to read in fresh data

% GOAL: find correlations between small subsections and the rest of the
% map 

%% read in data

ncfilename = 'arctic_seaice_climate_indicators_nh_v01r00_1979-2016.nc';
%ahra = ncread(ncfilename, 'AHRA'); % ahra melt onset          M
cfo = ncread(ncfilename, 'CFO'); % cont freeze onset          F
cmo = ncread(ncfilename, 'CMO'); % cont melt onset            M
doa = ncread(ncfilename, 'DOA'); %day of advance              F
doc = ncread(ncfilename, 'DOC'); %day of closing              F
doo = ncread(ncfilename, 'DOO'); %day of opening              M
dor = ncread(ncfilename, 'DOR'); %day of retreat              M
efo = ncread(ncfilename, 'EFO'); %early freeze onset          F
emo = ncread(ncfilename, 'EMO'); %early melt onset            M
iifp = ncread(ncfilename, 'IIFP'); %inner ice free period
oifp = ncread(ncfilename, 'OIFP'); %outer ice free period
sgip = ncread(ncfilename, 'SGIP'); %seasonal gain of ice period
siz = ncread(ncfilename, 'SIZ'); %seasonal ice zone
slip = ncread(ncfilename, 'SLIP'); %seasonal loss of ice period

lat = ncread(ncfilename, 'latitude'); %latitude
lon = ncread(ncfilename, 'longitude'); %longitude
year = ncread(ncfilename, 'year'); %year
land = hdfread('amsr_gsfc_25n.hdf', 'landmask');

%% cleaning
nodatamap= min(emo, [], 3) < 0;
yesdatalist = find(nodatamap < 1); % contains list of all indices with usable data

% orient land mask the correct way
land = flip(land);
land = rot90(land);
land = rot90(land);
land = rot90(land);


%% find new metrics for the whole field

% find meltDays: doo - cmo (days between date of opening and continuous melt onset)
dooSubCmo = nan(304, 448, 38); % preallocate
dorSubCmo = nan(304, 448, 38); % preallocate

% only use pixels in yesdatalist to avoid bad datapoints
for m = 1:length(yesdatalist)
     n = yesdatalist(m);
  
     % for each elem in yesdatalist, returns its row/column in a 304 x 448
     % matrix
     [r,c] = ind2sub([304, 448], n); 
    
     % for each elem in yesdatlist, find (doo - cmo)
     dooSubCmo(r,c,:) = doo(r,c,:) - cmo(r,c,:);
     dorSubCmo(r,c,:) = dor(r,c,:) - cmo(r,c,:);
   
end

%% find average data values over subsection

% create subsection(s): these element-by-element structures will contain
% all the relevant informatin about each subsection so code is easily
% reusable

% % comment out to use these locations instead
clear subsec;

subsec(1) = chukchi;
subsec(1).area = pick_area_from_start(subsec(1));

subsec(2).name = utqiagvik;
subsec(2).area = pick_area_from_start(subsec(2));

subsec(3).name = beaufort;
subsec(3).area = pick_area_from_start(subsec(3));

subsec(4).name = ebeaufort;
subsec(4).area = pick_area_from_start(subsec(4));

subsec(5).name = hudsonbay;
subsec(5).area = pick_area_from_start(subsec(5));

% UNUSED SUBSECS... 
subsec(6).name = barents;
subsec(6).area = pick_area_from_start(subsec(6));
subsec(7).name = kara;
subsec(7).area = pick_area_from_start(subsec(7));
subsec(8).name = laptev;
subsec(8).area = pick_area_from_start(subsec(8));
subsec(9).name = ESibS;
subsec(9).area = pick_area_from_start(subsec(9));
subsec(10).name = ESibN;
subsec(10).area = pick_area_from_start(subsec(10));



% calculate averages for each subsection and add that to structure
for i = 1:length(subsec)
    subsec(i).doaAvg = getAvg(doa, subsec(i).area);
    subsec(i).docAvg = getAvg(doc, subsec(i).area);

end


%% find Correlations between subsection and the whole field

% adds coeffs and pvalues to every element in the structure
for i = 1:length(subsec)   
    % doa vs. ___ (unused correlations)
%    [subsec(i).cc_doa_dor, subsec(i).pp_doa_dor] = findCorrelation(yesdatalist, subsec(i).doaAvg, dor, 1, 37, 2, 38);
%    [subsec(i).cc_doa_doo, subsec(i).pp_doa_doo] = findCorrelation(yesdatalist, subsec(i).doaAvg, doo, 1, 37, 2, 38);
%    [subsec(i).cc_doa_cmo, subsec(i).pp_doa_cmo] = findCorrelation(yesdatalist, subsec(i).doaAvg, cmo, 1, 37, 2, 38);
%     [subsec(i).cc_doa_dooSubCmo, subsec(i).pp_doa_dooSubCmo] = findCorrelation(yesdatalist, subsec(i).doaAvg, dooSubCmo, 1, 37, 2, 38);
%     [subsec(i).cc_doa_dorSubCmo, subsec(i).pp_doa_dorSubCmo] = findCorrelation(yesdatalist, subsec(i).doaAvg, dorSubCmo, 1, 37, 2, 38);
%     
    % doc vs. ___
   [subsec(i).cc_doc_dor, subsec(i).pp_doc_dor] = findCorrelation(yesdatalist, subsec(i).docAvg, dor, 1, 37, 2, 38);
   [subsec(i).cc_doc_doo, subsec(i).pp_doc_doo] = findCorrelation(yesdatalist, subsec(i).docAvg, doo, 1, 37, 2, 38);
   [subsec(i).cc_doc_cmo, subsec(i).pp_doc_cmo] = findCorrelation(yesdatalist, subsec(i).docAvg, cmo, 1, 37, 2, 38);
   [subsec(i).cc_doc_dooSubCmo, subsec(i).pp_doc_dooSubCmo] = findCorrelation(yesdatalist, subsec(i).docAvg, dooSubCmo, 1, 37, 2, 38);
   [subsec(i).cc_doc_dorSubCmo, subsec(i).pp_doc_dorSubCmo] = findCorrelation(yesdatalist, subsec(i).docAvg, dorSubCmo, 1, 37, 2, 38);
   
   % FUD Analysis
   [subsec(i).cc_doc_doc, subsec(i).pp_doc_doc] = findCorrelation(yesdatalist, subsec(i).docAvg, doc, 1, 37, 1, 37);
   [subsec(i).cc_doc_doa, subsec(i).pp_doc_doa] = findCorrelation(yesdatalist, subsec(i).docAvg, doa, 1, 37, 1, 37);
   [subsec(i).cc_doa_doc, subsec(i).pp_doa_doc] = findCorrelation(yesdatalist, subsec(i).doaAvg, doc, 1, 37, 1, 37);
   [subsec(i).cc_doa_doa, subsec(i).pp_doa_doa] = findCorrelation(yesdatalist, subsec(i).doaAvg, doa, 1, 37, 1, 37);
 
end


%% plot things

 plotSubsecVsWhole(nodatamap, land, subsec, ["doc"], ["cmo", "doo", "dor", "dooSubCmo"], 12, 0.05);


%% helper functions

% returns a 1d list of the average data values over that subsection for each year
function [avgs] = getAvg(data, subsec)
    avgs = nan(38,1);
    
    for year = 1:38
        subsecData = data(:, :, year); % slice off one year
        subsecData = subsecData(subsec); % get data only from the subsection
        subsecData = subsecData(subsecData>0); % keep only subsec data that's greater than 0
        tot = sum(subsecData, 'all'); % find the sum of all data points      
        numDataPoints = size(subsecData);
        avg = tot/numDataPoints(1); % find the average of those data points
        avgs(year) = avg; %add that year's average to the list of averages
    end
    
end

% returns matrices of correlation coeffs and p-vals between list of avgs
% and data over the subsection
function [cc, pp] = findCorrelation(yesdatalist, avgs, data, y1a, y1b, y2a, y2b)
    cc = nan(304,448);
    pp = nan(304,448);

    % find actual correlation values and populate nan matrices
    for m = 1:length(yesdatalist)
         n = yesdatalist(m);

         % for each elem in yesdatalist, returns its row/column in a 304 x 448
         % matrix
         [r,c] = ind2sub([304, 448], n); 

         % get the matrices of correlation coeffs/p vals for each elem in yes
         % data list & the subsection using custom method
         [x, p] = corrcoeffsubset(avgs, data, r, c, y1a, y1b, y2a, y2b); cc(n) = x(2); pp(n) = p(2);

    end
end

% finds corrcoeffs and p-vals for data over the specified number of years
function [x, p] = corrcoeffsubset(data1, data2, r, c, y1a, y1b, y2a, y2b) 

    y = data1(y1a:y1b); % list of average in subsection for each year
    d = data2(r, c, y2a:y2b); % values in pixel specified by [r,c] over same years
    d = d(:); % make d into 1d array so that it matches the shape of y

    % only find correlation on data where both datalists are >0 
    y1 = y(y>0 & d>0); 
    d1 = d(y>0 & d>0);


    % if there are enough data points, find matrix of correlation
    % coefficients and p values
    if length(y1) > 10 && length(d1) > 10
        [x, p] = corrcoef(y1, d1);
    else % otherwise return no signficant p val and no corrcoeff
        x = nan(1,2);
        p = [1 1];
    end

end

function mapgrid = pick_area_from_start(loc_struct)

rows = loc_struct.start_row_pixel : loc_struct.start_row_pixel+10;
cols = loc_struct.start_col_pixel : loc_struct.start_col_pixel+10;
background = zeros(304,448);
background(rows, cols) = 1;
mapgrid = background > 0;

end