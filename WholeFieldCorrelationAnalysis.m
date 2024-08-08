%addpath /Users/acb8/Documents/Research/MatlabUseful/
% uncomment to read in fresh data

% GOAL: 

%% read in data

ncfilename = 'arctic_seaice_climate_indicators_nh_v01r00_1979-2016.nc';
ahra = ncread(ncfilename, 'AHRA'); % ahra melt onset          M
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
yesdatalist = find(nodatamap < 1);

% orient land mask the correct way
land = flip(land);
land = rot90(land);
land = rot90(land);
land = rot90(land);

%% for each grid cell

% change these to get summer versus winter
y2a = 1;
y2b = 37;
y1a = 1;
y1b = 37;

% creates array of NaN values for each grid cell
% create for loop with top vals and side vals...
clear wholeField;
[wholeField(1).cc_doa_dor, wholeField(1).pp_doa_dor] = findCorrelation(yesdatalist, doa, dor, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doa_doo, wholeField(1).pp_doa_doo] = findCorrelation(yesdatalist, doa, doo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doa_ahra, wholeField(1).pp_doa_ahra] = findCorrelation(yesdatalist, doa, ahra, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doa_cmo, wholeField(1).pp_doa_cmo] = findCorrelation(yesdatalist, doa, cmo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doa_emo, wholeField(1).pp_doa_emo] = findCorrelation(yesdatalist, doa, emo, y1a, y1b, y2a, y2b);

[wholeField(1).cc_cfo_dor, wholeField(1).pp_cfo_dor] = findCorrelation(yesdatalist, cfo, dor, y1a, y1b, y2a, y2b);
[wholeField(1).cc_cfo_doo, wholeField(1).pp_cfo_doo] = findCorrelation(yesdatalist, cfo, doo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_cfo_ahra, wholeField(1).pp_cfo_ahra] = findCorrelation(yesdatalist, cfo, ahra, y1a, y1b, y2a, y2b);
[wholeField(1).cc_cfo_cmo, wholeField(1).pp_cfo_cmo] = findCorrelation(yesdatalist, cfo, cmo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_cfo_emo, wholeField(1).pp_cfo_emo] = findCorrelation(yesdatalist, cfo, emo, y1a, y1b, y2a, y2b);

[wholeField(1).cc_doc_dor, wholeField(1).pp_doc_dor] = findCorrelation(yesdatalist, doc, dor, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doc_doo, wholeField(1).pp_doc_doo] = findCorrelation(yesdatalist, doc, doo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doc_ahra, wholeField(1).pp_doc_ahra] = findCorrelation(yesdatalist, doc, ahra, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doc_cmo, wholeField(1).pp_doc_cmo] = findCorrelation(yesdatalist, doc, cmo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_doc_emo, wholeField(1).pp_doc_emo] = findCorrelation(yesdatalist, doc, emo, y1a, y1b, y2a, y2b);

[wholeField(1).cc_efo_dor, wholeField(1).pp_efo_dor] = findCorrelation(yesdatalist, efo, dor, y1a, y1b, y2a, y2b);
[wholeField(1).cc_efo_doo, wholeField(1).pp_efo_doo] = findCorrelation(yesdatalist, efo, doo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_efo_ahra, wholeField(1).pp_efo_ahra] = findCorrelation(yesdatalist, efo, ahra, y1a, y1b, y2a, y2b);
[wholeField(1).cc_efo_cmo, wholeField(1).pp_efo_cmo] = findCorrelation(yesdatalist, efo, cmo, y1a, y1b, y2a, y2b);
[wholeField(1).cc_efo_emo, wholeField(1).pp_efo_emo] = findCorrelation(yesdatalist, efo, emo, y1a, y1b, y2a, y2b);


%% plot things

%plotWholeFieldAnalysis(wholeField, land, ["doa", "doc", "cfo", "efo"],["emo", "cmo", "doo", "dor"], 0.05);

%% helper functions--------------------------------------------------------

function [cc, pp] = findCorrelation(yesdatalist, data1, data2, y1a, y1b, y2a, y2b)
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
         [x, p] = corrcoeffsubset(data1, data2, r, c, y1a, y1b, y2a, y2b); cc(n) = x(2); pp(n) = p(2);

    end
end

function [x, p] = corrcoeffsubset(data1, data2, r, c, y1a, y1b, y2a, y2b)

     y = data1(r, c, y1a:y1b);
     d = data2(r, c, y2a:y2b);
     
     % ????? y1 is all vals of y > 0 and d>0
     % ?? why is y1 different from d1 ??
     y1 = y(y>0 & d>0); 
     d1 = d(y>0 & d>0);
     
     % if there are enough data points, find matrix of correlation
     % coefficients and p values
     if length(y1) > 10 && length(d1)> 10
        [x, p] = corrcoef(y1, d1);
     else % otherwise return no correlation
         x = nan(1,2);
         p = [1 1];
     end

end