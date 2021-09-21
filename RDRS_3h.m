%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name         : RDRS_3h 
% Author       : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca) 
% Date         : Tue Sep 21 15:17:00 2021
% Description  : This code converts the Regional Deterministic Reforecast System (RDRS_v2; hourly dataset) netCDF files to 
%                Matlab format and temporally aggregates them to a 3h time step, and spatially averages them to the catchment 
%                scale.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
%%                                            THE ONLY PART TO MODIFY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declarations

% Directories
dataPath = 'Path\of\the\RDRS\data'  ; addpath(dataPath);  # To Modify
shpPath  = 'Shapefile\path'         ; addpath(shpPath);   # To Modify
OutPath = 'Output\path'             ; addpath(OutPath);   # To Modify

if ~exist(fullfile(OutPath), 'dir')
    mkdir(fullfile(OutPath));         addpath(OutPath);
end

% Catchments name
nameC = {'name 1';'name 2';'name 3';... 'name n'};    # To Modify
nBV   = numel(nameC);

% Period of the time series (.nc files)
dateStart = 'YYYY/MM/DD 12:00:00';  # To Modify
dateEnd   = 'YYYY/MM/DD 12:00:00';  # To Modify
dateRef   = datenum(dateStart):1:datenum(dateEnd);
nDays     = numel(dateRef);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                            DO NOT TOUCH FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: CREATE CATCHMENT MASK

% Import NetCDF coordinates
fileToRead=fullfile(dataPath,strcat(datestr(dateRef(1),'yyyymmddHH'),'.nc'));
ncid = netcdf.open(fileToRead,'NC_NOWRITE');
lat0 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'),'single');
lon0 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'),'single');

% Conversion of Lat and Lon to match with the .shp format
lon0b = lon0';
lat0b = lat0';

% close netCDF file
netcdf.close(ncid);

% Build catchment masks
for iCatch = 1:nBV
    % Import ctch shape
    [S]=shaperead(fullfile(shpPath,sprintf('%s.shp',nameC{iCatch})));
    % Get points inside the catchment
    inGrid_tmp = inpolygon(lon0b,lat0b,S.X,S.Y);
    % Transpose mask for NetCDF compatibility (y,x,T)
    inGrid_tmp = transpose(inGrid_tmp);
    % Convert into nan/value mask for data extraction
    inGrid_tmp = double(inGrid_tmp);
    inGrid_tmp(inGrid_tmp==0) = NaN;
    % Store for NetCDF extraction
    inGrid.(sprintf('C%s',nameC{iCatch})) = inGrid_tmp;
end

% Housekeeping
clear ncid S lat0  lon0  inGrid_tmp

%% Index to accumulate get the variables at 3h time step
nhr= 24; % time dimension: 24 hr 
index = 1:nhr;
elem  = [repmat(3,1,floor(nhr/3))];
endv  = nhr-sum(elem);
if(~endv)
    endv = [];
end
index = mat2cell(index,1,[elem,endv])';

%% Initialization
PttmpCatch = NaN(8,nDays,nBV);
TtmpCatch  = NaN(8,nDays,nBV);
TminCatch  = NaN(8,nDays,nBV);
TmaxCatch  = NaN(8,nDays,nBV);
Datehr     = NaN(nhr,nDays);

% ---------------------------------------------------------------------------------------------------------------------------   
%% Step 2: Retrieve NetCDF data
% ---------------------------------------------------------------------------------------------------------------------------  
for iDate = 1:nDays
    
    % Display process
    if rem( iDate,round(numel(dateRef)/50,0) ) == 0
        mntoc = round(toc/60,1);
        fprintf('%2.0f %% of files read - time elapsed %s minutes \n',iDate/numel(dateRef)*100, mntoc)
    end
    
    % Open NetCDF file   
    fileToRead = fullfile(dataPath,strcat(datestr(dateRef(iDate),'yyyymmddHH'),'.nc'));
    ncid       = netcdf.open(fileToRead,'NC_NOWRITE');
    
    % Retrieve RDRS variables
    % Precipitation
    dataP = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'RDRS_v2_P_PR0_SFC'),'single');
    dataP = (dataP .*1000); % m --> mm
    
     % Temperature
    dataT = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'RDRS_v2_P_TT_1.5m'),'single');
    
    % time
    date = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'),'single');
    Datehr(:,iDate) = double(date)./nhr + datenum('2000-01-01 12:00:00') + iDate;
    
    % Close NetCDF file
    netcdf.close(ncid);
    
    % Trick for computing cathcment mean at catchment scale
    Ptmp00 = arrayfun(@(ihr) squeeze(dataP(:,:,ihr)),1:nhr,'UniformOutput',0);
    Tmp00  = arrayfun(@(ihr) squeeze(dataT(:,:,ihr)),1:nhr,'UniformOutput',0);
    
    
    % -----------------------------------------------------------------------------------------------------------------------
    %% Step 3: Compute mean at the catchment scale - Catchment loop
    % -----------------------------------------------------------------------------------------------------------------------
   
    for iCatch = 1:nBV
        
        % Retrieve catchment mask
        inan = isnan(inGrid.(sprintf('C%s',nameC{iCatch})));
        
        Ptmp = arrayfun(@(ihr) mean(Ptmp00{ihr}(~inan)),1:nhr);
        Tmp  = arrayfun(@(ihr) mean(Tmp00{ihr}(~inan)),1:nhr);
        
        Ptmp = transpose(Ptmp);
        Ttmp = transpose(Tmp);
             
        PttmpCatch(:,iDate,iCatch) = cell2mat(cellfun(@(x) sum(Ptmp(x)),index,'un',0));
        TtmpCatch(:,iDate,iCatch)  = cell2mat(cellfun(@(x) mean(Ttmp(x)),index,'un',0));
        TminCatch(:,iDate,iCatch)  = cell2mat(cellfun(@(x) min(Ttmp(x)),index,'un',0));
        TmaxCatch(:,iDate,iCatch)  = cell2mat(cellfun(@(x) max(Ttmp(x)),index,'un',0));
    end
    
    
end

% ---------------------------------------------------------------------------------------------------------------------------
% Define output file name
% ---------------------------------------------------------------------------------------------------------------------------

outfile = sprintf('%s/RDRS_3hr.mat',OutPath);
% Export
save(outfile,'PttmpCatch', 'TtmpCatch', 'TminCatch','TmaxCatch','Datehr', '-v6');

% ---------------------------------------------------------------------------------------------------------------------------
% Save output - catchment-wise
% ---------------------------------------------------------------------------------------------------------------------------

% Date
dateRef = transpose(reshape(Datehr,1,[]));
dateRef = dateRef(3:3:end);
Date    = datevec(dateRef);

for iCatch=1:nBV
    
    % Extract catchment values
    Pt    = transpose(reshape(PttmpCatch(:,:,iCatch),1,[]));
    T     = transpose(reshape(TtmpCatch(:,:,iCatch),1,[]));
    Tmax  = transpose(reshape(TmaxCatch(:,:,iCatch),1,[]));
    Tmin  = transpose(reshape(TminCatch(:,:,iCatch),1,[]));
      
    % Export
    outfile = sprintf('%s/RDRS_3h_%s.mat',OutPath,nameC{iCatch});
    save(outfile,'Pt','T','Tmin','Tmax','Date','-v7.3');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                            END :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
