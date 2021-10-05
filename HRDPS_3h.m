%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name         : HRDPS_3h
% Author       : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
% Date         : Tue Sep 21 15:17:00 2021
% Description  : This code converts the High-Resolution Deterministic Prediction System (HRDPS) temperature and precipitation 
%                into Matlab variables and temporally aggregates them to a 3h time step, and spatially averages them to the 
%                catchment scale.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
%%                                            THE ONLY PART TO MODIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declarations

% Directories
Dir.dataPath = '\netCDF\files\Path'              ; addpath(Dir.dataPath);  % To Modify
Dir.shpPath  = '\Shapefiles\Path'                ; addpath(Dir.shpPath);   % To Modify
OutPath = '\output\directory\Patha';                                       % To Modify

if ~exist(fullfile(OutPath), 'dir')
    mkdir(fullfile(OutPath));         addpath(OutPath);
end

% Catchments name
nameC =  {'Catch1';'Catch2';'...CatchN'};   % To Modify
nBV   = numel(nameC);

% Forecasts setting
nbLT = 48; % Lead Time
ts   = 3;  % time step

% Getting the netCDF files name
cd(Dir.dataPath)
ncFiles  = dir('*.nc');

for ifile = 1: size(ncFiles,1)
     tmp = split(convertCharsToStrings(ncFiles(ifile).name), ".nc");
     nameFile(ifile) = tmp(1)   
end

nameFile = unique(nameFile);
dateRef = datenum(nameFile,'yyyymmddHH');
nDays    = numel(dateRef);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                            DO NOT TOUCH FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: CREATE CATCHMENT MASK

% Import NetCDF coordinates
fileToRead=fullfile(Dir.dataPath,nameFile(1));
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
    [S]=shaperead(fullfile(Dir.shpPath,sprintf('%s.shp',nameC{iCatch})));
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

%% Index to accumulate and get the variables at 3h time step
nhr= 48; % time dimension
index = 1:nhr;
elem  = [repmat(ts,1,floor(nhr/ts))];
endv  = nhr-sum(elem);
if(~endv)
    endv = [];
end
index = mat2cell(index,1,[elem,endv])';

%% Initialization
PttmpCatch = NaN(nDays,nhr/ts,nBV);
TtmpCatch  = NaN(nDays,nhr/ts,nBV);
TminCatch  = NaN(nDays,nhr/ts,nBV);
TmaxCatch  = NaN(nDays,nhr/ts,nBV);
Datehr     = NaN(nDays,nhr);

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
    fileToRead = fullfile(Dir.dataPath,nameFile(iDate));
    ncid       = netcdf.open(fileToRead,'NC_NOWRITE');
    
    % Retrieve RDRS variables
    % time
    date = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'),'single');
    NonZero = find(date>0); % Eliminate the lead time corresponding to the issue time.
    date = date(NonZero);
    Datehr(iDate,date) = double(date)./24 + dateRef(iDate);
    ntime = numel(date);
    
    % Precipitation
    dataP = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'HRDPS_P_PR_SFC'),'single');
    dataP = (dataP .*1000); % m --> mm
    dataP = dataP(:,:,NonZero);
    
    % Temperature
    dataT = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'HRDPS_P_TT_09950'),'single');
    dataT = dataT(:,:,NonZero);
         
    % Close NetCDF file
    netcdf.close(ncid);
    
    % Trick for computing cathcment mean at catchment scale
        Ptmp00 = arrayfun(@(iLT) squeeze(dataP(:,:,iLT)),1:ntime,'UniformOutput',0);
        Tmp00  = arrayfun(@(iLT) squeeze(dataT(:,:,iLT)),1:ntime,'UniformOutput',0);
    
   
    % -----------------------------------------------------------------------------------------------------------------------
    %% Step 3: Compute mean at the catchment scale - Catchment loop
    % -----------------------------------------------------------------------------------------------------------------------
    
    for iCatch = 1:nBV
        
        % Retrieve catchment mask
        inan = isnan(inGrid.(sprintf('C%s',nameC{iCatch})));
        
        Ptmp = arrayfun(@(iLT) mean(Ptmp00{iLT}(~inan)),1:ntime);
        Tmp  = arrayfun(@(iLT) mean(Tmp00{iLT}(~inan)),1:ntime);
        
        Pttmp = NaN(nbLT,1);
        Tttmp = NaN(nbLT,1);
        
        Pttmp(date,1) = Ptmp;         
        Pttmp(2:end) = diff(Pttmp); % Decummulate precipitation
       
        Tttmp(date,1) = Tmp;  
                      
        PttmpCatch(iDate,:,iCatch) = cell2mat(cellfun(@(x) sum(Pttmp(x)),index,'un',0));
        TtmpCatch(iDate,:,iCatch)  = cell2mat(cellfun(@(x) mean(Tttmp(x)),index,'un',0));
        TminCatch(iDate,:,iCatch)  = cell2mat(cellfun(@(x) min(Tttmp(x)),index,'un',0));
        TmaxCatch(iDate,:,iCatch)  = cell2mat(cellfun(@(x) max(Tttmp(x)),index,'un',0));
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                            END :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
