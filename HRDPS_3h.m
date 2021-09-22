%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name         : HRDPS_3h
% Author       : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
% Date         : Tue Sep 21 15:17:00 2021
% Description  : This code converts the High-Resolution Deterministic Prediction System (HRDPS) netCDF files to
%                Matlab format and temporally aggregates them to a 3h time step, and spatially averages them to the catchment
%                scale.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
%%                                            THE ONLY PART TO MODIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declarations

% Directories
Dir.dataPath = 'G:\THESE\PROJECT_PHD\Reanalys_Test\HRDPS'                               ; addpath(Dir.dataPath);  % To Modify
Dir.shpPath  = 'G:\THESE\GIS_COUNTRIES\QUEBEC'                                          ; addpath(Dir.shpPath);   % To Modify
OutPath = 'G:\THESE\PROJECT_PHD\Reanalys_Test\HOOPLA\HOOPLA-master\Data\3h\Det_HRDPS_met_fcast';% To Modify

if ~exist(fullfile(OutPath), 'dir')
    mkdir(fullfile(OutPath));         addpath(OutPath);
end

% Catchments name
nameC =  {'023402';'023422';'023429'};   % To Modify
nBV   = numel(nameC);

% Forecasts setting
nbLT = 48; % Lead Time
ts   = 3;  % time step

% Getting the netCDF files name
cd(Dir.dataPath)
ncFiles  = dir('*.nc');
nDays    = size(ncFiles,1);
dateRef  = NaN(nDays,1);

for ifile = 1: nDays
    nameFile(ifile) = convertCharsToStrings(ncFiles(ifile).name);
    tmp =  split(nameFile(ifile), ".nc");
    dateRef(ifile,1) = datenum(datetime(tmp(1),'InputFormat','yyyymmddHH','Format', 'yyyy-MM-dd HH'));
end
dateRef = datevec(dateRef);
dateRef(:,5) = 0;
dateRef = datenum(dateRef);

% Select files from 2019
dateRef  = dateRef(577:end);
nDays    = numel(dateRef);
nameFile = nameFile(577:end);

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
    [S]=shaperead(fullfile(Dir.shpPath,nameC{iCatch},sprintf('%s.shp',nameC{iCatch})));
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

% ---------------------------------------------------------------------------------------------------------------------------
% Save output - catchment-wise - HOOPLA
% ---------------------------------------------------------------------------------------------------------------------------

% Date
Date_LeadTimes = Datehr;
Date = datevec(dateRef);
leadTime = (1:48)./8;

for iCatch=1:nBV
    
    % Extract catchment values
    Pt = PttmpCatch(:,:,iCatch);
    T = TtmpCatch(:,:,iCatch);
    Tmax = TmaxCatch(:,:,iCatch);
    Tmin = TminCatch(:,:,iCatch);
    
   
    % Export
    outfile = sprintf('%s/Met_fcast_%s.mat',OutPath,nameC{iCatch});
    save(outfile,'leadTime','Pt','T','Tmin','Tmax','Date','Date_LeadTimes','-v7.3');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                            END :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
