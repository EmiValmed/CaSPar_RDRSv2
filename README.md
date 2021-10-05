# CaSPar_RDRSv2

This code converts the precipitation and temperature variables from the Regional Deterministic Reforecast System (RDRS_v2; hourly dataset) into Matlab variables at a desired time step and aggregated at the basin scale.

## RDRS-v2 variables in the code:

|**Variable** | **Variable long name**|	**Unit**|	**Level**|
|----------------|-----------------------------------|-----|----|
|RDRS_v2_P_PR0_SFC| Quantity of precipitation (model) | [m] | SFC|
|RDRS_v2_P_TT_1.5m | Air temperature | [°C]	|1.5m|

## Information required to use the code
The following information must be specified in the block comprised by the lines 9-31 (identified as "THE ONLY PART TO MODIFY") :point_down:

```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
%%                                            THE ONLY PART TO MODIFY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Declarations

% Directories
dataPath = 'Path\of\the\RDRS\data'  ; addpath(dataPath);  % To Modify
shpPath  = 'Shapefile\path'         ; addpath(shpPath);   % To Modify
OutPath = 'Output\path'             ;                     % To Modify

if ~exist(fullfile(OutPath), 'dir')
    mkdir(fullfile(OutPath));         addpath(OutPath);
end

% Catchments name
nameC = {'name1';'name2';'name3'; '...nameN'};    % To Modify
nBV   = numel(nameC);

% Defining time step
TS = 3; % time step in hr % To Modify

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
* Catchments' name (**nameC**).  
* Time step (**TS**)
* netCDF files path (**dataPath**)
* Shapefiles path (**shpPath**). Shapefiles must have the name of the catchment (e.g., *name1.shp, name2.shp,...nameN.shp*) 
* Output (**OutPath**)

> You can modify the code to work with other variables :smiley:.
> For more information on the available products visit the [CaSPar](https://github.com/julemai/CaSPAr/wiki/Available-products) website.

## Contact
Please, feel free to contact me if you find any bugs or have any questions. :smiley:
:e-mail:: emixi-sthefany.valdez-medina.1@ulaval.ca 
