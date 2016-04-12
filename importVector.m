function T = importVector(filename, startRow, endRow, nCols)
% IMPORTVECTOR Import numeric data from a text file as a matrix.
%
% USAGE:-------------------------------------------------------------------
% T = importVector(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
% T = importVector(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
%   Example:
%   T = importVector('DeploymentOne6_25to6_50hr.dat');
%
% DESCRIPTION:-------------------------------------------------------------
% Parse Nortek ADV '.dat' file into MATLAB Table format
%
% INPUTS:------------------------------------------------------------------
%
% REQUIRED:
% filename      string of full path of .dat file or name of file in path
% 
% OPTIONAL:
% startRow:     starting row number to read in
% endRow:       final row number to read in
% nCols:        number of data columns
%
% OUTPUTS:-----------------------------------------------------------------
% 
% T             Matlab Table of parsed data
%
% SEE ALSO:----------------------------------------------------------------
%
% motionCorrectVector.m
% textscan.m
%
% AUTHOR:------------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
%
% REFERENCE:---------------------------------------------------------------
% Long, M.H. and Nicholson, D.P (2016) Air-Sea Gas Exchange Determined from 
% an Aquatic Eddy Covariance Floating Platform. J. Geophys Res. Oceans
% ADD DOI HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end
if nargin < 4
    nCols = 38;
    varNames = {'BurstCount','Ens_Count65536','V_X','V_Y','V_Z','Amp_B1_count','Amp_B2_count','Amp_B3_count','SNR_B1','SNR_B2','SNR_B3','Corr_B1','CorrB2','Corr_B3','P_db','Alog1','Alog2','Chsum1','Ens_Count','Timers','Accel_X','Accel_Y','Accel_Z','Ang_X','Ang_Y','Ang_Z','Mag_X','Mag_Y','Mag_Z','OM_11','OM_12','OM_13','OM_21','OM_22','OM_23','OM_31','OM_32','OM_33'};
end
    

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.
formatSpec = [repmat('%f',[1,nCols]) '%[^\n\r]'];
%% Open the text file.
fileID = fopen(filename,'r');


%% Read header and assign variable names
% read first line
l = fgetl(fileID);
% if column names are not specified, they will be auto-generated from
% header line
if ~exist('varNames','var')
    varNames =  matlab.lang.makeValidName(strsplit(l,'\t'));
end
%

%% Read columns of data according to format string.

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', max(startRow(1)-2,0), 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
T = table(dataArray{1:end-1}, 'VariableNames', varNames);
