%% Import data from text file.
% Script for importing data from the following text file:
%
%    /home/vinith/ICP_SIMON/bun.conf
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2014/11/18 11:34:28

%% Initialize variables.
filename = '/home/vinith/ICP_SIMON/bun.conf';
delimiter = ' ';

%% Format string for each line of text:
%   column2: text (%s)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
%	column7: double (%f)
%   column8: double (%f)
%	column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
dataArray([2, 3, 4, 5, 6, 7, 8]) = cellfun(@(x) num2cell(x), dataArray([2, 3, 4, 5, 6, 7, 8]), 'UniformOutput', false);
bun = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;