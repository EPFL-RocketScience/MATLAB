clear all
close all

% open file
fid = fopen('TEST', 'rb');

% read file (PRECISION = 'real*4' : float 32bits, 4 bytes)
%           (SIZE = [n inf] : n*4 bytes per time step of measures, inf -> read until end of file)
data = fread(fid, [8 inf], 'real*4')'; 

% close file
fclose(fid);