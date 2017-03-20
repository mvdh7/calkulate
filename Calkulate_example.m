% This script provides an example set-up for using CalkCRM and Calkulate to
%  determine TA from VINDTA titration data.
% It requires the files "calk_fileCRM142.dat", "calk_file_1.dat" and
%  "calk_file_2.dat" to be on the MATLAB search path, along with the
%  functions CalkCRM, and Calkulate.
% If all goes well, then the values returned in the variable "ta" should
%  be: [2154.8; 2153.5].
% Written by Matthew P. Humphreys. Last updated 2016-06-28.

%% CalkCRM: inputs
% Filename from VINDTA:
datfile = 'calk_file_CRM142.dat';
% Approximate acid molarity / mol/l:
acidmolar_approx = 0.1;
% Acid density / kg/l:
acidrho = 1.0212; 
% TA pipette volume / ml:
pipvol = 99.5;
% Dickson CRM batch number:
CRMbatch = 142;
% Titrino pipette volume correction factor / ml-dispensed/ml-recorded:
burette_cxn = 1;
% Leave empty, or input structure specifying all constants' codes:
sc = [];
% Leave empty to use titration file temperature, otherwise provide an
%  over-ride temperature / degC:
tforce = [];

%% CalkCRM: run the program
acidmolar = CalkCRM(datfile,acidmolar_approx,acidrho,pipvol, ...
    CRMbatch,burette_cxn,sc,tforce);

%% Calkulate: inputs
% Note that inputs that are the same as used by CalkCRM above have been
%  commented out here. If CalkCRM is not being used, simply remove the
%  appropriate comment tags.
% Filenames from VINDTA:
datfile = {'calk_file_1.dat'; 'calk_file_2.dat'};
% % Acid molarity / mol/l (if not determined by CalkCRM above):
% acidmolar = 0.1035;
% % Acid density / kg/l:
% acidrho = 1.0212; 
% % TA pipette volume / ml:
% pipvol = 99.5;
% Sample practical salinity:
sal = 33.651 * [1;1];
% Sample DIC / mol/kg-sw
dic = 2026.91 * 1e-6 * [1;1];
% Sample phosphate / mol/kg-sw
phos = 0.42 * 1e-6 * [1;1];
% % Titrino pipette volume correction factor / ml-dispensed/ml-recorded:
% burette_cxn = 1;
% % Leave empty, or input structure specifying all constants' codes:
% sc = [];
% % Leave empty to use titration file temperature, otherwise provide an
% %  over-ride temperature / degC:
% tforce = [];

%% Calkulate: pre-allocate output
% Total alkalinity / mol/kg-sw:
ta = NaN(size(sal));

%% Calkulate: run the program
% To speed up calculations on large datasets by using parallel computing,
%  simply change the "for" in the following line to "parfor":
for M = 1:length(ta)
    ta(M) = Calkulate(datfile{M},acidmolar,acidrho,pipvol, ...
        sal(M),dic(M),phos(M),burette_cxn,sc,tforce);
end %for/parfor M