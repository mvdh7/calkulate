function [AT,EMF0,AT_RMS,AT_Npts] = calk_VINDTA(datfile,Vsamp,Cacid, ...
    S,CT,PT,SiT,burette_cx,Tk_force)
% Calculate total alkalinity from VINDTA titration data
% For Calkulate v2.0
% Written by Matthew P. Humphreys, last updated 2018-12-20
% 
% === INPUTS ===========
%    datfile: string, titration filename, 'path/to/datfile.dat'
%      Vsamp: sample volume in ml
%      Cacid: HCl acid titrant concentration in mol/kg
%          S: practical salinity
%         CT: dissolved inorganic carbon in mol/kg
%         PT: phosphate in mol/kg
%        SiT: silicate in mol/kg
% burette_cx: HCl burette volume correction factor (use 1 if not needed)
%   Tk_force: overrides .dat file temperature (use [] if not needed) in K
% 
% === OUTPUTS ==========
%         AT: best-fit total alkalinity in mol/kg
%       EMF0: best-fit EMF0 in mV
%     AT_RMS: root-mean-square of total alkalinity estimates in mol/kg
%    AT_Npts: number of total alkalinity estimates in fitted pH range
% 
% -------------------------------------------------------------------------

% Solve for total alkalinity
if ~isempty(Tk_force)
    Vfit = py.calkulate.VINDTA.MPH(datfile,Vsamp,Cacid,S,CT,PT,SiT, ...
        burette_cx,Tk_force);
    
else
    Vfit = py.calkulate.VINDTA.MPH(datfile,Vsamp,Cacid,S,CT,PT,SiT, ...
        burette_cx);
    
end %if else

% Convert Python to MATLAB - AT and EMF0
AT_EMF0 = double(py.array.array('d',Vfit{'x'}));

AT   = AT_EMF0(1);
EMF0 = AT_EMF0(2);

% Convert Python to MATLAB - AT diagnostics
AT_resid = double(py.array.array('d',Vfit{'fun'}));

AT_RMS  = rms  (AT_resid);
AT_Npts = numel(AT_resid);

end %function calk_VINDTA
