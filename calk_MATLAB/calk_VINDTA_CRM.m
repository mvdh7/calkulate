function Cacid = calk_VINDTA_CRM(datfile,Vsamp,AT_cert,S,CT,PT,SiT, ...
    burette_cx,Tk_force)
% Get HCl acid titrant concentration from CRM titration with known AT
% For Calkulate v2.0
% Written by Matthew P. Humphreys, last updated 2018-12-20
% 
% === INPUTS ===========
%    datfile: string, titration filename, 'path/to/datfile.dat'
%      Vsamp: sample volume in ml
%    AT_cert: certified total alkalinity in mol/kg
%          S: practical salinity
%         CT: dissolved inorganic carbon in mol/kg
%         PT: phosphate in mol/kg
%        SiT: silicate  in mol/kg
% burette_cx: HCl burette volume correction factor (use 1 if unsure)
%   Tk_force: overrides .dat file temperature (use [] if not needed) in K
% 
% === OUTPUT ===========
%      Cacid: best-fit HCl concentration in mol/kg
% 
% -------------------------------------------------------------------------

% Solve for acid concentration
if ~isempty(Tk_force)
    Cacid = py.calkulate.VINDTA.MPH_CRM(datfile,Vsamp,AT_cert, ...
        S,CT,PT,SiT,burette_cx,Tk_force);

else
    Cacid = py.calkulate.VINDTA.MPH_CRM(datfile,Vsamp,AT_cert, ...
        S,CT,PT,SiT,burette_cx);

end %if else

% Convert Python to MATLAB
Cacid = double(py.array.array('d',Cacid(1)));
    
end %function Cacid
