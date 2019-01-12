function [Macid,pH,Tk,Msamp,Cacid,S,XT,KX] = calk_Dickson1981
% Import simulated titration from Dickson (1981) Table 1,
%  doi:10.1016/0198-0149(81)90121-7
% 
% === OUTPUTS ==========
% Macid: acid titrant mass in kg
%    pH: Free scale pH
%    Tk: temperature in K
% Msamp: sample mass in kg
% Cacid: acid concentration in mol/kg
%     S: salinity
%    XT: total concentrations, in mol/kg-sw, of:
%        (1) AT  total alkalinity
%        (2) CT  dissolved inorganic carbon
%        (3) BT  borate
%        (4) ST  sulfate
%        (5) FT  fluoride
%        (6) PT  phosphate
%        (7) SiT silicate
%    KX: dissociation constants, on the Free pH scale, for:
%        ( 1) KC1   carbonic acid, first
%        ( 2) KC2   carbonic acid, second
%        ( 3) KB    boric acid
%        ( 4) KH2O  water
%        ( 5) KHSO4 bisulfate
%        ( 6) KHF   hydrogen fluoride
%        ( 7) KP1   phosphoric acid, first
%        ( 8) KP2   phosphoric acid, second
%        ( 9) KP3   phosphoric acid, third
%       [(10) KSi   silicic acid - not used by Dickson (1981)]

% Import dataset
Dickson1981 = cell(py.calkulate.gettit.Dickson1981());

% Convert from Python to MATLAB
Macid = double(py.array.array('d',Dickson1981{1}));
pH    = double(py.array.array('d',Dickson1981{2}));
Tk    = double(py.array.array('d',Dickson1981{3}));
Msamp = Dickson1981{4};
Cacid = Dickson1981{5};
S     = Dickson1981{6};

XT = NaN(7,1);
for T = 1:numel(XT)
    XT(T) = Dickson1981{7}{T};
end %for K

KX = NaN(10,1);
for K = 1:numel(KX)
    KKX = double(py.array.array('d',Dickson1981{8}{K}));
    KX(K) = KKX(1);
end %for K
KX(10) = NaN;

end %function calk_Dickson1981
