function [TA_final,E0_final,t,f2xint,meanT,pipmass,sc,ta,e0] ...
    = Calkulate(datfile,acidmolar,acidrho,pipvol, ...
        sal,dic,phos,burette_cxn,sc,tforce,suppress_warnings)
%Calkulate Calculates total alkalinity from VINDTA .dat file output
% Matthew P. Humphreys, 2015-01-20. Last updated 2016-07-07 [v1.0.2]
% This calculates total alkalinity (TA) from the .dat files output by the
%  VINDTA instrument (Marianda, Kiel, Germany). The method is based on the
%  Gran plot (G52) as modified by Hansson & Jagner (HJ73) and Bradshaw et
%  al. (BBSW81). However, DIC takes a value that is not permitted to change
%  during the iterative process.
% Use function CalkCRM to generate <acidmolar> input from CRM data.
% See end of main function for citation list & reference codes.
% === Please cite as: ===
% Humphreys, M. P. (2015). Calculating seawater total alkalinity from open-
%  cell titration data using a modified Gran plot technique, in: Measurements
%  and Concepts in Marine Carbonate Chemistry. PhD Thesis, Ocean and Earth
%  Science, University of Southampton, UK, pp. 25-44

%% === 0 === USER INPUTS and example values === 0 ===
% datfile = '0-0  1  (0)CRM-134-1085-A.dat'; % 'filename.dat' from VINDTA
% acidmolar = 0.10496; % Acid molarity / mol/L
% acidrho = 1.0212; % Acid density / kg/L
% pipvol = 99.995; % TA pipette volume / mL
% sal = 33.384; % Sample salinity
% dic = 2022.04 * 1e-6; % Sample DIC / mol/kg
% phos = 0.48 * 1e-6; % Sample phosphate / mol/kg
% burette_cxn = 4.9831/5; % Titrino pipette volume correction factor
% sc = []; % Leave empty or input structure containing all constants' codes
% tforce = []; % Leave empty to use titration file temperature
% suppress_warnings is an optional input used by CalkCRM. Ignore
max_iterations = 21; % Maximum number of iterations to run

%% === 1 === SELECT CONSTANTS === 1 ===
% The options are -
%  Ionic concentrations
%   Total sulfate, TSO4: MR66,KDCP67
%   Total fluoride, THF: W71,R65,KDCP67
%   Total boron, TB: LKB10,U74,KDCP67
%  Dissociation constants
%   Sulfate, kHSO4: D90b,WM13
%   Fluoride, kHF: PF87,DR79
%   Carbonic acid, k1k2: LDK00,GP89

% Sanitise <sc> inputs, and fill missing or bad fields with defaults
sc_fields   = {'TSO4' 'THF' 'TB' 'kHSO4' 'kHF' 'k1k2'};
sc_options  = {{'MR66' 'KDCP67'} % TSO4
               {'W71' 'R65' 'KDCP67'} % THF
               {'LKB10' 'U74' 'KDCP67'} % TB
               {'D90b' 'WM13'} % kHSO4
               {'PF87' 'DR79'} % kHF
               {'LDK00' 'GP89'}}; % k1k2
for F = 1:length(sc_fields) % F for Field
    if ~isfield(sc,sc_fields{F})
        sc.(sc_fields{F}) = sc_options{F}{1};
    end %if
    if ~ismember(sc.(sc_fields{F}),sc_options{F})
        disp(['-CALK- User input ''' sc.(sc_fields{F}) ''' for sc.' ...
            sc_fields{F} ' not recognised; reverting to default ''' ...
            sc_options{F}{1} '''']);
        sc.(sc_fields{F}) = sc_options{F}{1};
    end %if
end %for F

%% === 2 === READ VINDTA OUTPUT: .dat FILE === 2 ===
% This gives the option of inputting the titration data directly as a table
if ischar(datfile)
    t = calk_datfile(datfile);
else
    t = datfile;
end %if

% Correct cell temperature to <tforce>, if required
if ~isempty(tforce)
    t.t(:) = tforce;
end %if

% Initial .dat file calculations
if isempty(burette_cxn), burette_cxn = 1; end %if
t.vol = t.vol_raw * burette_cxn; % Volume acid added / mL
t.tk = t.t + 273.15; % Cell temperature / K
t.mass = t.vol*1e-3 * acidrho; % Acid added / kg

% Seawater sample analysis density (MP81) & associated calculations
piprho = calk_MP81(t.t(1),sal) * 1e-3; % Sample analysis density / kg/L
pipmass = pipvol * piprho * 1e-3; % Sample mass / kg
acidmolal = acidmolar/acidrho; % Acid molality / mol/kg

%% === 3 === NERNST EQUATION CONSTANTS === 3 ===
% Constants R & F are the 2010 CODATA recommended values
%  from http://physics.nist.gov/cuu/Constants
R =  8.3144621; % Gas constant / J/mol/K
F = 96.4853365; % Faraday constant / kC/mol
t.nernst = R * t.tk / F; % J/kC

%% === 4 === IONIC CONCENTRATIONS === 4 ===
% Total Sulfate / mol/kg
switch sc.TSO4
    case 'MR66',   TSO4 = calk_MR66(sal); % MR66 (recommended by DSC07)
    case 'KDCP67', TSO4 = 4.008 / 142.041; % KDCP67 (synthetic seawater)
end %switch
% Total Fluoride / mol/kg
switch sc.THF
    case 'W71',    THF = calk_W71(sal); % W71
    case 'R65',    THF = calk_R65(sal); % R65 (recommended by DSC07)
    case 'KDCP67', THF = 0.003 / 41.988; % KDCP67 (synthetic seawater)
end %switch
% Total Boron / mol/kg
switch sc.TB
    case 'LKB10',  TB = calk_LKB10(sal); % LKB10
    case 'U74',    TB = calk_U74(sal); % U74
    case 'KDCP67', TB = 0.026 / 61.8317; % KDCP67 (synthetic seawater)
end %switch

%% === 5 === DISSOCIATION CONSTANTS == 5 ===
% All converted to the Free pH scale
% Bisulfate kHSO4
switch sc.kHSO4
    case 'D90b', [~,pkHSO4] = calk_D90(t.tk,sal); % D90b (rec. DSC07)
    case 'WM13', pkHSO4 = calk_WM13(t.tk,sal); % WM13
end %switch
t.kHSO4 = 10.^-pkHSO4;
t.t2f = log10(1 + TSO4./t.kHSO4); % pH scale conversion: Total to Free
% Fluoride kHF
switch sc.kHF
    case 'PF87', pkHF_T = calk_PF87(t.tk,sal); % PF87 (rec. DSC07)
        t.kHF = 10.^-(pkHF_T + t.t2f);
    case 'DR79', pkHF = calk_DR79(t.tk,sal); % DR79
        t.kHF = 10.^-pkHF;
end %switch
t.s2f = log10(1 + TSO4./t.kHSO4 + THF./t.kHF); % pH conv.: SWS to Free
% Carbonic acid K1 and K2
switch sc.k1k2
    case 'LDK00', [pk1_T,pk2_T] = calk_LDK00(t.tk,sal); % LDK00 (rec.DSC07)
        t.k1 = 10.^-(pk1_T + t.t2f);
        t.k2 = 10.^-(pk2_T + t.t2f);
    case 'GP89', [pk1_SWS,pk2_SWS] = calk_GP89(t.tk,sal); % GP89
        t.k1 = 10.^-(pk1_SWS + t.s2f);
        t.k2 = 10.^-(pk2_SWS + t.s2f);
end %switch
% Boric acid (D90a)
pKB = calk_D90(t.tk,sal); % D90a (rec. DSC07)
t.kB = 10.^-(pKB + t.t2f);

% Water and phosphoric acid (KP67,DR79,JW79,M95,DSC07)
[pkw_T,pkP1_T,pkP2_T,pkP3_T,pkSi_T] = calk_DSC07(t.tk,sal);
t.kw = 10.^-(pkw_T + t.t2f);
t.kP1 = 10.^-(pkP1_T + t.t2f);
t.kP2 = 10.^-(pkP2_T + t.t2f);
t.kP3 = 10.^-(pkP3_T + t.t2f);

% Silicic acid (BM71,SNI81,M95,DSC07)
t.kSi = 10.^-(pkSi_T + t.t2f);

% Dilution correction factor
t.masscxn = pipmass ./ (pipmass + t.mass);

%% === 6 === INITIAL ESTIMATES for TA and E0 === 6 ===
% Originally from G52, used by HJ73 and BBSW81, written here from
%  DAA03 eq'ns (10) and (11). NB: DAA03's "F1" is the others' "F2". Here,
%  I use the label G.
t.G = (pipmass + t.mass) .* exp(t.emf ./ t.nernst); % DAA03 eq'n (10)
t.Glogic = t.G > 0.1*max(t.G); % This may need to be adjusted. It
% should catch the +ve linear section on the right hand side of 
% scatter(t.vol,t.G). This value is only used the first time, 
% not for subsequent iterations.

%% === 7 === ITERATION PRE-ALLOCATIONS === 7 ===
ta = NaN(max_iterations,1);
e0 = NaN(size(ta));
I = 1;
icomplete = false;

%% === 8 === ITERATION BEGINS === 8 ===
% for I = 1:length(ta)
while I <= length(ta)
    
f2fit = regstats(t.G(t.Glogic(:,I),I),t.mass(t.Glogic(:,I)), ...
    'linear','beta');
f2xint = -f2fit.beta(1)/f2fit.beta(2);
ta(I) = f2xint * acidmolal/pipmass; % TA / mol/kg
% DAA03 eq'n (11):
if I == 1
    t.e0(:,I) = t.emf - t.nernst .* log((t.mass-f2xint)*acidmolal ...
        ./ (pipmass+t.mass));
else
    t.e0(:,I) = t.emf - t.nernst .* log(((t.mass-f2xint)*acidmolal ...
        - pipmass*(t.hf(:,I-1) + t.hso4(:,I-1))) ./ (pipmass+t.mass));
end %if else

t.e0(~t.Glogic(:,I),I) = NaN;
e0(I) = nanmean(t.e0(:,I)); % E0 / mV

%% == 8.1 == CHECK DIFFERENCE from PREVIOUS ITERATION == 8.1 ==
PPC = 5e-3; % permitted % change
if I > 2
    if abs(100*(ta(I) - ta(I-1))/nanmean(ta)) < PPC ...
    || abs(100*(ta(I) - ta(I-2))/nanmean(ta)) < PPC
        icomplete = true;
        break
    end %if
end %if

%% == 8.2 == pH CALCULATION == 8.2 ==
% Modified Gran calculation from BBSW81's Appendix 2: based on HJ73, but 
%  with phosphate & silicate effects added in.
t.h(:,I) = 10.^((t.emf - e0(I))./(log(10)*t.nernst)); % HJ73 eq'n (19)
t.ph(:,I) = -log10(t.h(:,I)); % Free scale

%% == 8.3 == OTHER PROTON EQUATIONS == 8.3 ==
% As for pH calculation, based on HJ73 and BBSW81
% (Silicate only affects the DIC determination, so is not included)
t.dic(:,I) = dic;
t.bicarb(:,I) = t.masscxn.*t.dic(:,I)./(t.h(:,I)./t.k1 + 1); % HJ73 eq. 20
t.hso4(:,I) = t.masscxn.*TSO4./(1 + t.kHSO4./t.h(:,I));      % HJ73 eq. 21
t.hf(:,I) = t.masscxn.*THF./(1 + t.kHF./t.h(:,I));           % HJ73 eq. 22
t.borate(:,I) = t.masscxn.*TB./(t.h(:,I)./t.kB + 1);         % HJ73 eq. 23
t.oh(:,I) = t.kw ./ t.h(:,I);                                % HJ73 eq. 24
% From BBSW81 Appendix 2:
t.p_p2(:,I) = t.masscxn.*phos .* (1 - t.kP1.*t.kP2./(t.h(:,I).^2)) ...
    ./ (1 + t.kP1./t.h(:,I) ...
    + t.kP2.*t.kP3./t.h(:,I).^2 + t.kP1.*t.kP2.*t.kP3./t.h(:,I).^3);
% NB: [PO43-] is negligible in F2 pH range (BBSW81 Appendix 2)

%% == 8.4 == EXTEND TABLE on first loop == 8.4 ==
if I == 1
    itvars = {'G' 'Glogic' 'e0' 'h' 'ph' 'dic' 'bicarb' 'hso4' ...
        'hf' 'borate' 'oh' 'p_p2'};
    for V = 1:length(itvars)
        t.(itvars{V}) = [t.(itvars{V}) NaN(height(t),length(ta)-1)];
    end %for V
    t.Glogic(isnan(t.Glogic)) = 0;
    t.Glogic = logical(t.Glogic);
end %if
    
%% == 8.5 == CALCULATE F2, GET NEXT ESTIMATES of TA and E0 == 8.5 ==
% HJ73 eq'n (16), including phosphate like BBSW81 Appendix 2, but
%  remaining in the pH Free scale, and in mass not volume:
if I < length(ta)
    t.G(:,I+1) = (pipmass + t.mass).*(t.h(:,I) + t.hso4(:,I) ...
        + t.hf(:,I) - t.bicarb(:,I) + t.p_p2(:,I) - t.oh(:,I) ...
        - t.borate(:,I));
    t.Glogic(:,I+1) = t.ph(:,I) < 4 & t.ph(:,I) > 3;
end %if

%% == 8.6 == EVALUATION OF F1 and CALCULATION OF DIC == 8.6 ==
% We do not evaluate F1 or iterate a value for DIC, because:
%  (1) It should not affect the TA value
%  (2) Titration is open cell, so DIC result is not meaningful
%  (3) There should always be an independently-measured DIC value to input
% So DIC is instead held constant at this input value.
% Experimenting with input <dic> suggests that it does not significantly
%  affect the TA result: changing DIC from 2000 to 1000 micromol/kg changes
%  TA by only 2 micromol/kg. Permitting smaller changes in DIC during the
%  calculation therefore should make only negligible changes to TA (i.e. 
%  much less than the measurement precision).

I = I + 1;
% === ITERATION ENDS ===
end %while

%% === 9 === GET FINAL OUTPUTS === 9 ===
% Remove NaNs
FL = ~isnan(ta);
ta = ta(FL);
e0 = e0(FL);

for V = 1:length(itvars)
    t.(itvars{V}) = t.(itvars{V})(:,FL);
end %for V

% Get final values
TA_final = mean(abs(ta(end-1:end))) * 1e6; % micromol/kg
E0_final = mean(abs(e0(end-1:end))); % mV
meanT = mean(t.t); % degC

% Check for non-convergence
if nargin == 10
    suppress_warnings = false;
end %if

if ~icomplete && ~suppress_warnings
    disp(['Calkulate: WARNING! TA result did not converge for sample "' ...
        datfile '"']);
end %if

%% === 10 === REFERENCES === 10 ===
%    G52: Gran, 1952, Analyst 77. doi:10.1039/AN9527700661
%    R65: Riley, 1965, Deep-Sea Res 12(2). doi:10.1016/0011-7471(65)90027-6
%   MR66: Morris & Riley, 1966, Deep-Sea Res 13(4). doi:10.1016/0011-7471(66)90601-2
% KDCP67: Kester et al., 1967, Limnol Oceanogr 12(1). doi:10.4319/lo.1967.12.1.0176
%   KP67: Kester & Pytkowicz, 1967, Limnol Oceanogr 12(2). doi:10.4319/lo.1967.12.2.0243
%  WLD69: Wooster et al., 1969, Limnol Oceanogr 14(3). doi:10.4319/lo.1969.14.3.0437
%    W71: Warner, 1971, Deep-Sea Res 18(12). doi:10.1016/0011-7471(71)90030-1
%   HJ73: Hansson & Jagner, 1973, Anal Chim Acta 65(2). doi:10.1016/S0003-2670(01)82503-4
%    U74: Uppström, 1974, Deep-Sea Res 21(2). doi:10.1016/0011-7471(74)90074-6
%   BM76: Baes & Mesmer, 1976, The Hydrolysis of Cations. Wiley
%   DR79: Dickson & Riley, 1979, Mar Chem 7(2). doi:10.1016/0304-4203(79)90002-1
%   JW79: Johansson & Wedborg, 1979, Mar Chem 8(1). doi:10.1016/0304-4203(79)90032-X
% BBSW81: Bradshaw et al., 1981, Earth Planet Sci Lett 55(1). doi:10.1016/0012-821X(81)90090-X
%   MP81: Millero & Poisson, 1981, Deep-Sea Res Pt A 28(6). doi:10.1016/0198-0149(81)90122-9
%  SNI81: Sjöberg et al., 1981, Mar Chem 10(6). doi:10.1016/0304-4203(81)90005-0
%   PF87: Perez & Fraga, 1987, Mar Chem 21(2). doi:10.1016/0304-4203(87)90036-3
%   GP89: Goyet & Poisson, 1989, Deep-Sea Res Pt A 36(11). doi:10.1016/0198-0149(89)90064-2
%   D90a: Dickson, 1990, Deep-Sea Res Pt A 37(5). doi:10.1016/0198-0149(90)90004-F
%   D90b: Dickson, 1990, J Chem Thermodyn 22(2). doi:10.1016/0021-9614(90)90074-Z
%    B92: Butler, 1992, Mar Chem 38(3-4). doi:10.1016/0304-4203(92)90037-B
%  CRP94: Clegg et al., 1994, J Chem Soc F Trans 90(13). doi:10.1039/FT9949001875
%    M95: Millero, 1995, Geochim Cosmochim Acta 59(4). doi:10.1016/0016-7037(94)00354-O
%  LDK00: Lueker et al., 2000, Mar Chem 70(1-3). doi:10.1016/S0304-4203(00)00022-0
%  DAA03: Dickson et al., 2003, Mar Chem 80(2-3). doi:10.1016/S0304-4203(02)00133-0
%  DSC07: Dickson et al., 2007, Guide to best practices for CO2 measurements. PICES Special Publication 3
%  LKB10: Lee et al., 2010, Geochim Cosmochim Acta 74(6). doi:10.1016/j.gca.2009.12.027
%   WM13: Waters & Millero, 2013, Mar Chem 149. doi:10.1016/j.marchem.2012.11.003

%% === 11 === VERSION INFORMATION === 11 ===
% v1 [2016-06-28]
% v1.0.1 [2016-07-07]
% ------ Added <suppress_warnings> input (optional, used by CalkCRM)
% ------ Minor corrections to CalkCRM
% v1.0.2 [2016-07-07]
% ------ Technical adjustments to convergence testing mechanisms in both
%        Calkulate and CalkCRM

end %function Calkulate


function t = calk_datfile(filename)
%calk_datfile Load VINDTA titration data from .dat file
% Input <filename> = 'filename.dat'; output <t> = titration data

% Load titration data from .dat file
fileid = fopen(filename,'r');
tdata = textscan(fileid,'%f%f%f', 'headerlines',2);
fclose(fileid);

% Extract titration variables
t.vol_raw  = tdata{1}; % Volume of acid added / mL
t.emf      = tdata{2}; % Electrode EMF        / mV
t.t        = tdata{3}; % Temperature          / deg C
t = struct2table(t);

end %function calk_datfile


function TSO4 = calk_MR66(sal)
%calk_MR66 Total sulfate in seawater / mol/kg
% Input <sal> = salinity; output <TSO4> = total sulfate / mol/kg
% Total sulfate = [sulfate(2-)] + [bisulfate(-)]

% Sulfate:chlorinity (MR66)
SO4_Cl = 0.14000;

% Atomic mass (DSC07)
SO4_mass = 32.065 + 4*15.999; % g/mol

% Salinity:chlorinity (WLD69,DSC07)
TSO4 = (SO4_Cl/SO4_mass) * sal/1.80655; % mol/kg

end %function calk_MR66

function TF = calk_R65(sal)
%calk_R65 Total fluorine <TF> in seawater / mol/kg
% Input <sal> = salinity; output <TF> = total fluorine / mol/kg
% Total fluorine = [fluoride(-)] + [hydrogen fluoride]

% Fluorine:chlorinity (R65)
F_Cl = 6.7e-5;

% Atomic mass (DSC07)
F_mass = 18.998; % g/mol

% Salinity:chlorinity (WLD69,DSC07)
TF = (F_Cl/F_mass) * sal/1.80655; % mol/kg

end %function calk_W71

function TF = calk_W71(sal)
%calk_W71 Total fluorine <TF> in seawater / mol/kg
% Input <sal> = salinity; output <TF> = total fluorine / mol/kg
% Total fluorine = [fluoride(-)] + [hydrogen fluoride]

% Fluorine:chlorinity (W71)
F_Cl = 6.75e-5;

% Atomic mass (DSC07)
F_mass = 18.998; % g/mol

% Salinity:chlorinity (WLD69,DSC07)
TF = (F_Cl/F_mass) * sal/1.80655; % mol/kg

end %function calk_W71

function TB = calk_U74(sal)
%calk_U74 Total boron <TB> in seawater / mol/kg
% Input <sal> = salinity; output <TB> = total boron / mol/kg
% Total boron = [boric acid] + [borate(-)]

% Boron:chlorinity (U74)
B_Cl = 0.232;

% Atomic mass (DSC07)
B_mass = 10.811; % g/mol

% Salinity:chlorinity (WLD69,DSC07)
TB = (1e-3*B_Cl/B_mass) * sal/1.80655; % mol/kg

end %function calk_U74

function pkHF = calk_DR79(tk,sal)
%calk_DR79 HF stoichiometric dissociation constant <pkHF>, pH-Free
% Reaction: HF <=> H{+} + F{-}
% Inputs <tk> = temperature / K, <sal> = salinity
% Output <pkHF> = stoich. diss. constant, Free pH scale

% Ionic strength (DSC07)
ionstr = 19.924 * sal ./ (1000 - 1.005 * sal);

% Calculate pkHF (DR79)
ln_kF = 1590.2./tk - 12.641 + 1.525*sqrt(ionstr) + log(1 - 0.001005*sal);
pkHF = -log10(exp(ln_kF));

end %function calk_DR79

function rho = calk_MP81(t,sal)
%calk_MP81 Seawater density <rho> / kg/m^3 at pressure = 1 atmosphere
% Inputs <t> = temperature / deg C, <sal> = salinity
% Output <rho> = density / kg/m^3

% Calculate density (MP81)
rho  = 999.842594                        ...
     +   6.793952e-2 * t                 ...
     -   9.095290e-3 * t.^2              ...
     +   1.001685e-4 * t.^3              ...
     -   1.120083e-6 * t.^4              ...
     +   6.536336e-9 * t.^5              ...
 + (     0.824493                        ...
     -   4.0899e-3   * t                 ...
     +   7.6438e-5   * t.^2              ...
     -   8.2467e-7   * t.^3              ...
     +   5.3875e-9   * t.^4 ) .*sal      ...
 + ( -   5.72466e-3                      ...
     +   1.0227e-4   * t                 ...
     -   1.6546e-6   * t.^2 ) .*sal.^1.5 ...
 +       4.8314e-4            .*sal.^2   ;

end %function calk_MP81

function pkHF_T = calk_PF87(tk,sal)
%calk_PF87 HF stoichiometric dissociation constant <pkHF>, pH-Total
% Reaction: HF <=> H{+} + F{-}
% Inputs <tk> = temperature / K, <sal> = salinity
% Output <pkHF> = stoich. diss. constant, Total pH scale

% Calculate pkHF (PF87)
ln_bHF = -874./tk - 0.111.*sal.^0.5 + 9.68;
ln_kHF = -ln_bHF;
pkHF_T = -log10(exp(ln_kHF));

end %function calk_PF87

function [pk1_SWS,pk2_SWS] = calk_GP89(tk,sal)
%calk_GP89 Carbonic acid stoichiometric dissociation constants, pH-SWS
% Reaction: CO2 + H2O <=> HCO3{-} + H{+} <=> CO3{2-} + 2H{+}
% Inputs <tk> = temperature / K, <sal> = salinity
% Outputs <pk1_SWS>, <pk2_SWS> = stoich. diss. constants, SWS pH scale

% Calculate constants (GP89)
pk1_SWS =  812.27./tk + 3.356 - 0.00171*sal.*log(tk) + 0.000091*sal.^2;
pk2_SWS = 1450.87./tk + 4.604 - 0.00385*sal.*log(tk) + 0.000182*sal.^2;

end %function calk_GP89

function [pkB_T,pkHSO4] = calk_D90(tk,sal)
%calk_D90 Boric acid (pH-Total) & bisulfate (pH-Free) eq'm constants
% Inputs <tk> = temperature / K, <sal> = salinity
% Outputs <pkB_T>, <pkHSO4> = equilibrium constants

%% BORIC ACID - Total pH scale (D90a)
% Reaction: B(OH)3 + H2O <=> B(OH)4{-} + H{+}
% D90a, equation 23
ln_KB = (-8966.90   -  2890.53  *sal.^0.5 - 77.942   *sal ...
        + 1.728*sal.^1.5 - 0.0996*sal.^2) ./ tk ...
     +    148.0248 +   137.1942*sal.^0.5 +   1.62142*sal ...
     - (   24.4344 +    25.085 *sal.^0.5 +   0.2474 *sal) .* log(tk) ...
     +      0.053105*sal.^0.5 .* tk;
pkB_T = -log10(exp(ln_KB));

%% BISULFATE - Free pH scale (D90b)
% Reaction: HSO4{-} <=> SO4{2-} + H{+}
% D90b is paywalled, so equations here are as reprinted by DSC07
ionstr = 19.942*sal ./ (1000 - 1.005*sal); % Ionic strength
ln_KS = -4276.1./tk + 141.328 - 23.093*log(tk) ...
    + ((-13856./tk) + 324.57 - 47.986*log(tk)) .* ionstr.^0.5 ...
    + ((35474./tk) - 771.54 + 114.723*log(tk)) .* ionstr ...
    - (2698./tk) .* ionstr.^1.5 ...
    + (1776./tk) .* ionstr.^2 ...
    + log(1 - 0.001005*sal);
pkHSO4 = -log10(exp(ln_KS));

end %function calk_D90

function [pk1_T,pk2_T] = calk_LDK00(tk,sal)
%calk_LDK00 Carbonic acid stoichiometric dissociation constants, pH-Total
% Reaction: CO2 + H2O <=> HCO3{-} + H{+} <=> CO3{2-} + 2H{+}
% Inputs <tk> = temperature / K, <sal> = salinity
% Outputs <pk1_T>, <pk2_T> = stoich. diss. constants, Total pH scale

% Calculate constants (LDK00)
pk1_T = 3633.86./tk - 61.2172 + 9.6777.*log(tk) ...
    - 0.011555*sal + 0.0001152*sal.^2;
pk2_T = 471.78./tk + 25.929 - 3.16967.*log(tk) ...
    - 0.01781*sal + 0.0001122*sal.^2;

end %function calk_LDK00

function [pkw_T,pkP1_T,pkP2_T,pkP3_T,pkSi_T] = calk_DSC07(tk,sal)
%calk_DSC07 Water & phosphoric acid dissociation constants, pH-Total
% Inputs <tk> = temperature / K, <sal> = salinity
% Outputs <pkw_T>, <pkP1_T>, <pkP2_T>, <pkP3_T>, <pkSi_T>
%   = stoich. diss. constants, Total pH scale

% These equations are all taken from DSC07, but they are originally from
%  M95; DSC07 have subtracted 0.015 from the constant term in each case,
%  to approximately convert pH-SWS to pH-Total.

%% WATER (DSC07)
% kw = [H+][OH-]
ln_kw = (-13847.26./tk) + 148.9652 - 23.6521*log(tk) ...
    + ((118.67./tk) - 5.977 + 1.0495*log(tk)) .* sal.^0.5 ...
    - 0.01615 * sal;
pkw_T = -log10(exp(ln_kw));

%% PHOSPHORIC ACID
% These equations, first presented by M95, are based on a composite of data
%  from KP67, DR79 and JW79. NB: the equations below are from DSC07.
% kP1 = [H+][H2PO4-]/[H3PO4]
ln_kP1 = -4576.752./tk + 115.525 - 18.453*log(tk) ...
    + (-106.736./tk + 0.69171) .* sal.^0.5 ...
    + (-0.65643./tk - 0.01844) .* sal;
pkP1_T = -log10(exp(ln_kP1));
% kP2 = [H+][HPO42-]/[H2PO4-]
ln_kP2 = -8814.715./tk + 172.0883 - 27.927*log(tk) ...
    + (-160.34./tk + 1.3566) .* sal.^0.5 ...
    + (0.37335./tk - 0.05778) .* sal;
pkP2_T = -log10(exp(ln_kP2));
% kP3 = [H+][PO43-]/[HPO42-]
ln_kP3 = -3070.75./tk - 18.141 ...
    + (17.27039./tk + 2.81197) .* sal.^0.5 ...
    + (-44.99486./tk - 0.09984) .* sal;
pkP3_T = -log10(exp(ln_kP3));

%% SILICIC ACID
% Ionic strength <I> is first calculated following DSC07:
I = 19.924*sal / (1000 - 1.005*sal);
% The following equation, first presented by M95, is then based on a 
%  composite of data from BM76 and SNI81. NB: the equation below is taken
%  from DSC07.
% kSi = [H+][SiO(OH)3-]/[Si(OH)4]
ln_kSi = -8904.2./tk + 117.385 - 19.334*log(tk) ...
    + (-458.79./tk + 3.5913) * I^0.5 ...
    + (188.74./tk - 1.5998) * I ...
    + (-12.1652./tk + 0.07871) * I^2 ...
    + log(1 - 0.001005*sal);
pkSi_T = -log10(exp(ln_kSi));

end %function calk_DSC07

function TB = calk_LKB10(sal)
%calk_LKB10 Total boron <TB> in seawater / mol/kg
% Input <sal> = salinity; output <TB> = total boron / mol/kg
% Total boron = [boric acid] + [borate(-)]

B_sal = 0.1336;     % = boron/salinity / (mg/kg)/ (LKB10)
B = B_sal * sal;    % = boron / mg/kg
B_RAM = 10.811e3;   % = relative atomic mass of boron / mg/mol (DSC07)
TB = B/B_RAM;       % = total boron / mol/kg

end %function calk_LKB10

function pkHSO4 = calk_WM13(tk,sal)
%calk_WM13 Bisulfate dissociation constant, pH-Free
% Reaction: HSO4{-} <=> SO4{2-} + H{+}
% Inputs <tk> = temperature / K, <sal> = salinity
% Outputs <pkB_T>, <pkHSO4> = equilibrium constants

%% Equation 29
% Constants (from Corrigendum to WM13, Table 6)
c1  =  4.24666         ;
c2  = -0.152671        ;
c3  =  2.67059   * 1e-2;
c4  = -4.2128    * 1e-5;
c5  =  0.2542181       ;
c6  = -5.09534   * 1e-3;
c7  =  7.1589    * 1e-4;
c8  = -2.91179   * 1e-3;
c9  =  2.09968   * 1e-5;
c10 = -4.03724   * 1e-5;
% Equation
log10_KK = (c1 + c2*tk + c3*tk.*log(tk) + c4*tk.^2) .* sal.^0.5 ...
        + (c5 + c6*tk + c7*tk.*log(tk))             .* sal      ...
        + (c8 + c9*tk)                              .* sal.^1.5 ...
        +  c10                                       * sal.^2   ;

%% Equation (30) 
% Constants (taken directly from CRP94)
a1  =    562.69486        ;
a2  = -  102.5154         ;
a3  = -    1.117033 * 1e-4;
a4  =      0.2477538      ;
a5  = -13273.75           ;
% Equation
log10_Ko = a1 + a2*log(tk) + a3*tk.^2 + a4*tk + a5./tk;

%% Calculate pkHSO4 (WM13)
pkHSO4 = -(log10_KK + log10_Ko);

end %function calk_WM13
