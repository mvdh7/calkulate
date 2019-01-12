%% This section returns virtually identical results to Python Gran version

% Differences are that MATLAB Gran convergence checks against last 2
%  results, in case result is bouncing between 2 values, and then also
%  reports the mean of the last two results at the end.
% CalkCRM also is a less accurate fitting algorithm than the Python
%  version, taking that into account the results are the same.

% datfile           = '0-0  0  (0)CRM-151-0494-2.dat';
datfile           = '0-0  0  (0)CRM-151-0855-3.dat';
acidmolar         = 0.1 * rho_acid(298.15); % mol/l
acidrho           = rho_acid(298.15); % kg/l
pipvol            = 100; % ml
sal               = 33.345;  % practical
dic               = 2033.83 * 1e-6; % mol/kg
phos              =    0.56 * 1e-6; % mol/kg
si                =    3.5  * 1e-6; % mol/kg
burette_cxn       = 1;
sc                = [];
tforce            = [];
suppress_warnings = false;

TA_cert = 2225.56;
tic
[TA_final,E0_final,t,f2xint,meanT,pipmass,sc,ta,e0] ...
    = Calkulate(datfile,acidmolar,acidrho,pipvol, ...
        sal,dic,phos,burette_cxn,sc,tforce,suppress_warnings);
    
acidmolar_cal = CalkCRM(datfile, ...
    acidmolar,acidrho,pipvol,151,burette_cxn,sc,tforce);

[TA_final_cal,E0_final_cal,t_cal,~,~,~,~,ta_cal,e0_cal] ...
    = Calkulate(datfile,acidmolar_cal,acidrho,pipvol, ...
        sal,dic,phos,burette_cxn,sc,tforce,suppress_warnings);
    
TA_diff = TA_final_cal - TA_cert;
toc

%% Test MPH methods in Python
initpy('C:\Users\yau17reu\anaconda\Anaconda3\envs\calkenv\pythonw.exe')
tic
pygran = py.calkulate.VINDTA.halfGran(datfile,pipvol, ...
    acidmolar_cal*acidrho,sal,dic,phos);

pyacid = numpy2double(py.calkulate.VINDTA.MPH_CRM(datfile,pipvol, ...
    TA_cert*1e-6,sal,dic,phos,si));

pymph = py.calkulate.VINDTA.MPH(datfile,pipvol, ...
    pyacid(1),sal,dic,phos,si);

pymph_TA_EMF0 = numpy2double(pymph{'x'});

pymph_TA = pymph_TA_EMF0(1) * 1e6;
toc
