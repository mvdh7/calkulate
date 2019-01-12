%% Load Dickson (1981) simulated titration data
clear d
d.mass_acid = (0:0.05:2.5)';
d.ph = [8.065650; 7.925895; 7.752062; 7.539922; 7.312923; 7.111723;
    6.948715; 6.816116; 6.704825; 6.608421; 6.522665; 6.444704;
    6.372551; 6.304758; 6.240231; 6.178101; 6.117655; 6.058275;
    5.999403; 5.940504; 5.881044; 5.820455; 5.758107; 5.693259;
    5.625006; 5.552183; 5.473224; 5.385920; 5.286975; 5.171175;
    5.029724; 4.847252; 4.601818; 4.304988; 4.046597; 3.860034;
    3.722943; 3.616542; 3.530163; 3.457664; 3.395285; 3.340587;
    3.291906; 3.248062; 3.208187; 3.171626; 3.137874; 3.106531;
    3.077278; 3.049854; 3.024045]; % -log[H+]
d = struct2table(d);

d_pure = d;

m_sample = 200; % g

TA = 0.00245; % mol/kg-sw
TB = 0.00042; % mol/kg-sw
TC = 0.00220; % mol/kg-sw
TS = 0.02824; % mol/kg-sw
TF = 0.00007; % mol/kg-sw

% Constants on Free H scale
Kw    = 4.32e-14;
K1    = 1.00e-06;
K1K2  = 8.20e-16;
KB    = 1.78e-09;
bHSO4 = 1.23e+01;
bHF   = 4.08e+02;

% Run Calkulate
sc.TSO4  = 'D81';
sc.THF   = 'D81';
sc.TB    = 'D81';
sc.kHSO4 = 'D81';
sc.kHF   = 'D81';
sc.k1k2  = 'D81';

d.vol_raw  = d.mass_acid; % Volume of acid added / mL

EMF0 = 630;

d.h = 10.^-d.ph;
d.emf = EMF0 + log(d.h) * 8.3144621 * 298.15 / 96.4853365;

d.t = 25*ones(height(d),1); % Temperature          / deg C

sal = 35;

vol_sample = m_sample * 1000 / cMP81(25,35);
acidmolar = 0.3;

disp('MATLAB Gran solver:')
tic
[TA_final,E0_final,t,f2xint,meanT,pipmass,sc,ta,e0] ...
    = Calkulate(d,acidmolar,1,vol_sample,sal,TC,0,1,sc,25);
toc

%% Test Python
python_exe = ['C:\\Users\\yau17reu\\anaconda\\Anaconda3\\envs\\' ...
    'calkulate\\python.exe'];
[~,~,pyloaded] = pyversion;
if ~pyloaded
    pyversion(python_exe)
end %if

%%
disp('Python Gran solver:')

tic
py_acid_mass   = vec2numpy(d.mass_acid * 1e-3);
py_emf         = vec2numpy(d.emf             );
py_Tk          = vec2numpy(d.t + 273.15      );

Gran_out = cell(py.calkulate.solve.Gran( ...
    m_sample/1000,py_acid_mass,py_emf,py_Tk,acidmolar, ...
    TC,K1,TS,1/bHSO4,TF,1/bHF,TB,KB,Kw));

TA_final_py =              Gran_out{1} ;
E0_final_py =              Gran_out{2} ;
i_py        =              Gran_out{3} ;
TA_py       = numpy2double(Gran_out{4});
toc

%% Plot it
figure(1); clf;

plot(d.m_acid,d.ph)

%% Save .csv
writetable(d_pure,'dickson81tit.csv');
