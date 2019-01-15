calk_initpy( ... Windows
    'C:\Users\yau17reu\anaconda\Anaconda3\envs\spritzer\pythonw.exe')

% ---------------------------------------------------- Inputs to vary -----

clear t
t.Vacid = (0:0.15:4.2)'; % acid volume in ml
t = struct2table(t);

t.Tk = 298.15 * ones(size(t.Vacid)); % temperature in K

Vsamp = 100; % sample volume in ml
Cacid = 0.1; % acid concentration in mol/kg

S = 35; % practical salinity

AT  = 2250 * 1e-6; % total alkalinity           in mol/kg-sw
CT  = 2100 * 1e-6; % dissolved inorganic carbon in mol/kg-sw
PT  =    1 * 1e-6; % total phosphate            in mol/kg-sw
SiT =   10 * 1e-6; % total silicate             in mol/kg-sw

EMF0 = 660; % EMF0 in mV

% ------------------------------------------------------ Calculations -----

% Acid and sample volume --> mass
t.Macid = t.Vacid .* calk_dens_acid(t.Tk) * 1e-3; % acid mass in kg
Msamp = Vsamp * calk_dens_sw(t.Tk,S) * 1e-3;
Msamp = Msamp(1); % sample mass in kg

% Solve for [H+] and pH
[t.H,t.pH] = calk_simH(t.Macid,t.Tk,Msamp,Cacid,S,AT,CT,PT,SiT);

% Convert H to EMF
t.EMF = calk_H2EMF(t.H,EMF0,t.Tk);

% ------------------------------------------------- Save as .dat file -----

fid = fopen('datfiles/calk_simtit.dat','w');

fprintf(fid,'\r\n\r\n');
fprintf(fid,'%.3f\t%.3f\t%.3f\r\n',[t.Vacid t.EMF t.Tk-273.15]');

fclose(fid);
