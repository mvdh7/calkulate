% Initialise Python
% calk_initpy('//anaconda/envs/spritzer/bin/python') % Mac
calk_initpy( ... Windows
    'C:\Users\yau17reu\anaconda\Anaconda3\envs\spritzer\pythonw.exe')

% Input settings
datpath = 'datfiles/';
% datfile = '0-0  0  (0)CRM-144-0435-4.dat';
datfile = 'Ruth/9999_20190205_7_0_0_2.dat';
Vsamp = 100.0123; % ml
[CT,AT_cert,S,PT,SiT] = dicksonCRM(144);
CT      = CT      * 1e-6;
AT_cert = AT_cert * 1e-6;
PT      = PT      * 1e-6;
SiT     = SiT     * 1e-6;
burette_cx = 1;
Tk_force = [];
printpath = 'figures/';

% % Either guess Cacid...
% Cacid = 0.1;

% ... or get proper Cacid from a CRM:
Cacid = calk_VINDTA_CRM([datpath datfile],Vsamp,AT_cert,S,CT,PT,SiT, ...
    burette_cx,Tk_force);

% Plot the figure!
t = calk_ptl(datpath,datfile,Vsamp,Cacid,CT,S,PT,SiT, ...
    burette_cx,Tk_force,[]);
