% Initialise Python
calk_initpy('//anaconda/envs/spritzer/bin/python') % Mac
% calk_initpy( ... Windows
%     'C:\Users\yau17reu\anaconda\Anaconda3\envs\spritzer\pythonw.exe')

% Settings
datpath = 'datfiles/';
datfile = '0-0  0  (0)CRM-144-0435-4.dat';
Vsamp = 100;
% Cacid = 0.1;
[CT,AT_cert,S,PT,SiT] = dicksonCRM(144);
CT      = CT      * 1e-6;
AT_cert = AT_cert * 1e-6;
PT      = PT      * 1e-6;
SiT     = SiT     * 1e-6;
burette_cx = 1;
Tk_force = [];
printpath = 'figures/';

% Get proper Cacid
Cacid = calk_VINDTA_CRM([datpath datfile],Vsamp,AT_cert,S,CT,PT,SiT, ...
    burette_cx,Tk_force);

% Plot the figure!
t = calk_ptl(datpath,datfile,Vsamp,Cacid,CT,S,PT,SiT, ...
    burette_cx,Tk_force,printpath);
