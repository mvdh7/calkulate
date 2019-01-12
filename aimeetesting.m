initpy('C:\Users\yau17reu\anaconda\Anaconda3\envs\calkenv\pythonw.exe')

fpath = 'E:/Dropbox/MATLAB/My toolkits/Calkulate/calkulate-uea/';

calk_VINDTA_CRM([fpath '142-0  0  (0)crm_0788_01_01.dat'],110.71, ...
    2.2276e-3,33.389,2.0381e-3,0.29e-6,3.3e-6,1,[])

testgft = py.numpy.genfromtxt([fpath 'aimeedat.dat']);
