pyversion('/anaconda/envs/calkenv/bin/python');

%%
datfile = '0-0  0  (0)CRM-151-0494-2.dat';
vdata = py.calkulate.gettit.VINDTA(datfile);

XT = py.calkulate.conc.XT(35);
KX = py.calkulate.dissoc.KX(298.15,35);

test = cell(py.calkulate.solve.Gran_VINDTA(datfile,0.1,0.1, ...
    35,2000e-6,0));
