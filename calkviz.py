import calkulate as calk
datfile = 'datfiles/CRM-144-0435-4.dat'
Vsamp = 100. # ml
psal = 33.571
CT = 0.00203153
PT = 3.1e-7
SiT = 2.5e-6
Cacid = 0.1
calk.plot.everything(datfile, Vsamp, psal, CT, PT, SiT, Cacid)
test = calk.io.vindta(datfile)
