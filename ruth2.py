import calkulate as calk
import numpy as np

Macid = np.linspace(0,4.2,29)
Tk    = np.full_like(Macid,298.15)

XT = calk.conc.XT(0.,CT=00e-6)
XT[0] = 00e-6
KX = calk.dissoc.KX_F(Tk,35.,XT[3],XT[4])

H = calk.sim.H(Macid,100.,0.1,XT,KX)
pH = -np.log10(H)

from matplotlib import pyplot as plt

plt.scatter(Macid,pH)
