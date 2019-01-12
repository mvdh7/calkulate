import calkulate as calk
import numpy as np

# Import Dickson (1981) Table 1
Macid, pH, Tk, Msamp, Cacid, sal, XT, KX = calk.gettit.Dickson1981()
AT_cert = XT[0]
CT_cert = XT[1]

# Re-simulate pH
H  = calk.sim.H(Macid,Msamp,Cacid,XT,KX)
pH = -np.log10(H)

# Re-simulate pH with no DIC for DAA03 method
XT_DAA03 = np.copy(XT)
XT_DAA03[1] = 0

H_DAA03  = calk.sim.H(Macid,Msamp,Cacid,XT_DAA03,KX)
pH_DAA03 = -np.log10(H_DAA03)

# Convert pH to EMF using arbitrary EMF0
#H = 10**-pH # if not re-simulating

EMF0 = 660.
EMF_pure       = calk.solve.H2EMF(H      ,EMF0,Tk)
EMF_pure_DAA03 = calk.solve.H2EMF(H_DAA03,EMF0,Tk)

Ureps = int(500)

METHODS = ['MPH', 'hG', 'D81', 'DAA03a', 'DAA03b']

AT = {}
E0 = {}
for METHOD in METHODS:
    AT[METHOD] = np.full(Ureps,np.nan)
    E0[METHOD] = np.full(Ureps,np.nan)

for U in range(Ureps):

    # Add some NOISE
    Unoise = np.random.normal(size=np.size(EMF_pure), loc=0, scale=0.1)
    EMF       = EMF_pure       + Unoise
    EMF_DAA03 = EMF_pure_DAA03 + Unoise
    
    # Solve for AT etc. following Dickson (1981)
    AT_CT_f_D81 = calk.solve.Dickson1981(Macid,EMF,Tk,Msamp,Cacid,XT,KX)
    
    AT['D81'][U], CT_D81, f_D81 = AT_CT_f_D81[0]['x']
    E0['D81'][U] = AT_CT_f_D81[1]
    
    # Solve for AT etc. following Humphreys
    AT['MPH'][U], E0['MPH'][U] = calk.solve.MPH(Macid,EMF,Tk,Msamp,Cacid,
                                            *XT,KX)['x']
    
    # Solve for AT etc. with half-Gran
    AT['hG'][U], E0['hG'][U] = calk.solve.halfGran(Macid,EMF,Tk,Msamp,Cacid,
                                               *XT,*KX)
    
    # Solve for AT etc. with DAA03 - DOESN'T REALLY WORK
    AT_f_DAA03a = calk.solve.DAA03(Macid,EMF,Tk,Msamp,Cacid,*XT,KX)
    
    AT['DAA03a'][U], f_DAA03a = AT_f_DAA03a[0]['x']
    E0['DAA03a'][U] = AT_f_DAA03a[1]
    
    # Solve for AT with DAA03 with no DIC in the titration - WORKS
    AT_f_DAA03b = calk.solve.DAA03(Macid,EMF_DAA03,Tk,Msamp,Cacid,*XT,KX)
    
    AT['DAA03b'][U], f_DAA03b = AT_f_DAA03b[0]['x']
    E0['DAA03b'][U] = AT_f_DAA03b[1]
    
for METHOD in METHODS:
    print('%6s: %.1f +- %.1f ; %.2f +- %.2f' % (METHOD,
                                                np.nanmean(AT[METHOD] * 1e6),
                                                np.nanstd (AT[METHOD] * 1e6),
                                                np.nanmean(E0[METHOD]),
                                                np.nanstd (E0[METHOD])))
