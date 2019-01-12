import calkulate as calk
import numpy as np

#%% This section returns virtually identical results to Python Gran version

# Differences are that MATLAB Gran convergence checks against last 2
#  results, in case result is bouncing between 2 values, and then also
#  reports the mean of the last two results at the end.
# CalkCRM also is a less accurate fitting algorithm than the Python
#  version, taking that into account the results are the same.

# ===== A real .dat file =====

datfile = 'datfiles/0-0  0  (0)CRM-151-0494-2.dat'
#datfile = '0-0  0  (0)CRM-151-0855-3.dat'

Vacid,EMF,Tk = calk.gettit.VINDTA(datfile)

Macid = Vacid * calk.dens.acid(Tk) * 1e-3

Cacid = 0.1 # mol/kg

S   =   33.345
CT  = 2033.83  * 1e-6
PT  =    0.56  * 1e-6
SiT =    3.5   * 1e-6

AT_cert = 2225.56 * 1e-6

Vsamp = 100 # ml
Msamp = Vsamp * calk.dens.sw(Tk[0],S) * 1e-3 # kg

XT = calk.conc.XT(S,CT,PT,SiT)
KX = calk.dissoc.KX_F(Tk,S,XT[3],XT[4])

# ===== Dickson's table =====

#Macid, pH, Tk, Msamp, Cacid, sal, XT, KX = calk.gettit.Dickson1981()
#AT_cert = XT[0]
#
#H = 10**-pH                             #         Use Dickson's table
##H = calk.sim.H(Macid,Msamp,Cacid,XT,KX) # Re-simulate Dickson's table
#
#EMF = calk.solve.H2EMF(H,660,Tk)
#
#KXL = KX

# ===== same for both below here =====

F1 = calk.solve.F1(Macid,EMF,Tk,Msamp)

LF1 = np.logical_and(F1 > 0.1*np.max(F1), F1 < 0.9 * np.max(F1))

ATg  = calk.solve.Gran_guess_AT(Macid[LF1],F1[LF1],Msamp,Cacid)
EMF0g_vec = calk.solve.Gran_EMF0(Macid[LF1],EMF[LF1],Tk[LF1],Msamp,Cacid,ATg)
EMF0g = np.mean(EMF0g_vec)

Hg = calk.solve.EMF2H(EMF,EMF0g,Tk)
pHg = -np.log10(Hg)

LF = np.logical_and(pHg > 3, pHg < 4)

mu = Msamp / (Msamp + Macid)

[_, CT, BT, ST, FT, PT, SiT] = XT
[KC1,KC2, KB, KH2O, KHSO4, KHF, KP1,KP2,KP3, KSi] = KX

Gran_bicarb = mu * CT / (Hg / KC1 + 1)
Gran_HSO4   = mu * ST / (1 + KHSO4 / Hg)
Gran_HF     = mu * FT / (1 + KHF   / Hg)
Gran_borate = mu * BT / (1 + Hg / KB)
Gran_OH     = KH2O / Hg
Gran_P_P2   = mu * PT * (1 - KP1*KP2/(Hg**2)) \
    / (1 + KP1/Hg + KP2*KP3/Hg**2  \
         + KP1*KP2*KP3/Hg**3)

F1n = (  Hg  \
       + Gran_HSO4       \
       + Gran_HF         \
       - Gran_bicarb     \
       - Gran_OH         \
       - Gran_borate     \
       + Gran_P_P2     ) * (Msamp + Macid)

Gran_AT_final,Gran_E0_final = calk.solve.halfGran(Macid,
    EMF,Tk,Msamp,Cacid,*XT,*KX, pHrange=[3.,4.], suppress_warnings=False)

Cacid_cal = calk.calib.halfGran(Macid,EMF,Tk,Msamp,AT_cert,XT,KX)['x'][0]

Gran_AT_final_cal,Gran_E0_final_cal = calk.solve.halfGran(Macid,
    EMF,Tk,Msamp,Cacid_cal,*XT,*KX, pHrange=[3.,4.], suppress_warnings=False)

AT_diff = (Gran_AT_final_cal - AT_cert) * 1e6

#%% Test MPH method

Cacid_MPH_E0 = calk.calib.MPH(Macid,EMF,Tk,Msamp,AT_cert,XT,KX)['x'][0]

MPHFULL = calk.solve.MPH(Macid,EMF,Tk,Msamp,Cacid_MPH_E0,*XT,KX)

AT_MPH_E0_cal,EMF0_MPH_cal = MPHFULL['x']

AT_MPH_E0_diff = AT_MPH_E0_cal - AT_cert

# Test VINDTA function
MPHvin_Cacid = calk.VINDTA.MPH_CRM(datfile,Vsamp,AT_cert,S,CT,PT,SiT)
MPHvin = calk.VINDTA.MPH(datfile,Vsamp,Cacid_MPH_E0,S,CT,PT,SiT)

##%%
##
#XT[0] = AT_MPH_E0_cal
#calk.viz.tit(Macid,EMF,Tk,EMF0_MPH_cal,Msamp,Cacid_MPH_E0,XT,KX,LF)