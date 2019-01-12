import calkulate as calk

# Do a real VINDTA .dat file
datfile = 'datfiles/0-0  0  (0)CRM-144-0435-4.dat'
Vacid,EMF,Tk = calk.gettit.VINDTA(datfile)

Dacid = calk.dens.acid(Tk[0]) # kg/l

Macid = Vacid * Dacid * 1e-3 # kg

Vsamp = 100. # ml

S = 33.571

Msamp = Vsamp * calk.dens.sw(Tk[0],S) / 1e6 # kg

# Get *XT and *KX
AT_cert = 0.00223860
# ^^^ solver works with other values but not with the real one???
CT  = 0.00203153
PT  = 3.1e-7
SiT = 2.5e-6

XT = calk.conc.XT(S,CT,PT,SiT)

KX = calk.dissoc.KX_F(Tk,S,XT[3],XT[4])

# Solve

Cacid = 0.1

Macid,EMF,Tk,Msamp, F1g,Lg, ATg,EMF0g,pHg,L = calk.VINDTA.guessGran(datfile,
    Vsamp,Cacid,S)

#test = calk.VINDTA.MPH(datfile,Vsamp,Cacid,S,CT,PT,SiT,1,None)

##%% Solve for acid molarity
#Cacid = calk.solve.Gran_CRM(Macid,EMF,Tk,Msamp,AT_cert,XT,KX)
#
##Vacid, EMF, Tk = gettit.VINDTA(datfile)
##    
##Macid = Vacid * dens.acid(Tk[0]) * 1e-3 # kg
##
##XT = conc.XT(S,CT,PT)
##KX = dissoc.KX(Tk,S)
##
##return Gran(Macid,EMF,Tk,Msamp,Cacid,*XT,*KX)
#
#Gran_AT_final,Gran_E0_final,i,Gran_E0,Gran_pH = calk.solve.Gran(Macid,EMF,Tk,
#    Msamp,Cacid,*XT,*KX)
#Gran_AT_final = Gran_AT_final * 1e6
#
#Gran_AT_final_VINDTA,Gran_E0_final_VINDTA,_,_,_ = calk.solve.Gran_VINDTA( \
#    datfile,Msamp,Cacid,S,CT_cert,PT)
#Gran_AT_final_VINDTA = Gran_AT_final_VINDTA * 1e6

##%% Solve for alkalinity
#import numpy as np
##
##Cacids = np.linspace(0.09,0.11,20)
##
##Gran_ATs = np.full_like(Cacids,np.nan)
##Gran_E0s = np.full_like(Cacids,np.nan)
##
##for i,Ca in enumerate(Cacids):
##    
##    Gran_ATs[i],Gran_E0s[i],_,_ = calk.solve.Gran(Macid,EMF,Tk,
##        Msamp,Ca,*XT,*KX)
##
##Gran_ATs = Gran_ATs * 1e6
#
#H = calk.solve.EMF2H(EMF,Gran_E0_final,Tk[0])
#pH = -np.log10(H)
#
#L = np.logical_and(pH > 3,pH < 4.5)
#KXL = [KXi[L] for KXi in KX]
#
#from scipy.optimize import least_squares as olsq
#
#MPH_Cacid = olsq(lambda Cacid: \
#    calk.solve.MPH(Macid[L],H[L],Msamp,Cacid,*XT,KXL)[0] - AT_cert,0.1)['x']
#
#MPH_AT,f_MPH = calk.solve.MPH(Macid[L],H[L],Msamp,MPH_Cacid,*XT,KXL)
#
#MPH_AT_pH = calk.solve.MPH_pH(Macid[L],H[L],Msamp,MPH_Cacid,f_MPH,*XT,KXL)
#
#dE0_MPH = calk.solve.f2dEMF0(Tk,f_MPH)
#
## Sticky issue - we're solving for EMF0 _and_ TA, but when fitting to find
##   the acid concentration, EMF0 isn't held steady.
## Should EMF0 be constant? Should we pick and set a constant value for it,
##   rather than using a different value per sample?
## Investigate this with the Darwin-Guam dataset ASAP.
#
#f_test = 1.#1.033027453
#
#testCacid = calk.solve.MPH_CRM_pH(Macid[L],H[L],Msamp,AT_cert,f_test,*XT,KXL)
#testCRMTA = calk.solve.MPH_pH(Macid[L],H[L],Msamp,testCacid,f_test,*XT,KXL)
#testCRMTAf,fCRMTA = calk.solve.MPH(Macid[L],H[L],Msamp,testCacid,*XT,KXL)
#
#testCRMx = calk.solve.MPH_CRM(Macid[L],H[L],Msamp,AT_cert,*XT,KXL)
#
#
#
#
##%%
#Vacid, EMF, Tk = calk.gettit.VINDTA(datfile)
#    
#Macid = Vacid * calk.dens.acid(Tk)      * 1e-3 # kg
#Msamp = Vsamp * calk.dens.sw  (Tk[0],S) * 1e-6 # kg
#
#XT = calk.conc.XT(S,CT_cert,PT)
#KX = calk.dissoc.KX_F(Tk,S,XT[3],XT[4])
#
#F1 = calk.solve.F1(Macid,EMF,Tk,Msamp)
#LF = F1 > 0.1*np.max(F1)
#AT0  = calk.solve.Gran_guess_AT(Macid[LF],F1[LF],0.1,Msamp)
#EMF0 = calk.solve.Gran_EMF0(Macid[LF],EMF[LF],Tk[LF],Msamp,0.1,AT0)
#EMF0 = np.mean(EMF0)
#
#H0 = calk.solve.EMF2H(EMF,EMF0,Tk)
#pH0 = -np.log10(H0)
#L = np.logical_and(pH0 > 3.,pH0 < 4.)
#
#KXL = [KXi[L] for KXi in KX]
#
#testVINDTA0 = calk.solve.MPH_CRM(Macid[L],H0[L],Msamp,AT_cert,*XT,KXL)
#testVINDTA1 = calk.solve.MPH_CRM_VINDTA(datfile,Vsamp,AT_cert,S,CT_cert,PT)
#
#testgetV = calk.solve.MPH_VINDTA(datfile,Vsamp,testVINDTA1[0],S,CT_cert,PT)
