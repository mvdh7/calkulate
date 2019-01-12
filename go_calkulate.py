import calkulate as calk
import numpy as np

# Import Dickson (1981) test titration dataset (Table 1, no phosphate)
Macid,pH,Tk,Msamp,Cacid,sal,XT,KX = calk.gettit.Dickson1981()
H = 10**-pH

#[KC1,KC2,KB,Kw,KHSO4,KHF] = KX
#[AT,CT,BT,ST,FT,PT] = XT

# Test with phosphate too (vs D81 Table 4)
XT[-1] = 0.00001 # PT / mol/kg

# Re-simulate Dickson's test titration dataset (Table 1, no phosphate)
Hsim = calk.sim.H(Macid,Msamp,Cacid,*XT,*KX)

# Add random measurement errors
Hsim = Hsim + np.random.normal(loc=0,scale=Hsim*0.01,size=np.size(Hsim))

# Convert to pH
pHsim = -np.log10(Hsim)

# Select range for fitting
L = np.logical_and(pHsim > 3,pHsim < 4)

# Solve following DAA03
AT_DAA03,f_DAA03 = calk.solve.DAA03(Macid[L],Hsim[L],Msamp,Cacid,*XT,*KX[:-3])

# Solve following DAA03 but with f fixed at 1
AT_DAA03_pH = calk.solve.DAA03_pH(Macid[L],Hsim[L],Msamp,Cacid,*XT,*KX[:-3])

# Solve my own way
AT_MPH,f_MPH = calk.solve.MPH(Macid[L],Hsim[L],Msamp,Cacid,*XT,*KX)

dE0_MPH = calk.solve.f2dEMF0(Tk,f_MPH)

#%% Gran plot approach
Gran_EMF0 = 630.
EMF = calk.solve.H2EMF(Hsim,Gran_EMF0,Tk)

Gran_AT_final,Gran_E0_final,i,Gran_AT = calk.solve.Gran(Macid,EMF,Tk,
    Msamp,Cacid,*XT,*KX)
