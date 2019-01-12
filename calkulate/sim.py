from scipy.optimize import least_squares as olsq
from numpy import full_like, inf, nan

# Simulate total alkalinity from known pH and total concentrations

def AT(H,mu,xATx,CT,BT,ST,FT,PT,SiT,KC1,KC2,KB,Kw,KHSO4,KHF,KP1,KP2,KP3,KSi):
       
    # XT are values at start of titration
    
    # Carbonates
    CO2aq  = mu * CT / (1 + KC1/H + KC1*KC2/H**2)
    bicarb = KC1 * CO2aq  / H
    carb   = KC2 * bicarb / H
    
    # Borate
    B4 = mu * BT * KB / (H + KB)
    
    # Hydroxide
    OH = Kw / H
    
    # Bisulfate
    HSO4 = mu * ST * H / (KHSO4 + H)
    
    # Hydrogen fluoride
    HF   = mu * FT * H / (KHF   + H)
    
    # Phosphates
    P0 = mu * PT / (1 + KP1/H + KP1*KP2/H**2 + KP1*KP2*KP3/H**3)     # [H3PO4]
    P2 = mu * PT / (H**2/(KP1*KP2) + H/KP2 + 1 + KP3/H)              # [HPO42-]
    P3 = mu * PT / (H**3/(KP1*KP2*KP3) + H**2/(KP2*KP3) + H/KP3 + 1) # [PO43-]
    
    # Silicate
    SiOOH3 = mu * SiT * KSi / (H + KSi)
    
    # Total alkalinity
    AT = bicarb + 2*carb + B4 + OH - H - HSO4 - HF - P0 + P2 + 2*P3 + SiOOH3
    
    return AT, [bicarb,2*carb, B4, OH, -H, -HSO4, -HF, -P0,P2,2*P3, SiOOH3]


# Simulate pH from known total alkalinity and total concentrations

def H(Macid,Msamp,Cacid,XT,KX):
    
    H = full_like(Macid,nan)
    
    mu = Msamp / (Msamp + Macid)
    
    for i,Macid_i in enumerate(Macid):
    
        H[i] = olsq(lambda H: \
            AT(H,mu[i],*XT,*KX)[0] - mu[i] * XT[0] + Macid_i*Cacid                       \
                / (Macid_i + Msamp),
            1e-8, method='trf', bounds=(1e-15,inf))['x']
        
    return H
