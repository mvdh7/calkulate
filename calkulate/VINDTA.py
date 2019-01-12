from . import calib, conc, dens, dissoc, gettit, solve

# ===== INPUTS AND THEIR UNITS ================================================


# Vsamp    = sample volume                    / ml
# Cacid    = acid molality                    / mol/kg
# S        = practical salinity
# AT_cert  = certified total alkalinity       / mol/kg-sw
# CT       = dissolved inorganic carbon       / mol/kg-sw
# PT       = phosphate                        / mol/kg-sw
# SiT      = silicate                         / mol/kg-sw
# Tk_force = titration temperature (optional) / K


# ===== PREPARATORY FUNCTION ==================================================


def prep(datfile,Vsamp,S,CT,PT,SiT, burette_cx=1, Tk_force=None):
    
    Vacid, EMF, Tk = gettit.VINDTA(datfile)

    if Tk_force is not None:
        Tk[:] = Tk_force

    Msamp =              Vsamp * dens.sw(Tk[0],S) * 1e-3 # sample mass / kg    
    Macid = burette_cx * Vacid * dens.acid(Tk)    * 1e-3 # acid mass   / kg
    
    XT = conc.XT(S,CT,PT,SiT)          # total concentrations  / mol/kg-sw
    KX = dissoc.KX_F(Tk,S,XT[3],XT[4]) # dissociation constants
    
    return Macid, EMF, Tk, Msamp, XT, KX


# ===== HALF-GRAN FUNCTIONS ===================================================


def halfGran(datfile,Vsamp,Cacid,S,CT,PT, burette_cx=1, Tk_force=None):
    
    Macid,EMF,Tk,Msamp,XT,KX = prep(datfile,Vsamp,S,CT,PT,0,
                                    burette_cx,Tk_force)
    
    return solve.halfGran(Macid,EMF,Tk,Msamp,Cacid,*XT,*KX)


def halfGran_CRM(datfile,Vsamp,AT_cert,S,CT,PT, burette_cx=1, Tk_force=None):
    
    Macid,EMF,Tk,Msamp,XT,KX = prep(datfile,Vsamp,S,CT,PT,0,
                                    burette_cx,Tk_force)
    
    Cacid = calib.halfGran(Macid,EMF,Tk,Msamp,AT_cert,XT,KX)['x'][0]
    
    AT,EMF0,_,_,_,_ = solve.halfGran(Macid,EMF,Tk,Msamp,Cacid,*XT,*KX)
    
    return Cacid, AT, EMF0


# ===== MPH FUNCTIONS =========================================================


def MPH(datfile,Vsamp,Cacid,S,CT,PT,SiT, burette_cx=1, Tk_force=None):
    
    Macid,EMF,Tk,Msamp,XT,KX = prep(datfile,Vsamp,S,CT,PT,SiT,
                                    burette_cx,Tk_force)
        
    return solve.MPH(Macid,EMF,Tk,Msamp,Cacid,*XT,KX)


def MPH_CRM(datfile,Vsamp,AT_cert,S,CT,PT,SiT, burette_cx=1, Tk_force=None):
    
    Macid,EMF,Tk,Msamp,XT,KX = prep(datfile,Vsamp,S,CT,PT,SiT,
                                    burette_cx,Tk_force)
    
    Cacid = calib.MPH(Macid,EMF,Tk,Msamp,AT_cert,XT,KX)['x'][0]
    
    AT,EMF0 = solve.MPH(Macid,EMF,Tk,Msamp,Cacid,*XT,KX)['x']
    
    return Cacid, AT, EMF0
