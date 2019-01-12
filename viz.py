from matplotlib import pyplot as plt
from numpy import array, log10, std
from numpy import abs as np_abs
from numpy import min as np_min
from numpy import max as np_max
from . import sim, solve

def tit(Macid,EMF,Tk,EMF0,Msamp,Cacid,XT,KX,LF):
    
    F1 = solve.F1(Macid,EMF,Tk,Msamp)
    
    
    mu = Msamp / (Msamp + Macid)
    
    lsqH  = solve.EMF2H(EMF,EMF0,Tk)
    lsqAT,lsqconcs = sim.AT(lsqH,mu,*XT,*KX)
    sumconc = lsqAT - XT[0]*mu + Macid*Cacid/(Macid+Msamp)
    
    fclrs = ['xkcd:navy', 'xkcd:navy', 'xkcd:blue', 'xkcd:teal', 'xkcd:teal',
             'xkcd:red', 'xkcd:tangerine', 'xkcd:purple', 'xkcd:purple',
             'xkcd:purple', 'xkcd:green']
    
    fstyl = ['-','--','-','-','--','-','-','-','--',':','-']
    
    fxl = array([np_min(Macid),np_max(Macid)])
    
    fig,ax = plt.subplots(2,2)
    
    # (a) lsq function
    ax[0,0].plot(Macid,sumconc*1e6, c='k', linewidth=2)
    ax[0,0].scatter(Macid[LF],sumconc[LF]*1e6, c='r')
    ax[0,0].set_xlim(fxl)
    
    ax[0,0].set_title(std(sumconc[LF]*1e6))
    
    ax[0,0].grid(alpha=0.5)
    
    ax[0,0].set_xlabel('Acid mass / kg')
    ax[0,0].set_ylabel('AT diff / umol/kg')
    
    # (b) F1 function
    ax[0,1].plot(Macid,F1, c='k')
    ax[0,1].scatter(Macid[LF],F1[LF], c='r')
    ax[0,1].grid(alpha=0.5)
    
    ax[0,1].set_xlim(fxl)
    
    ax[0,1].set_xlabel('Acid mass / kg')
    ax[0,1].set_ylabel('F1')
    
    # (c) alkalinity components
    ax[1,0].set_xlim(fxl)
    
    for i,conc in enumerate(lsqconcs):
        ax[1,0].plot(Macid,log10(np_abs(conc)), c=fclrs[i], linestyle=fstyl[i])
    
    ax[1,0].plot(Macid,log10(XT[0] + sumconc),'k', linewidth=2)
    
    ax[1,0].legend(('bicarb','2*carb', 'B4', 'OH', 'H', 'HSO4', 'HF', 
                  'P0','P2','2*P3', 'SiOOH3', 'total'))
    
    ax[1,0].set_xlabel('Acid mass / kg')
    ax[1,0].set_ylabel('log10 concentration')
    
    # (d) pH range used
    ax[1,1].plot(-log10(solve.EMF2H(EMF,EMF0,Tk)),F1, c='k')
    ax[1,1].scatter(-log10(solve.EMF2H(EMF[LF],EMF0,Tk[LF])),F1[LF], c='r')
    
    ax[1,1].grid(alpha=0.5)
    
    ax[1,1].set_xlabel('pH')
    ax[1,1].set_ylabel('F1')
