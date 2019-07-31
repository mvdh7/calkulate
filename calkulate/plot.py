# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
from numpy import full_like, linspace, log10, logical_and, mean, nan, sqrt
from numpy import max as np_max
from matplotlib.pyplot import (figure, gca, rcParams, subplots_adjust, title,
    xlabel, ylabel)
from matplotlib.cm import tab20
from . import solve, vindta

def _graphprep(datfile, Vsamp, psal, CT, PT, SiT, Cacid):
    """Preparatory calculations for plotting."""
    Macid, EMF, tempK, Msamp, XT, KXF = vindta.prep(datfile, Vsamp, psal, CT, 
        PT, SiT)
    f1g = solve.f1(Macid, EMF, tempK, Msamp)
    Lg = logical_and(f1g > 0.1*np_max(f1g), f1g < 0.9*np_max(f1g))
    ATg, emf0g, _, pHg = solve.guessGran(Macid, EMF, tempK, Msamp, Cacid)
    gemf0 = solve.Gran_EMF0(Macid[Lg], EMF[Lg], tempK[Lg], Msamp, Cacid, ATg)
    L = logical_and(pHg > 3, pHg < 4)
    AT_emf0 = solve.complete(Macid, EMF, tempK, Msamp, Cacid, XT, KXF)
    H = solve.emf2h(EMF, AT_emf0['x'][1], tempK)
    solution = vindta.simAT(Macid, tempK, H, Msamp, psal, CT, PT, SiT)
    mu = Msamp / (Msamp + Macid)
    ATpts = (solution[0] + Macid*Cacid/(Macid + Msamp))/mu
    RMS = sqrt(mean(AT_emf0['fun']**2))
    Npts = len(AT_emf0['fun'])
    rgb = full_like(EMF, nan, dtype=str)
    for i in range(len(rgb)):
        if Lg[i] and L[i]:
            rgb[i] = 'g'
        elif Lg[i] and not L[i]:
            rgb[i] = 'y'
        elif L[i] and not Lg[i]:
            rgb[i] = 'b'
        else:
            rgb[i] = 'w'
    return (Macid, EMF, Msamp, f1g, Lg, ATg, emf0g, gemf0, AT_emf0, solution,
        ATpts, RMS, Npts, rgb)

def graphax1(fig, gs, Macid, EMF, rgb, AT_emf0, Msamp, Cacid, ATg):
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(Macid*1e3, EMF, c=rgb, edgecolors='k')
    ax1.axvline(1e3*AT_emf0['x'][0]*Msamp/Cacid, color='b', linestyle='--')
    ax1.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ylabel('EMF (mV)')
    title('a) Final EMF$^{0}$ = %f mV' % AT_emf0['x'][1], fontsize=10)
    return ax1

def graphax2(fig, gs, Macid, f1g, rgb, ATg, Msamp, Cacid):
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.scatter(Macid*1e3, f1g, c=rgb, edgecolors='k')
    ax2.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ylabel('F$_{1}$')
    title('b) First guess AT = %f μmol kg$^{-1}$' % (ATg*1e6), fontsize=10)
    return ax2
    
def graphax3(fig, gs, ax1, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid):
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
    ax3.scatter(Macid[Lg]*1e3, gemf0, c=rgb[Lg], edgecolors='k')
    ax3.axhline(emf0g, color='y')
    ax3.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ylabel('First guess \n EMF$^{0}$ (mV)')
    title('c) First guess EMF$^{0}$ = %f mV' % emf0g, fontsize=10)
    return ax3

def graphax4(fig, gs, Macid, solution):
    ax4 = fig.add_subplot(gs[3,0])
    ax4.scatter(Macid*1e3, solution[0]*1e6)
    ylabel('Estimated titration \n alkalinity (μmol kg$^{-1}$)')
    xlabel('Mass of acid (g)')
    title('d) Linear decrease in AT during titration', fontsize=10)
    return ax4

def graphax5(fig, gs, Macid, ATpts, rgb, AT_emf0, RMS, Npts):
    ax5 = fig.add_subplot(gs[0,1])
    ax5.scatter(Macid*1e3, ATpts*1e6, c=rgb, edgecolors='k')
    ax5.axhline(AT_emf0['x'][0]*1e6, color='b')
    ylabel('AT from pH \n (μmol kg$^{-1}$)')
    title('e) Final AT = %f μmol kg$^{-1} \pm$ %f (n = %i)' % \
         ((AT_emf0['x'][0]*1e6), (RMS*1e6), Npts), fontsize=10)
    return ax5

def graphax6(fig, gs, Macid, ATpts, solution):
    ax6 = fig.add_subplot(gs[1:4, 1])
    ax6.plot(Macid*1e3, -log10(ATpts), label='AT', marker='o', markersize=3, \
             c='k')
    names = ['+[HCO$_{3}^{-}$]', '+2[CO$_{3}^{2-}$]', '+[B(OH)$_{4}$]', \
             '+[OH$^{-}$]', '-[H$^{+}$]', '-[HSO$_{4}^{-}$]', '-[HF]', \
             '-[H$_{3}$PO$_{4}$]', '+[HPO$_{4}^{2-}$]', '+2[PO$_{4}^{3-}$]', \
             '+[SiO(OH)$_{3}$]']
    c = iter(tab20(linspace(0, 0.9, len(names))))
    for k, name in enumerate(names):
        if name.startswith('-'):
            solution[1][k] = solution[1][k]*-1
    for j, name in enumerate(names):
        ax6.plot(Macid*1e3, -log10(solution[1][j]), label=name, marker='x', \
                 markersize=3, c=next(c))
    ax = gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ax6.legend(bbox_to_anchor=(1.05, 1))
    xlabel('Mass of acid (g)')
    ylabel('-log$_{10}$(Concentration from pH)(mol kg$^{-1}$)')
    title('f) AT components', fontsize=10)
    return ax6

def graph(datfile, Vsamp, psal, CT, PT, SiT, Cacid):
    """Plot everything for a single titration."""
    (Macid, EMF, Msamp, f1g, Lg, ATg, emf0g, gemf0, AT_emf0, solution, ATpts,
        RMS, Npts, rgb) = _graphprep(datfile, Vsamp, psal, CT, PT, SiT, Cacid)
    fig = figure(figsize=[17, 10])
    rcParams.update({'font.size': 10})
    gs = fig.add_gridspec(4, 2)
    subplots_adjust(wspace=0.3, hspace=0.8)
    ax1 = graphax1(fig, gs, Macid, EMF, rgb, AT_emf0, Msamp, Cacid, ATg)
    ax2 = graphax2(fig, gs, Macid, f1g, rgb, ATg, Msamp, Cacid)
    ax3 = graphax3(fig, gs, ax1, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid)
    ax4 = graphax4(fig, gs, Macid, solution)
    ax5 = graphax5(fig, gs, Macid, ATpts, rgb, AT_emf0, RMS, Npts)
    ax6 = graphax6(fig, gs, Macid, ATpts, solution)
    return fig
