# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
from numpy import full_like, linspace, log10, logical_and, mean, nan, sqrt
from numpy import min as np_min
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

def emf_titration(ax, Macid, emf, rgb, AT_emf0, Msamp, Cacid, ATg, sublabel):
    """EMF change as acid is added throughout a titration."""
    ax.scatter(Macid*1e3, emf, c=rgb, edgecolors='k', clip_on=False)
    ax.axvline(1e3*AT_emf0['x'][0]*Msamp/Cacid, color='b', linestyle='--')
    ax.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ax.set_xlim([0, np_max(Macid)*1e3])
    yrange = np_max(emf) - np_min(emf)
    ax.set_ylim([np_min(emf)-yrange*0.05, np_max(emf)+yrange*0.05])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('EMF / mV')
    ax.set_title('{}Final EMF$^\circ$ = {:.2f} mV'.format(sublabel,
        AT_emf0['x'][1]), fontsize=10)
    return ax

def f1_Gran(ax, Macid, f1g, rgb, ATg, Msamp, Cacid, sublabel):
    """F1 Gran plot function for the first alkalinity estimate."""
    ax.scatter(Macid*1e3, f1g*1e-7, c=rgb, edgecolors='k', clip_on=False)
    ax.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylim([0, np_max(f1g*1.05e-7)])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$F_1 \cdot 10^{-7}$')
    ax.set_title('{}First-guess alkalinity = {:.1f} μmol/kg'.format(sublabel,
        ATg*1e6), fontsize=10)
    return ax
    
def emf0_estimate(ax, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid, sublabel):
    """First estimate of EMF0."""
    ax.scatter(Macid[Lg]*1e3, gemf0, c=rgb[Lg], edgecolors='k', clip_on=False)
    ax.axhline(emf0g, color='y')
    ax.axvline(1e3*ATg*Msamp/Cacid, color='y', linestyle='--')
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('First-guess EMF$^\circ$ / mV')
    ax.set_title('{}First-guess EMF$^\circ$ = {:.2f} mV'.format(sublabel,
        emf0g), fontsize=10)
    return ax

def AT_titration(ax, Macid, solution, sublabel):
    """Linear decrease in alkalinity as acid is added during a titration."""
    ax.scatter(Macid*1e3, solution[0]*1e6, clip_on=False)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylabel('Alkalinity / μmol kg$^{-1}$')
    ax.set_xlabel('Acid mass / g')
    title(sublabel, fontsize=10)
    return ax

def AT_estimates(ax, Macid, ATpts, rgb, AT_emf0, RMS, Npts, sublabel):
    """Original sample alkalinity estimated from each titration point."""
    ax.scatter(Macid*1e3, ATpts*1e6, c=rgb, edgecolors='k', clip_on=False)
    ax.axhline(AT_emf0['x'][0]*1e6, color='b')
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylabel('AT from pH / μmol kg$^{-1}$')
    ax.set_title('{}Final alkalinity = ({:.1f} $\pm$ {:.1f}) μmol/kg ($n$ = {})'.format(
        sublabel, AT_emf0['x'][0]*1e6, RMS*1e6, Npts), fontsize=10)
    return ax

def components(ax, Macid, ATpts, solution, sublabel):
    """Changes in every component of alkalinity throughout a titration."""
    ax.plot(Macid*1e3, -log10(ATpts), label='AT', marker='o', markersize=3,
             c='k', clip_on=False)
    names = [
        '+[HCO$_{3}^{-}$]',
        '+2[CO$_{3}^{2-}$]',
        '+[B(OH)$_{4}$]',
        '+[OH$^{-}$]',
        '-[H$^{+}$]',
        '-[HSO$_{4}^{-}$]',
        '-[HF]',
        '-[H$_{3}$PO$_{4}$]',
        '+[HPO$_{4}^{2-}$]',
        '+2[PO$_{4}^{3-}$]',
        '+[SiO(OH)$_{3}$]',
    ]
    c = iter(tab20(linspace(0, 0.9, len(names))))
    for k, name in enumerate(names):
        if name.startswith('-'):
            solution[1][k] = solution[1][k]*-1
    for j, name in enumerate(names):
        ax.plot(Macid*1e3, -log10(solution[1][j]), label=name, marker='x',
                 markersize=3, c=next(c), clip_on=False)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_xlabel('Mass of acid / g')
    ax.set_ylabel('-log$_{10}$(concentration from pH / mol kg$^{-1}$)')
    ax.set_title(sublabel, fontsize=10)
    return ax

def graph(datfile, Vsamp, psal, CT, PT, SiT, Cacid):
    """Plot everything for a single titration."""
    (Macid, EMF, Msamp, f1g, Lg, ATg, emf0g, gemf0, AT_emf0, solution, ATpts,
        RMS, Npts, rgb) = _graphprep(datfile, Vsamp, psal, CT, PT, SiT, Cacid)
    fig = figure(figsize=[17, 10])
    rcParams.update({'font.size': 10})
    gs = fig.add_gridspec(4, 2)
    subplots_adjust(wspace=0.3, hspace=0.8)
    ax1 = fig.add_subplot(gs[0, 0])
    emf_titration(ax1, Macid, EMF, rgb, AT_emf0, Msamp, Cacid, ATg, '(a) ')
    ax2 = fig.add_subplot(gs[1, 0])
    f1_Gran(ax2, Macid, f1g, rgb, ATg, Msamp, Cacid, '(b) ')
    ax3 = fig.add_subplot(gs[2, 0])
    emf0_estimate(ax3, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid, '(c) ')
    ax4 = fig.add_subplot(gs[3, 0])
    AT_titration(ax4, Macid, solution, '(d) ')
    ax5 = fig.add_subplot(gs[0, 1])
    AT_estimates(ax5, Macid, ATpts, rgb, AT_emf0, RMS, Npts, '(e) ')
    ax6 = fig.add_subplot(gs[1:4, 1])
    components(ax6, Macid, ATpts, solution, '(f)')
    return fig
