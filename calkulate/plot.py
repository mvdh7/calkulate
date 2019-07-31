# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Visualise Calkulate calculations."""
from numpy import array, log10, logical_and, mean, sqrt, zeros
from numpy import min as np_min
from numpy import max as np_max
from matplotlib.pyplot import figure, rcParams, subplots_adjust
from . import solve, vindta

_rgb_guess = array([0.96, 0.86, 0.04])
_rgb_final = array([0.21, 0.46, 1])
_rgb_both = array([0.27, 0.8, 0.54])

def prep(datfile, Vsamp, psal, CT, PT, SiT, Cacid):
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
    rgb = zeros((len(EMF), 3))
    for i in range(len(rgb)):
        if Lg[i] and L[i]:
            rgb[i] = _rgb_both
        elif Lg[i] and not L[i]:
            rgb[i] = _rgb_guess
        elif L[i] and not Lg[i]:
            rgb[i] = _rgb_final
        else:
            rgb[i] = array([1, 1, 1])
    return (Macid, EMF, Msamp, f1g, Lg, ATg, emf0g, gemf0, AT_emf0, solution,
        ATpts, RMS, Npts, rgb)

def emf_titration(ax, Macid, emf, rgb, AT_emf0, Msamp, Cacid, ATg, sublabel):
    """EMF change as acid is added throughout a titration."""
    ax.axvline(1e3*AT_emf0['x'][0]*Msamp/Cacid, color=_rgb_final,
        linestyle='--', zorder=1)
    ax.axvline(1e3*ATg*Msamp/Cacid, color=_rgb_guess, linestyle='--', zorder=1)
    ax.scatter(Macid*1e3, emf, c=rgb, edgecolors='k', clip_on=False, zorder=2)
    ax.set_xlim([0, np_max(Macid)*1e3])
    yrange = np_max(emf) - np_min(emf)
    ax.set_ylim([np_min(emf) - yrange*0.05, np_max(emf) + yrange*0.05])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('EMF / mV')
    ax.set_title('{} Final EMF$^\circ$ = {:.2f} mV'.format(sublabel,
        AT_emf0['x'][1]), fontsize=10)
    return ax

def f1_Gran(ax, Macid, f1g, rgb, ATg, Msamp, Cacid, sublabel):
    """F1 Gran plot function for the first alkalinity estimate."""
    ax.axvline(1e3*ATg*Msamp/Cacid, color=_rgb_guess, linestyle='--', zorder=1)
    ax.scatter(Macid*1e3, f1g*1e-7, c=rgb, edgecolors='k', clip_on=False,
        zorder=2)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylim([0, np_max(f1g*1.05e-7)])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$F_1 \cdot 10^{-7}$')
    ax.set_title('{} First-guess alkalinity = {:.1f} μmol/kg'.format(sublabel,
        ATg*1e6), fontsize=10)
    return ax
    
def emf0_estimate(ax, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid, sublabel):
    """First estimate of EMF0."""
    ax.axhline(emf0g, color=_rgb_guess, zorder=1)
    ax.axvline(1e3*ATg*Msamp/Cacid, color=_rgb_guess, linestyle='--', zorder=1)
    ax.scatter(Macid[Lg]*1e3, gemf0, c=rgb[Lg], edgecolors='k', clip_on=False,
        zorder=2)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('First-guess EMF$^\circ$ / mV')
    ax.set_title('{} First-guess EMF$^\circ$ = {:.1f} mV'.format(sublabel,
        emf0g), fontsize=10)
    return ax

def AT_titration(ax, rgb, Macid, solution, sublabel):
    """Linear decrease in alkalinity as acid is added during a titration."""
    ax.scatter(Macid*1e3, solution[0]*1e6, clip_on=False, edgecolors='k',
        c=rgb)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylabel('Alkalinity / μmol kg$^{-1}$')
    ax.set_xlabel('Acid mass / g')
    ax.set_title(sublabel, fontsize=10)
    return ax

def AT_estimates(ax, Macid, ATpts, rgb, AT_emf0, RMS, Npts, sublabel):
    """Original sample alkalinity estimated from each titration point."""
    ax.axhline(AT_emf0['x'][0]*1e6, color=_rgb_final, zorder=1)
    ax.scatter(Macid*1e3, ATpts*1e6, c=rgb, edgecolors='k', clip_on=False,
        zorder=2)
    ax.set_xlim([0, np_max(Macid)*1e3])
    yrange = (np_max(ATpts) - np_min(ATpts))*1e6
    ax.set_ylim([np_min(ATpts*1e6) - yrange*0.05,
        np_max(ATpts*1e6 + yrange*0.05)])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('AT from pH / μmol kg$^{-1}$')
    ax.set_title(('{} Final alkalinity = ({:.1f} $\pm$ {:.1f}) μmol/kg' + 
        ' ($n$ = {})').format(sublabel, AT_emf0['x'][0]*1e6, RMS*1e6, Npts),
        fontsize=10)
    return ax

def components(ax, Macid, ATpts, solution, sublabel):
    """Every component of alkalinity throughout a titration."""
    ax.plot(Macid*1e3, -log10(ATpts), label='AT', marker='o', markersize=3,
             c='k', clip_on=False)
    rgbs = array([
        [0.4, 0.4, 0.4],
        [0.4, 0.4, 0.4],
        [0.19, 0.31, 0.97],
        [1, 0.05, 0.05],
        [1, 0.05, 0.05],
        [0.95, 0.8, 0.14],
        [0.12, 0.74, 0.12],
        [1, 0.5, 0],
        [1, 0.5, 0],
        [1, 0.5, 0],
        [0.94, 0.56, 0.63],
    ])
    names = [
        '+[HCO$_3^-$]',
        '+2[CO$_3^{2-}$]',
        '+[B(OH)$_4^-$]',
        '+[OH$^-$]',
        '$-$[H$^+$]',
        '$-$[HSO$_4^-$]',
        '$-$[HF]',
        '$-$[H$_3$PO$_4$]',
        '+[HPO$_4^{2-}$]',
        '+2[PO$_4^{3-}$]',
        '+[SiO(OH)$_3$]',
    ]
    for k, name in enumerate(names):
        if name.startswith('$-$'):
            solution[1][k] = -solution[1][k]
    for j, name in enumerate(names):
        ax.plot(Macid*1e3, -log10(solution[1][j]), label=name, marker='x',
            markersize=3, c=rgbs[j], clip_on=False)
    ax.set_xlim([0, np_max(Macid)*1e3])
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$-$log$_{10}$(concentration from pH / mol kg$^{-1}$)')
    ax.set_title(sublabel, fontsize=10)
    return ax

def titration(datfile, Vsamp, psal, CT, PT, SiT, Cacid):
    """Plot everything for a single titration with a VINDTA-style .dat file."""
    (Macid, EMF, Msamp, f1g, Lg, ATg, emf0g, gemf0, AT_emf0, solution, ATpts,
        RMS, Npts, rgb) = prep(datfile, Vsamp, psal, CT, PT, SiT, Cacid)
    fig = figure(figsize=[17, 10])
    rcParams.update({'font.size': 10})
    gs = fig.add_gridspec(4, 2)
    subplots_adjust(wspace=0.3, hspace=0.8)
    ax1 = fig.add_subplot(gs[0, 0])
    emf_titration(ax1, Macid, EMF, rgb, AT_emf0, Msamp, Cacid, ATg, '(a)')
    ax2 = fig.add_subplot(gs[1, 0])
    f1_Gran(ax2, Macid, f1g, rgb, ATg, Msamp, Cacid, '(b)')
    ax3 = fig.add_subplot(gs[2, 0])
    emf0_estimate(ax3, Macid, Lg, gemf0, rgb, emf0g, ATg, Msamp, Cacid, '(c)')
    ax4 = fig.add_subplot(gs[3, 0])
    AT_titration(ax4, rgb, Macid, solution, '(d)')
    ax5 = fig.add_subplot(gs[0, 1])
    AT_estimates(ax5, Macid, ATpts, rgb, AT_emf0, RMS, Npts, '(e)')
    ax6 = fig.add_subplot(gs[1:4, 1])
    components(ax6, Macid, ATpts, solution, '(f)')
    return fig
