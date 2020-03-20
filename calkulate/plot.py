# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Visualise Calkulate calculations."""
from numpy import array, log10, logical_and, mean, sqrt, zeros
from numpy import any as np_any
from numpy import min as np_min
from numpy import max as np_max
from matplotlib.pyplot import figure, rcParams, subplots_adjust
from . import datfile, simulate, solve

_rgb_guess = array([0.96, 0.86, 0.04])
_rgb_final = array([0.21, 0.46, 1])
_rgb_both = array([0.27, 0.8, 0.54])

def prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, concAcid, WhichKs=10, WhoseKSO4=1, WhoseKF=1, WhoseTB=2,
        totalAmmonia=0, totalH2Sulfide=0):
    """Preparatory calculations for plotting."""
    massAcid, emf, tempK, massSample, concTotals, eqConstants = \
        datfile.prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF,
        WhoseTB=WhoseTB, totalAmmonia=totalAmmonia,
        totalH2Sulfide=totalH2Sulfide)
    f1Guess = solve.f1(massAcid, emf, tempK, massSample)
    LGuess = logical_and(
        f1Guess > 0.1*np_max(f1Guess),
        f1Guess < 0.9*np_max(f1Guess),
    )
    alkGuess, emf0Guess, _, pHGuess = solve.guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    granEmf0 = solve.granEmf0Guess(massAcid[LGuess], emf[LGuess],
        tempK[LGuess], massSample, concAcid, alkGuess)
    L = logical_and(pHGuess > 3, pHGuess < 4)
    alk_emf0 = solve.complete(massAcid, emf, tempK, massSample, concAcid,
        concTotals, eqConstants)
    h = solve.emf2h(emf, alk_emf0['x'][1], tempK)
    mu = solve.mu(massAcid, massSample)
    alkSim = simulate.alk(h, mu, concTotals, eqConstants)
    alk0Sim = (alkSim[0] + massAcid*concAcid/(massAcid + massSample))/mu
    RMS = sqrt(mean(alk_emf0['fun']**2))
    Npts = len(alk_emf0['fun'])
    rgb = zeros((len(emf), 3))
    for i in range(len(rgb)):
        if LGuess[i] and L[i]:
            rgb[i] = _rgb_both
        elif LGuess[i] and not L[i]:
            rgb[i] = _rgb_guess
        elif L[i] and not LGuess[i]:
            rgb[i] = _rgb_final
        else:
            rgb[i] = array([1, 1, 1])
    return (massAcid, emf, massSample, f1Guess, LGuess, alkGuess, emf0Guess,
        granEmf0, alk_emf0, alkSim, alk0Sim, RMS, Npts, rgb)

def emfTitration(ax, massAcid, emf, massSample, concAcid, alk_emf0, alkGuess,
        rgb, sublabel):
    """EMF change as acid is added throughout a titration."""
    ax.axvline(1e3*alk_emf0['x'][0]*massSample/concAcid, color=_rgb_final,
        linestyle='--', zorder=1)
    ax.axvline(1e3*alkGuess*massSample/concAcid, color=_rgb_guess,
        linestyle='--', zorder=1)
    ax.scatter(massAcid*1e3, emf, c=rgb, edgecolors='k', clip_on=False,
        zorder=2)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    yrange = np_max(emf) - np_min(emf)
    ax.set_ylim([np_min(emf) - yrange*0.05, np_max(emf) + yrange*0.05])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('EMF / mV')
    ax.set_title('{} Final EMF$^\circ$ = {:.2f} mV'.format(sublabel,
        alk_emf0['x'][1]), fontsize=10)
    return ax

def f1Gran(ax, massAcid, massSample, concAcid, f1Guess, alkGuess, rgb,
        sublabel):
    """F1 Gran plot function for the first alkalinity estimate."""
    ax.axvline(1e3*alkGuess*massSample/concAcid, color=_rgb_guess,
        linestyle='--', zorder=1)
    ax.scatter(massAcid*1e3, f1Guess*1e-7, c=rgb, edgecolors='k',
        clip_on=False, zorder=2)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    ax.set_ylim([0, np_max(f1Guess*1.05e-7)])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$F_1 \cdot 10^{-7}$')
    ax.set_title('{} First-guess alkalinity = {:.1f} μmol/kg'.format(sublabel,
        alkGuess*1e6), fontsize=10)
    return ax

def emf0Estimates(ax, massAcid, massSample, concAcid, granEmf0, emf0Guess,
         alkGuess, LGuess, rgb, sublabel):
    """First estimate of EMF0."""
    ax.axhline(emf0Guess, color=_rgb_guess, zorder=1)
    ax.axvline(1e3*alkGuess*massSample/concAcid, color=_rgb_guess,
        linestyle='--', zorder=1)
    ax.scatter(massAcid[LGuess]*1e3, granEmf0, c=rgb[LGuess], edgecolors='k',
        clip_on=False, zorder=2)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('First-guess EMF$^\circ$ / mV')
    ax.set_title('{} First-guess EMF$^\circ$ = {:.1f} mV'.format(sublabel,
        emf0Guess), fontsize=10)
    return ax

def alkTitration(ax, massAcid, alkSim, rgb, sublabel):
    """Linear decrease in alkalinity as acid is added during a titration."""
    ax.scatter(massAcid*1e3, alkSim[0]*1e6, clip_on=False, edgecolors='k',
        c=rgb)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    ax.set_ylabel('Alkalinity / μmol kg$^{-1}$')
    ax.set_xlabel('Acid mass / g')
    ax.set_title(sublabel, fontsize=10)
    return ax

def alkEstimates(ax, massAcid, alk0Sim, rgb, alk_emf0, RMS, Npts, sublabel):
    """Original sample alkalinity estimated from each titration point."""
    ax.axhline(alk_emf0['x'][0]*1e6, color=_rgb_final, zorder=1)
    ax.scatter(massAcid*1e3, alk0Sim*1e6, c=rgb, edgecolors='k', clip_on=False,
        zorder=2)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    yrange = (np_max(alk0Sim) - np_min(alk0Sim))*1e6
    ax.set_ylim([np_min(alk0Sim*1e6) - yrange*0.05,
        np_max(alk0Sim*1e6 + yrange*0.05)])
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$A_\mathrm{T}$ from pH / μmol$\cdot$kg$^{-1}$')
    ax.set_title(('{} Final alkalinity = ({:.1f} $\pm$ {:.1f}) μmol/kg' +
        ' ($n$ = {})').format(sublabel, alk_emf0['x'][0]*1e6, RMS*1e6, Npts),
        fontsize=10)
    return ax

def components(ax, massAcid, alk0Sim, alkSim, sublabel):
    """Every component of alkalinity throughout a titration."""
    ax.plot(massAcid*1e3, -log10(alk0Sim), label='Total alk.',
        marker='o', markersize=3, c='k', clip_on=False)
    # Keys in rgbs_names match those in simulate.alk's components dict
    rgbs_names = {
        '+HCO3': ['xkcd:grey', '+[HCO$_3^-$]'],
        '+2*CO3': ['xkcd:grey', '+2[CO$_3^{2-}$]'],
        '+B(OH)4': ['xkcd:blue', '+[B(OH)$_4^-$]'],
        '+OH': ['xkcd:red', '+[OH$^-$]'],
        '-H': ['xkcd:red', '$-$[H$^+$]'],
        '-HSO4': ['xkcd:golden yellow', '$-$[HSO$_4^-$]'],
        '-HF': ['xkcd:green', '$-$[HF]'],
        '-H3PO4': ['xkcd:orange', '$-$[H$_3$PO$_4$]'],
        '+HPO4': ['xkcd:orange', '+[HPO$_4^{2-}$]'],
        '+2*PO4': ['xkcd:orange', '+2[PO$_4^{3-}$]'],
        '+SiO(OH)3': ['xkcd:pink', '+[SiO(OH)$_3$]'],
        '+NH3': ['xkcd:aquamarine', '+[NH$_3$]'],
        '+HS': ['xkcd:indigo', '+[HS$^-$]']
    }
    for component in alkSim[1].keys():
        if component.startswith('-'): # this is a bit sketchy
            yVar = -alkSim[1][component]
        else:
            yVar = alkSim[1][component]
        if np_any(yVar > 0):
            ax.plot(massAcid*1e3, -log10(yVar), label=rgbs_names[component][1],
                marker='x', markersize=3, c=rgbs_names[component][0],
                clip_on=False)
    ax.set_xlim([0, np_max(massAcid)*1e3])
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_xlabel('Acid mass / g')
    ax.set_ylabel('$-$log$_{10}$(concentration from pH / mol$\cdot$kg$^{-1}$)')
    ax.set_title(sublabel, fontsize=10)
    return ax

def everything(datfile, volSample, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, concAcid, WhichKs=10, WhoseKSO4=1, WhoseKF=1, WhoseTB=2,
        totalAmmonia=0, totalH2Sulfide=0):
    """Plot everything for a single titration from a VINDTA-style .dat file."""
    (massAcid, emf, massSample, f1Guess, LGuess, alkGuess, emf0Guess, granEmf0,
        alk_emf0, alkSim, alk0Sim, RMS, Npts, rgb) = prep(datfile, volSample,
        pSal, totalCarbonate, totalPhosphate, totalSilicate, concAcid,
        WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF, WhoseTB=WhoseTB,
        totalAmmonia=totalAmmonia, totalH2Sulfide=totalH2Sulfide)
    fig = figure(figsize=[17, 10])
    rcParams.update({'font.size': 10})
    gs = fig.add_gridspec(4, 2)
    subplots_adjust(wspace=0.3, hspace=0.8)
    ax1 = fig.add_subplot(gs[0, 0])
    emfTitration(ax1, massAcid, emf, massSample, concAcid, alk_emf0, alkGuess,
        rgb, '(a)')
    ax2 = fig.add_subplot(gs[1, 0])
    f1Gran(ax2, massAcid, massSample, concAcid, f1Guess, alkGuess, rgb, '(b)')
    ax3 = fig.add_subplot(gs[2, 0])
    emf0Estimates(ax3, massAcid, massSample, concAcid, granEmf0, emf0Guess,
        alkGuess, LGuess, rgb, '(c)')
    ax4 = fig.add_subplot(gs[3, 0])
    alkTitration(ax4, massAcid, alkSim, rgb, '(d)')
    ax5 = fig.add_subplot(gs[0, 1])
    alkEstimates(ax5, massAcid, alk0Sim, rgb, alk_emf0, RMS, Npts, '(e)')
    ax6 = fig.add_subplot(gs[1:4, 1])
    components(ax6, massAcid, alk0Sim, alkSim, '(f)')
    return fig
