"""
Example script for generating emission line profiles

Shows for bith simple xwindline, as well as specific Fe-Kalpha case
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../')

from pyxwind import xw_line


def main():

    xlax, feax = define_figure()


    ##################
    # Defining xwind oibject
    ##################
    xw = xw_line(log_mdot_w=-2,
                 r_in=1000,
                 r_out=2*1000,
                 d_foci=np.sqrt(3)*1000,
                 fcov=0.8,
                 log_vinf=-3,
                 rv=0.5*1000,
                 beta=1,
                 kappa=-1,
                 Afe=1)


    #Defining input energy grid
    #Note, these are bin edges! Output grid will have length len(Ebns)-1
    Ebins = np.linspace(5, 8, 10000)
    
    #simple xwindline fist
    ph_xwl = xw.calc_xwindline(inc=60, E0=6, Ebins=Ebins, vturb=100)

    #and now xwindfe
    #also setting to give equivalent width of Fe Kalpha. If this is False, then only returns one value
    ph_fe, Fe_EW = xw.calc_xwindfe(inc=60, N0=1e-2, gamma=1.9, Ebins=Ebins,
                            vturb=100, give_EW=True)


    ###########
    # Plotting
    ##########
    Emid = 0.5*(Ebins[1:] + Ebins[:-1]) #Using linear center since currently linearly spaced bins

    #For xwindline only shoing line prof
    xlax.plot(Emid, ph_xwl, color='seagreen')

    #For xindfe also showing incident continuum
    inspec = 1e-2*Emid**(-1.9)
    feax.plot(Emid, ph_fe, color='dodgerblue', ls='dashed') #line profile
    feax.plot(Emid, inspec, ls='dotted', color='k') #incident spec
    feax.plot(Emid, ph_fe+inspec, color='dodgerblue') #total
    

    feax.set_yscale('log')
    feax.set_ylim(1e-7, 1e-3)

    xlax.set_xlim(5.5, 6.5)
    feax.set_xlim(6, 7)

    xlax.text(s='xwindline', x=0.01, y=0.95, transform=xlax.transAxes)
    feax.text(s='xwindfe', x=0.95, y=0.95, transform=feax.transAxes,
              horizontalalignment='right')
    feax.text(s=f'EW = {Fe_EW:.3f} keV', x=0.95, y=0.9, transform=feax.transAxes,
              horizontalalignment='right')


    xlax.set_ylabel('N(E)   Aribitrary Units')
    xlax.set_xlabel('Energy   (keV)')

    feax.set_xlabel('Energy   (keV)')
    feax.yaxis.tick_right()
    feax.yaxis.set_label_position('right')
    feax.tick_params(axis='y', which='both', left=True, right=True)
    feax.set_ylabel(r'N(E)   photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$')

   
    plt.show()


    return



def define_figure():
    """
    Does figure definitions

    Returns
    -------
    xlax : pyplot axis
        Axis for standard xwindline
    feax : pyplot axis
        Axis for windfe 
    """

    plt.rcParams.update({'font.family':'Times New Roman',
                         'font.size':12,
                         'lines.linewidth':1.5,
                         'mathtext.fontset':'stix',
                         'xtick.direction':'in',
                         'ytick.direction':'in',
                         'xtick.bottom':True,
                         'xtick.top':True,
                         'ytick.left':True,
                         'ytick.right':True})

    fig = plt.figure(figsize=(12, 6))
    grd = plt.GridSpec(1, 2, wspace=0.05)

    xlax = fig.add_subplot(grd[:, 0])
    feax = fig.add_subplot(grd[:, 1])

    return xlax, feax



if __name__ == '__main__':
    main()
