"""
Example script for taking some input spectrum and convolving with XWIND
using xw_conv

As an illustration this script uses an input of two separated narrow gaussians
Gaussian components. But note that xw_conv can in theory be used for
continuum as well (under certain assumptions...)
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
from pyxwind import xw_conv


def main():


    ########################
    #Initiating xwind object
    ########################
    xw = xw_conv(log_mdot_w=-2,
                 r_in=1000,
                 r_out=1.5*1000,
                 d_foci=np.sqrt(3)*1000,
                 fcov=0.8,
                 log_vinf=-2.5,
                 rv=1000,
                 beta=1,
                 kappa=-1)

    
    ########################
    #Generating input spec
    ########################
    ebins = np.geomspace(0.1, 10, 10000) #bin edges, keV
    emids = 10**(np.log10(ebins[1:]*ebins[:-1])/2) #Geometric centre

    ph_input = np.zeros_like(emids)
    ph_input += gau(emids, 1, 1e-3, 0.5) + gau(emids, 1, 1e-3, 5) 

    xw.set_input_spectrum(ear=ebins, ph=ph_input, eunit='keV', phunit='photon_keV') #adding input spectrum to xwind object


    ########################
    #Doing convolution
    ########################
    ph_out = xw.do_convolve(inc=45, vturb=100)


    ########################
    #Plotting
    ########################
    #Note, for a very quick plot, just use: xw.plot_spec()
    #Doing manualy here for slightly more controll

    plt.rcParams.update({'font.family':'Times New Roman',
                         'font.size':12,
                         'lines.linewidth':1,
                         'mathtext.fontset':'stix',
                         'xtick.direction':'in',
                         'ytick.direction':'in',
                         'xtick.bottom':True,
                         'xtick.top':True,
                         'ytick.left':True,
                         'ytick.right':True})

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)


    ax.loglog(emids, ph_input, color='k', ls='dotted') #input spec
    ax.loglog(emids, ph_out, color='blue') #convolved spec

    #limits
    ax.set_ylim(1e-3, 2)
    ax.set_xlim(0.3, 7)

    #labels
    ax.set_xlabel('Energy   (keV')
    ax.set_ylabel(r'photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$')

    #Making zoomed insets to highlight convolution
    #low energy gaussian first
    axin1 = ax.inset_axes([0.2, 0.5, 0.45, 0.45], xlim=(0.48,0.52), ylim=(0.05, 1.2))
    axin1.loglog(emids, ph_input, color='k', ls='dotted') #input spec
    axin1.loglog(emids, ph_out, color='blue') #convolved spec
    axin1.tick_params(axis='both', which='both', labelsize=0, labelcolor='white')

    ax.indicate_inset_zoom(axin1, edgecolor='k')

    #higher energy gaussian
    axin2 = ax.inset_axes([0.3, 0.02, 0.45, 0.45], xlim=(4.8,5.2), ylim=(8e-3, 1.2))
    axin2.loglog(emids, ph_input, color='k', ls='dotted') #input spec
    axin2.loglog(emids, ph_out, color='blue') #convolved spec
    axin2.tick_params(axis='both', which='both', labelsize=0, labelcolor='white')

    ax.indicate_inset_zoom(axin2, edgecolor='k')

    plt.show()


    return



def gau(emids, A, s0, E0):
    """
    Simple Gaussian

    Parameters
    ----------
    emids : array
        Energy bin midpoints 
        Units : keV
    A : float
        Amplitude
    s0 : float
        Standard deviation (width)
        Units : keV
    E0: float
        Central energy
        Units : keV 

    Returns
    -------
    ph : array
        Gaussian photon flux array
    """

    ph = A*np.exp(-0.5*((emids-E0)/s0)**2)
    return ph




if __name__ == '__main__':
    main()
