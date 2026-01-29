"""
Example script for initating an xwind model object, and then extracting the corresponding
density, fluoresence, and velocity profile

For details on XWIND, see Hagen et al. (2026, submitted)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_cols

#adding xwind to pythonpath
import os
import sys
sys.path.append('../')

from pyxwind import xw_line #Using xw_line in this example. But would be identical for xw_conv


def minimal_example():
    """
    Shows how to quickly plot the wind profiles, simply leaving the code defulats 
    for normalisations etc.

    Note how this has minimal information (colourbar units, etc)
    """



    dax, eax, vlax, vpax = define_figure()
    

    #######################
    # Defining XWIND object
    #######################
    xw = xw_line(Mbh=1e8,     #Mass is only required to give density profile! Irrelevant for line
                 log_mdot_w=-2,
                 r_in=1000,
                 r_out=2*1000,
                 d_foci=np.sqrt(3)*1000,
                 fcov=0.8,
                 log_vinf=-2,
                 rv=0.5*1000,
                 beta=1,
                 kappa=-1,
                 Afe=1)

    
    #######################
    # Plotting wind profiles
    #######################
    #density
    xw.plot_wind('ndens', ax=dax)    
    #relative emissivity
    xw.plot_wind('rel_vol_emiss', ax=eax)
    #Streamline velocity
    xw.plot_wind('vl', ax=vlax)
    #Azimuthal velocity
    xw.plot_wind('vphi', ax=vpax)

    plt.show()




    return


def detailed_example():
    """
    Still using XWIND methods, does a slightly better version of the same plot as in 
    minimal_example()
    """

    dax, eax, vlax, vpax = define_figure()

    #######################
    # Defining XWIND object
    #######################
    xw = xw_line(Mbh=1e8,     #Mass is only required to give density profile! Irrelevant for line
                 log_mdot_w=-2,
                 r_in=1000,
                 r_out=2*1000,
                 d_foci=np.sqrt(3)*1000,
                 fcov=0.8,
                 log_vinf=-2,
                 rv=4*1000,
                 beta=1,
                 kappa=-1,
                 Afe=1)

    
    #######################
    # Plotting wind profiles
    #######################
    #density
    #Can define extermal normalisation (e.g. power law rather than log) by using cscale argument
    xw._calc_ndens() #Getting density solution, so know where to set limits
    dens_norm = mpl_cols.PowerNorm(vmin=np.log10(np.amin(xw.ndens)), vmax=np.log10(np.amax(xw.ndens)), gamma=0.4)
    xw.plot_wind('ndens', ax=dax, cmap='jet', cscale=dens_norm, inc_cbar=True,
                  cbar_loc='left', show_axlabel=False)   
 
    #relative emissivity 
    xw.plot_wind('rel_vol_emiss', ax=eax, cmap='hot', cnorm='log', inc_cbar=True,
                 cbar_loc='right', show_axlabel=False, cbar_label='Normalised Volume Emissivity') #can also set own cbar labels

    #Streamline velocity
    xw.plot_wind('vl', ax=vlax, cmap='viridis', cnorm='log', inc_cbar=True,
                  cbar_loc='left', show_axlabel=False)
    #Azimuthal velocity
    xw.plot_wind('vphi', ax=vpax, cmap='viridis', cnorm='log', inc_cbar=True,
                  cbar_loc='right', show_axlabel=False, cbar_label=r'$v_{\phi}/c$')


    plt.show()


def define_figure():
    """
    Stays constant for all examples

    Returns
    -------
    dax : pyplot axis
        Density axis
    eax : pyplot axis
        Relative emissivity axis
    vlax : pyplot axis
        Streamline velocity axis
    vpax : pyplot axis
        Azimuthal velocity axis
    """
    
    plt.rcParams.update({'font.family':'Times New Roman',
                         'font.size':12,
                         'mathtext.fontset':'stix',
                         'xtick.direction':'in',
                         'ytick.direction':'in',
                         'xtick.bottom':True,
                         'xtick.top':True,
                         'ytick.left':True,
                         'ytick.right':True})

    fig = plt.figure(figsize=(10, 8)) #doing 2x2 panels
    grd = plt.GridSpec(2, 2, hspace=0, wspace=0)

    dax = fig.add_subplot(grd[0, 0]) #density
    eax = fig.add_subplot(grd[0, 1], sharex=dax, sharey=dax) #emissivity
    vlax = fig.add_subplot(grd[1, 0], sharex=dax, sharey=dax) #Velocity along the streamline
    vpax = fig.add_subplot(grd[1, 1], sharex=dax, sharey=dax) #Azimuthal velocity

    ###
    #labels etx
    ###
    dax.tick_params(axis='both', which='both', labelsize=0, labelcolor='white')
    eax.tick_params(axis='both', which='both', labelsize=0, labelcolor='white')
    vlax.tick_params(axis='y', which='both', labelsize=0, labelcolor='white')
    vpax.tick_params(axis='y', which='both', labelsize=0, labelcolor='white')

    vlax.set_xlabel(r'x   ($R_{G}$)')  
    vpax.set_xlabel(r'x   ($R_{G}$)')   


    return dax, eax, vlax, vpax



if __name__ == '__main__':
    #minimal_example()
    detailed_example()


