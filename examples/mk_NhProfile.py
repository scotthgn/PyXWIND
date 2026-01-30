"""
Demonstrates how to generate an Nh profile for a gien wind
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
from pyxwind import xw_line


def main():


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


    #getting column-density profile
    NHarr, theta_arr = xw.calc_NH_profile()

    #if you want just a single line of sight, use calc_los_NH. e.g.
    NH_los = xw.calc_los_NH(inc=60)
    print(f'Nh at inc=60 deg: {NH_los} cm^-2')


    #############################
    #Plotting
    #############################
    plt.rcParams.update({'font.family':'Times New Roman',
                         'font.size':12,
                         'mathtext.fontset':'stix',
                         'xtick.direction':'in',
                         'ytick.direction':'in',
                         'xtick.bottom':True,
                         'xtick.top':True,
                         'ytick.left':True,
                         'ytick.right':True})
    
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)

    ax.plot(theta_arr, NHarr, color='dodgerblue', lw=1.5)
    ax.set_yscale('log')

    ax.set_xlabel('Line of sight   (deg)')
    ax.set_ylabel(r'Column-Density, $N_H$   (cm$^{-2}$)')
    
    plt.show()


    return




if __name__ == '__main__':
    main()
