#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 10:04:37 2025

@author: Scott Hagen

Python wrapper for the windline subroutine
"""

import numpy as np
import matplotlib.pyplot as plt

from numpy.typing import ArrayLike

from _pyxwind_ObjectHandler import PyXWIND_object
from utils.xwind import xwindline, xwindfe


class xw_line(PyXWIND_object):
    
    
    def __init__(self,
                 Mbh: float = 1e8,
                 log_mdot_w: float = -3,
                 r_in: float = 500,
                 r_out: float = 1000,
                 d_foci: float = 500,
                 fcov: float = 0.5,
                 log_vinf: float = -3,
                 rv: float = 100,
                 beta: float = 1,
                 kappa: float = -1,
                 Afe: float = 1):
        """
        Parameters
        ----------
        Mbh : float
            Black hole mass.
            NOTE - does not affect line-profile. Only needed for wind property
            profiles in physical units
            UNITS: Msol
            DEFAULT: 1e8.
        log_mdot_w : float
            log10 Eddington scaled mass-outflow rate
            UNITS: Mdot/Mdot_edd
            DEFAULT: 0.1
        r_in : float
            Inner launch radius
            UNITS: Rg
            DEFAULT: 500
        r_out : float
            Outer launch radius
            UNITS: Rg
            DEFAULT: 1000
        d_foci : float
            Distance below the origin to wind foci
            UNITS: Rg
            DEFAULT: 500
        fcov : float
            Covering fraction of wind, including both sides of the disc
            UNITS: Omega/4pi
            DEFAULT: 0.5
        log_vinf : float
            log10 Outflow velcotiy at infinity
            UNITS: c
            DEFAULT: 1e-2
        rv : flaot
            Velocity scale length. Distance along streamline where the velocity
            reaches 0.5 vinf
            UNITS: Rg
            DEFAULT: 100
        beta : float
            Velocity law exponent. Determines acceleration along streamline
            UNITS: dimensionless
            DEFAULT: 1
        kappa : float
            Wind density gradient exponent at surface.
            Defined as dMdot/dR propto R^(-kappa)
            IF negative, then treated as log number density for a constant 
            density wind
            UNITS: dimensionless IF >= 0
                   log cm^-3     IF < 0
            Default: 1
        Afe : float
            Iron abundance, relative to solar
            Units : [Fe/Fe_sol]
            Default: 1
        """
        
        super().__init__(Mbh, log_mdot_w, r_in, r_out, d_foci, fcov, log_vinf, rv,
                         beta, kappa, Afe)


    def calc_xwindline(self, inc: float, E0: float, Ebins: ArrayLike,
                       vturb: float = 0, edens_unit: bool = True) -> ArrayLike:
        """
        Calculates windline profile.
        Normalisation set s.t the total photon flux integrates to unity
        
        Parameters
        ----------
        inc : float
            Observer inclination, w.r.t to the z-axis
            (i.e i=0 implies face on, i=90 implies edge on)
            Units : deg
        E0 : float
            Line rest frame energy
            Units : keV
        Ebins : ArrayLike
            Output energy bin edges
            (note internal energy reolution is set to E/E0 = 5e-4)
            Units : keV
        vturb : float
            Turbulent velocity. Gives smoothing of line profile via 
            convolution with a Gaussian
            Units : km/s
        edens_unit : bool
            If True, converts output array to photons/s/cm^2/keV
            DEFAULT: True
            

        Returns
        -------
        ph : ArrayLike, shape=len(Ebins)-1
            Output line spectrum
            Units : photons/s/cm^2/keV if edens_units==True (DEFAULT)
                    photons/s/cm^2/bin if edens_units==False

        """
        
        wpars = [self.mdot_w, self.r_in, self.r_out, self.d_foci,
                 self.fcov, self.vinf, self.rv, self.beta, vturb, self.kappa, inc, E0]
        
        ph = xwindline(wpars, Ebins, len(Ebins)-1) #ph/s/cm^2/bin
        
        if edens_unit:
            ph = ph/(Ebins[1:] - Ebins[:-1]) #ph/s/cm^2/keV
        
        return ph
    
    
    def calc_xwindfe(self, inc: float, N0: float, gamma: float, Ebins: ArrayLike,
                     vturb: float = 0, edens_unit: bool = True, give_EW: bool = False) -> ArrayLike:
        """
        Windline profile specifically for Fe-Kalpha
        Normalisation calculated from the wind density/fluoresence profile

        Parameters
        ----------
        inc : float
            Observer inclination, w.r.t to the z-axis
            (i.e i=0 implies face on, i=90 implies edge on)
            Units : deg
        N0 : float
            Normalisation of illuminating X-ray spectrum
            Defined as: photons/s/cm^2 at 1 keV
        gamma : float
            Spectral index for incident power-law spectrum
        Ebins : ArrayLike
            Output energy bin edges
            (note internal energy reolution is set to E/E0 = 5e-4)
            Units : keV
        vturb : float
            Turbulent velocity. Gives smoothing of line profile via 
            convolution with a Gaussian
            Units : km/s
        edens_units : bool, optional
            If True, converts output array to photons/s/cm^2/keV
            DEFAULT: True
        give_EW : bool, optional
            If True, then evaluates and also returns the equivalent width of the 
            line in units of keV

        Returns
        -------
        ph : ArrayLike, shape=len(Ebins)-1
            Output line spectrum
            Units : photons/s/cm^2/keV if edens_units==True (DEFAULT)
                    photons/s/cm^2/bin if edens_units==False

        """
        
        wpars = [self.mdot_w, self.r_in, self.r_out, self.d_foci,
                 self.fcov, self.vinf, self.rv, self.beta, vturb, self.kappa, inc, 
                 self.Afe, N0, gamma]
        
        ph = xwindfe(wpars, Ebins, len(Ebins)-1) #ph/s/cm^2/bin
        
        if give_EW:
            ph_tot = np.sum(ph)
            cont = N0 * 6.4**(-gamma)
            EW = ph_tot/cont #keV
        
        if edens_unit:
            ph = ph/(Ebins[1:] - Ebins[:-1]) #ph/s/cm^2/keV
        
        if give_EW:
            return ph, EW
        else:
            return ph
    
        
        

if __name__ == '__main__':
    
    #testing
    #wnd = xw_line(1e8, 0.1, 1000, 10000, 2000, 0.8, 1e-2, 500, 1, -1)
    wnd = xw_line(2e7, -4, 2022.22, 5055.55, 3502.59, 0.4998, -4,
                  200, 1, 1)
    
    ear = np.linspace(5, 7, 1000)
    
    
    #line profile
    #ph = wnd.calc_xwindline(45, 6.0, Ebins=ear, edens_unit=True)
    ph = wnd.calc_xwindfe(45, 8.496e-2, 1.8, Ebins=ear, edens_unit=True)
    
    fig_line = plt.figure(figsize=(7, 7))
    ax_line = fig_line.add_subplot(111)
    ax_line.plot(0.5*(ear[1:] + ear[:-1]), ph)
    ax_line.set_yscale('log')
    
    """
    
    #wind profiles
    fig_prof = plt.figure(figsize=(8, 8))
    ax_wnd = fig_prof.add_subplot(111)
    wnd.plot_wind('ndens', show=False, ax=ax_wnd)
    
    fig_lum = plt.figure(figsize=(8, 8))
    ax_lum = fig_lum.add_subplot(111)
    wnd.plot_wind('vol_emiss', show=False, ax=ax_lum, cmap='hot')
    

    #Nh
    nh = wnd.calc_los_NH(90)
    NHprof, thetas = wnd.calc_NH_profile()
    
    fig_nh = plt.figure(figsize=(7, 7))
    ax_nh = fig_nh.add_subplot(111)
    ax_nh.plot(thetas, NHprof)
    ax_nh.set_yscale('log')
    
    """
    plt.show()
    
    
    
    #print(ph)
    #print(len(ph), len(ear))
    
    
    