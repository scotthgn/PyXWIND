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
from utils.xwind import windline


class pywindline(PyXWIND_object):
    
    
    def __init__(self,
                 Mbh: float = 1e8,
                 mdot_w: float = 0.1,
                 r_in: float = 500,
                 r_out: float = 1000,
                 d_foci: float = 500,
                 fcov: float = 0.5,
                 vinf: float = 1e-2,
                 rv: float = 100,
                 beta: float = 1,
                 kappa: float = -1,
                 Afe: float = 1):
        """
        Parameters
        ----------
        Mbh : float
            Black hole mass.
            UNITS: Msol
            DEFAULT: 1e8.
        mdot_w : float
            Eddington scaled mass-outflow rate
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
        vinf : float
            Outflow velcotiy at infinity
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
        
        super().__init__(Mbh, mdot_w, r_in, r_out, d_foci, fcov, vinf, rv,
                         beta, kappa, Afe)


    def calc_windline(self, inc: float, E0: float, Ebins: ArrayLike,
                      edens_unit: bool = True) -> ArrayLike:
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
        
        wpars = [self.Mbh, self.mdot_w, self.r_in, self.r_out, self.d_foci,
                 self.fcov, self.vinf, self.rv, self.beta, self.kappa, inc, 
                 self.Afe, E0, 1, 2.0]
        
        ph = windline(wpars, Ebins, True, len(Ebins)-1) #ph/s/cm^2/bin
        
        if edens_unit:
            ph = ph/(Ebins[1:] - Ebins[:-1]) #ph/s/cm^2/keV
        
        return ph
    
    
    def calc_windfe(self, inc: float, N0: float, gamma: float, Ebins: ArrayLike,
                    edens_unit: bool = True) -> ArrayLike:
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
        edens_units : bool, optional
            If True, converts output array to photons/s/cm^2/keV
            DEFAULT: True

        Returns
        -------
        ph : ArrayLike, shape=len(Ebins)-1
            Output line spectrum
            Units : photons/s/cm^2/keV if edens_units==True (DEFAULT)
                    photons/s/cm^2/bin if edens_units==False

        """
        
        wpars = [self.Mbh, self.mdot_w, self.r_in, self.r_out, self.d_foci,
                 self.fcov, self.vinf, self.rv, self.beta, self.kappa, inc, 
                 self.Afe, 6.4, N0, gamma]
        
        ph = windline(wpars, Ebins, False, len(Ebins)-1) #ph/s/cm^2/bin
        
        if edens_unit:
            ph = ph/(Ebins[1:] - Ebins[:-1]) #ph/s/cm^2/keV
        
        return ph
        

if __name__ == '__main__':
    
    #xw = pywindline(r_in=40, kappa=1)
    #xw.calc_windline(45, 6.4, np.linspace(4.5, 7.5, 1000))#
    
    #looking for the seg fautl!
    wnd = pywindline(1e8, 0.1, 10, 1500, 1000, 0.6, 1e-2, 500, 1, -1)
    ear = np.linspace(4, 7, 1000)
    ph = wnd.calc_windline(45, 6.0, Ebins=ear, edens_unit=False)
    
    plt.plot(0.5*(ear[1:] + ear[:-1]), ph)
    plt.show()
    
    #print(ph)
    #print(len(ph), len(ear))
    
    
    