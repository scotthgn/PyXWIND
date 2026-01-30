#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 12:13:09 2025

@author: Scott Hagen

Python wrapper for windconv subroutine
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from numpy.typing import ArrayLike

from ._pyxwind_ObjectHandler import PyXWIND_object
from .utils.xwind import xwindconv


class xw_conv(PyXWIND_object):
    """
    Convolution kernel for windline.
    
    For a given wind object and input spectrum, does the convoluiton between
    the windline profile and the input spectrum
    """
    
    def __init__(self,
                 Mbh: float = 1e8,
                 log_mdot_w: float = -3,
                 r_in: float = 1000,
                 r_out: float = 1500,
                 d_foci: float = 1700,
                 fcov: float = 0.5,
                 log_vinf: float = -3,
                 rv: float = 1000,
                 beta: float = 1,
                 kappa: float = -1):
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
            DEFAULT: 1000
        r_out : float
            Outer launch radius
            UNITS: Rg
            DEFAULT: 1500
        d_foci : float
            Distance below the origin to wind foci
            UNITS: Rg
            DEFAULT: 1700
        fcov : float
            Covering fraction of wind, including both sides of the disc
            UNITS: Omega/4pi
            DEFAULT: 0.5
        log_vinf : float
            Outflow velcotiy at infinity
            UNITS: c
            DEFAULT: 1e-2
        rv : flaot
            Velocity scale length. Distance along streamline where the velocity
            reaches 0.5 vinf
            UNITS: Rg
            DEFAULT: 1000
        beta : float
            Velocity law exponent. Determines acceleration along streamline
            UNITS: dimensionless
            DEFAULT: 1
        kappa : float
            Wind density gradient exponent at surface.
            Defined as dMdot/dR propto R^(-kappa)
            Default: 1
        """
        
        super().__init__(Mbh, log_mdot_w, r_in, r_out, d_foci, fcov, log_vinf, rv,
                         beta, kappa)
        
        return
    
    
    def set_input_spectrum(self, ear: ArrayLike, ph: ArrayLike, 
                           eunit: str = 'keV', phunit: str = 'photon_keV',
                           check_unit: bool = True):
        """
        Defines the input spectrum that is convolved with the line profile.
        
        Internally the code uses units of photons/s/cm^2/bin. For a given 
        eunit and phunit this will automaatically convert to internal units.

        NOTE: ONLY reasonable spectral units are accepted. A list of these is
        given in the deiscription of eunit and phunit respectively (below)        

        Parameters
        ----------
        ear : ArrayLike
            Energy bin edges
            UNIT : given by eunit param
        ph : ArrayLike
            Photon/energy flux
            Must have length len(ear) - 1
            UNIT : given by phunit par
        eunit : str, {'keV', 'Hz', 'AA'}
            Unit of input energy bins. Here AA refers to Angstrøm
            DEFUALT : keV
        phunit : str, {'photon_keV', 'photon_Hz', 'photon_AA',
                       'cgs_keV', 'cgs_Hz', 'cgs_AA',
                       'SI_keV', 'SI_Hz', 'SI_AA'}
            Spectral unit. MUST be in a form of flux density. The input string 
            is specified in two parts, seperated by an underscore (_).
            The first part details the flux-energy unit, and are assumed to be
            in their typical combinations. ie:
                - photon = photon/s/cm^2
                - cgs = erg/s/cm^2
                - SI = W/m^2
            
            The second part indicates the spectral-energy unit:
                - keV
                - Hz
                - Angstrom
            
            Thus photon_keV = photon/s/cm^2/keV
                 cgs_Hz = erg/s/cm^2/Hz
                 SI_AA = W/m^2/AA
                 etc.
            
            DEFAULT: photon_keV
        check_unit : bool, optional
            If True then checks input units and converts to internal units.
            This should ALWAYS be set to True, unless you are passing input
            arrays in units of keV (for ear) and photons/s/cm^2/bin (for ph).
            Given as a parameter s.t long intensive fit routines can skip an
            unecessary unit check by the user handling this externally
            DEFAULT : True

        Returns
        -------
        None.

        """
        
        self.ear_input, self.ph_input = ear, ph
        self._convert_inputSpec(ear, ph, eunit, phunit)
        
        #plt.plot(0.5*(self.ear[1:]+self.ear[:-1]), self.ph)
        #plt.show()
        
        return
    
    
    def do_convolve(self, inc: float, vturb: float = 0, conv_unit: bool = True):
        """
        

        Parameters
        ----------
        inc : float
            Observer inclination, measured w.r.t the z-axis
            UNITS : deg
        conv_unit : bool, optional
            If True, will convert output spectrum back to the user input units
            If False, output is given in photons/s/cm^2/bin

        Returns
        -------
        None.

        """
        
        #checking input spectrum is defined
        if hasattr(self, 'ph'):
            pass
        else:
            raise AttributeError('Input spectrum does not exist! Define this'
                                 ' using set_input_spectrum(ear, ph, eunit, phunit) method')
        
        
        wpars = [self.mdot_w, self.r_in, self.r_out, self.d_foci,
                 self.fcov, self.vinf, self.rv, self.beta, vturb, self.kappa, inc]
        
        
        #print()
        #photon/s/cm^2/bin
        self.ph_conv = xwindconv(wpars, self.ear, self.ph, len(self.ear)-1)
        
        #print()
        #print(self.ph_conv)
        #print()
        
        if conv_unit:
            self.ph_conv_usr = self._convert_outputSpec(self.ph_conv)
            return self.ph_conv_usr
        else:
            return self.ph_conv
        
    ###################################################################
    #---- PLOTTING
    ###################################################################    
    
    def plot_spec(self, 
                  ax: None|type(plt.axis) = None,
                  as_log: bool = True,
                  show_inspec: bool = True,
                  show: bool = False):
        """
        Plots convolved spectrum (solid) and input spec (if asked, dotted)
        Does plot in user defined units

        Returns
        -------
        None.

        """
        
        if hasattr(self, 'ph_conv_usr'):
            pass
        else:
            if hasattr(self, 'ph_conv'):
                self.ph_conv_usr = self._convert_outputSpec(self.ph_conv)
            else:
                raise AttributeError('Convolved model spectrum does not exist!'
                                     ' Run do_convolve(inc) method first!')
        
        
        if ax is None:
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111)
        
        if show_inspec:
            ax.plot(0.5*(self.ear_input[1:] + self.ear_input[:-1]), self.ph_input,
                    ls='dotted', color='k')
        
        ax.plot(0.5*(self.ear_input[1:] + self.ear_input[:-1]), self.ph_conv_usr,
                color='blue')
        
        if as_log:
            ax.set_xscale('log')
            ax.set_yscale('log')
            
        ax.set_xlabel(self.eunit)
        if self.phunit_flx.lower() == 'photon':
            ax.set_ylabel(r'photon/s/cm$^{}$/{}'.format(2, self.phunit_en))
        elif self.phunit_flx.lower() == 'cgs':
            ax.set_ylabel(r'erg/s/cm$^{}$/{}'.format(2, self.phunit_en))
        elif self.phunit_flx.lower() == 'si':
            ax.set_ylabel(r'W/m$^{}$/{}'.format(2, self.phunit_en))
        else:
            raise ValueError
    
        if show:
            plt.show()
    
        
    ###################################################################
    #---- UNIT HANDLING
    ###################################################################
    
    def _convert_outputSpec(self, ph_out: ArrayLike) -> ArrayLike:
        """
        Converts output (i.e convolved) model spectrum from internal units
        (photons/s/cm^2/bin) to whatever the user defined upon input

        Parameters
        ----------
        ph_out : ArrayLike
            Model photon spectrum

        Returns
        -------
        ph_usr : ArrayLike
            Model specturm in user defined units
        """
        
        ph_usr = ph_out / (self.ear[1:] - self.ear[:-1]) #ph/s/cm^2/keV
        
        if self.phunit_en.lower() == 'kev':
            pass
        elif self.phunit_en.lower() == 'hz':
            ph_usr = (ph_usr*u.ph/u.s/u.cm**2/u.keV).to(
                u.ph/u.s/u.cm**2/u.Hz, equivalencies=u.spectral()) #ph/s/cm^2/Hz
        elif self.phunit_en.lower() == 'aa':
            ph_usr * (ph_usr*u.ph/u.s/u.cm**2/u.keV).to(u.ph/u.s/u.cm**2/u.AA, 
                equivalencies=u.spectral_density(0.5*(self.ear[1:]+self.ear[:-1])*u.keV)) #ph/s/cm^2/AA
        else:
            raise ValueError
            
        
        if self.phunit_flx.lower() == 'photon':
            pass
        elif self.phunit_flx.lower() == 'cgs':
            keverg = (1.0*u.erg).to(u.keV).value #keV/erg
            ph_usr = ph_usr * 0.5*(self.ear[1:]+self.ear[:-1]) #keV/s/cm^2/phunit_en
            ph_usr = ph_usr / keverg #erg/s/cm^2/phunit_en
        elif self.phunit_flx.lower == 'si':
            kevJ = (1.0 * u.J/u.m**2).to(u.keV/u.cm**2).value #(keV/cm^2)/(J/m^2)
            ph_usr = ph_usr * 0.5*(self.ear[1:]+self.ear[:-1]) #keV/s/cm^2/phunit_en
            ph_usr = ph_usr / kevJ #J/s/m^2/phunit_en
        else:
            raise ValueError
        
        
        if self.isflipped:
            ph_usr = np.flip(ph_usr)
            
        return ph_usr
    
    
    def _convert_inputSpec(self, ear_in: ArrayLike, ph_in: ArrayLike,
                           eunit: str, phunit: str) -> ArrayLike:
        """
        Converts input energy grid and spectrum to internal units of keV for
        ear, and photons/s/cm^2/bin for ph

        Parameters
        ----------
        ear_in : ArrayLike
            Input energy bins
        ph_in : ArrayLike
            Input spectrum
        eunit : str
            Input energy unit
        phunit : str
            Input flux-density unit

        """
        
        #first splitting ph_unit parts
        phunit_flx, phunit_en = phunit.split('_')
        
        #checkin valid units
        self._check_units(eunit, phunit_flx, phunit_en)
        #if valid strogin unit strings for later use
        self.eunit = eunit
        self.phunit_flx, self.phunit_en = phunit_flx, phunit_en
        
        #converting energy bins
        self.ear = self._convert_earInUnit(ear_in, eunit)
        
        #making sure arrays are organised in rising energy
        self.isflipped = False
        if self.ear[1] - self.ear[0] < 0:
            self.ear = np.flip(self.ear)
            ph_in = np.flip(ph_in)
            self.isflipped = True
        
        #converting spectrum units
        self.ph = self._convert_PhInUnit(ph_in, phunit_flx, phunit_en)
        


     
    def _convert_PhInUnit(self, ph_in: ArrayLike, phunit_flx: str,
                                phunit_en: str) -> ArrayLike:
        """
        Converts spectral photon flux to photon/s/cm^2/bin

        Parameters
        ----------
        ph_in : ArrayLike
            Input spectrum
        phunit_flx : str
            flux unit
        phunit_en : str
            energy density unit

        Returns
        -------
        ph_out : ArrayLike
            Input spectrum in units of photons/s/cm^2/bin

        """
        
        if phunit_flx.lower() == 'photon':
            ph_out = ph_in #photon/s/cm^2/phunit_en
        elif phunit_flx.lower() == 'cgs':
            keverg = (1.0*u.erg).to(u.keV).value #keV/erg
            ph_out = ph_in * keverg #keV/cm^2/s/phunit_en
            ph_out = ph_out/(0.5*(self.ear[1:]+self.ear[:-1])) #photon/s/cm^2/phunit_en
        elif phunit_flx.lower() == 'si':
            kevJ = (1.0 * u.J/u.m**2).to(u.keV/u.cm**2).value #(keV/cm^2)/(J/m^2)
            ph_out = ph_in * kevJ #keV/s/cm^2/phunit_en
            ph_out = ph_out/(0.5*(self.ear[1:]+self.ear[:-1])) #photon/s/cm^2/phunit_en
        else:
            raise ValueError(f'{phunit_flx} invalid')
            
        
        if phunit_en.lower() == 'kev':
            pass
        elif phunit_en.lower() == 'hz':
            ph_out = (ph_out*u.ph/u.s/u.cm**2/u.Hz).to(
                u.ph/u.s/u.cm**2/u.keV, equivalencies=u.spectral()).value #ph/s/cm^2/keV
        elif phunit_en.lower() == 'aa':
            ph_out = (ph_out*u.ph/u.s/u.cm**2/u.AA).to(u.ph/u.s/u.cm**2/u.AA, 
                    equivalencies=u.spectral_density(0.5*(self.ear[1:]+self.ear[:-1])*u.keV)).value 
        else:
            raise ValueError(f'{phunit_en} not valid')
        
        #conversion factor from input unit to photons/s/cm^2/bin
        ph_out = ph_out * (self.ear[1:] - self.ear[:-1]) #(photon/s/cm^2/bin) / phunit
        
        return ph_out
            
    
    def _convert_earInUnit(self, ear_in: ArrayLike, eunit: str) -> ArrayLike:
        """
        Converts unit of input energy bin to keV

        Parameters
        ----------
        ear_in : ArrayLike
            Input energy bins

        Returns
        -------
        ear_out : ArrayLike
            Energy bins in keV
        """
        
        if eunit.lower() == 'kev':
            ear_out = ear_in
        elif eunit.lower() == 'hz':
            ear_out = (ear_in*u.Hz).to(u.keV, equivalencies=u.spectral()).value #keV/Hz
        elif eunit.lower() == 'aa':
            ear_out = (ear_in*u.AA).to(u.keV, equivalencies=u.spectral()).value #keV/AA
        else:
            raise ValueError(f'{eunit} invalid')

        return ear_out
    

    def _check_units(self, eunit: str, phunit_flx: str, phunit_en:str):
        """
        Chacks validity of input unit strings
        Raises error if invalid input

        Parameters
        ----------
        eunit : str
            DESCRIPTION.
        phunit_flx : str
            DESCRIPTION.
        phunit_en : str
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if eunit.lower() in ['kev', 'hz', 'aa']:
            pass
        else:
            print(f'{eunit} not valid for ear')
            print('Must be in [keV, Hz, AA]')
            raise ValueError
            
        if phunit_en.lower() in ['kev', 'hz', 'aa']:
            pass
        else:
            print(f'{phunit_en} not valid for second part of phunit')
            print('Must be in [keV, Hz, AA]')
            raise ValueError

        if phunit_flx.lower() in ['photon', 'cgs', 'si']:
            pass
        else:
            print(f'{phunit_flx} not valid for first part of phunit')
            print('Must be in [photon, cgs, SI]')
            raise ValueError


if __name__ == '__main__':
    
    phc = xw_conv(r_in=40, log_vinf=-1)
    def gau(eas, s0, x0):
        g = np.exp(-0.5 * ((eas - x0)/s0)**2)
        return g
    
    ear_in = np.geomspace(0.1, 10, 1000)
    #print(ear_in[1] - ear_in[0])
    #print(ear_in)
    ph_in1 = gau(0.5*(ear_in[1:]+ear_in[:-1]), 1e-2, 5)
    ph_in2 = gau(0.5*(ear_in[1:]+ear_in[:-1]), 1e-2, 1)
    #print(max(ph_in))
    
    #plt.plot(0.5*(ear_in[1:] + ear_in[:-1]), ph_in)
    #plt.show()
    
    phc.set_input_spectrum(ear_in, ph_in1+ph_in2, eunit='kev', phunit='photon_kev')
    phc.do_convolve(45)
    phc.plot_spec(show=True, as_log=False)
    
    
