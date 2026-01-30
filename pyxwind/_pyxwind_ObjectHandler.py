#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 15:13:00 2025

@author: Scott Hagen

Python wrapper for xwind. The main line-profile calculations are contained in
the utils directory. As these are written in fortran (for performance reasons)
they first need to be compiled with numpy f2py bindings in order to create a
python callable library. Scripts are provided to do this (assuming the user has
an adequate fortran compilor installed, e.g, gfortran)

This module also contains functionality for handling plot grids, s.t output
wind structures/parameters can be visualised.
"""

import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt

from numpy.typing import ArrayLike


class PyXWIND_object:
    """
    Main XWIND object, which all subsequent models inherit from.
    
    Contains methods for calculating wind properties:
        - Density profile
            - Corresponding line-of-sight column-density
        - Velocity profile
        - Absorption and Transmition
        - Fluoresence/Emissivity profile
        
    Also contains helper methods for plotting wind properties
    """
    
    
    dcos_theta = 0.002
    dlog_r = 0.01
    dphi = 0.001
    
    Emin_in = 6.9
    Emax_in = 40
    Nbins_E = 100
    
    def __init__(self,
                 Mbh: float = 1e8,
                 log_mdot_w: float = -3.0,
                 r_in: float = 500,
                 r_out: float = 1000,
                 d_foci: float = 500,
                 fcov: float = 0.5,
                 log_vinf: float = -3.0,
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
        log_mdot_w : float
            Log10 Eddington scaled mass-outflow rate
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

        #Read pars
        self.Mbh = float(Mbh)
        self.mdot_w, self.log_mdot_w = float(10**log_mdot_w), log_mdot_w
        self.r_in = float(r_in)
        self.r_out = float(r_out)
        self.d_foci = float(d_foci)
        self.fcov = float(fcov)
        self.vinf, self.log_vinf = float(10**log_vinf), log_vinf
        self.rv = float(rv)        
        self.beta = float(beta)
        self.kappa = float(kappa)
        self.Afe = float(Afe)
        
        
        #geometric checks
        self._check_fcov()
        
        
        #physical scales
        self._set_constants()
        self._set_MdotEdd()
        self.Rg = (self.G*self.Mbh)/self.c**2.0
        

        return
    
    def _init_geometry(self):
        """
        Initiates geometry arrays
        Only required if looking at wind properties (e.g density/velocity profile, etc).
        For users only interested in the line shape this is ignored, as 
        it is handled internally in the fortran codes
        
        Calculated automatically if required

        """
        
        #generating bins in cos theta and log r, defining geometry
        self._calc_rmax()
        self.cos_th_bins = np.arange(0, self.fcov+self.dcos_theta, self.dcos_theta)
        self.log_rbins = np.arange(0, np.log10(self.rmax)+self.dlog_r, self.dlog_r)

        self.log_rmids, self.cos_th_mids = np.meshgrid(self.log_rbins[:-1]+0.5*self.dlog_r,
                                                       self.cos_th_bins[:-1]+0.5*self.dcos_theta,
                                                       indexing='ij')
        self.rmids = 10**(self.log_rmids)
        
        self._gene_grid_mask()
        self._calc_rlaunch()
        self._calc_lws()
    
    
    def _init_incidentEgrd(self):
        """
        Initiates energy grid and absorption coefficient grid for the incident
        spectrum

        Calculates automatically if required
        """
    
        self.Egrd_in = np.geomspace(self.Emin_in, self.Emax_in, self.Nbins_E)
        self._gene_crossSection()
        self._calc_absCoeff()
    
    def _set_constants(self):
        """
        Sets physical constants in sensible units (cgs)
        Stores as object attributes

        Returns
        -------
        None.

        """
        
        self.G = (const.G * const.M_sun).to(u.cm**3 / u.s**2).value #cm^3/Msol/s^2
        self.sigma_sb = const.sigma_sb.to(u.erg/u.s/u.K**4/u.cm**2).value 
        self.c = const.c.to(u.cm/u.s).value #speed of light, cm/s
        self.h = const.h.to(u.erg*u.s).value #Plank constant, erg s
        self.k_B = const.k_B.to(u.erg/u.K).value #Boltzmann constant, #erg/K
        self.m_e = const.m_e.to(u.g).value  #electron mass, g
        self.m_p = const.m_p.to(u.g).value #proton mass, g


    def _set_MdotEdd(self):
        """
        Sets Eddington mass-accretion rate (in g/s)
        Used to scale mass outflow rate, which should be given in units of
        Eddington
        Assumes Ledd = eta Mdot_edd c^2
        
        Since BH spin is not part of the model, will assume a radiative 
        efficiency for a non-spinning black hole (i.e eta = 0.057)

        Returns
        -------
        None.

        """
        
        Ledd = 1.39e38 * self.Mbh #ergs/s
        self.Mdot_edd = Ledd/(0.057 * self.c**2) #g/s
    
    ###########################################################################
    #### GEOMETRIC CONSIDERATIONS
    ###########################################################################
    
    def _check_fcov(self):
        """
        Checks if covering fraction is reached, given the input parameters!

        If fcov too big, throws a warning, and then updates to new value        

        Returns
        -------
        None.

        """
    
        if np.arccos(self.fcov) <= np.arctan(self.r_in/self.d_foci):
            print('WARNING! Covering fraction fcov is never reached for input wind!')
            print(self.fcov, self.r_in, self.d_foci)
            th_min = np.arctan(self.r_in/self.d_foci) + 0.1
            self.fcov = np.cos(th_min)
            print(f'Updating covering fraction to: {self.fcov}')
    
        else:
            pass
    
    
    def _calc_los_intercept(self, theta: float, r_launch: float) -> float:
        """
        Calculates the intercept between a streamline launched at r_launch
        and a line-of-sight given by theta

        Parameters
        ----------
        theta : float
            Line-of-sight angle 
            Units : rad
        r_launch : float
            Streamline launch radius
            Units : Rg

        Returns
        -------
        r_int : float
            Radius (in spherical polars) of intercept
            Units : Rg
        """
        
        r_int = r_launch/(np.sin(theta) - (r_launch/self.d_foci) * np.cos(theta))
        return r_int
        
    
    def _calc_rmax(self):
        """
        Calculates the boundary radius 

        Returns
        -------
        None.

        """
        
        th_min = np.arccos(self.fcov)
        self.rmax = self._calc_los_intercept(th_min, self.r_in)
    
    def _gene_grid_mask(self):
        """
        Generates the mask to be applied to the wind grid
        
        Since unifrom grid in theta and r - mask out values where grid midpoint
        outside wind boundaries

        Returns
        -------
        None.

        """
        
        #start with all masked
        self.grid_mask = np.full_like(self.rmids, True, dtype=bool)
        #iterating over theta
        for i in range(self.grid_mask.shape[1]):
            th_mid = np.arccos(self.cos_th_mids[0, i]) #single l.o.s, so th const
            r_inner = self._calc_los_intercept(th_mid, self.r_in)
            
            #finding outer radius (either set by rmax or outer wind boundary)
            if th_mid > np.arctan(self.r_out/self.d_foci):
                rw_bound = self._calc_los_intercept(th_mid, self.r_out)
                if rw_bound < self.rmax:
                    rgrd_max = rw_bound
                else:
                    rgrd_max = self.rmax
            else:
                rgrd_max = self.rmax
            
            idx_in = np.argwhere(np.logical_and(self.rmids[:, i]>=r_inner,
                                                self.rmids[:, i]<=rgrd_max))[:, 0]
            self.grid_mask[idx_in, i] = False
        
        self.rmids = np.ma.array(self.rmids, mask=self.grid_mask)
        self.log_rmids = np.ma.array(self.log_rmids, mask=self.grid_mask)
        self.cos_th_mids = np.ma.array(self.cos_th_mids, mask=self.grid_mask)
    
    
    def _calc_rlaunch(self):
        """
        Calculates launch radii for streamllines corresponding to grid-points

        Returns
        -------
        r_ls

        """

        self.cos_beta = self.d_foci + self.rmids*self.cos_th_mids
        self.cos_beta /= np.sqrt(self.rmids**2 + 2*self.rmids*self.d_foci*self.cos_th_mids + self.d_foci**2)

        self.r_ls = self.d_foci * np.ma.tan(np.ma.arccos(self.cos_beta))

    
    def _calc_lws(self):
        """
        Calculates distance along streamlines for each grid-point

        Returns
        -------
        None.

        """

        th_mids = np.ma.arccos(self.cos_th_mids)
        self.lws = np.sqrt(self.rmids**2 + self.r_ls**2 - 2*self.rmids*self.r_ls*np.sin(th_mids))



    ###########################################################################
    #### VELOCITY PROFILE
    ###########################################################################
    def calc_velocity_profile(self):
        """
        Calculates velocity profile across wind grid
        Stores as class attribute

        Returns
        -------
        None.

        """
        
        if hasattr(self, 'log_rmids'):
            pass
        else:
            self._init_geometry()
            
        self.vl_grd = self.calc_vl(self.lws)
        self.vphi_grd = self.calc_vphi(self.rmids, self.r_ls)
        
    
    def calc_vl(self, lw: float|ArrayLike) -> float|ArrayLike:
        """
        Calculates velocity along the wind streamline at a given point along 
        the streamline
        
        (Note - takes lw parameter to allow for seperate plotting of velocity
         along a streamline. For line-profile calculations handled internally)
        
        Parameters
        ----------
        lw : float or array
            Length along streamline (as measured from base of the wind)
            Units : Rg
            
        Returns
        -------
        v_l : float or array
            Wind velocity along the streamline
            Units : c
        """

        v_l = self.vinf * (1 - (self.rv/(self.rv + lw)))**self.beta
        return v_l
    
    
    def calc_vphi(self, 
                  r: float|ArrayLike,
                  rl: float|ArrayLike) -> float|ArrayLike:
        """
        Calculates azimuthal velocity for a given point in the streamline
        Assumes Keplerian at the base of the wind, then conserves angular momentum

        Parameters
        ----------
        r : float or array
            Distance from black hole (in spherical polars)
            Units : Rg
        rl : float or array
            Launch radius of streamline
            Units : Rg
            
        Returns
        -------
        v_phi : float or array
            Tangental wind velocity
            Units: c
        """
    
        v_phi = (rl/r) * np.sqrt(1/rl)
        return v_phi    
    
    
    ###########################################################################
    #### DENSITY AND COLUMN-DENSITIES
    ###########################################################################
    
    def _calc_ndens(self):
        """
        Calculates the Hydrogen number density within each grid cell
        Assumes conservation of mass and a constant mass outflow rate

        Returns
        -------
        None.

        """
        if hasattr(self, 'vl_grd'):
            pass
        else:
            self.calc_velocity_profile()

        #first getting mass-outflow per unit area at base (for each grid-point)
        dmdot_dA = self.mdot_w*self.Mdot_edd*(self.kappa + 2) * self.r_ls**(self.kappa)
        dmdot_dA /= (4*np.pi * (self.r_out**(self.kappa+2) - self.r_in**(self.kappa+2)))
        dmdot_dA *= self.Rg**(-2)
        
        #now getting number dens
        self.ndens = (1/(1.23*self.m_p*self.vl_grd*self.c)) * dmdot_dA
    
    
    def calc_los_NH(self, inc: float) -> float:
        """
        Calculates the column-density along a line-of sight

        Parameters
        ----------
        inc : float
            Line of sight inclination
            Units : deg

        Returns
        -------
        Nh : float
            Line of sight column-density through the wind
            Units : cm^-2

        """
        
        if hasattr(self, 'ndens'):
            pass
        else:
            self._calc_ndens()
        
        
        cos_inc = np.cos(np.deg2rad(inc))
        #index of low bin edge
        idx_low = int(len(self.cos_th_bins[self.cos_th_bins <= cos_inc]) - 1)
        drs = 10**(self.log_rbins[1:]) - 10**(self.log_rbins[:-1]) #radial bin widths in Rg
        if self.cos_th_bins[idx_low] == self.cos_th_bins[-1]:
            return 0
        else:
            return np.ma.sum(drs*self.Rg * self.ndens[:, idx_low])
    
    
    def calc_NH_profile(self) -> ArrayLike:
        """
        Calculates NH profile of wind as a function of theta.
        Uses internal cos_theta grid
        

        Returns
        -------
        NHarr : ArrayLike
            Column-density through the wind as a function of theta
            Units : cm^-2
        thetas : ArrayLike
            Also gives the relevant theta array (evaluated in centre of cos_th bins)
            This is to avoid the user having to deal with internal array structure
            Units : deg

        """
        
        if hasattr(self, 'ndens'):
            pass
        else:
            self._calc_ndens()
    
        #generating output theta arr
        thetas = np.acos(self.cos_th_bins[:-1] + 0.5*self.dcos_theta)
        thetas = np.rad2deg(thetas)
        
        #generating outpt Nh arr
        NHarr = np.empty_like(thetas)
        drs = 10**(self.log_rbins[1:]) - 10**(self.log_rbins[:-1]) #radial bin widths in Rg
        for i in range(len(thetas)):
            NHarr[i] = np.ma.sum(drs*self.Rg * self.ndens[:, i])
        
        return NHarr, thetas
        
    
    ###########################################################################
    #### ABSORPTION AND EMISSIVITY
    ###########################################################################
    def define_inputSpec(self, Gamma: float, N0: float):
        """
        Defines input spectrum
        Currently just a power-law to be consistent with lineprofile calculations

        Parameters
        ----------
        Gamma : float
            Photon index
        N0 : float
            Normalisation at 1keV
            Units : photons s^-1 cm^-2 keV^-1

        Returns
        -------
        None.

        """
        
        if hasattr(self, 'Egrd_in'):
            pass
        else:
            self._init_incidentEgrd()
        
        self.inspec = N0 * self.Egrd_in**(-Gamma)
        
    
    
    def _gene_crossSection(self):
        """
        Initiates array for Fe-Kalpha cross section
        
        Defined as 0 if E < Ek
                sigma_0 * exp(E/Ek)^(-alpha) if E > Ek

        Returns
        -------
        None.

        """
        
        Ek = 7.1 #Edge energy - fixed for now
        sigma0 = 3.37e-20
        alpha = 2.67
        self.cross_sec = np.zeros_like(self.Egrd_in)
        self.cross_sec[self.Egrd_in >= Ek] = sigma0 * (
                                        self.Egrd_in[self.Egrd_in>Ek]/Ek)**(-alpha)

    
    
    def _calc_absCoeff(self):
        """
        Calculates the absorption coefficient needed to calculated
        the spectrum transmitted through each grid point, such that for a given
        input spectrum Nx(E) 
        The spectrum transmitted through grid-point i is given by:
            Ni(E) = Nx(E) * absCoeff
        
        This is multiplicative along a line of sight, such that for grid-point
        k this is (in math notation):
            absCoeff = Pi_{i}^k exp(-sigma(E) * Afe * ni * delta r)
            
        where Pi^{i} indicates a multiplication over all elements from 0 to k,
        Afe is the iron abundance, ni is the number density for grid-point i
        

        Returns
        -------
        None.

        """

        if hasattr(self, 'ndens'):
            pass
        else:
            self._calc_ndens()


        Afe = 4.68e-5 #abundance from Anderson \& Grevesse 1989 
    
        ln_absCoeff = Afe * self.rmids * self.dlog_r * self.Rg * np.cumsum(self.ndens, axis=0)
        ln_absCoeff = ln_absCoeff[:, :, np.newaxis] * self.cross_sec
        self.abs_coeff = np.exp(-ln_absCoeff)

        
    
    def calc_transmittedSpec(self):
        """
        For a given input spectrum - calculates the transmitted spectrum
        in each grid-point
        
        If no input spectrum is defined gives a warning, and defulats to a 
        flat power-law

        Returns
        -------
        None.

        """
        
        if hasattr(self, 'inspec'):
            pass
        else:
            print('No input spectrum given!!!')
            print('Defualting to a flat power-law (Gamma=2)')
            print('Call .define_inputSpec() first to use own input spectrum')
            self.define_inputSpec(2, 10)
            
        self.transmitted_spec = self.inspec * self.abs_coeff
    
    
    def calc_fluoresence(self):
        """
        Calculates the fluorescent yield within each wind cell

        Returns
        -------
        None.

        """
        
        
        if hasattr(self, 'transmitted_spec'):
            pass
        else:
            self.calc_transmittedSpec()
        
        self.fluoresence = np.empty_like(self.rmids)
        #for each polar angle
        for i in range(self.transmitted_spec.shape[1]):
            idx_non_mask = np.argwhere(self.grid_mask[:, i] == False)[:, 0]     
            t_theta_spec = self.transmitted_spec[idx_non_mask, i, :]
            
            inspec = np.row_stack((self.inspec, t_theta_spec[:-1]))
            diffspec = inspec - t_theta_spec
            
            self.fluoresence[idx_non_mask, i] = np.trapz(diffspec, self.Egrd_in, axis=1)
        
        self.fluoresence *= self.dcos_theta*self.dphi*0.3 #0.3 is Fe-Kalpha yield
        #flatteing array to only valid inputs - for performace boost
        #self.fluoresence_valid = self.fluoresence[self.grid_mask==False] 


    def calc_volumeEmissivity(self):
        """
        Calculates volumetric emissivity of each wind cell
        
        Simply divides thr fluoresence by the volume-element of each cell
        (i.e dV = R^2 sin(theta) dR dtheta dphi)

        Returns
        -------
        None
        """

        if hasattr(self, 'fluoresence'):    
            pass
        else:
            self.calc_fluoresence()
        
        drs = 10**(self.log_rbins[1:]) - 10**(self.log_rbins[:-1])
        dV = self.rmids**2 * drs[:, np.newaxis] * self.Rg**3 * self.dcos_theta * self.dphi

        self.vol_emiss = self.fluoresence/dV


    ###########################################################################
    #---- PLOTTING METHODS
    ###########################################################################
    
    def plot_wind(self,
                  profile_type: str,
                  ax: None|type(plt.axis) = None,
                  cmap: str = 'jet',
                  cnorm: str = 'log',
                  cscale = None,
                  inc_cbar: bool = False,
                  cbar_label: str = 'none',
                  cbar_loc: str = 'right',
                  show_axlabel: bool = True,
                  show: bool = False):
        """
        Plots a 2D map of the wind geometry (in x, z).
        The colour gradient is set by the profile_type parameter
        
        Parameters
        ----------
        profile_type : {'ndens', 'vl', 'vphi', 'fluoresence', 'rel_emiss', 'vol_emiss'}
            What wind parameter to plot in profile
            Shown as colour gradient
        ax : None|plt.Axis
            Pyplot axis to plot on
            If None, generates own figure + axis
            DEFAULT : None
        cmap : str, optional
            Pyplot colourmap to use. 
            DEFAULT: 'jet'.
        cnorm : {'log', 'linear'}
            Colourmap normalisation to use
            Either a string, or a user defined normalisation scale
            DEFAULT: 'log'
        cscale : pyplot normalisation scale
            Optional scale to be passed if the user want more control over 
            colour normalisation
        inc_cbar : bool
            if True, then adds a colourbar to the axis
            DEFAULT: False
        cbar_label : str
            Label for colourbar
            If none then automatically defined via profile type
        cbar_loc : {'bottom', 'left', 'right', 'top'}
            Colourbar location
            DEFAULT : right
        show_axlabel : bool
            Whether to write axis labels
        show : bool
            Whether to show plot. 
            DEFAULT: False.

        Returns
        -------
        None.

        """
        
        if ax is None:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
        
        rel_attr = False
        if profile_type == 'vl' or profile_type == 'vphi':
            attr = f'{profile_type}_grd'
            cpar = f'{profile_type}/c'
            cunit = ''
        elif profile_type == 'ndens':
            attr = profile_type
            cpar = 'n'
            cunit = r'cm$^{-3}$'
        elif profile_type == 'fluoresence':
            attr = profile_type
            cpar = ''
            cunit = r'photons s$^{-1}$ cm$^{-2}$'
        elif profile_type == 'rel_emiss':
            attr = 'fluoresence'
            rel_attr = True
            cpar = 'Emissivity'
            cunit = ''
        elif profile_type == 'vol_emiss':
            attr = 'vol_emiss'
            cpar = 'Volumetric Emissivity'
            cunit = ''
        elif profile_type == 'rel_vol_emiss':
            attr = 'vol_emiss'
            rel_attr = True
            cpar = 'Volumetric emissivity'
            cunit = r'$N_{\gamma}/(dV N_{\gamma, max}$'
        else:
            raise NameError('f{profile_type} invalid. profile_type must be:'
                            ' ndens, vlm vphi, or fluoresence')
        
        pc = self._do_pcolormapWind(ax, attr, cmap, cnorm, cscale, rel_attr)
        
        if inc_cbar == True:
            if cbar_label == 'none':
                if cnorm == 'log':
                    cbar_label = r'$\log_{{{}}} {}$   {}'.format(10, cpar, cunit)
                else:
                    cbar_label = f'{cpar}   {cunit}'
            
            plt.colorbar(pc, ax=ax, location=cbar_loc, label=cbar_label)
        
        if show_axlabel:
            ax.set_xlabel(r'x   $R_{G}$')
            ax.set_ylabel(r'z   $R_{G}$')
        
        if show:
            plt.show()
            
        
        return ax
    
    
    def _do_pcolormapWind(self,
                          ax: plt.axis,
                          attr: str,
                          cmap: str,
                          cnorm: str,
                          cscale,
                          rel_attr: bool) -> plt.pcolor:
        """
        Does the colourmapping of the wind attribute

        Parameters
        ----------
        ax : plt.Axis
            DESCRIPTION.
        attr : str
            DESCRIPTION.
        cmap : str
            DESCRIPTION.
        cnorm : str
            DESCRIPTION.

        Returns
        -------
        pc : plt.pcolor
            colourmesh object

        """
        
        if hasattr(self, 'x_plt_mesh'):
            pass
        else:
            self._define_plotGrid()
        
        if hasattr(self, attr):
            pass
        else:
            if attr == 'ndens':
                self._calc_ndens()
            elif attr == 'vl_grd' or attr == 'vphi_grd':
                self.calc_velocity_profile()
            elif attr == 'fluoresence':
                self.calc_fluoresence()
            elif attr == 'vol_emiss':
                self.calc_volumeEmissivity()
            else:
                raise NameError(f'{attr} not valid!')
            
        cm_attr = getattr(self, attr)
        
        if rel_attr:
            cm_attr = cm_attr/np.amax(cm_attr)
        
        if cscale is None:
            pc = ax.pcolor(self.x_plt_mesh, self.z_plt_mesh, cm_attr,
                           vmin=np.amin(cm_attr), vmax=np.amax(cm_attr), 
                           norm=cnorm, cmap=cmap)
        else:
            if cnorm == 'log':
                pc = ax.pcolor(self.x_plt_mesh, self.z_plt_mesh, np.log10(cm_attr), 
                               norm=cscale, cmap=cmap)
            else:
                pc = ax.pcolor(self.x_plt_mesh, self.z_plt_mesh, cm_attr,
                               norm=cscale, cmap=cmap)
            
        return pc
    
    def _define_plotGrid(self):
        """
        Defines the 2D grid for plotting - usies cartesian coords

        Returns
        -------
        None.

        """
        
        if hasattr(self, 'log_rbins'):
            pass
        else:
            self._init_geometry()
        
        rbin_mesh, th_bin_mesh = np.meshgrid(10**(self.log_rbins), np.arccos(self.cos_th_bins),
                                             indexing='ij')
        
        #extedning mask to work with grid
        plt_mask = np.column_stack((self.grid_mask, self.grid_mask[:, -1]))
        plt_mask = np.row_stack((plt_mask, np.append(self.grid_mask[-1, :], self.grid_mask[-1, -1])))
        
        self.x_plt_mesh = np.ma.array(rbin_mesh*np.sin(th_bin_mesh),
                                      mask=plt_mask)
        self.z_plt_mesh = np.ma.array(rbin_mesh*np.cos(th_bin_mesh),
                                      mask=plt_mask)





if __name__ == '__main__':
    xw = PyXWIND_object(r_out=1000, mdot_w=0.01)
    #xw.calc_fluoresence()
    
    xw.plot_wind('ndens', inc_cbar=True, cbar_loc='right', show=True)
