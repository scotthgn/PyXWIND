# PyXWIND

This is a python wrapper for the [XWIND](https://github.com/scotthgn/XWIND/tree/main?tab=readme-ov-file) model. This is an analytic model for calculating the line profiles of outflowing winds in AGN. Specifically, this is designed for slow(ish) BLR scale winds.  It is originally designed with the Fe-K $\alpha$ complex in mind, but also allows for the convolution of the same transfer functions for any line species (under the strong assumption that their emissivity traces density in the same way as netral Fe-K $\alpha$). In addition, this python version includes functionality for extracting the information on the wind structure (i.e density, emissivity, velcity profiles, etc), which allows for the user to plot the wind structures corresponding to their model fits. If this is used in your research, please cite: **Hagen et al. (2026, submitted to A&A)**.

Internally, the routines used to calculate the line profiles are written in fortran, due to significant perfomance benefits. This is done using the `numpy` `f2py` bindings. Beyond compiling the fortran code (see Installation and the included `compile.sh` script), the user does not need to worry about this, as everythin on the user end is handled by the `PyXWIND_object` class, with the relevant subclasses being `xw_line` and `xw_conv` (see below for more details).


## Requirements

1. Python 3.10 or above (tested with Python v.3.12.4 and v.3.13.5)
2. `numpy` (tested on version 1.26.4 and 2.3.3)
3. `matplotlib` (tested on version 3.9.1 and 3.10.6)
4. `astropy` (tested on version 7.1.0)
5. fortran compilor/python binder: </br>
&emsp; &emsp; &#9656; If **python 11 or lower**: Any fortran complior should do </br>
&emsp; &emsp; &#9656; If **python 12 and above**: `meson` is required. Can be found at: https://mesonbuild.com/ </br>
&emsp; &emsp; &emsp; &emsp; This requirement comes from changes in `numpy` `f2py` bindings with python12 </br>



## Installation

Currently:
1. Clone the repository
2. From within the main directory run the compile script compile.sh

```
./compile.sh
```

## Usage

Once the compile script is run, and `PyXWIND` is added to the `pythonpath`, it can be imported as:

```
import pyxwind
````

If you do not want to edit the `pythonpath` manually, you can also do:
```
import sys
sys.path.append('/path/to/PyXWIND')
import pyxwind
```

Within `PyXWIND` there are two main classes for the user to use: `xw_line` and `xw_conv`. The main difference between these is how the line profiles are calculated. In `xw_line` they are treated as additive components, whereas in `xw_conv` you need to supply an input spectrum which is then convolved with the `xwind` line shape. These are imported as:

```
from pyxwind import xw_line, xw_conv
```

For examples on how to use these see the scripts in the `examples` directory.

## Class parameters
This is an overview of the parameters needed when initiating either `xw_line` or `xw_conv`. These are the same for both classes.

&#9656; M_BH: float, optional </br>
&emsp; &emsp; Central black hole mass </br>
&emsp; &emsp; NOTE: **only plays a role if you need the physical density profile** </br>
&emsp; &emsp; For line profiles this **does not matter**, since those work in dimensionless units </br>
&emsp; &emsp; Units: $M_{\odot}$ </br>
&emsp; &emsp; Default: 1e8 </br>

</br>

&#9656; log_mdot_w: float, optional </br>
&emsp; &emsp; Wind mass outflow rate, scaled by the Eddington mass accretion rate </br>
&emsp; &emsp; Units: $\log{10} (\dot{M}\_{w} / \dot{M}\_{\rm{Edd}}) $ </br>
&emsp; &emsp; Default: -3 </br>

</br>

&#9656; r_in: float, optional </br>
&emsp; &emsp; Inner wind launch radius </br>
&emsp; &emsp; Units: $R\_{G}$ </br>
&emsp; &emsp; Defualt: 1000 </br>

</br>

&#9656; r_out: float </br>
&emsp; &emsp; Outer launch radius </br>
&emsp; &emsp; Units: $R\_{G}$ </br>
&emsp; &emsp; Default: 1500 </br>

</br>

&#9656; d_foci </br>
&emsp; &emsp; Distance below the origi of the wind focus </br>
&emsp; &emsp; Units: $R\_{G}$ </br>
&emsp; &emsp; Default: 1700 </br>

</br>

&#9656; fcov: float, optional </br>
&emsp; &emsp; Wind covering fraction, as seen from the central source </br>
&emsp; &emsp; Units: $\frac{\Omega}{4 \pi}$ </br>
&emsp; &emsp; Default: 0.5 </br>

</br>

&#9656; log_vinf: float, optional </br>
&emsp; &emsp; Outflow velocity at infinity </br>
&emsp; &emsp; Units: $v/c$ </br>
&emsp; &emsp; Default: -3 </br>

</br>

&#9656; rv: float, optional </br>
&emsp; &emsp; Velocity scale length (i.e. distance along streamline where $v\_{l} = 0.5 v\_{\infty}$ </br>
&emsp; &emsp; Units: $R\_{G}$ </br>
&emsp; &emsp; Default: 1000 </br>

</br>

&#9656; beta: float, optional </br>
&emsp; &emsp; Velocity law exponent. Determines acceleration along streamline </br>
&emsp; &emsp; Units: Dimensionless </br>
&emsp; &emsp; Default: 1 </br>

</br>

&#9656; kappa: float, optional </br>
&emsp; &emsp; Wind density gradient exponent. Determines launch efficiency with radius </br>
&emsp; &emsp; Units: Dimensionless </br>
&emsp; &emsp; Default: -1 </br>

</br>

&#9656; Afe: float, optional </br>
&emsp; &emsp; Iron abundance, relative to solar </br>
&emsp; &emsp; Units: [Fe/Fe_sol] </br>
&emsp; &emsp; Default: 1 </br>

## Key class methods

Below there is an overview of the key class methods. First giving the ones unique to `xw_line` and `xw_conv`, before finishing with the ones that exist in both classes.

### `xw_line` specific methods:

&#9656; `calc_xwindline(inc, E0, Ebins, vturb=0, edens_unit=True)` </br>
&emsp; &emsp; Calculates a sinlge line profile, assuming a delta function in the rest frame. Can be used for any rest frame energy </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ----------- </br>
&emsp; &emsp; &#9656; inc: float </br>
&emsp; &emsp; &emsp; Observer inclination </br>
&emsp; &emsp; &emsp; Units: degrees </br>
&emsp; &emsp; &#9656; E0: float </br>
&emsp; &emsp; &emsp; Rest frame energy </br>
&emsp; &emsp; &emsp; Units: keV </br>
&emsp; &emsp; &#9656; Ebins: ArrayLike </br>
&emsp; &emsp; &emsp; Energy bin edges (output will lave size len(Ebins)-1) </br>
&emsp; &emsp; &emsp; Units: keV </br>
&emsp; &emsp; &#9656; vturb: float, optional </br>
&emsp; &emsp; &emsp; Turbulent velocity (assumed constant throughout wind) </br>
&emsp; &emsp; &emsp; Units: km/s </br>
&emsp; &emsp; &emsp; Default: 0 </br>
&emsp; &emsp; &#9656; edens_unit: bool, optional </br>
&emsp; &emsp; &emsp; If True convets output to phhotons/s/cm^2/keV </br>
&emsp; &emsp; &emsp; If False output is photons/s/cm^2/bin </br>
&emsp; &emsp; &emsp; Default: True </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; -------- </br>
&emsp; &emsp; &#9656; ph: ArrayLike </br>
&emsp; &emsp; &emsp; Output photon flux array </br>
&emsp; &emsp; &emsp; Shape: len(Ebins)-1 </br>
&emsp; &emsp; &emsp; Units: set by edens_unit par

</br>

&#9656; `calc_xwindfe(inc, N0, gamma, Ebins, vturb=0, edens_unit=True, give_EW=False)` </br>
&emsp; &emsp; xwindline profile specific to Fe-K $\alpha$. </br>
&emsp; &emsp; Normalisation calculated self-consistently from density profile and input spectrum </br>
&emsp; &emsp; Uses 7-Lorentzian Holzer et al. (1997) profile for the rest frame emission </br>
&emsp; &emsp; This gives both Fe-K $\alpha_{1}$ and Fe-K $\alpha_{2}$ </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ------------ </br>
&emsp; &emsp; &#9656; inc: float </br>
&emsp; &emsp; &emsp; Observer inclination </br>
&emsp; &emsp; &emsp; Units: degrees </br>
&emsp; &emsp; &#9656; N0: float </br>
&emsp; &emsp; &emsp; Normalisation of illuminating X-ray spectrum (i.e central source) </br>
&emsp; &emsp; &emsp; Units: photons/s/cm^2 at 1 keV </br>
&emsp; &emsp; &#9656; gamma: float </br>
&emsp; &emsp; &emsp; Photon index of illuminating X-ray spectrum </br>
&emsp; &emsp; &emsp; Units: Dimensionless </br>
&emsp; &emsp; &#9656; Ebins: ArrayLike </br>
&emsp; &emsp; &emsp; Energy bin edges </br>
&emsp; &emsp; &emsp; Units: keV </br>
&emsp; &emsp; &#9656; vturb: float, optional </br>
&emsp; &emsp; &emsp; Turbulent velocity </br>
&emsp; &emsp; &emsp; Units: km/s </br>
&emsp; &emsp; &emsp; Default: 0 </br>
&emsp; &emsp; &#9656; edens_units: bool, optional </br>
&emsp; &emsp; &emsp; If True, output has units: photons/s/cm^2/keV </br>
&emsp; &emsp; &emsp; If False, output has units: photons/s/cm^2/bin </br>
&emsp; &emsp; &emsp; Default: True </br>
&emsp; &emsp; &#9656; give_EW: bool, optional </br>
&emsp; &emsp; &emsp; If True, also returnes line equivalent width (in units keV) </br>
&emsp; &emsp; &emsp; If False, only returnes photon spectrum </br>
&emsp; &emsp; &emsp; Default: False </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; --------- </br>
&emsp; &emsp; &#9656; ph: ArrayLike </br>
&emsp; &emsp; &emsp; Output photon flux array </br>
&emsp; &emsp; &emsp; Size: len(Ebins)-1 </br>
&emsp; &emsp; &emsp; Units: set by edens_units par </br>
&emsp; &emsp; &#9656; EW: float </br>
&emsp; &emsp; &emsp; ONLY returned IF give_EW=True </br>
&emsp; &emsp; &emsp; Line equivalent width </br>
&emsp; &emsp; &emsp; Units: keV </br>

</br>

### `xw_conv` specific methods

&#9656; `set_input_spectrum(ear, ph, eunit='keV', phunit='photon_keV', check_unit=True)` </br>
&emsp; &emsp; Defines the input spectrum to be convolved </br>
&emsp; &emsp; Internally the code uses units of photons/s/cm^2/bin </br>
&emsp; &emsp; For a given eunit and phunit, this method automatically handles unit conversions </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ------------ </br>
&emsp; &emsp; &#9656; ear: ArrayLike </br>
&emsp; &emsp; &emsp; Energy bin edges </br>
&emsp; &emsp; &#9656; ph: ArrayLike </br>
&emsp; &emsp; &emsp; Photon flux array </br>
&emsp; &emsp; &emsp; Must be in a form of flux density (see phunit) </br>
&emsp; &emsp; &emsp; Size: len(ear)-1 </br>
&emsp; &emsp; &#9656; eunit: str, optional </br>
&emsp; &emsp; &emsp; Sets the input units of ear </br>
&emsp; &emsp; &emsp; Options: {'keV', 'Hz', 'AA'} </br>
&emsp; &emsp; &emsp; Note, AA=Ångstrom </br>
&emsp; &emsp; &emsp; Default: keV </br>
&emsp; &emsp; &#9656; phunit: str, optional </br>
&emsp; &emsp; &emsp; Sets the input unit of ph (note, must be a form of flux density) </br>
&emsp; &emsp; &emsp; Default: 'photon_keV'
&emsp; &emsp; &emsp; Options: {'photon_keV', 'photon_Hz', 'photon_AA', </br>
&emsp; &emsp; &emsp; &emsp; 'cgs_keV', 'cgs_Hz', 'cgs_AA', </br>
&emsp; &emsp; &emsp; &emsp; 'SI_keV', 'SI_Hz', 'SI_AA'} </br>
&emsp; &emsp; &emsp; Note: photon=photons/s/cm^2 </br>
&emsp; &emsp; &emsp;       cgs=erg/s/cm^2 </br>
&emsp; &emsp; &emsp;       SI=W/m^2 </br>
&emsp; &emsp; &emsp;       Such that photon_keV=photon/s/cm^2/keV </br>
&emsp; &emsp; &#9656; check_unit: bool, optional </br>
&emsp; &emsp; &emsp; If True, then checks the input unit and converts to internal units </br>
&emsp; &emsp; &emsp; NOTE: This should **always** be set to True unless you are specifically </br>
&emsp; &emsp; &emsp; dealing with things in units of photons/s/cm^2/bin </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; --------- </br>
&emsp; &emsp; &#9656; None </br>

</br>

&#9656; `do_convolve(inc, vturb=0, conv_unit=True)` </br>
&emsp; &emsp; Does the model convolution </br>
&emsp; &emsp; Assumes `set_input_spectrum()` has been defined frist </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ------------
&emsp; &emsp; &#9656; inc: float </br>
&emsp; &emsp; &emsp; Observer inclination </br>
&emsp; &emsp; &emsp; Units: degrees </br>
&emsp; &emsp; &#9656; vturb: float, optional </br>
&emsp; &emsp; &emsp; Turbulent velocity </br>
&emsp; &emsp; &emsp; Units: km/s </br>
&emsp; &emsp; &emsp; Default: 0 </br>
&emsp; &emsp; &#9656; conv_unit: bool, optional </br>
&emsp; &emsp; &emsp; If True, converts output back to units given by user in `set_input_spectrum()` </br>
&emsp; &emsp; &emsp; If False, output is given in photons/s/cm^2/bin </br>
&emsp; &emsp; &emsp; Default: True </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; --------- </br>
&emsp; &emsp; &#9656; ph_conv: ArrayLike </br>
&emsp; &emsp; &emsp; Convolved spectrum </br>
</br>

### Universal methods
&#9656; `calc_los_NH(inc)` </br>
&emsp; &emsp; Calculates column-density along a given line-of-sight </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ------------ </br>
&emsp; &emsp; &#9656; inc: float </br>
&emsp; &emsp; &emsp; Line of sight inclination (measured from z-axis) </br>
&emsp; &emsp; &emsp; Units: degrees </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; --------- </br>
&emsp; &emsp; &#9656; Nh: float </br>
&emsp; &emsp; &emsp; Line of sight column-density </br>
&emsp; &emsp; &emsp; Units: cm^{-2} </br>

</br>

&#9656; `calc_NH_profile()` </br>
&emsp; &emsp; Calculates Nh profile of the wind as a function of theta (line-of-sight) </br>
&emsp; &emsp; Uses internal grid in cos_theta </br>
</br>
&emsp; &emsp; Returns </br>
&emsp; &emsp; --------- </br>
&emsp; &emsp; &#9656; NHarr: ArrayLike </br>
&emsp; &emsp; &emsp; Array of line-of-sight column-densities </br>
&emsp; &emsp; &emsp; Units: cm^{-2} </br>
&emsp; &emsp; &#9656; thetas: ArrayLike </br>
&emsp; &emsp; &emsp; Array of lines-of-sight (corresponding to NHarr) </br>
&emsp; &emsp; &emsp; Units: degrees </br>

</br>

&#9656; `plot_wind(profile_type, ax=None, cmap='jet', cnorm='log', cscale=None, inc_cbar=False, cbar_label='none', cbar_loc='right', show_axlabel=True, show=False)` </br>
&emsp; &emsp; Plots a 2D map (in x,z plane) of the wind geometry </br>
&emsp; &emsp; The colour gradient, set by profile_type, corresponds to a physical wind property (i.e density, emissivity, etc) </br>
</br>
&emsp; &emsp; Parameters </br>
&emsp; &emsp; ------------- </br>
&emsp; &emsp; &#9656; profile_type: str </br>
&emsp; &emsp; &emsp; What wind parameter to plot, shown as the colour gradient </br>
&emsp; &emsp; &emsp; Options: {'ndens', 'vl', 'vphi', 'fluoresence', 'rel_emiss', 'vol_emiss'} </br>
&emsp; &emsp; &#9656; ax: None|pyplot axis, optional </br>
&emsp; &emsp; &emsp; Pyplot axis toplot on </br>
&emsp; &emsp; &emsp; If None, then generates own figure+axis </br>
&emsp; &emsp; &emsp; Default: None </br>
&emsp; &emsp; &#9656; cmap: str, optional </br>
&emsp; &emsp; &emsp; Colourmap to use (must correspond to a valid pyplot colourmap string) </br>
&emsp; &emsp; &emsp; Default: 'jet' </br>
&emsp; &emsp; &#9656; cnorm: str, optional </br>
&emsp; &emsp; &emsp; Colourmap normalisation scale </br>
&emsp; &emsp; &emsp; Options: {'log', 'linear'} </br>
&emsp; &emsp; &emsp; Default: 'log' </br>
&emsp; &emsp; &#9656; cscale: pyplot normalisation object, optional </br>
&emsp; &emsp; &emsp; User defined normalisation scale (overwrites cnorm) </br>
&emsp; &emsp; &emsp; Must correspond to a matplotlib.colors normalisation object </br>
&emsp; &emsp; &emsp; Default: None </br>
&emsp; &emsp; &#9656; inc_cbar: bool, optional </br>
&emsp; &emsp; &emsp; If True, includes colourbar </br>
&emsp; &emsp; &emsp; Default: False </br>
&emsp; &emsp; &#9656; cbar_label: str, optional </br>
&emsp; &emsp; &emsp; Label for the colourbar </br>
&emsp; &emsp; &emsp; If 'none' then automatically defined by the profile_type </br>
&emsp; &emsp; &emsp; Default: 'none' </br>
&emsp; &emsp; &#9656; cbar_loc: str, optional </br>
&emsp; &emsp; &emsp; Where to place the colourbar </br>
&emsp; &emsp; &emsp; Options: {'bottom', 'top', 'left', 'right'} </br>
&emsp; &emsp; &emsp; Default: 'right' </br>
&emsp; &emsp; &#9656; show_axlabel: bool, optional </br>
&emsp; &emsp; &emsp; If True, writes labels on x and z axis </br>
&emsp; &emsp; &emsp; Default: True </br>
&emsp; &emsp; &#9656; show: bool, optional </br>
&emsp; &emsp; &emsp; If True, calls plt.show() to display plot </br>
&emsp; &emsp; &emsp; Default: False </br>

