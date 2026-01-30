Here you will find a set of scripts giving examples of basic PyXWIND usage. These are:

### mk_lineProfile.py
This demonstrates initiating a XWIND object, using the `xw_line` class. It then shows how to create a line profile, using both `xwindline`
and `xwindfe`. The results are then plotted. In the latter case, the plot also shows the input spectrum, and resulting line Equivalent Width.
The relevant parameters are aviablable for the user to explore.

### convolve_spec.py
This demonstrates the use of PyXWIND as a convolution routine, using the `xw_conv` class. This starts by defining an input spectrum, 
consisting of two seperated narrow Gaussian components (for demonstrative purposes. This could easily be extended to work on any
input spectrum). Then defining an XWIND object, it applies the convolution routine, to then give the output spectrum. Both input
and outputs are plotted.

### mk_windProfile.py
This demonstrates one of the main advantages of PyXWIND over standard Xspec XWIND. The ability to plot profiles of the wind structure
(i.e density, emissivity, velocity, etc). The demonstration is done using `xw_line`, but it would be identical in `xw_conv` (since they
are both subclasses of the `PyXWIND_object` class). It starts by defining an XWIND object, which sets the wind structure.
It then plots the density, volumetric emissiity, and velocity profiles, using the internal `plot_wind()` method.

### mk_NhProfile.py
Demonstrated how to calculate a column-density profile (i.e $N_{H}$ as a function of inclination) for a given wind structure.
As usual, starts by definig an XWIND object. It then calls the `calc_Nh_profile()` method to give the result $N_{H}$ profile. 
