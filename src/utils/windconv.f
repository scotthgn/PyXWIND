c      Routine to smooth an input spectrum
c      Uses windline to calculate a line-profile at some arbitrary
c      rest energy. Then convolves in log energy, since
c      \Delta log E = \Delta log epsilon, where epsilon = E/E0
c
c      CAVEAT! The line-profile is re-normalised to conserve photon
c      number. Hence Mbh and mdot are ONLY used to set the density
c      profile, required for calculating the emissivity


       subroutine windconv(param,ear,photar,ne)
       implicit none

       integer npars
       parameter(npars=11)

       integer ne
       real ear(0:ne), photar(ne), param(npars)

       real photar_wnd(ne), wnd_pars(15)
       logical donorm

       integer i

cf2py  intent(in) param
cf2py  intent(in) ear
cf2py  intent(in) photar
cf2py  intent(in) ne

cf2py  intent(out) photar       

c      ARGS
c      -----
c      param : Array containing model parameters
c      ear : Array containing energy bin edges (keV)
c      photar : Array containing input spectrum to be convolved
c      ne : Number of energy bins

c      PARS
c      -----
c      param(1):    Mbh, Black hole mass, Msol
c      param(2):    mdot_w, wind mass outflow rate, Mdot/Mdot_edd
c      param(3):    r_in, Inner Launch radius, Rg
c      param(4):    r_out, Outer launch radius, Rg
c      param(5):    d_foci, Distance to wind foci, Rg
c      param(6):    fcov, wind covering fraction, Omega/4pi
c      param(7):    vinf, outflow velcoity at infinity, c
c      param(8):    rv, wind velocity scale length, Rg
c      param(9):    vexp, Wind velocity exponent
c      param(10):   kappa, radial density law exponent
c      param(11):   inc, observer inclination, deg


c      Uses windline code to calculate initial line profile
c      Start by filling parameter array for windline
       do i=1, npars, 1
          wnd_pars(i) = param(i)
       end do

       wnd_pars(12) = 1.0 !Abundance, irrelevant since re-normalises
       wnd_pars(13) = 0.5*(ear(ne) + ear(0)) !E0 - arbitrary so set to middle of energy range
       wnd_pars(14) = 1.0 !incident spec norm - arbitrary since re-normalises
       wnd_pars(15) = 2.0 !Incident spectrum photon index, arbitrary since re-normalises
       
c      Calculating windline and applying convolution
       donorm = .true.
       call windline(wnd_pars,ear,photar_wnd,donorm,ne)
       call convolve_line(ear,photar,ear/wnd_pars(13),photar_wnd,ne)

       return
       end
