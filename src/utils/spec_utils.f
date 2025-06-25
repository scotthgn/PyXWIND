c     Contains usfule utility routines
c     Currently includes re-normalising emission line and checking spacing in energy grid

      
       subroutine re_norm(ph,ph_tot,ne)
c      Re-bins an input spectrum
c      Input assumed to be in units of photons/(s/cm^2)/bin
c      Hence summing rather than integrating
c      (Converting to photons/s/cm^2/keV will naturally fold in the integration factor) 
       implicit none

       integer ne
       real ph(ne)
       real ph_tot

       real ph_sum
       integer i
       
cf2py  intent(in) ne
cf2py  intent(in) ph_tot
cf2py  intent(in) ph
c
cf2py  intent(out) ph
       
c      PARS
c      -----
c      ph : Input photon spectrum to be re-normalised
c      ph_tot : New normalisation (total photon flux)
c      ne : Number of energy bins
       
       ph_sum = 0.0
       do i=1, ne, 1
          ph_sum = ph_sum + ph(i)
       end do
       ph = ph*(ph_tot/ph_sum)
       
       return
       end
       


       subroutine check_ebinsLog(ear,is_log,ne)
c      Checks if energy bins are evenly distributed in log
       implicit none

       integer ne
       real ear(0:ne)
       logical is_log

       real precision, erat, crat, fchange
       integer i
       
cf2py  intent(in) ear
cf2py  intent(in) ne
cf2py  intent(out) is_log

c      PARS
c      -----
c      ear : Input energy bin edges
c      is_log : Boolean flag - updated by subroutine
c      ne : Number of energy bins

c      Initialising
       is_log = .true.
       precision = 1e-6 !since using real
       erat = ear(1)/ear(0)
c      Iterating over bins
       aitr: do i=1, ne, 1
          crat = ear(i)/ear(i-1)
          fchange = abs(erat-crat)/erat
          if (fchange.gt.precision) then
             is_log = .false.
             exit aitr
          else
             continue
          end if
       end do aitr

       return
       end



       subroutine check_ebinsLin(ear,is_lin,ne)
c      Checks if energy bins are linearly spaced
       implicit none

       integer ne
       real ear(0:ne)
       logical is_lin

       real dini, di, dchange, precision
       integer i

cf2py  intent(in) ear
cf2py  intent(in) ne
cf2py  intent(out) is_lin

c      PARS
c      -----
c      ear : Input energy bin edges
c      is_log : Boolean flag - updated by subroutine
c      ne : Number of energy bins

       is_lin = .true.
       precision = 1e-6
       dini = ear(1) - ear(0)
       litr: do i=1, ne, 1
          di = ear(i) - ear(i-1)
          dchange = abs(dini - di)
          if (dchange.gt.precision) then
             is_lin = .false.
             exit litr
          else
             continue
          end if
       end do litr

       return
       end
       
