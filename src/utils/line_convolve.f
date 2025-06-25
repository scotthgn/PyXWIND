c      Contains subroutines for doing convolution between an emission
c      line profile and an input spectrum
c
c      Since line width in energy units depends on central energy
c      (i.e line profile naturally in E/E0), this routine accounts for
c      the change in absolute line width throughout the convolution

       subroutine convolve_line(ear,ph,e_kern,ph_kern,ne)
c      Does the convolution - only routine any user should be calling here
       implicit none

       integer ne
       real ear(0:ne), e_kern(0:ne)
       real ph(ne), ph_in(ne), ph_kern(ne)

       real emid_i, emid_j, ee_ij
       real ee_min, ee_max
       integer idx_tran

       logical is_log, is_lin
       real ldE, deps

       integer i, j

cf2py  intent(in) ear
cf2py  intent(in) ee_kern
cf2py  intent(in) ph_kern
cf2py  intent(in) ne
cf2py  intent(in) ph
 
cf2py  intent(out) ph

c      PARS
c      -----
c      ear : Input spectrum energy bins
c      ph : Input/Output spectrum (in photons/s/cm^2/bin) - to be convolved
c      ee_kern : Convolution kernel (line profile) fractional energy bin (E/E0)
c      ph_kern : Convlution kernel (i.e line profile) - normalised to unity
c      ne : Number of energy bins


c      Copying input array and zeroing output
       do i=1, ne, 1
          ph_in(i) = ph(i)
          ph(i) = 0.0
       end do


c      Checking if energy grid is log or lin
       call check_ebinsLog(e_kern,is_log,ne)
       if (is_log) then
          ldE = log10(e_kern(1)) - log10(e_kern(0))
       else
          call check_ebinsLin(e_kern,is_lin,ne)
          if (is_lin) then
             ldE = e_kern(1) - e_kern(0)
          end if
       end if
       
c      Performing convolution
       ee_min = e_kern(0)
       ee_max = e_kern(ne)
       do i=1, ne, 1
          emid_i = 0.5*(ear(i) + ear(i-1))
          do j=1, ne, 1
             emid_j = 0.5*(ear(j) + ear(j-1))
             ee_ij = emid_i/emid_j

             if ((ee_ij.ge.ee_max).or.(ee_ij.le.ee_min)) then
                continue
             else   
                if (is_log) then
                   idx_tran = ceiling((log10(ee_ij)-log10(ee_min))/ldE)
                   ph(i) = ph(i) + ph_in(j)*ph_kern(idx_tran)

                else if (is_lin) then
                   deps = ldE/emid_j
                   idx_tran = ceiling((ee_ij-ee_min)/ldE)
                   ph(i) = ph(i) + ph_in(j) * ph_kern(idx_tran)

                else
                   call direct_bin_search(e_kern,ph,ph_in,ph_kern,ee_ij,
     $                  i,j,ne)
                end if
             end if
          end do
       end do

       
       return
       end



       subroutine direct_bin_search(ee_ar,ph,ph_in,ph_kern,ee_ij,i,j,ne)
c      Subroutine for doing a direct seach through each energy bin
c      Only needs to be called if input energy grid is neither linear nor log
c      as in this case cannot directly calculate bin idx.
c
c      This is slow for large arrays (when combined with the main convolution routine)
c      so recommendation is to use either linearly or geometrically spaced bins
c
c      Is called automatically by convolve_line if needed
       implicit none

       integer ne
       real ee_ar(0:ne), ph(ne), ph_in(ne), ph_kern(ne)
       real ee_ij
       integer i, j, k

cf2py  intent(in) ee_ar
cf2py  intent(in) ph
cf2py  intent(in) ph_in
cf2py  intent(in) ph_kern
cf2py  intent(in) ee_ij
cf2py  intent(in) i
cf2py  intent(in) j
cf2py  intent(in) ne

cf2py  intent(out) ph

c      PARS
c      -----
c      ee_ar : kernel energy array (i.e in E/E0)
c      ph : Photon array being convolved
c      ph_in : Original (copy of) photon array
c      ph_kern : Convolution kernel (i.e line profile)
c      ee_ij : Current energy bin in E/E0
c      i : index of upper bin loop
c      j : index of lower bin loop
c      ne : Number of energy bins

       citr: do k=1, ne, 1
          if ((ee_ij.lt.ee_ar(k)).and.(ee_ij.gt.ee_ar(k-1))) then
             ph(i) = ph(i) + ph_in(j)*ph_kern(k)
             exit citr
          else
             continue
          end if
       end do citr

       return
       end

