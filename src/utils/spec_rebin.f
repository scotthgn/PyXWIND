c      subroutines for re-binning an spectrum onto a new energy grid
c      Works in photons/s/cm^2 - so simply re-distributes photon counts
c      among new bins - conserving photon number
c     -----------------------------------------------------------------
      
       subroutine re_bin(ear,ear_new,ph,ph_new,ne,Nnew)
       implicit none

       integer ne, Nnew
       real ear(0:ne), ph(ne)
       real ear_new(0:Nnew), ph_new(Nnew)

       logical is_log, is_lin
       real ldE

       integer i

cf2py  intent(in) ear  
cf2py  intent(in) ear_new
cf2py  intent(in) ph
cf2py  intent(in) ne
cf2py  intent(in) Nnew
       
cf2py  intent(out) ph_new


       
       !initialising output array
       do i=1, Nnew, 1
          ph_new(i) = 0.0
       end do

       !Checking if input grid is log or lin
       call check_ebinsLog(ear_new,is_log,Nnew)
       if (is_log) then
          ldE = log10(ear_new(1)) - log10(ear_new(0))
       else
          call check_ebinsLin(ear_new,is_lin,Nnew)
          if (is_lin) then
             ldE = ear_new(1) - ear_new(0)
          end if
       end if
       
       if ((is_log).or.(is_lin)) then
          call rebin_evenE(ear,ne,ear_new,Nnew,ldE,ph,ph_new,is_log)
       else
          call rebin_noeven(ear,ear_new,ph,ph_new,ne,Nnew)
       end if
       
       return
       end

      
       subroutine rebin_evenE(ear,ne,ear_new,Nnew,dE,ph,ph_new,is_log)
c      If evenly spaced bins (in either lin or log) then can calculate
c      bin indexes analytically - no need to iterate
       implicit none

       integer ne, Nnew
       real ear(0:ne), ph(ne)
       real ear_new(0:Nnew), ph_new(Nnew)
       real dE
       
       integer enew_idx1, enew_idx2
       logical is_log

       real dE_old, fb1, fb2, fbi
       integer ndiff, eni_idx
       
       integer i, j

cf2py  intent(in) ear
cf2py  intent(in) ne       
cf2py  intent(in) ear_new
cf2py  intent(in) Nnew
cf2py  intent(in) dE
cf2py  intent(in) ph
cf2py  intent(in) ph_new
cf2py  intent(in) is_log

cf2py  intent(out) ph_new       
       
       do i=1, ne, 1
          if (is_log) then
             enew_idx1 = ceiling((log10(ear(i-1))-log10(ear_new(0)))/dE)
             enew_idx2 = ceiling((log10(ear(i))-log10(ear_new(0)))/dE)
          else
             enew_idx1 = ceiling((ear(i-1)-ear_new(0))/dE)
             enew_idx2 = ceiling((ear(i)-ear_new(0))/dE)
          end if
          dE_old = ear(i) - ear(i-1) !needed for fractional overlap
          
          if ((enew_idx1.ge.Nnew).or.(enew_idx2.le.0)) then !upper edge case
             continue

          else if (enew_idx1.eq.enew_idx2) then !new bin fully containes old
             ph_new(enew_idx1) = ph_new(enew_idx1) + ph(i)
             
          else if ((enew_idx2 - enew_idx1).gt.1) then !new bin smaller than old bins
             ndiff = enew_idx2 - enew_idx1
             do j=0, ndiff, 1
                eni_idx = enew_idx1 + j
                if (j.eq.0) then
                   fbi = (ear_new(enew_idx1) - ear(i-1))/dE_old
                else if (j.eq.ndiff) then
                   fbi = (ear(i) - ear_new(enew_idx2-1))/dE_old
                else
                   fbi = (ear_new(eni_idx) - ear_new(eni_idx-1))/dE_old 
                end if
                ph_new(eni_idx) = ph_new(eni_idx) + fbi*ph(i)
             end do
             
          else !new bins bigger, but not aligned
             fb1 = (ear_new(enew_idx1) - ear(i-1))/dE_old !fractional overlap
             fb2 = 1.0 - fb1

             ph_new(enew_idx1) = ph_new(enew_idx1) + fb1*ph(i)
             ph_new(enew_idx2) = ph_new(enew_idx2) + fb2*ph(i)
          end if
       end do
       
       return
       end



       subroutine rebin_noeven(ear,ear_new,ph,ph_new,ne,Nnew)
c      If energy bins not evenly spaced (in either lin or log)
c      then this does the binning iteratively
       implicit none

       integer ne, Nnew
       real ear(0:ne), ph(ne)
       real ear_new(0:Nnew), ph_new(Nnew)

       real fbi, dE_old
       
       integer i, j

       aitr: do i=1, ne, 1
          if (ear(i).le.ear_new(0)) then
             continue
          else if (ear(i-1).ge.ear_new(Nnew)) then
             exit aitr
          else
             dE_old = ear(i) - ear(i-1)
             bitr: do j=1, Nnew, 1
                if (ear_new(j).le.ear(i-1)) then
                   continue
                else if (ear_new(j-1).ge.ear(i)) then
                   exit bitr
                else
                   if (ear_new(j-1).le.ear(i-1)) then
                      if (ear_new(j).ge.ear(i)) then !total overlap
                         fbi = 1.0
                      else
                         fbi = (ear_new(j)-ear(i-1))/dE_old !partial, lower edge
                      end if
                   else
                      if (ear_new(j).le.ear(i)) then !new bin smaller than old
                         fbi = (ear_new(j)-ear_new(j-1))/dE_old
                      else  !partial, upper edge
                         fbi = (ear(i)-ear_new(j-1))/dE_old
                      end if
                   end if
                   ph_new(j) = ph_new(j) + fbi*ph(i)
                end if
              end do bitr
           end if
       end do aitr

       return
       end
