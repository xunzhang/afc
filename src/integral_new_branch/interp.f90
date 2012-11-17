!----------------------------------------------------------------------
!
!
!
!   interp.f90
!
!
!
!----------------------------------------------------------------------


      subroutine interp(ia, local_coords, charx)
      
      use parameter_def
      use common_data    
      implicit none

      ! argument
      integer, intent(in) :: ia
      real(double), intent(in) :: local_coords(3)
      real(double), intent(out) :: charx
      
      ! local variables 
      integer, parameter :: nord1 = 2
      integer, parameter :: nord2 = 1
      integer :: ib
      integer :: kpl, jmn, jmx, iip, jjp
      real(double) :: rnu
      real(double) :: pintch, term, norm, ted

      charx = 0.0d0
      ib = basistag(ia)
      rnu = sqrt((local_coords(1) - coord(1, ia)) ** 2 + &
                 (local_coords(2) - coord(2, ia)) ** 2 + &
                 (local_coords(3) - coord(3, ia)) ** 2)
      rnu = max1(minrnu12, rnu)
      
      pintch = 0.0d0
      kpl = basis(ib)%ngrid + basis(ib)%hinv * log(rnu / basis(ib)%radius)
      jmn = max0(1, min0(basis(ib)%ngrid - nord1, kpl - nord2))
      jmx = jmn + nord1
      do jjp = jmn, jmx
         
         term = 1.0d0
         norm = 1.0d0
         do iip = jmn, jmx
            if(iip /= jjp) then
               term = term * (rnu - basis(ib)%rgrid(iip))
               norm = norm * (basis(ib)%rgrid(jjp) - basis(ib)%rgrid(iip))
            end if
         end do
         
         ted = term / norm
         pintch = pintch + ted * rho(jjp, ia)
      end do
      
      charx = charx + pintch

      end subroutine interp
