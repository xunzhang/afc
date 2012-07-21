!--------------------------------------------------------------------------------------------------
!
!
!	pfunc.f90
!
!	Dependency: parameter_def.f90
!
!
!
!--------------------------------------------------------------------------------------------------

      function dist(a, b) result(distance)

      use parameter_def
      !implicit none  

      ! argument
      real(double), intent(in) :: a(3), b(3) 
      real(double) :: distance  
      ! local variables
      integer :: i

      distance = 0.0d0
      do i = 1, 3
         distance = distance + (a(i) - b(i)) ** 2            
      end do

      distance = sqrt(distance)

      end function dist


      function correlation(coord, center1, center2) result(mu)

      use parameter_def
      !implicit none

      ! argument
      real(double), intent(in) :: coord(3), center1(3), center2(3)
      real(double) :: mu

      mu = (dist(coord, center1) - dist(coord, center2)) / dist(center1, center2)

      end function correlation


      function func_p(mu) result(p)

      use parameter_def
      !implicit none

      ! argument
      real(double), intent(in) :: mu
      real(double) :: p

      p = 1.5 * mu - 0.5 * mu ** 3;

      end function func_p 


      function func_s(mu) result(s)

      use parameter_def
      !implicit none

      ! argument
      real(double), intent(in) :: mu
      real(double) :: s
      ! local variables
      real(double) :: tmp

      tmp = 0.0d0
      tmp = func_p(mu)
      tmp = func_p(tmp)
      tmp = func_p(tmp)
      s = 0.5 * (1 - tmp)   

      end function func_s


      subroutine normalize(pwgts, na, np)

      use parameter_def
      implicit none
      ! argument
      real(double), intent(inout) :: pwgts(na, np)
      integer, intent(in) :: na, np
      ! local variables
      integer :: i, j, k
      real(double) :: psum

      do i = 1, np 
         do j = 1, na
            psum = 0.0d0
            do k = 1, na
            psum = psum + pwgts(k, i)
            end do
            pwgts(j, i) = pwgts(j, i) / psum
         end do
      end do

      end subroutine normalize


      subroutine cal_patition(acenter_coords, na, coords, np, pwgts)

      use parameter_def
      !implicit none

      ! argument
      integer, intent(in) :: na, np
      real(double), intent(in) :: acenter_coords(3, na), coords(3, np)
      real(double), intent(out) :: pwgts(na, np)
      ! local variables
      integer :: ia, ip, ja
      real(double) :: mu, sres

      do ip = 1, np
         do ia = 1, na
            pwgts(ia, ip) = 1.0d0
            do ja = 1, na
               if(ia /= ja) then

                  mu = 0.0d0
                  sres = 0.0d0
                  mu = correlation(coords(:, ip), acenter_coords(:, ia), acenter_coords(:, ja))
                  sres = func_s(mu)
                  pwgts(ia, ip) = pwgts(ia, ip) * sres

               end if
            end do
         end do
      end do

      call normalize(pwgts, na, np) 

      end subroutine cal_patition

