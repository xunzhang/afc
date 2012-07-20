!------------------------------------------------------------------------------
!
!
!	radint.f90
!	
!	Calcalute the integral pattern as "[0,+] f(R)*R^2 dR" using second guass-chebyshev method 
!
!	Dependency: parameter_def.f90
!
!
!------------------------------------------------------------------------------
 
      module radint    
 
      subroutine Sec_Gauss_Chebyshev(n, p, r, wgts)

      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      integer, intent(in) :: p
      real(double), intent(out) :: r(n), wgts(n)
      ! local variables
      integer :: ip
      real(double) :: coff1, coff2
      real(double) :: theta
      real(double) :: tmp_sum, tmp_diff
      real(double), allocatable :: x(:)
      
      allocate(x(n))

      if(n.lt.1) then
         write(*, *)
         write(*, '("Error: n < 1")')
         write(*, *)
         stop
      end if
      
      coff1 = 0.0d0
      coff2 = 0.0d0
      coff1 = pi / (n + 1)
      coff2 = 2 * p ** 3 
      coff2 = coff2 * coff1

      do ip = 1, n

         theta = 0.0d0
         tmp_sum = 0.0d0
         tmp_diff = 0.0d0
         
         x(ip) = 0.0d0
         r(ip) = 0.0d0
         wgts(ip) = 0.0d0

         theta = ip * coff1
         x(ip) = cos(theta) ! 'xk = cos(k*pi/(n+1))'
         tmp_sum = 1 + x(ip)
         tmp_diff = 1 - x(ip)
         !write(*, *) tmp_sum
         !print *, "hel"
         !write(*, *) tmp_diff
         r(ip) = tmp_sum / tmp_diff
         r(ip) = p * r(ip) ! 'rk = p*(1+xk)/(1-xk)'
         wgts(ip) = (tmp_sum ** 2.5) / (tmp_diff ** 3.5)
         !print *, "del"
         !write(*, *) wgts(ip)
         wgts(ip) = coff2 * wgts(ip)
         
      end do
      end subroutine Sec_Gauss_Chebyshev 


      subroutine Dirichlet_func(x, y)

      use parameter_def
      implicit none
      real(double), intent(in) :: x
      real(double), intent(out) :: y
      
      y = sin(x) / x ** 3

      end subroutine Dirichlet_func


      program test_Sec_Gauss_Chebyshev
      !-----------------------------------------------------------------------------------------------------------------------------
      !	
      !	use func: f(R) = sinR/R^3 to test the precision the Dirichlet result of int(f(R)R^2, [0,+]) is pi/2=1.570796
      !
      !-----------------------------------------------------------------------------------------------------------------------------

         use parameter_def
         implicit none
         integer :: i
         integer :: num = 50, p = 1
         real(double) r(50), wgts(50), value(50)
         real(double) res
         
         call Sec_Gauss_Chebyshev(num, p, r, wgts)
          
         res = 0.0d0
      
         do i = 1, num
            value(i) = 0.0d0
            call Dirichlet_func(r(i), value(i))
            res = res + value(i) * wgts(i)
            write(*, *)
            print *, "r(i) is "
            write(*, *) r(i)
            print *, "value(i) is "
            write(*, *) value(i)
            print *, "wgts(i) is "
            write(*, *) wgts(i)
            print *, "res is "
            write(*, *) res
            write(*, *)
         end do
         
         print *, "Sec_Gauss_CHebyshev int Test Finished!"
      end

      end module radint
