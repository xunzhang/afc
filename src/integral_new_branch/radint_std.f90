!-------------------------------------------------------------------------------------------------------------
!	
!
!	radint.f90
!
!	Calculate the integral pattern as "[-1,1] sqrt(1-x**2)*f(x) dx", using Second-Guass-Chebyshev method 
!	
!	Dependency: parameter_def.f90
!
!
!-------------------------------------------------------------------------------------------------------------

      
      subroutine Sec_Gauss_Chebyshev_std(n, x, wgts)
      
      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      real(double), intent(out) :: x(n), wgts(n)
      ! local variables
      integer :: ip
      real(double) :: coff
      real(double) :: theta
      real(double) :: tmp_sum, tmp_diff

      if(n.lt.1) then
         write(*, *)
         write(*, '("Error: n < 1")')
         write(*, *)
         stop
      end if
      
      coff = 0.0d0
      coff = pi / (n + 1)

      do ip = 1, n

         theta = 0.0d0
         tmp_sum = 0.0d0
         tmp_diff = 0.0d0
         
         theta = ip * coff
         x(ip) = cos(theta) ! 'xk = cos(k*pi/(n+1))'
         
         tmp_sum = 1 + x(ip)
         tmp_diff = 1 - x(ip)
         
         wgts(ip) = tmp_sum * tmp_diff
         wgts(ip) = coff * wgts(ip)
         
      end do
      end subroutine Sec_Gauss_Chebyshev_std


      subroutine test_func(x, y) ! test_func = 1/sqrt((1-x^2))
      
      use parameter_def
      implicit none
      real(double), intent(in) :: x
      real(double), intent(out) :: y
      
      y = 1 / sqrt(1 - x ** 2)

      end subroutine test_func

      
      subroutine test_func2(x, y) ! test_func2 = 2 * sin((1+x)/(1-x))/((1+x)/(1-x))^3 * ((1+x)^1.5 / (1-x)^4.5)
      use parameter_def
      implicit none
      real(double), intent(in) :: x
      real(double), intent(out) :: y
      real(double) :: tmp_add, tmp_diff, tmp
      
      tmp_add = 1 + x
      tmp_diff = 1 - x
      tmp = tmp_add / tmp_diff
      
      y = 2 * (sin(tmp) / (tmp ** 3)) * ((tmp_add ** 1.5) / (tmp_diff**4.5))
      
      end subroutine test_func2


      program test_Sec_Gauss_Chebyshev_std ! result is 2 for test_func & pi/2 for test_func2
         
         use parameter_def
         implicit none
         integer :: i
         integer :: num = 5000, p = 1
         real(double) x(5000), wgts(5000), value(5000)
         real(double) res
         
         call Sec_Gauss_Chebyshev_std(num, x, wgts)
         
         res = 0.0d0
      
         do i = 1, num
            
            value(i) = 0.0d0
            
            call test_func(x(i), value(i))
            
            res = res + value(i) * wgts(i)
            
            write(*, *)
            print *, "x(i) is"
            write(*, *) x(i)
            print *, "value(i) is"
            write(*, *) value(i)
            print *, "wgts(i) is"
            write(*, *) wgts(i)
            write(*, *)
         end do
         
         write(*, *) res
         write(*, *)
         print *, "Sec_Gauss_CHebyshev int Test Finished!"
      
      end
