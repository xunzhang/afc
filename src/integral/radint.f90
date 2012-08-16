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
 
      subroutine Sec_Gauss_Chebyshev(n, p, r, wgts)

      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      real(double), intent(in) :: p
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
         write(*, '("Error in radint.f90: n < 1")')
         write(*, *)
         stop
      end if
     
      write(*, *) "ppppppppppppppppppppppppp is", p

      coff1 = 0.0d0
      coff2 = 0.0d0
      coff1 = pi / (n + 1)
      coff2 = 2.0d0 * (p ** 3)
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
