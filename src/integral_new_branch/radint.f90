!------------------------------------------------------------------------------
!
!	radint.f90
!	
!	Calcalute the integral pattern as "[0,+] f(R)*R^2 dR" using second guass-chebyshev method 
!
!	Dependency: parameter_def.f90
!
!------------------------------------------------------------------------------
 

      subroutine Uniform(n, r_bnd, r, wgts)

      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      real(double), intent(in) :: r_bnd
      real(double), intent(out) :: r(n), wgts(n)
      ! local variables
      integer :: ip
      real(double) :: interval

      interval = r_bnd / n
      do ip = 1, n
         r(ip) = ip * interval
         wgts(ip) = interval
      end do
      end subroutine Uniform 


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
         r(ip) = tmp_sum / tmp_diff
         r(ip) = p * r(ip) ! 'rk = p*(1+xk)/(1-xk)'
         wgts(ip) = (tmp_sum ** 2.5) / (tmp_diff ** 3.5)
         wgts(ip) = coff2 * wgts(ip)
         
      end do
      deallocate(x)
      end subroutine Sec_Gauss_Chebyshev 


      subroutine LOG_MK1996(n, p, r, wgts)

      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      real(double), intent(in) :: p
      real(double), intent(out) :: r(n), wgts(n)
      ! local variables
      integer :: ip
      real(double), allocatable :: x(:)
      real(double) temp

      allocate(x(n))
      do ip = 1, n
         x(ip) = ip * 1.0d0 / (n + 1)
         temp = log(1 - x(ip)**3)
         r(ip) = -p * temp
         wgts(ip) = p**3 * 3 * x(ip)**2 * temp**2
         wgts(ip) = wgts(ip) / ((n + 1) * (1 - x(ip)**3))
      end do
      deallocate(x)
      end subroutine LOG_MK1996 


      subroutine MultiExp(n, p, r, wgts)

      use parameter_def
      implicit none
      ! argument
      integer, intent(in) :: n
      real(double), intent(in) :: p
      real(double), intent(out) :: r(n), wgts(n)
      ! local variables
      integer :: ip
      real(double), allocatable :: x(:)

      allocate(x(n))

      select case(n)
      case(1)
         call RWLN1(x, wgts)
      case(2)
         call RWLN2(x, wgts)
      case(3)
         call RWLN3(x, wgts)
      case(4)
         call RWLN4(x, wgts)
      case(5)
         call RWLN5(x, wgts)
      case(6)
         call RWLN6(x, wgts)    
      case(8)
         call RWLN8(x, wgts)
      case(10)
         call RWLN10(x, wgts)     
      case(15)
         call RWLN15(x, wgts)
      case(20)
         call RWLN20(x, wgts)
      case default
         write(*, *)
         write(*, '("Error: invalid number of covering points used!")')
         write(*, *)    
      end select
      
      do ip = 1, n
         r(ip) = -p * log(x(ip))
         wgts(ip) = p**3 * wgts(ip) / x(ip)
      end do
         
      deallocate(x)
       
      end subroutine MultiExp 
