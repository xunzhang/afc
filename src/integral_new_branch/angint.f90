!-------------------------------------------------------------------------
! 	angint.f90
! ------------------------------------------------------------------------
!	Dependency:
! 		
!-------------------------------------------------------------------------


      subroutine lebsam(n, coords, wgts)
      use parameter_def
      !use parameter_def

      implicit none

! argument
      integer, intent(in) :: n
      real(double), intent(out) :: coords(3, n), wgts(n)

! local variables
      integer :: m
      real(double), allocatable :: x(:), y(:), z(:)
      integer :: i

      if(n.lt.1) then
         write(*, *)
         write(*, '("Error in angint.f90: n < 1")')
         write(*, *)
         stop
      end if
      
      m = n
      
      allocate(x(n), y(n), z(n))
      
      select case(n)
      case(6)
         call LD0006(x, y, z, wgts, m)
      case(14)
         call LD0014(x, y, z, wgts, m)
      case(26)
         call LD0026(x, y, z, wgts, m)
      case(38)
         call LD0038(x, y, z, wgts, m)
      case(50)
         call LD0050(x, y, z, wgts, m)
      case(74)
         call LD0074(x, y, z, wgts, m)
      case(86)
         call LD0086(x, y, z, wgts, m)
      case(110)
         call LD0110(x, y, z, wgts, m)
      case(146)
         call LD0146(x, y, z, wgts, m)
      case(170)
         call LD0170(x, y, z, wgts, m)
      case(194)
         call LD0194(x, y, z, wgts, m)
      case(230)
         call LD0230(x, y, z, wgts, m)
      case(266)
         call LD0266(x, y, z, wgts, m)
      case(302)
         call LD0302(x, y, z, wgts, m)
      case(350)
         call LD0350(x, y, z, wgts, m)
      case(434)
         call LD0434(x, y, z, wgts, m)
      case(590)
         call LD0590(x, y, z, wgts, m)
      case(770)
         call LD0770(x, y, z, wgts, m)
      case(974)
         call LD0974(x, y, z, wgts, m)
      case(1202)
         call LD1202(x, y, z, wgts, m)
      case(1454)
         call LD1454(x, y, z, wgts, m)
      case(1730)
         call LD1730(x, y, z, wgts, m)
      case(2030)
         call LD2030(x, y, z, wgts, m)
      case(2354)
         call LD2354(x, y, z, wgts, m)
      case(2702)
         call LD2702(x, y, z, wgts, m)
      case(3074)
         call LD3074(x, y, z, wgts, m)
      case(3470)
         call LD3470(x, y, z, wgts, m)
      case(3890)
         call LD3890(x, y, z, wgts, m)
      case(4334)
         call LD4334(x, y, z, wgts, m)
      case(4802)
         call LD4802(x, y, z, wgts, m)
      case(5294)
         call LD5294(x, y, z, wgts, m)
      case(5810)
         call LD5810(x, y, z, wgts, m)
      case default
         write(*, *)
         write(*, '("Error: invalid number of covering points!")')
         write(*, *)
      end select
     
      do i = 1, n
         wgts(i) = wgts(i) * 4 * pi
      end do

      coords(1, :) = x(:)
      coords(2, :) = y(:)
      coords(3, :) = z(:)
      
      deallocate(x, y, z)
      
      end subroutine lebsam 
