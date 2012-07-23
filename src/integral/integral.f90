!--------------------------------------------------------------------------
!
!
!	integral.f90
!
!	Explaination:
!
!	Dependency: parameter_def.f90, radint.f90, angint.f90, pfunc.f90
!
!
!
!---------------------------------------------------------------------------
      
      
      subroutine integral3d(na, ns, nr, int_res)
      
      use parameter_def
      ! argument
      integer, intent(in) :: na, ns, nr
      real(double), intent(out) :: int_res
      ! local variables
      integer :: ka, js, ir, np
      real(double) :: int_rval(nr), int_sval(ns)
      real(double) :: rcoords(nr), rwgts(nr)
      real(double) :: scoords(3, ns), swgts(ns)
      real(double) :: accords(3, na)
      real(double) :: pcoords(3, ns * nr), pwgts(na, ns * nr)
      real(double) :: sine
      real(double) :: coord(3)
      
      np = ns * nr

      call Sec_Gauss_Chebyshev(nr, 1, rcoords, rwgts) 
      
      call lebsam(ns, scoords, swgts) 
      
      do ka = 1, na
         accords(1, ka) = grid3d(ka)%x
         accords(2, ka) = grid3d(ka)%y
         accords(3, ka) = grid3d(ka)%z
      end do

      do ir = 1, nr
         do is = 1, ns
            pcoords(1, (ir - 1) * ns + is) = rcoords(ir) * scoords(1, is) 
            pcoords(2, (ir - 1) * ns + is) = rcoords(ir) * scoords(2, is) 
            pcoords(3, (ir - 1) * ns + is) = rcoords(ir) * scoords(3, is) 
         end do
      end do
      
      call cal_patition(acoords, na, pcoords, np, pwgts)
       
      int_res = 0.0d0
      do ka = 1, na
         int_sval(ka) = 0.0d0
         do js = 1, ns
            int_rval(js) = 0.0d0
            do ir = 1, nr
               sine = 0.0d0
               coord(1) = 0.0d0
               coord(2) = 0.0d0
               coord(3) = 0.0d0
               coord(1) = rcoords(ir) * scoords(1, is) 
               coord(2) = rcoords(ir) * scoords(2, is)
               coord(3) = rcoords(ir) * scoords(3, is)
               sine = coord(3) / rcoords(ir)
               sine = sqrt(1 - tmp ** 2)
               ! or F(p(coord))
               int_rval(js) = int_rval(js) + rwgts(ir) * ( pwgts(ka, (ir - 1) * ns + is) * F(p(ir, is)) * sine ) 
            end do
            int_sval(ka) = int_sval(ka) + swgts(js) * int_rval(js)
         end do
         int_res = int_res + int_sval(ka)
      end do 
      
      end subroutine integral3d
