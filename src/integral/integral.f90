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
      
      function exp_func(x) result(y)
         use parameter_def
 
         ! argument
         real(double), intent(in) :: x
         real(double) :: y
         y = exp(10 - x)

      end function exp_func
     
      subroutine integral3d(na, ns, nr, int_res)
      
      use parameter_def
      use common_data
      !implicit none

      ! argument
      integer, intent(in) :: na, ns, nr
      real(double), intent(out) :: int_res
      ! local variables
      integer :: is, ka, js, ir, np
      integer :: ib, ja
      real(double) :: int_rval(ns), int_sval(na)
      ! real(double) :: int_rval(nr), int_sval(ns)
      real(double) :: rcoords(nr), rwgts(nr)
      real(double) :: scoords(3, ns), swgts(ns)
      real(double) :: acoords(3, na)
      !real(double) :: pcoords(3, ns * nr), pwgts(na, ns * nr)
      real(double) :: pcoords(3, 17000), pwgts(18, 17000)
      real(double) :: sine
      real(double) :: charx
      real(double) :: tmp_coord(3)
      !real(double) :: charx(ns * nr)
      real(double) :: zn_atom
      real(double) :: para

      print *, "I am here"
      np = ns * nr
      para = 0.6614040960026232
      
      ! generate points in radial direction, the nr rcoords is relative coordinates to atom
      call Sec_Gauss_Chebyshev(nr, para, rcoords, rwgts) 
      
      ! generate points in angular direction, the ns scoords is relative coordinates to atom
      call lebsam(ns, scoords, swgts) 
     
      ! get the coordinates of na atoms, notice raw * rbuff
      do ka = 1, na
         ja = atomtag(ka)
         acoords(1, ka) = coord(1, ja)
         acoords(2, ka) = coord(2, ja)
         acoords(3, ka) = coord(3, ja)
      end do
      
      ! to each atom, still the relative coordinates, but combine radial and angular
      do ir = 1, nr
         do is = 1, ns
            pcoords(1, (ir - 1) * ns + is) = rcoords(ir) * scoords(1, is) 
            pcoords(2, (ir - 1) * ns + is) = rcoords(ir) * scoords(2, is) 
            pcoords(3, (ir - 1) * ns + is) = rcoords(ir) * scoords(3, is) 
         end do
      end do
      
      call cal_patition(acoords, na, pcoords, np, pwgts)
      
      int_res = 0.0d0
      do ka = 9, 9
         int_sval(ka) = 0.0d0
         do js = 1, ns
            int_rval(js) = 0.0d0
            do ir = 1, nr
               sine = 0.0d0
               tmp_coord(1) = 0.0d0
               tmp_coord(2) = 0.0d0
               tmp_coord(3) = 0.0d0
               tmp_coord(1) = rcoords(ir) * scoords(1, js) 
               tmp_coord(2) = rcoords(ir) * scoords(2, js)
               tmp_coord(3) = rcoords(ir) * scoords(3, js)
               sine = tmp_coord(3) / rcoords(ir)
               sine = sqrt(1 - sine ** 2)
               
               call interp(ka, tmp_coord, charx)
               ! or F(p(coord))
               ! int_rval(js) = int_rval(js) + rwgts(ir) * ( pwgts(ka, (ir - 1) * ns + is) * F(p(ir, is)) * sine ) 
               ! int_rval(js) = int_rval(js) + rwgts(ir) * (pwgts(ka, (ir - 1) * ns + js) * charx * sine) 
               int_rval(js) = int_rval(js) + rwgts(ir) * charx * sine
               ! charx = exp_func(rcoords(ir))
               ! int_rval(js) = int_rval(js) + rwgts(ir) * charx
               write(*, *) "rcoords", rcoords(ir)
               write(*, *) "rwgts(ir) is", rwgts(ir)
               write(*, *) "swgts is", swgts(js)
               write(*, *) "charx", charx
               write(*, *)
            end do
            int_sval(ka) = int_sval(ka) + swgts(js) * int_rval(js)  
            write(*, *) "int result", int_sval(ka) * 4 * pi
         end do
         int_sval(ka) = int_sval(ka) * 4 * pi
         int_res = int_res + int_sval(ka)
      end do 
      
      do ka = 1, 1
         ib = basistag(ka)
         zn_atom = basis(ib) % zn
         print *, "zn print"
         write(*, *) zn_atom
      end do

      end subroutine integral3d
