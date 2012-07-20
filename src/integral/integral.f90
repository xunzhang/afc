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
      
      function integrand(pcoords, )
      end function integrand
      
      
      subroutine 3d_integral(na, ns, nr, int_res)
      
      use parameter_def
      ! argument
      integer, intent(in) :: na, ns, nr
      real(double), intent(out) :: int_res
      ! local variables
      integer :: ka, js, ir, np
      real(double) :: int_rval(nr), int_sval(ns)
      real(double) :: rcoords(nr), rwgts(nr)
      real(double) :: scoords(3, ns), swgts(ns)
      real(double) :: accords(na, 3)
      np = ns * nr
      real(double) :: pcoords(np, 3), pwgts(np)

      call Sec_Gauss_Chebyshev(nr, 1, rcoords, rwgts) 
      
      call lebsam(ns, scoords, swgts) 
      
      accords = 
      
      do ir = 1, nr
         do is = 1, ns
            pcoords(,) = rcoords(ir) * 
         end do
      end do
      
      call cal_patition(acoords, na, pcoords, np, pwgts)
       
      int_res = 0.0d0
      do ka = 1, na
         int_sval(ka) = 0.0d0
         do js = 1, ns
            int_rval(js) = 0.0d0
            do ir = 1, nr
               int_rval(js) = int_rval(js) + rwgts(ir) * &
               ( pwgts(ka, ir * ns + is) * F(p(ir, is)) * sin(ssss) )
            end do
            int_sval(ka) = int_sval(ka) + swgts(js) * int_rval(js)
         end do
         int_res = int_res + int_sval(ka)
      end do 
      
      end subroutine 3d_integral
