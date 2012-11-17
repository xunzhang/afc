!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	ifort Lebedev-Laikov.f90 radint.f90 angint.f90 pfunc.f90 test_pfunc.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      program test_cal_patition
      
         use parameter_def
      
         integer :: i, j, k
         integer :: na, np, nr, ns
         real(double) :: atom_coords(3, 2), pcoords(3, 5500)
         real(double) :: pwgts(2, 5500) 
         
         real(double) :: rcoords(50), rwgts(50) 
         real(double) :: scoords(3, 110), swgts(110)
         
         na = 2 
         nr = 50
         ns = 110
         np = 5500
        
         atom_coords(1, 1) = 1.0d0
         atom_coords(2, 1) = 1.0d0
         atom_coords(3, 1) = 1.0d0
 
         atom_coords(1, 2) = 2.0d0
         atom_coords(2, 2) = 2.0d0
         atom_coords(3, 2) = 2.0d0
        
         call Sec_Gauss_Chebyshev(nr, 1, rcoords, rwgts)

         call lebsam(ns, scoords, swgts)

         do i = 1, nr
            do j = 1, ns
               do k = 1, 3
                  pcoords(k, (i - 1) * ns + j) = rcoords(i) * scoords(k, j)
               end do
            end do
         end do

         call cal_patition(atom_coords, na, pcoords, np, pwgts)
         
         do i = 1, na
            do j = 1, np
               write(*, *)
               print *, "atom i: "
               write(*, *) i
               print *, "int point j: "
               write(*, *) j
               print *, "its weight is:"
               write(*, *) pwgts(i, j)
               write(*, *)
            end do
         end do
         print *, "Test pfunc.f90 finished!"

      end
