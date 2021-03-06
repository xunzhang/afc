!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	ifort angint.f90 Lebedev-Laikov.f90 test_angint.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      program test_lebsam

         use parameter_def
         implicit none
         
         integer :: i
         integer :: num = 110
         real(double) scoords(3, 110), swgts(110)
         real(double) area

         call lebsam(num, scoords, swgts)

         area = 0.0d0
         do i = 1, num
            write(*, *)
            write(*, *) scoords(1, i)
            write(*, *) scoords(2, i)
            write(*, *) scoords(3, i)
            write(*, *) swgts(i)
            area = area + swgts(i)
            write(*, *)
         end do

         print *, "Test angint.f90 finished!"
         write(*, *) area
      end
