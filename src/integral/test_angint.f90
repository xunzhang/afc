

      program test_lebsam

         use parameter_def
         implicit none
         
         integer :: i
         integer :: num = 110
         real(double) scoords(3, 110), swgts(110)

         call lebsam(num, scoords, swgts)

         do i = 1, num
            write(*, *)
            write(*, *) scoords(1, i)
            write(*, *) scoords(2, i)
            write(*, *) scoords(3, i)
            write(*, *) swgts(i)
            write(*, *)
         end do

         print *, "Test Finished!"
      end
