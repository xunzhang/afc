      function unit_func(r, r0) result(y)
      
      use parameter_def
      implicit none

      ! argument
      real(double), intent(in) :: r, r0
      real(double) :: y
      
      if(r<=r0) then
        y = 1
      else 
        y = 0.0d0
      endif
      end function unit_func 
      
      
      function exp_func(x) result(y)
      
      use parameter_def
      implicit none

      ! argument
      real(double), intent(in) :: x
      real(double) :: y
      print *, x 
      y = exp(10 - x)

      end function exp_func



      program test_single_atom_integral
      
         use parameter_def
          
         integer :: i
         real(double) :: para = 7.0 !0.6614040960026232 
         integer :: nr, ns = 5810 
         integer :: ir, js
         real(double), allocatable :: rcoords(:), rwgts(:)
         !real(double) :: rcoords(50000), rwgts(50000)
         real(double) :: scoords(3, 5810), swgts(5810)
         real(double) :: func_value, sine, volumn
         real(double) :: temp, temp2, r0, pr
         integer :: numbers
         character*10 :: argv(2)
         
         call getarg(1, argv(1))
         read(argv(1)(:), '(i)') numbers

         call getarg(2, argv(2))
         read(argv(2), *) temp2

         call getarg(3, argv(3))
         read(argv(3), *) temp

         nr = numbers 
         para = temp2 / temp
         r0 = temp2
         pr = temp
         
         allocate(rcoords(nr))
         allocate(rwgts(nr))
         
         !call Uniform(nr, 10.0d0, rcoords, rwgts)
         call LOG_MK1996(nr, para, rcoords, rwgts)
         !call Sec_Gauss_Chebyshev(nr, para, rcoords, rwgts)
         !call MultiExp(nr, para, rcoords, rwgts)
         
         call lebsam(ns, scoords, swgts) 
         volumn = 0.0d0
         do ir = 1, nr
            temp = 0.0d0
            do js = 1, ns
               !sine = sqrt(1 - scoords(3, js) ** 2)
               temp = temp + swgts(js) * unit_func(rcoords(ir), r0) !* sine
            end do
            !print *, "check volumn"
            !write (*, *) temp
            !write (*, *) rcoords(ir)
            volumn = volumn + rwgts(ir) * temp !* rcoords(ir) ** 2
         end do
         !print *, "volumn result is"
         !write (*, *) 'nr ', nr, 'r0 ', r0, 'para ', 'radio', pr, 'volumn ', volumn
         write (*, *) volumn
         deallocate(rcoords)
         deallocate(rwgts)
      end
