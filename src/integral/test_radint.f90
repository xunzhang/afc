


      function Dirichlet_func(x) result(y)
      
      use parameter_def
      implicit none

      ! argument
      real(double), intent(in) :: x
      real(double) :: y
      
      y = sin(x) / x ** 3

      end function Dirichlet_func
      
      
      function exp_func(x) result(y)
      
      use parameter_def
      implicit none

      ! argument
      real(double), intent(in) :: x
      real(double) :: y
      print *, x 
      y = exp(10 - x)

      end function exp_func


      !	use 'f(R) = sinR/R^3' to test precision of int(f(R)*R^2, [0,+]) is 1.570796(pi/2)  --wrong test func, but why?
      ! use 'f(R) = e^(10-R)' to test precision of int(f(R)*R^2, [0,+]) is 44052.93155
      program test_Sec_Gauss_Chebyshev
      
         use parameter_def
      
         integer :: i
         integer :: num = 5000, p = 1
         real(double) r(5000), rwgts(5000), func_value(5000)
         real(double) res

         call Sec_Gauss_Chebyshev(num, p, r, rwgts)
      
         res = 0.0d0
         
         do i = 1, num
            func_value(i) = 0.0d0
            ! func_value(i) = Dirichlet_func(r(i))
            func_value(i) = exp_func(r(i))
            res = res + func_value(i) * rwgts(i)
            write(*, *)
            print *, "r(i) is"
            write(*, *) r(i)
            print *, "func_value(i) is"
            write(*, *) func_value(i)
            print *, "rwgts(i) is"
            write(*, *) rwgts(i)
            print *, "res is"
            write(*, *) res
            write(*, *)
         end do
      print *, "Test radint.f90 finished!"   
      end
