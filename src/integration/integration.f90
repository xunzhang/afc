      
      function interp(x, y, z) result(val)
         
         use parameter_def
         use common_data
         implicit none
         ! arguments
         real(double), intent(in) :: x, y, z
         real(double) :: val
         ! local variables
         integer, parameter :: nord1 = 2
         integer, parameter :: nord2 = 1
         integer :: ia, ib, ip
         real(double) :: rnu
         integer :: kpl, jmn, jmx, jjp, iip
         real(double) :: pintch, term, norm, ted
         integer :: i, j, k, l
          
         val = 0.0d0
         do ia = 1, natom
            ib = basistag(ia)
            rnu = dsqrt((x - coord(1, ia)) ** 2 + &
                        (y - coord(2, ia)) ** 2 + &
                        (z - coord(3, ia)) ** 2   )
            rnu = dmax1(minrnu12, rnu)
            pintch = 0.0d0
            if(rnu <= basis(ib) % radius) then
               kpl = basis(ib)%ngrid + basis(ib) % hinv * dlog(rnu / basis(ib)%radius)
               jmn = max0(1, min0(basis(ib)%ngrid-nord1, kpl - nord2))
               jmx = jmn + nord1
               do jjp = jmn, jmx
                  term = 1.0d0
                  norm = 1.0d0
                  do iip = jmn, jmx
                     if(iip /= jjp) then
                        term = term * (rnu - basis(ib) % rgrid(iip))
                        norm = norm * (basis(ib) %rgrid(jjp) - basis(ib) % rgrid(iip))
                     end if
                  end do
                  ted = term / norm
                  pintch = pintch + ted * rho(jjp, ia)
               end do
            end if
            val = val + pintch
         end do
          
      end function interp 
      

      function unit_func(r, r0) result(y)
  
         use parameter_def
         implicit none
         real(double), intent(in) :: r, r0
         real(double) :: y

         if(r <= r0) then
            y = 1.0d0
         else
            y = 0.0d0
         endif

      end function unit_func
 
 
      function exp_func(x) result(y)
 
         use parameter_def
         implicit none 
         real(double), intent(in) :: x
         real(double) :: y

         y = exp(10 - x)
 
      end function exp_func

       
      function dist3d(a, b) result(distance)
  
         use parameter_def
         ! argument
         real(double), intent(in) :: a(3), b(3)
         real(double) :: distance
         ! local variables
         integer :: i

         distance = 0.0d0
         do i = 1, 3
            distance = distance + (a(i) - b(i)) ** 2
         end do  
         distance = sqrt(distance)
  
      end function dist3d
 

      ! coord and charx is global visible
      subroutine integral3d(ncenter, nsphpt, nradpt, int_value)
      
      use parameter_def
      use common_data
      ! argument 
      integer, intent(in) :: ncenter, nsphpt, nradpt
      real(double), intent(out) :: int_value
      ! local vars
      integer :: iatom, i, j, k
      real(double) :: parm, radr, rwgt, miu, pwgt, single_int_value, sample_coord(3), r0, ra, rb, rab, tmp
      type(intpt), allocatable :: sample_point(:)
      real(double), allocatable :: scoords(:, :), swgts(:), sval(:, :), pdt_pwgt(:), func_val(:) 
      
      allocate(sample_point(nradpt * nsphpt), scoords(3, nsphpt), swgts(nsphpt), sval(ncenter, ncenter), pdt_pwgt(ncenter), func_val(nradpt * nsphpt))
      
      parm = 1.0d0
      
      call lebsam(nsphpt, scoords, swgts) 
      
      int_value = 0.0d0 
      do iatom = 1, ncenter
         ! call lebsam(nsphpt, scoords, swgts)
         do i = 1, nradpt
            !call Sec_Gaussian_Chebyshev(i, nradpt, parm, radr, rwgt)
            call LOG_MK1996(i, nradpt, parm, radr, rwgt)
            sample_point( (i - 1) * nsphpt + 1 : i * nsphpt) % x = radr * scoords(1, :)
            sample_point( (i - 1) * nsphpt + 1 : i * nsphpt) % y = radr * scoords(2, :)
            sample_point( (i - 1) * nsphpt + 1 : i * nsphpt) % z = radr * scoords(3, :)
            sample_point( (i - 1) * nsphpt + 1 : i * nsphpt) % wgt = 4 * pi * rwgt * swgts
         end do
         
         !do i = 1, nradpt * nsphpt
         !   write(*, *) "dim1 is: ", (sample_point(i) % x)
         !   write(*, *) "dim2 is: ", (sample_point(i) % y)
         !   write(*, *) "dim3 is: ", (sample_point(i) % z)
         !end do

         ! transfer to absolute coords
         sample_point % x = sample_point % x + coord(1, iatom)
         sample_point % y = sample_point % y + coord(2, iatom)
         sample_point % z = sample_point % z + coord(3, iatom)
         !write(*, *) "coord dim1 is ", coord(1, iatom)
         !write(*, *) "coord dim2 is: ", coord(2, iatom)
         !write(*, *) "coord dim3 is: ", coord(3, iatom)
         
         ! do interp to get charx in sampled position
         do i = 1, nradpt * nsphpt
            func_val(i) = interp(sample_point(i) % x, sample_point(i) % y, sample_point(i) % z)
            !write(*, *) "func_val is: ", func_val(i)
         end do 
         
         single_int_value = 0.0d0 
         ! single atom integration, traverse every sampling points
         do i = 1, nradpt * nsphpt
            sval = 1.0d0
            sample_coord(1) = sample_point(i) % x
            sample_coord(2) = sample_point(i) % y
            sample_coord(3) = sample_point(i) % z
            do j = 1, ncenter ! for every sampling point:i of atom:iatom, cal its pwgt relative to j
               ra = dsqrt( (sample_coord(1)-coord(1, j))**2 + (sample_coord(2)-coord(2,j))**2 + (sample_coord(3)-coord(3,j))**2 )
               do k = 1, ncenter ! to cal its pwgt relative to j, traverse k which is not equal to j
                  if(k == j) cycle
                  rb = dsqrt( (sample_coord(1)-coord(1, k))**2 + (sample_coord(2)-coord(2,k))**2 + (sample_coord(3)-coord(3,k))**2 )
                  rab = dsqrt( (coord(1, j)-coord(1, k))**2 + (coord(2, j)-coord(2,k))**2 + (coord(3, j)-coord(3,k))**2 )
                  miu = (ra - rb) / rab
                  tmp = 1.5d0 * (miu) - 0.5d0 * miu **3
                  tmp = 1.5d0 * (tmp) - 0.5d0 * tmp **3
                  tmp = 1.5d0 * (tmp) - 0.5d0 * tmp **3
                  sval(j, k) = 0.5d0 * (1 - tmp)
               end do
            end do 
            ! accumulate pweight 
            pdt_pwgt = 1.0d0
            do k = 1, ncenter 
               pdt_pwgt = pdt_pwgt * sval(:, k)
            end do
            ! Normalization pweight
            if(sum(pdt_pwgt) /= 0) then
               pwgt = pdt_pwgt(iatom) / sum(pdt_pwgt)
            end if
            single_int_value = single_int_value + pwgt * sample_point(i) % wgt * func_val(i)
         end do 
         int_value = int_value + single_int_value
         write(*, *) "single integration is: ", single_int_value
      end do
      write(*, *) "integration is: ", int_value
      end subroutine integral3d
