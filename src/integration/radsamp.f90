      
      subroutine Uniform(indx, nradpt, r_bnd, parm, radr, rwgt)
         use parameter_def
         implicit none
         ! argument
         integer, intent(in) :: indx, nradpt
         real(double), intent(in) :: r_bnd, parm
         real(double), intent(out) :: radr, rwgt
         ! local variables
         real(double) :: interval
         interval = 0.0d0
         radr = 0.0d0
         rwgt = 0.0d0
         interval = r_bnd / nradpt
         radr = indx * interval
         rwgt = interval
      end subroutine Uniform
       

      subroutine Sec_Gaussian_Chebyshev(indx, nradpt, parm, radr, rwgt)
         use parameter_def
         implicit none
         ! argument
         integer, intent(in) :: indx, nradpt
         real(double), intent(in) :: parm
         real(double), intent(out) :: radr, rwgt
         ! local variables
         real(double) :: radx
         radx = 0.0d0
         radr = 0.0d0
         rwgt = 1.0d0
         radx = cos(real(indx) * pi / real(nradpt + 1))
         radr = parm * (1 + radx) / (1 - radx)
         rwgt = 2.0d0 * pi * parm ** 3 / (nradpt + 1)
         rwgt = rwgt * (1 + radx) ** 2.5 / (1 - radx) ** 3.5
      end subroutine Sec_Gaussian_Chebyshev


      subroutine LOG_MK1996(indx, nradpt, parm, radr, rwgt)
         use parameter_def
         implicit none
         ! argument
         integer, intent(in) :: indx, nradpt
         real(double), intent(in) :: parm
         real(double), intent(out) :: radr, rwgt
         ! local variables
         real(double) :: radx
         radx = 0.0d0
         radr = 0.0d0
         rwgt = 1.0d0
         radx = real(indx) / real(nradpt + 1)
         radr = -parm * log(1 - radx ** 3)
         rwgt = 3 * radx ** 2 * log(1 - radx ** 3) ** 2
         rwgt = parm ** 3 * rwgt / ((nradpt + 1) * (1 - radx ** 3))
      end subroutine LOG_MK1996


      !subroutine MultiExp2003()
      !end subroutine MultiExp2003
