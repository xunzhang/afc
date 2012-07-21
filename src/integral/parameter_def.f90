!-------------------------------------------------------------------------------  
!
!   Parameter definition. The "KIND" parameters are refered to the 
!   famous book "Numerical Recipes".
!
!      integer, parameter :: i4b = selected_int_kind(9)
!      integer, parameter :: i2b = selected_int_kind(4)
!      integer, parameter :: i1b = selected_int_kind(2)
!      integer, parameter :: sp  = kind(1.0)
!      integer, parameter :: dp  = kind(1.0d0)
!      integer, parameter :: spc = kind((1.0,1.0))
!      integer, parameter :: dpc = kind((1.0d0,1.0d0))
!      integer, parameter :: lgt = kind(.true.)
!
!------------------------------------------------------------------------------- 
       
      module parameter_def

      implicit none

! kind parameter of double precision 
      integer, parameter :: double = kind(1.0d0)

      real(double), parameter :: pi = 3.14159265358979d0       
      real(double), parameter :: au_ev = 27.21138344d0       
      real(double), parameter :: au_ai = 0.5291772083d0 

! the minimum distance between a point and a nuclear     
      real(double), parameter :: minrnu12 = 1.0d-12 
      real(double), parameter :: minrnu15 = 1.0d-15         

! the minimum value of atomic distance (in au)     
      real(double), parameter :: dguard = 1.0d0   
            
! maximum number of atomic radial wave function (Rnl) in each basis 
      integer, parameter :: mxnls = 8
      
! maximum number of radial points (logarithmic distribution) 
      integer, parameter :: mxrpt = 300

! maximum l value 
      integer, parameter :: mxlva = 6

! about wave function file split
      integer, parameter :: iostar = 300
      integer, parameter :: mxfdiv = 100



      end module parameter_def

!-------------------------------------------------------------------------------        