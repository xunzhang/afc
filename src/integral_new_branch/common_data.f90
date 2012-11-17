!------------------------------------------------------------------------------- 
!
!   Use "common_data" module for the mass data share and exchange. Most 
!   variables in the common blocks in DVM code are transfered into this 
!   "common_data" module for easy reference, and aviod the old feature 
!   "common" in Fortran77. Most variables are renamed.
!
!   This module should be "use" at the top of the calling programs, and 
!   all the "allocatable" variables must be "allocate" and "dellocate" 
!   in the main program, or it will cause troubles.
!
!    
!-------------------------------------------------------------------------------
!   
!   Dependency:
!
!     parameter_def
!
!-------------------------------------------------------------------------------

      module common_data

      use parameter_def

      implicit none
      
      type basis_type
         character(len=80) :: title
         integer           :: ngrid
         integer           :: nshell
         real(double)      :: radius
         real(double)      :: hinv         
         real(double)      :: zn
         real(double)      :: xion
         real(double)      :: rbar1
         real(double)      :: vbar
         real(double)      :: rbar2         
         real(double)      :: rnuc         
         integer           :: naov
         integer           :: naoc
         integer           :: lmax
         integer           :: nfrz
         integer           :: ngrid3d
         real(double)      :: alpha
         real(double)      :: r0
         real(double)      :: beta
         real(double)      :: gamma
         real(double)      :: w1tab  
         real(double), dimension(mxrpt)        :: rgrid  
         real(double), dimension(mxrpt)        :: vatom
         integer,      dimension(mxnls)        :: nls
         real(double), dimension(mxnls)        :: occ
         real(double), dimension(mxnls)        :: eig
         integer,      dimension(mxnls)        :: kgrid         
         integer,      dimension(mxnls)        :: l
         real(double), dimension(mxrpt, mxnls) :: rnl
         real(double), dimension(mxrpt)        :: rho
         real(double), dimension(mxrpt)        :: rhov
         real(double), dimension(mxrpt)        :: vcoul
         integer,      dimension(mxnls)        :: nlsfrz         
         logical,      dimension(mxnls)        :: ifcore
      end type basis_type      

      type system_type
         integer                              :: natom
         integer                              :: nbuff
         real(double)                         :: rbuff
         integer                              :: naov
         integer                              :: naoc         
         real(double)                         :: nelect
         integer                              :: atomtagindx
         integer                              :: bufftagindx                  
         integer                              :: aovlistindx
         integer                              :: aoclistindx
         integer                              :: grid3dindx 
         integer                              :: mtxindx        
         integer                              :: eigindx        
         integer                              :: ngrid3d
         integer                              :: nfiledivtot
         integer                              :: ngriddivtot
         integer                              :: nfiledivval
         integer                              :: ngriddivval
         character(len=32), dimension(mxfdiv) :: filedivval 
      end type system_type

      type aovlist_type
         integer                        :: tick 
         integer                        :: atom
         integer                        :: n
         integer                        :: l
         integer                        :: m        
      end type aovlist_type

      type aoclist_type
         integer                        :: atom
         integer                        :: n
         integer                        :: l
         integer                        :: m
      end type aoclist_type

      type molist_type
         real(double)                   :: eig
         real(double)                   :: occ
         real(double)                   :: pfphi2         
         integer                        :: system
         integer                        :: eorder
         character(len=1)               :: spin
      end type molist_type

      type grid3d_type
         real(double)                   :: x
         real(double)                   :: y
         real(double)                   :: z
         real(double)                   :: w 
         real(double)                   :: pf 
         real(double)                   :: wpf                                                   
      end type grid3d_type 

      type potgrid3d_type
         real(double)                   :: charx
         real(double)                   :: spnsx
         real(double)                   :: vcoux
         real(double)                   :: vxcup 
         real(double)                   :: vxcdw
         real(double)                   :: veffup
         real(double)                   :: veffdw
         real(double)                   :: eexc
         real(double)                   :: vne  
      end type potgrid3d_type 

      type plot_type
         integer                        :: dim
         integer                        :: gridindx
         integer                        :: ngrid
         integer,      dimension(3)     :: kgrid
         real(double), dimension(3, 4)  :: vector      
      end type plot_type 

      type plotgrid_type
         real(double)                   :: x
         real(double)                   :: y
         real(double)                   :: z
         real(double)                   :: charx
         real(double)                   :: spnsx
         real(double)                   :: vcoul
         real(double)                   :: vxcup
         real(double)                   :: vxcdw
         real(double)                   :: vefld         
         real(double)                   :: diff
      end type plotgrid_type 

      type smplmixsav_type
         real(double)                   :: bc0    
         real(double)                   :: bs0   
         real(double)                   :: dis0  
         real(double)                   :: dss0  
         real(double)                   :: dis2  
         real(double)                   :: dss2  
         real(double)                   :: dis   
         real(double)                   :: dss   
         real(double)                   :: cmix  
         real(double)                   :: smix  
         real(double)                   :: pmixc 
         real(double)                   :: pmixs 
         real(double)                   :: bc    
         real(double)                   :: bs    
         real(double)                   :: pminc 
         real(double)                   :: pmins 
         integer                        :: ncyg  
         integer                        :: ncys  
         logical                        :: extrapc 
         logical                        :: extraps
         real(double)                   :: theta
      end type smplmixsav_type
      
      type(basis_type),     dimension(:), allocatable :: basis
      type(system_type),    dimension(:), allocatable :: system 
      type(aovlist_type),   dimension(:), allocatable :: aovlist
      type(aoclist_type),   dimension(:), allocatable :: aoclist      
      type(molist_type),    dimension(:), allocatable :: molist
      type(grid3d_type),    dimension(:), allocatable :: grid3d
      type(potgrid3d_type), dimension(:), allocatable :: potgrid3d
      type(plot_type),      dimension(:), allocatable :: plot
      type(plotgrid_type),  dimension(:), allocatable :: plotgrid
      type(plotgrid_type),  dimension(:), allocatable :: rgaussgrid


      real(double),    dimension(:, :), allocatable :: coord
      integer,            dimension(:), allocatable :: basistag
      integer,            dimension(:), allocatable :: atomtag
      integer,            dimension(:), allocatable :: bufftag      
      real(double),    dimension(:, :), allocatable :: rho
      real(double),    dimension(:, :), allocatable :: spn      
      real(double),    dimension(:, :), allocatable :: rhov
      real(double),    dimension(:, :), allocatable :: vcoul
      real(double),       dimension(:), allocatable :: hmtxup
      real(double),       dimension(:), allocatable :: hmtxdw
      real(double),       dimension(:), allocatable :: tmtx      
      real(double),       dimension(:), allocatable :: omtx
      real(double),       dimension(:), allocatable :: pfomtx      
      real(double),       dimension(:), allocatable :: eigup
      real(double),       dimension(:), allocatable :: eigdw      
      real(double),       dimension(:), allocatable :: vecup
      real(double),       dimension(:), allocatable :: vecdw
      real(double),       dimension(:), allocatable :: pfphi2up
      real(double),       dimension(:), allocatable :: pfphi2dw  
      real(double),       dimension(:), allocatable :: occup
      real(double),       dimension(:), allocatable :: occdw            
      
      
      real(double)                                  :: eps
      integer                                       :: ncmdarg
      character(len=256)                            :: cmdlarg
      character(len=25)                             :: prefix
      real(double)                                  :: unitscale
      character(len=80)                             :: jobid      
      character(len=32)                             :: basisfile
      integer                                       :: nbasis
      integer                                       :: natom
      integer                                       :: nbuff
      integer                                       :: nsystem
      integer                                       :: functional
      integer                                       :: spin
      integer                                       :: frozen    
      real(double)                                  :: nelect
      real(double)                                  :: charge      
      real(double)                                  :: ef  
      real(double)                                  :: beta 
      integer                                       :: pftype
      integer                                       :: ipop
      integer                                       :: isolve
      integer                                       :: ifefld
      real(double), dimension(1:3)                  :: efld
      integer                                       :: naov
      integer                                       :: naoc
      integer                                       :: nmo
      integer                                       :: itscf
      integer                                       :: scfstep 
      integer                                       :: mixtype
      real(double)                                  :: occshift
      integer                                       :: ifscfdiis 
      integer                                       :: nscfdiis
      real(double)                                  :: scfdiistol
      real(double)                                  :: scfrhotol
      real(double)                                  :: scfspntol      
      real(double)                                  :: scfrhomix
      real(double)                                  :: scfspnmix
      integer, dimension(1:3)                       :: grid3dseed
      integer                                       :: ngrid3d
      integer                                       :: nbytediv
      integer                                       :: ifbsfn 
      integer                                       :: smplmixload    
      integer                                       :: ifetot          
      integer                                       :: ifbypscfc            
      integer                                       :: printbasis 
      integer                                       :: printscf
      integer                                       :: ifpopu
      integer                                       :: ifpdos
      integer                                       :: nplot
      integer                                       :: plotinfo      
      integer                                       :: ifdebug
      integer                                       :: mtxdim
      real(double)                                  :: tac
      real(double)                                  :: scfrhodis
      logical                                       :: scfconvg
      integer                                       :: nplotgrid
      integer                                       :: nrgaussgrid


! saved files
      character(len=32)                             :: filebfnbak                  
      character(len=32)                             :: filescfbak                  
      character(len=32)                             :: filemixbak  
      character(len=32)                             :: filepdos
      character(len=32)                             :: filepopu

      
! variational parameters for simple mixing, here are mulliken populations
      integer                                       :: nrhofcn
      integer,            dimension(:), allocatable :: rhofcnindx            
      real(double),       dimension(:), allocatable :: rhofcn
      real(double),       dimension(:), allocatable :: spnfcn
      integer                                       :: nsmplmix 
      type(smplmixsav_type)                         :: smplmixsav   
      integer,            dimension(:), allocatable :: smplmixindx    
      real(double),    dimension(: ,:), allocatable :: smplmixrho
      real(double),    dimension(: ,:), allocatable :: smplmixspn


! temporary variables
      character(len=10000) :: tmpstring
      integer :: mxaov
      integer :: mxaoc
      integer :: mxvco
      integer :: mxatm      
      real(double), dimension(:), allocatable :: basv
      real(double), dimension(:), allocatable :: tbasv
      real(double), dimension(:), allocatable :: basc
      real(double), dimension(:), allocatable :: tbasc
      real(double), dimension(:), allocatable :: bascnor      
      real(double), dimension(:), allocatable :: vcocoe      
      real,         dimension(:), allocatable :: basv_sngl
      real,         dimension(:), allocatable :: tbasv_sngl
      real,         dimension(:), allocatable :: basc_sngl
      real,         dimension(:), allocatable :: tbasc_sngl      


! exclusive variables for Lapack and eigen module
! within the module, a(,) b(,) z(,) w(), work() are shared for LUDCP 
      integer :: lpk_n      
      integer :: lpk_lda    
      integer :: lpk_ldb    
      integer :: lpk_ldz    
      integer :: lpk_liwork 
      integer :: lpk_lwork  
      integer,      dimension(:),    allocatable :: lpk_iwork 
      integer,      dimension(:),    allocatable :: lpk_ifail      
      real(double), dimension(:, :), allocatable :: lpk_a
      real(double), dimension(:, :), allocatable :: lpk_b  
      real(double), dimension(:, :), allocatable :: lpk_z
      real(double), dimension(:),    allocatable :: lpk_w
      real(double), dimension(:),    allocatable :: lpk_work
      integer :: lpk_m   
      integer :: lpk_itype
      integer :: lpk_il
      integer :: lpk_iu
      integer :: lpk_info
      real(double) :: lpk_vl
      real(double) :: lpk_vu
      real(double) :: lpk_abstol
      character (len=1) :: lpk_jobz
      character (len=1) :: lpk_uplo
      character (len=1) :: lpk_range   

! exclusive variables for Mulliken population and PDOS analysis
      real(double), dimension(:), allocatable :: popu
      real(double), dimension(:), allocatable :: bondorder
      real(double), dimension(:), allocatable :: nccs_ij
      real(double), dimension(:), allocatable :: nccs_ii
      integer,      dimension(:), allocatable :: iatomtick

! exclusive variables for Broyden variational process
      real(double), dimension(:),    allocatable :: bryd_f
      real(double), dimension(:),    allocatable :: bryd_dumvi
      real(double), dimension(:),    allocatable :: bryd_ui      
      real(double), dimension(:),    allocatable :: bryd_vti
      real(double), dimension(:),    allocatable :: bryd_t1
      real(double), dimension(:),    allocatable :: bryd_df            
      real(double), dimension(:, :), allocatable :: bryd_vector
      real(double), dimension(:),    allocatable :: bryd_f_old
      real(double), dimension(:),    allocatable :: bryd_dumvi_old
      real(double), dimension(:),    allocatable :: bryd_ui_old
      real(double), dimension(:),    allocatable :: bryd_vti_old
      real(double), dimension(:, :), allocatable :: bryd_a
      real(double), dimension(:, :), allocatable :: bryd_b
      real(double), dimension(:, :), allocatable :: bryd_d
      real(double), dimension(:),    allocatable :: bryd_cm
      real(double), dimension(:),    allocatable :: bryd_w
      real(double)                               :: bryd_amix_old
      integer                                    :: bryd_index32
      integer                                    :: bryd_lastit_old
      
! exclusive variables for total energy calculation
      real(double) :: et_tsu     
      real(double) :: et_tmol    
      real(double) :: et_tcr     
      real(double) :: et_sc      
      real(double) :: et_ce1     
      real(double) :: et_ce2     
      real(double) :: et_ce3     
      real(double) :: et_emxc    
      real(double) :: et_etxc    
      real(double) :: et_etot    
      real(double) :: et_eb     


      
      end module common_data        

!------------------------------------------------------------------------------- 


