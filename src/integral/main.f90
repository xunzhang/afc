!-------------------------------------------------------------------------------  
!                                                                        
!   This is a Divide and Conquer (DAC) method for solving Kohn-Sham     
!   equations within Local Density Functional Theory (LDFT) scheme.     
!   Discrete Variational Method (DVM), a handy real-space Linear        
!   Combination of Atomic Orbitals for Molecular Orbitals (LCAO-MO)     
!   method is introduced for solving the subsystem problem in DAC.      
!                                                                        
!   Presently, only a simple way called Self Consistent Charge (SCC)    
!   is available for realizing the whole variational process, which     
!   is time saving but precision lost. For more details please refer    
!   to references on DVM.                                               
!                                                                        
!   Group Theory is not applied for Hamiltonian diagonalization of      
!   subsystem.                                                          
!                                                                        
!                                                                        
!   Prof. Chong-Yu Wang, 
!   Department of Physics, Tsinghua University, Beijing 100084, P. R. China 
!   
!   Dr. Shan-Ying Wang
!   Department of Physics, Tsinghua University, Beijing 100084, P. R. China 
!                                                         
!   26 May 2005                                                        
!                          
!                                              
!-------------------------------------------------------------------------------  
!   
!   Dependency:
!
!     parameter_def, common_data, simplemix, print_error, maceps, 
!     get_cmdlarg, print_hello, get_input, get_basis, initialization,
!     get_grid3d, get_pfunction, print_initial, print_timer, simplemix, 
!     eneprp, bsfngn0, bsfngn1, bsfngn2, modelrhovp, potgn2, hmtxgn, 
!     secula, get_pfphi2, sorteig, findef, popuana, updtscfc, potgn3, 
!     energy, scfctrl, print_scfc, popdump, pdosgn, plotgn 
!
!-------------------------------------------------------------------------------
      
      program main

      use parameter_def
      use common_data
      use simplemix
      use popuana
      use pdosgn      
      use plotgn
      
      implicit none

! local variables 
      integer :: inunit
      integer :: otunit
      character(len=32) :: inmast
      character(len=32) :: incoor
      character(len=32) :: otmast      

! local temporary variables
      integer :: status0 
      integer :: is      
      integer :: ii, jj, kk, ll
      real(double) :: int_result 

      call maceps (eps)
      
      call get_cmdlarg (ncmdarg, cmdlarg, prefix, status0)
      if (status0 /= 0) then
         call print_error (" an illegal prefix (not fortran name) for the job")
         stop 
      end if   

      inunit = 15
      otunit = 16
            
      inmast = trim(prefix) // ".input"
      incoor = trim(prefix) // ".incar"
      otmast = trim(prefix) // ".otput"

      open (unit = otunit, file = otmast, status = "unknown", iostat = status0)
      if (status0 /= 0) then
         call print_error (" main: open " // trim(otmast) // " file error")
         stop
      end if  

      call print_hello (otunit)
      call print_timer ("         ",      6, 0)      
      call print_timer ("Job begin",      6, 3)
      call print_timer ("Job begin", otunit, 3)


!===> begin scf initialization 

      call get_input (inunit, inmast, 1)
      
      if (.not. allocated(basis)) allocate(basis(nbasis))       

      call get_basis (inunit, 1)
      
      if (nplot > 0) then
         if (.not. allocated(plot)) allocate(plot(nplot))       
      end if   
              
      call get_input (inunit, inmast, 2)
      call get_basis (inunit, 2)
      
      if (.not. allocated(system))   allocate(system(nsystem))        
      if (.not. allocated(coord))    allocate(coord(3, natom))        
      if (.not. allocated(basistag)) allocate(basistag(natom))        
      if (.not. allocated(atomtag))  allocate(atomtag(natom))        
           
      call get_input (inunit, incoor, 3)
      call initialization (1)

      if (.not. allocated(bufftag))     allocate(bufftag(nbuff))        
      if (.not. allocated(rhofcnindx))  allocate(rhofcnindx(natom))        
      if (.not. allocated(smplmixindx)) allocate(smplmixindx(natom))        

      call initialization (2)
      
      if (.not. allocated(aovlist))    allocate(aovlist(naov)) 
      if (.not. allocated(aoclist))    allocate(aoclist(naoc))    
      if (.not. allocated(molist))     allocate(molist(nmo))
      if (.not. allocated(grid3d))     allocate(grid3d(ngrid3d))
      if (.not. allocated(potgrid3d))  allocate(potgrid3d(ngrid3d))      
      if (.not. allocated(rho))        allocate(rho(mxrpt, natom))
      if (.not. allocated(rhov))       allocate(rhov(mxrpt, natom))
      if (.not. allocated(vcoul))      allocate(vcoul(mxrpt, natom))         
      if (.not. allocated(hmtxup))     allocate(hmtxup(mtxdim))
      if (.not. allocated(tmtx))       allocate(tmtx(mtxdim))
      if (.not. allocated(omtx))       allocate(omtx(mtxdim))
      if (.not. allocated(pfomtx))     allocate(pfomtx(mtxdim))      
      if (.not. allocated(eigup))      allocate(eigup(naov))
      if (.not. allocated(vecup))      allocate(vecup(mtxdim))
      if (.not. allocated(pfphi2up))   allocate(pfphi2up(naov))
      if (.not. allocated(occup))      allocate(occup(naov))      
      if (.not. allocated(rhofcn))     allocate(rhofcn(nrhofcn))
      if (.not. allocated(smplmixrho)) allocate(smplmixrho(nsmplmix, 4)) 
      if (.not. allocated(basv))       allocate(basv(mxaov)) 
      if (.not. allocated(tbasv))      allocate(tbasv(mxaov)) 
      if (.not. allocated(basc))       allocate(basc(mxaoc)) 
      if (.not. allocated(tbasc))      allocate(tbasc(mxaoc)) 
      if (.not. allocated(bascnor))    allocate(bascnor(mxaoc)) 
      if (.not. allocated(vcocoe))     allocate(vcocoe(mxvco)) 
      if (.not. allocated(basv_sngl))  allocate(basv_sngl(mxaov)) 
      if (.not. allocated(tbasv_sngl)) allocate(tbasv_sngl(mxaov)) 
      if (.not. allocated(basc_sngl))  allocate(basc_sngl(mxaoc)) 
      if (.not. allocated(tbasc_sngl)) allocate(tbasc_sngl(mxaoc)) 

      if (spin /= 0) then
         if (.not. allocated(spn))        allocate(spn(mxrpt, natom))      
         if (.not. allocated(hmtxdw))     allocate(hmtxdw(mtxdim))
         if (.not. allocated(eigdw))      allocate(eigdw(naov)) 
         if (.not. allocated(vecdw))      allocate(vecdw(mtxdim))
         if (.not. allocated(pfphi2dw))   allocate(pfphi2dw(naov))
         if (.not. allocated(occdw))      allocate(occdw(naov)) 
         if (.not. allocated(spnfcn))     allocate(spnfcn(nrhofcn))
         if (.not. allocated(smplmixspn)) allocate(smplmixspn(nsmplmix, 4)) 
      end if

      if (.not. allocated(lpk_iwork))  allocate(lpk_iwork(lpk_liwork))
      if (.not. allocated(lpk_ifail))  allocate(lpk_ifail(lpk_n))
      if (.not. allocated(lpk_a))      allocate(lpk_a(lpk_lda,lpk_lda))
      if (.not. allocated(lpk_b))      allocate(lpk_b(lpk_ldb,lpk_ldb))
      if (.not. allocated(lpk_z))      allocate(lpk_z(lpk_ldz,lpk_ldz))
      if (.not. allocated(lpk_w))      allocate(lpk_w(lpk_n))
      if (.not. allocated(lpk_work))   allocate(lpk_work(lpk_lwork))      
      
      if (.not. allocated(popu))       allocate(popu(mxaov))
      if (.not. allocated(bondorder))  allocate(bondorder(mxatm*(mxatm+1)/2))
      if (.not. allocated(nccs_ij))    allocate(nccs_ij(mxaov*(mxaov+1)/2))
      if (.not. allocated(nccs_ii))    allocate(nccs_ii(mxaov))      
      if (.not. allocated(iatomtick))  allocate(iatomtick(mxaov))      
      
      if (ifscfdiis /= 0 ) then
         if (spin == 0) then
            if (.not. allocated(bryd_f))         allocate(bryd_f(nsmplmix))
            if (.not. allocated(bryd_dumvi))     allocate(bryd_dumvi(nsmplmix))
            if (.not. allocated(bryd_ui))        allocate(bryd_ui(nsmplmix))
            if (.not. allocated(bryd_vti))       allocate(bryd_vti(nsmplmix))
            if (.not. allocated(bryd_t1))        allocate(bryd_t1(nsmplmix))
            if (.not. allocated(bryd_df))        allocate(bryd_df(nsmplmix))
            if (.not. allocated(bryd_vector))    allocate(bryd_vector(nsmplmix, 2))
            if (.not. allocated(bryd_f_old))     allocate(bryd_f_old(nsmplmix))
            if (.not. allocated(bryd_dumvi_old)) allocate(bryd_dumvi_old(nsmplmix))
            if (.not. allocated(bryd_ui_old))    allocate(bryd_ui_old(nsmplmix*nscfdiis))
            if (.not. allocated(bryd_vti_old))   allocate(bryd_vti_old(nsmplmix*nscfdiis))
            if (.not. allocated(bryd_a))         allocate(bryd_a(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_b))         allocate(bryd_b(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_d))         allocate(bryd_d(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_cm))        allocate(bryd_cm(nscfdiis))
            if (.not. allocated(bryd_w))         allocate(bryd_w(nscfdiis))              
         else
            if (.not. allocated(bryd_f))         allocate(bryd_f(nsmplmix*2))
            if (.not. allocated(bryd_dumvi))     allocate(bryd_dumvi(nsmplmix*2))
            if (.not. allocated(bryd_ui))        allocate(bryd_ui(nsmplmix*2))
            if (.not. allocated(bryd_vti))       allocate(bryd_vti(nsmplmix*2))
            if (.not. allocated(bryd_t1))        allocate(bryd_t1(nsmplmix*2))
            if (.not. allocated(bryd_df))        allocate(bryd_df(nsmplmix*2))
            if (.not. allocated(bryd_vector))    allocate(bryd_vector(nsmplmix*2, 2))
            if (.not. allocated(bryd_f_old))     allocate(bryd_f_old(nsmplmix*2))
            if (.not. allocated(bryd_dumvi_old)) allocate(bryd_dumvi_old(nsmplmix*2))
            if (.not. allocated(bryd_ui_old))    allocate(bryd_ui_old(nsmplmix*2*nscfdiis))
            if (.not. allocated(bryd_vti_old))   allocate(bryd_vti_old(nsmplmix*2*nscfdiis))
            if (.not. allocated(bryd_a))         allocate(bryd_a(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_b))         allocate(bryd_b(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_d))         allocate(bryd_d(nscfdiis, nscfdiis))
            if (.not. allocated(bryd_cm))        allocate(bryd_cm(nscfdiis))
            if (.not. allocated(bryd_w))         allocate(bryd_w(nscfdiis))            
         end if
      end if    

      if (nplot > 0) then
         if (.not. allocated(plotgrid)) allocate(plotgrid(nplotgrid)) 
      end if   
      call print_timer ("initialization start", otunit, 3)
      call initialization (3)      
      call print_timer ("initialization End", otunit, 3)
      call get_grid3d
      call print_timer ("get_grid3d End", otunit, 3)
      call get_pfunction
      call print_timer ("get_pfunction End", otunit, 3)
      call print_initial (otunit)
      call print_timer ("print_initial End", otunit, 3)
      
      if (ifdebug > 0) then
         write(     6, "(' main: ifdebug make a stop here for test purpose')")
         write(otunit, "(' main: ifdebug make a stop here for test purpose')")
         stop
      end if   


!===> load variational parameters

      call smplmix0 
      call print_timer ("smplmix0 End", otunit, 3)


!===> enter SCFC calculation

      if (ifbypscfc == 0) then

         if (ifetot > 0) call eneprp
         
         select case (ifbsfn)
         case (0)
            do is = 1, nsystem
               call bsfngn0 (is, otunit)
            end do
      call print_timer ("bsfngn0 End", otunit, 3)
         case (1) 
            do is = 1, nsystem
               call bsfngn1 (is)
            end do
      call print_timer ("bsfngn1 End", otunit, 3)
         case default
            do is = 1, nsystem
               call bsfngn2 (is)
            end do
      call print_timer ("bsfngn2 End", otunit, 3)
         end select

         itscf = 1
         do while (itscf <= scfstep)
            call modelrhovp
      call print_timer ("modelrhovp End", otunit, 3)
            do is = 1, nsystem
      call print_timer ("===========", otunit, 3)
               call potgn2 (is)
      call print_timer ("potgn2 End", otunit, 3)
               call hmtxgn (is)
      call print_timer ("hmtxgn End", otunit, 3)
               call secula (is)
      call print_timer ("secula End", otunit, 3)
               call get_pfphi2 (is)
      call print_timer ("get_pfphi2 End", otunit, 3)
      call print_timer ("===========", otunit, 3)
            end do !end do is
            call sorteig
      call print_timer ("sorteig End", otunit, 3)
            call findef
      call print_timer ("findef End", otunit, 3)
            call getpopu
      call print_timer ("getpopu End", otunit, 3)
            call updtscfc
      call print_timer ("updtscfc  End", otunit, 3)
            if (ifetot > 1) then
               call potgn3
               call energy
            elseif (ifetot == 1) then
               if ((itscf+2) >= scfstep) then
                  call potgn3
                  call energy
               end if   
            end if
            call scfctrl
            call print_scfc (     6)
            call print_scfc (otunit)
         end do 
      call print_timer ("scfstep while-loop  End", otunit, 3)
      end if
      
!===> population, partial DOS and isosurface or contour plot data
 
      write(*, *) coord
      !if (ifpopu > 0) call popdump 
      !if (ifpdos > 0) call getpdos 
      !if ((nplot > 0) .and. (plotinfo > 0) .and. (plotinfo < 1000)) call getplot
      print *, "the intintint result before is"
      call integral3d (natom, 302, 80, int_result)
      print *, "the intintint result is"
      write(*, *) int_result
!===> deallocate all dynamic arrays even this seems unnecessary

      if (allocated(basis))       deallocate(basis)       
      if (allocated(plot))        deallocate(plot)       
      if (allocated(system))      deallocate(system)        
      if (allocated(coord))       deallocate(coord)        
      if (allocated(basistag))    deallocate(basistag)        
      if (allocated(atomtag))     deallocate(atomtag)        
      if (allocated(bufftag))     deallocate(bufftag)        
      if (allocated(rhofcnindx))  deallocate(rhofcnindx)        
      if (allocated(smplmixindx)) deallocate(smplmixindx)        
      if (allocated(aovlist))     deallocate(aovlist) 
      if (allocated(aoclist))     deallocate(aoclist)    
      if (allocated(molist))      deallocate(molist)
      if (allocated(grid3d))      deallocate(grid3d)
      if (allocated(potgrid3d))   deallocate(potgrid3d)      
      if (allocated(rho))         deallocate(rho)
      if (allocated(rhov))        deallocate(rhov)
      if (allocated(vcoul))       deallocate(vcoul)         
      if (allocated(hmtxup))      deallocate(hmtxup)
      if (allocated(tmtx))        deallocate(tmtx)
      if (allocated(omtx))        deallocate(omtx)
      if (allocated(pfomtx))      deallocate(pfomtx)      
      if (allocated(eigup))       deallocate(eigup)
      if (allocated(vecup))       deallocate(vecup)
      if (allocated(pfphi2up))    deallocate(pfphi2up)
      if (allocated(occup))       deallocate(occup)      
      if (allocated(rhofcn))      deallocate(rhofcn)
      if (allocated(smplmixrho))  deallocate(smplmixrho) 
      if (allocated(basv))        deallocate(basv) 
      if (allocated(tbasv))       deallocate(tbasv) 
      if (allocated(basc))        deallocate(basc) 
      if (allocated(tbasc))       deallocate(tbasc) 
      if (allocated(bascnor))     deallocate(bascnor) 
      if (allocated(vcocoe))      deallocate(vcocoe) 
      if (allocated(basv_sngl))   deallocate(basv_sngl) 
      if (allocated(tbasv_sngl))  deallocate(tbasv_sngl) 
      if (allocated(basc_sngl))   deallocate(basc_sngl) 
      if (allocated(tbasc_sngl))  deallocate(tbasc_sngl) 
      if (allocated(spn))         deallocate(spn)      
      if (allocated(hmtxdw))      deallocate(hmtxdw)
      if (allocated(eigdw))       deallocate(eigdw) 
      if (allocated(vecdw))       deallocate(vecdw)
      if (allocated(pfphi2dw))    deallocate(pfphi2dw)
      if (allocated(occdw))       deallocate(occdw) 
      if (allocated(spnfcn))      deallocate(spnfcn)
      if (allocated(smplmixspn))  deallocate(smplmixspn) 
      if (allocated(lpk_iwork))   deallocate(lpk_iwork)
      if (allocated(lpk_ifail))   deallocate(lpk_ifail)
      if (allocated(lpk_a))       deallocate(lpk_a)
      if (allocated(lpk_b))       deallocate(lpk_b)
      if (allocated(lpk_z))       deallocate(lpk_z)
      if (allocated(lpk_w))       deallocate(lpk_w)
      if (allocated(lpk_work))    deallocate(lpk_work)      
      if (allocated(popu))        deallocate(popu)
      if (allocated(bondorder))   deallocate(bondorder)
      if (allocated(nccs_ij))     deallocate(nccs_ij)
      if (allocated(nccs_ii))     deallocate(nccs_ii)      
      if (allocated(iatomtick))   deallocate(iatomtick)      
      if (allocated(bryd_f))         deallocate(bryd_f)
      if (allocated(bryd_dumvi))     deallocate(bryd_dumvi)
      if (allocated(bryd_ui))        deallocate(bryd_ui)
      if (allocated(bryd_vti))       deallocate(bryd_vti)
      if (allocated(bryd_t1))        deallocate(bryd_t1)
      if (allocated(bryd_df))        deallocate(bryd_df)
      if (allocated(bryd_vector))    deallocate(bryd_vector)
      if (allocated(bryd_f_old))     deallocate(bryd_f_old)
      if (allocated(bryd_dumvi_old)) deallocate(bryd_dumvi_old)
      if (allocated(bryd_ui_old))    deallocate(bryd_ui_old)
      if (allocated(bryd_vti_old))   deallocate(bryd_vti_old)
      if (allocated(bryd_a))         deallocate(bryd_a)
      if (allocated(bryd_b))         deallocate(bryd_b)
      if (allocated(bryd_d))         deallocate(bryd_d)
      if (allocated(bryd_cm))        deallocate(bryd_cm)
      if (allocated(bryd_w))         deallocate(bryd_w)              
      if (allocated(plotgrid))       deallocate(plotgrid) 

      call print_timer ("Job end",      6, 3)
      call print_timer ("Job end", otunit, 3)      
      
      
      end program main

!-------------------------------------------------------------------------------



