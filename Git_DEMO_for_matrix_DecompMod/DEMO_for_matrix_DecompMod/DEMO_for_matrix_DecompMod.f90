!  DEMO_for_matrix_DecompMod.f90 
!
!  FUNCTIONS:
!  DEMO_for_matrix_DecompMod - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: DEMO_for_matrix_DecompMod
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program DEMO_for_matrix_DecompMod


    implicit none

    ! Variables
    ! column level
   real t_soisno(10)    ! soil temperature (Kelvin)  (-nlevsno+110nlevgrnd)
   real psisat(10)      ! soil water potential at saturation for CN code (MPa)
   real soilpsi(10)     ! soil water potential in each soil layer (MPa)
   real dz(10)          ! soil layer thickness (m)
   real cwdc(10)          ! (gC/m2) coarse woody debris C
   real litr1c(10)        ! (gC/m2) litter labile C
   real litr2c(10)        ! (gC/m2) litter cellulose C
   real litr3c(10)        ! (gC/m2) litter lignin C
   real soil1c(10)        ! (gC/m2) soil organic matter C (fast pool)
   real soil2c(10)        ! (gC/m2) soil organic matter C (medium pool)
   real soil3c(10)        ! (gC/m2) soil organic matter C (slow pool)
   real soil4c(10)        ! (gC/m2) soil organic matter C (slowest pool)
   real cwdn(10)          ! (gN/m2) coarse woody debris N
   real litr1n(10)        ! (gN/m2) litter labile N
   real litr2n(10)        ! (gN/m2) litter cellulose N
   real litr3n(10)        ! (gN/m2) litter lignin N
   integer clandunit(10)      ! index into landunit level quantities
   integer itypelun(10)      ! landunit type
   ! pft level
   real rootfr(10)      ! fraction of roots in each soil layer  (nlevgrnd)

   real  fpi(10)           ! fraction of potential immobilization (no units)
   real  cwdc_to_litr2c(10)
   real  cwdc_to_litr3c(10)
   real  litr1_hr(10)
   real  litr1c_to_soil1c(10)
   real  litr2_hr(10)
   real  litr2c_to_soil2c(10)
   real  litr3_hr(10)
   real  litr3c_to_soil3c(10)
   real  soil1_hr(10)
   real  soil1c_to_soil2c(10)
   real  soil2_hr(10)
   real  soil2c_to_soil3c(10)
   real  soil3_hr(10)
   real  soil3c_to_soil4c(10)
   real  soil4_hr(10)
   real  cwdn_to_litr2n(10)
   real  cwdn_to_litr3n(10)
   real  potential_immob(10)
   real  litr1n_to_soil1n(10)
   real  sminn_to_soil1n_l1(10)
   real  litr2n_to_soil2n(10)
   real  sminn_to_soil2n_l2(10)
   real  litr3n_to_soil3n(10)
   real  sminn_to_soil3n_l3(10)
   real  soil1n_to_soil2n(10)
   real  sminn_to_soil2n_s1(10)
   real  soil2n_to_soil3n(10)
   real  sminn_to_soil3n_s2(10)
   real  soil3n_to_soil4n(10)
   real  sminn_to_soil4n_s3(10)
   real  soil4n_to_sminn(10)
   real  sminn_to_denit_l1s1(10)
   real  sminn_to_denit_l2s2(10)
   real  sminn_to_denit_l3s3(10)
   real  sminn_to_denit_s1s2(10)
   real  sminn_to_denit_s2s3(10)
   real  sminn_to_denit_s3s4(10)
   real  sminn_to_denit_s4(10)
   real  sminn_to_denit_excess(10)
   real  gross_nmin(10)
   real  net_nmin(10)

!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer c,j          !indices
   integer fc           !lake filter column index
   real dt           !decomp timestep (seconds)
   real dtd          !decomp timestep (days)
   real fr(10)   !column-level rooting fraction by soil depth
   real frw(10)          !rooting fraction weight
   real t_scalar(10)     !soil temperature scalar for decomp
   real minpsi, maxpsi        !limits for soil water scalar for decomp
   real psi                   !temporary soilpsi for water scalar
   real w_scalar(10)     !soil water scalar for decomp
   real rate_scalar  !combined rate scalar for decomp
   real cn_l1(10)        !C:N for litter 1
   real cn_l2(10)        !C:N for litter 2
   real cn_l3(10)        !C:N for litter 3
   real cn_s1        !C:N for SOM 1
   real cn_s2        !C:N for SOM 2
   real cn_s3        !C:N for SOM 3
   real cn_s4        !C:N for SOM 4
   real rf_l1s1      !respiration fraction litter 1 -> SOM 1
   real rf_l2s2      !respiration fraction litter 2 -> SOM 2
   real rf_l3s3      !respiration fraction litter 3 -> SOM 3
   real rf_s1s2      !respiration fraction SOM 1 -> SOM 2
   real rf_s2s3      !respiration fraction SOM 2 -> SOM 3
   real rf_s3s4      !respiration fraction SOM 3 -> SOM 4
   real k_l1         !decomposition rate constant litter 1
   real k_l2         !decomposition rate constant litter 2
   real k_l3         !decomposition rate constant litter 3
   real k_s1         !decomposition rate constant SOM 1
   real k_s2         !decomposition rate constant SOM 2
   real k_s3         !decomposition rate constant SOM 3
   real k_s4         !decomposition rate constant SOM 3
   real k_frag       !fragmentation rate constant CWD
   real ck_l1        !corrected decomposition rate constant litter 1
   real ck_l2        !corrected decomposition rate constant litter 2
   real ck_l3        !corrected decomposition rate constant litter 3
   real ck_s1        !corrected decomposition rate constant SOM 1
   real ck_s2        !corrected decomposition rate constant SOM 2
   real ck_s3        !corrected decomposition rate constant SOM 3
   real ck_s4        !corrected decomposition rate constant SOM 3
   real ck_frag      !corrected fragmentation rate constant CWD
   real cwd_fcel     !cellulose fraction of coarse woody debris
   real cwd_flig     !lignin fraction of coarse woody debris
   real cwdc_loss    !fragmentation rate for CWD carbon (gC/m2/s)
   real cwdn_loss    !fragmentation rate for CWD nitrogen (gN/m2/s)
   real plitr1c_loss(10) !potential C loss from litter 1
   real plitr2c_loss(10) !potential C loss from litter 2
   real plitr3c_loss(10) !potential C loss from litter 3
   real psoil1c_loss(10) !potential C loss from SOM 1
   real psoil2c_loss(10) !potential C loss from SOM 2
   real psoil3c_loss(10) !potential C loss from SOM 3
   real psoil4c_loss(10) !potential C loss from SOM 4
   real pmnf_l1s1(10)    !potential mineral N flux, litter 1 -> SOM 1
   real pmnf_l2s2(10)    !potential mineral N flux, litter 2 -> SOM 2
   real pmnf_l3s3(10)    !potential mineral N flux, litter 3 -> SOM 3
   real pmnf_s1s2(10)    !potential mineral N flux, SOM 1 -> SOM 2
   real pmnf_s2s3(10)    !potential mineral N flux, SOM 2 -> SOM 3
   real pmnf_s3s4(10)    !potential mineral N flux, SOM 3 -> SOM 4
   real pmnf_s4(10)      !potential mineral N flux, SOM 4
   real immob(10)        !potential N immobilization
   real ratio        !temporary variable
   real dnp          !denitrification proportion
   integer nlevdecomp ! bottom layer to consider for decomp controls
   real spinup_scalar         !multiplier for AD_SPINUP algorithm
   
   !===============================matrix variables========================================
    
   real,dimension(8,8)       ::a_matrix
   real,dimension(8,8)       ::kk_matrix
   real,dimension(8,8)       ::na_matrix

   real,dimension(8,8)       ::nk_matrix
   real,dimension(8,1)       ::c_matrix_current_1
   real,dimension(8,1)       ::c_matrix_next_1
   real,dimension(8,1)       ::c_matrix_in_1
   real cn_cwd,rate1,rate2,rate3
   
   open(55,file="C_outputs.csv")
   write(55,'(a500)')",runtime,cwdc,litr1c, litr2c,litr3c,soil1c,soil2c,soil3c,soil4c,&
                     &cwdc_to_litr2c,cwdc_to_litr3c,litr1_hr,litr1c_to_soil1c,&
                     &litr2_hr,litr2c_to_soil2c,litr3_hr,litr3c_to_soil3c,&
                     &soil1_hr,soil1c_to_soil2c,soil2_hr,soil2c_to_soil3c,&
                     &soil3_hr,soil3c_to_soil4c,soil4_hr"

   !=================================================chen revision 1 设置输入值=====================================================
    cwdc(1) = 0.          ! (gC/m2) coarse woody debris C
    litr1c(1) = 50.        ! (gC/m2) litter labile C
    litr2c(1) = 30.        ! (gC/m2) litter cellulose C
    litr3c(1) = 20.       ! (gC/m2) litter lignin C
    soil1c(1) = 500.       ! (gC/m2) soil organic matter C (fast pool)
    soil2c(1) = 1000.       ! (gC/m2) soil organic matter C (medium pool)
    soil3c(1) = 2000.       ! (gC/m2) soil organic matter C (slow pool)
    soil4c(1) = 4000.       ! (gC/m2) soil organic matter C (slowest pool)
    cwdn(1)   = 30.       ! (gN/m2) coarse woody debris N
    litr1n(1) = 2.       ! (gN/m2) litter labile N
    litr2n(1) = 2.       ! (gN/m2) litter cellulose N
    litr3n(1) = 2.       ! (gN/m2) litter lignin N
   !================================================================================================================================  
    dt = 1.0
   dtd = 1.0    !设置测试时间步长
      
    ! Body of test1
      ! set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0)   !chen:设置土壤碳氮比
   cn_s1 = 12.0
   cn_s2 = 12.0
   cn_s3 = 10.0
   cn_s4 = 10.0

   ! set respiration fractions for fluxes between compartments                !chen:设置库之间碳交换中呼吸的比例 1-CUE
   ! (from Biome-BGC v4.2.0)
   rf_l1s1 = 0.39
   rf_l2s2 = 0.55
   rf_l3s3 = 0.29
   rf_s1s2 = 0.28
   rf_s2s3 = 0.46
   rf_s3s4 = 0.55

   ! set the cellulose and lignin fractions for coarse woody debris          !chen：设置凋落物中木质素和纤维素的比例
   cwd_fcel = 0.76
   cwd_flig = 0.24

   ! set initial base rates for decomposition mass loss (1/day)
   ! (from Biome-BGC v4.2.0, using three SOM pools)
   ! Value inside log function is the discrete-time values for a
   ! daily time step model, and the result of the log function is
   ! the corresponding continuous-time decay rate (1/day), following
   ! Olson, 1963.                                                             !chen:设置初始的基准降解率
   k_l1 = -log(1.0-0.7)                                              
   k_l2 = -log(1.0-0.07)
   k_l3 = -log(1.0-0.014)
   k_s1 = -log(1.0-0.07)
   k_s2 = -log(1.0-0.014)
   k_s3 = -log(1.0-0.0014)
   k_s4 = -log(1.0-0.0001)
   k_frag = -log(1.0-0.001)

   ! calculate the new discrete-time decay rate for model timestep           !chen:根据模型计算秒时间尺度的降解率    
   k_l1 = 1.0-exp(-k_l1*dtd)
   k_l2 = 1.0-exp(-k_l2*dtd)
   k_l3 = 1.0-exp(-k_l3*dtd)
   k_s1 = 1.0-exp(-k_s1*dtd)
   k_s2 = 1.0-exp(-k_s2*dtd)
   k_s3 = 1.0-exp(-k_s3*dtd)
   k_s4 = 1.0-exp(-k_s4*dtd)
   k_frag = 1.0-exp(-k_frag*dtd)
   
   ! The following code implements the acceleration part of the AD spinup
   ! algorithm, by multiplying all of the SOM decomposition base rates by 10.0.

!   if (use_ad_spinup) then
!      spinup_scalar = 20.
!      k_s1 = k_s1 * spinup_scalar
!      k_s2 = k_s2 * spinup_scalar
!      k_s3 = k_s3 * spinup_scalar
!      k_s4 = k_s4 * spinup_scalar
!   end if
    
      ! set initial values for potential C and N fluxes                            !chen:设置初始0值
   plitr1c_loss(:) = 0.
   plitr2c_loss(:) = 0.
   plitr3c_loss(:) = 0.
   psoil1c_loss(:) = 0.
   psoil2c_loss(:) = 0.
   psoil3c_loss(:) = 0.
   psoil4c_loss(:) = 0.
   pmnf_l1s1(:) = 0.
   pmnf_l2s2(:) = 0.
   pmnf_l3s3(:) = 0.
   pmnf_s1s2(:) = 0.
   pmnf_s2s3(:) = 0.
   pmnf_s3s4(:) = 0.
   pmnf_s4(:) = 0.

   ! column loop to calculate potential decomp rates and total immobilization
   ! demand.
 
   c=1
      
      ! calculate litter compartment C:N ratios                                        !chen:计算凋落物CN比
      if (litr1n(c) > 0.) cn_l1(c) = litr1c(c)/litr1n(c)
      if (litr2n(c) > 0.) cn_l2(c) = litr2c(c)/litr2n(c)
      if (litr3n(c) > 0.) cn_l3(c) = litr3c(c)/litr3n(c)
      
      

      ! calculate the final rate scalar as the product of temperature and water       
      ! rate scalars, and correct the base decomp rates                                !chen:k乘以温度湿度系数

      rate_scalar = t_scalar(c) * w_scalar(c)
      rate_scalar = 0.8
      ck_l1 = k_l1 * rate_scalar
      ck_l2 = k_l2 * rate_scalar
      ck_l3 = k_l3 * rate_scalar
      ck_s1 = k_s1 * rate_scalar
      ck_s2 = k_s2 * rate_scalar
      ck_s3 = k_s3 * rate_scalar
      ck_s4 = k_s4 * rate_scalar
      ck_frag = k_frag * rate_scalar
!=======================================================================================================================      

!=================================================chen revision 1 matrix初始化==========================================
      kk_matrix   = 0.0 
      c_matrix_in_1 = 0.0 
      a_matrix    = 0.0 
      na_matrix   = 0.0
      c_matrix_current_1 = 0.0
      nk_matrix  = 0.0
!=======================================================================================================================      
            do j=1,8
                a_matrix(j,j) = -1.0 
                na_matrix(j,j) = -1.0
            enddo 
            
             

              
           do j =1,100
           !C matrix    
              a_matrix(3,1) = cwd_fcel
              a_matrix(4,1) = cwd_flig
              a_matrix(5,2) = 1.0 - rf_l1s1
              a_matrix(6,3) = 1.0 - rf_l2s2
              a_matrix(7,4) = 1.0 - rf_l3s3
              a_matrix(6,5) = 1.0 - rf_s1s2
              a_matrix(7,6) = 1.0 - rf_s2s3
              a_matrix(8,7) = 1.0 - rf_s3s4 
              
              kk_matrix(1,1) = ck_frag /dt
              kk_matrix(2,2) = ck_l1 /dt
              kk_matrix(3,3) = ck_l2 /dt 
              kk_matrix(4,4) = ck_l3 /dt 
              kk_matrix(5,5) = ck_s1 /dt
              kk_matrix(6,6) = ck_s2 /dt
              kk_matrix(7,7) = ck_s3 /dt
              kk_matrix(8,8) = ck_s4 /dt
              
          !N matrix not used yet      
              na_matrix(3,1) = cwd_fcel
              na_matrix(4,1) = cwd_flig
              if (litr1n(c) > 0.0) rate1 = cn_s1/cn_l1(c)
              na_matrix(5,2) = 1.0 - rf_l1s1 - rate1
              if (litr2n(c) > 0.0) rate2 = cn_s2/cn_l2(c)
              na_matrix(6,3) = 1.0 - rf_l2s2 - rate2
              if (litr3n(c) > 0.0) rate3 = cn_s3/cn_l3(c)
              na_matrix(7,4) = 1.0 - rf_l3s3 - rate3
              na_matrix(6,5) = 1.0 - rf_s1s2 - (cn_s2/cn_s1)
              na_matrix(7,6) = 1.0 - rf_s2s3 - (cn_s3/cn_s2)
              na_matrix(8,7) = 1.0 - rf_s3s4 - (cn_s4/cn_s3)
               
              nk_matrix(1,1) = ck_frag /cn_cwd /dt 
              nk_matrix(2,2) = ck_l1 /cn_l1(1) /dt
              nk_matrix(3,3) = ck_l2 /cn_l2(1) /dt 
              nk_matrix(4,4) = ck_l3 /cn_l3(1) /dt 
              nk_matrix(5,5) = ck_s1 /cn_s1 /dt
              nk_matrix(6,6) = ck_s2 /cn_s2 /dt
              nk_matrix(7,7) = ck_s3 /cn_s3 /dt
              nk_matrix(8,8) = ck_s4 /cn_s4 /dt
    !==================================chen: revision1, N库的潜在转移==================================================
              
              cwdn_to_litr2n(c) = cwdc(c) * nk_matrix(1,1) * na_matrix(3,1)
              cwdn_to_litr3n(c) = cwdc(c) * nk_matrix(1,1) * na_matrix(4,1)
              if (litr1n(c) > 0.0) then
              litr1n_to_soil1n(c) = litr1c(c) * nk_matrix(2,2) * na_matrix(5,2)! * cn_l1(1)/cn_s1
              end if
              if (litr2n(c) > 0.0) then
              litr2n_to_soil2n(c) = litr2c(c) * nk_matrix(3,3) * na_matrix(6,3)
              end if
              if (litr3n(c) > 0.0) then
              litr3n_to_soil3n(c) = litr3c(c) * nk_matrix(4,4) * na_matrix(7,4)
              end if 
             ! if (soil1n(c) > 0.0) then
              soil1n_to_soil2n(c) = soil1c(c) * nk_matrix(5,5) * na_matrix(6,5)
             ! end if
             ! if (soil2n(c) > 0.0) then
              soil2n_to_soil3n(c) = soil2c(c) * nk_matrix(6,6) * na_matrix(7,6)
             ! end if
             ! if (soil3n(c) > 0.0) then
              soil3n_to_soil4n(c) = soil3c(c) * nk_matrix(7,7) * na_matrix(8,7)
              !end if
        
              pmnf_l1s1(c) = litr1n_to_soil1n(c) * cn_l1(c)/cn_s1
              pmnf_l2s2(c) = litr2n_to_soil2n(c) * cn_l2(c)/cn_s2
              pmnf_l3s3(c) = litr3n_to_soil3n(c) * cn_l3(c)/cn_s3
              pmnf_s1s2(c) = soil1n_to_soil2n(c) * cn_s1/cn_s2
              pmnf_s2s3(c) = soil2n_to_soil3n(c) * cn_s2/cn_s3
              pmnf_s3s4(c) = soil3n_to_soil4n(c) * cn_s3/cn_s4
              pmnf_s4(c)   = -soil4c(c) * nk_matrix(8,8)
        !      print *,pmnf_l1s1(c),pmnf_l2s2(c),pmnf_l3s3(c),&
        !      pmnf_s1s2(c),pmnf_s2s3(c),pmnf_s3s4(c),pmnf_s4(c)
              
      immob(c) = 0.0
      ! litter 1 -> SOM 1
      if (pmnf_l1s1(c) > 0.0) then
         immob(c) = immob(c) + pmnf_l1s1(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 -> SOM 2
      if (pmnf_l2s2(c) > 0.0) then
         immob(c) = immob(c) + pmnf_l2s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 -> SOM 3
      if (pmnf_l3s3(c) > 0.0) then
         immob(c) = immob(c) + pmnf_l3s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 -> SOM 2
      if (pmnf_s1s2(c) > 0.0) then
         immob(c) = immob(c) + pmnf_s1s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 -> SOM 3
      if (pmnf_s2s3(c) > 0.0) then
         immob(c) = immob(c) + pmnf_s2s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 -> SOM 4
      if (pmnf_s3s4(c) > 0.0) then
         immob(c) = immob(c) + pmnf_s3s4(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4
      gross_nmin(c) = gross_nmin(c) - pmnf_s4(c)

      potential_immob(c) = immob(c)

  ! end do ! end column loop

   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.
   
   !predefined the N scalar
   fpi(1) = 0.5      
   c = 1
   dnp = 0.01

       !====================================================chen: revision1 更新输入输出变量=============================================================      
            !库间碳转移
            !异养呼吸
              cwdc_to_litr2c(c) = cwdc(c) * kk_matrix(1,1) * a_matrix(3,1) 
              cwdc_to_litr3c(c) = cwdc(c) * kk_matrix(1,1) * a_matrix(4,1) 
              !print *,pmnf_l1s1(1)
              if (pmnf_l1s1(c) > 0.) then
                litr1c_to_soil1c = litr1c(c) * kk_matrix(2,2) * a_matrix(5,2) * fpi(c)
                litr1_hr(c) = litr1c(c) * kk_matrix(2,2) * (1.0 - a_matrix(5,2))* fpi(c)
                kk_matrix(2,2) = ck_l1 /dt * fpi(c)
              else
                litr1c_to_soil1c = litr1c(c) * kk_matrix(2,2) * a_matrix(5,2)
                litr1_hr(c) = litr1c(c) * kk_matrix(2,2) * (1.0 - a_matrix(5,2))
              end if
              
              if (pmnf_l2s2(c) > 0.) then
                litr2c_to_soil2c = litr2c(c) * kk_matrix(3,3) * a_matrix(6,3) * fpi(c)
                litr2_hr(c) = litr2c(c) * kk_matrix(3,3) * (1.0 - a_matrix(6,3)) * fpi(c)
                kk_matrix(3,3) = ck_l2 /dt * fpi(c)
              else
                litr2c_to_soil2c = litr2c(c) * kk_matrix(3,3) * a_matrix(6,3)
                litr2_hr(c) = litr2c(c) * kk_matrix(3,3) * (1.0 - a_matrix(6,3))
              end if
              
              if (pmnf_l3s3(c) > 0.) then
                litr3c_to_soil3c = litr3c(c) * kk_matrix(4,4) * a_matrix(7,4)* fpi(c)
                litr3_hr(c) = litr3c(c) * kk_matrix(4,4) * (1.0 - a_matrix(7,4))*fpi(c)
                kk_matrix(4,4) = ck_l3 /dt * fpi(c)
              else
                litr3c_to_soil3c = litr3c(c) * kk_matrix(4,4) * a_matrix(7,4)
                litr3_hr(c) = litr3c(c) * kk_matrix(4,4) * (1.0 - a_matrix(7,4))
                
              end if
              
              if (pmnf_s1s2(c) > 0.) then
                soil1c_to_soil2c = soil1c(c) * kk_matrix(5,5) * a_matrix(6,5) * fpi(c)
                soil1_hr =  soil1c(c) * kk_matrix(5,5) * (1.0 - a_matrix(6,5)) * fpi(c)
                kk_matrix(5,5) = ck_s1 /dt * fpi(c)
              else
                soil1c_to_soil2c = soil1c(c) * kk_matrix(5,5) * a_matrix(6,5)   
                soil1_hr =  soil1c(c) * kk_matrix(5,5) * (1.0 - a_matrix(6,5))
              end if
              
              if (pmnf_s2s3(c) > 0.) then
                soil2c_to_soil3c = soil2c(c) * kk_matrix(6,6) * a_matrix(7,6) * fpi(c)
                soil2_hr =  soil2c(c) * kk_matrix(6,6) * (1.0 - a_matrix(7,6)) * fpi(c)
                kk_matrix(6,6) = ck_s2 /dt * fpi(c)
              
              else
                soil2c_to_soil3c = soil2c(c) * kk_matrix(6,6) * a_matrix(7,6)
                soil2_hr =  soil2c(c) * kk_matrix(6,6) * (1.0 - a_matrix(7,6))  
              end if
              
              if (pmnf_s3s4(c) > 0.) then
                soil3c_to_soil4c = soil3c(c) * kk_matrix(7,7) * a_matrix(8,7) * fpi(c)
                soil3_hr =  soil3c(c) * kk_matrix(7,7) * (1.0 - a_matrix(8,7)) * fpi(c)
                kk_matrix(7,7) = ck_s3 /dt * fpi(c)
              else
                soil3c_to_soil4c = soil3c(c) * kk_matrix(7,7) * a_matrix(8,7)
                soil3_hr =  soil3c(c) * kk_matrix(7,7) * (1.0 - a_matrix(8,7))    
              end if
                soil4_hr =  soil4c(c) * kk_matrix(8,8)

                
                
!==============================================================chen: revision 1 碳库更新=========================================================== 
              c_matrix_current_1(1,1) = cwdc(c) 
              c_matrix_current_1(2,1) = litr1c(c)
              c_matrix_current_1(3,1) = litr2c(c)
              c_matrix_current_1(4,1) = litr3c(c)
              c_matrix_current_1(5,1) = soil1c(c)
              c_matrix_current_1(6,1) = soil2c(c)
              c_matrix_current_1(7,1) = soil3c(c)
              c_matrix_current_1(8,1) = soil4c(c) 
              
              c_matrix_in_1(1,1)    = 0.
              c_matrix_in_1(2,1)    = 10 
              c_matrix_in_1(3,1)    = 5
              c_matrix_in_1(4,1)    = 5
              c_matrix_next_1(:,:) = c_matrix_in_1 + c_matrix_current_1 + &
                 matmul(matmul(a_matrix,kk_matrix),c_matrix_current_1)*dt      
              
              cwdc(c) = c_matrix_next_1(1,1)  
              litr1c(c) =  c_matrix_next_1(2,1)
              litr2c(c) = c_matrix_next_1(3,1)
              litr3c(c) = c_matrix_next_1(4,1)
              soil1c(c) = c_matrix_next_1(5,1)
              soil2c(c) = c_matrix_next_1(6,1)
              soil3c(c) = c_matrix_next_1(7,1)
              soil4c(c) =  c_matrix_next_1(8,1)
!================================================================================================================================================     
            !print *,c_matrix_next_1

   print *, cwdc(1),'cwdc'          ! (gC/m2) coarse woody debris C
   print *, litr1c(1),'litr1c'        ! (gC/m2) litter labile C
   print *, litr2c(1),'litr2c'        ! (gC/m2) litter cellulose C
   print *, litr3c(1),'litr3c'        ! (gC/m2) litter lignin C
   print *, soil1c(1),'soil1c'        ! (gC/m2) soil organic matter C (fast pool)
   print *, soil2c(1),'soil2c'        ! (gC/m2) soil organic matter C (medium pool)
   print *, soil3c(1),'soil3c'        ! (gC/m2) soil organic matter C (slow pool)
   print *, soil4c(1),'soil4c'        ! (gC/m2) soil organic matter C (slowest pool)
!   print *, cwdn(1),'cwdn'          ! (gN/m2) coarse woody debris N
!   print *, litr1n(1),'litr1n'        ! (gN/m2) litter labile N
!   print *, litr2n(1),'litr2n'        ! (gN/m2) litter cellulose N
!   print *, litr3n(1),'litr3n'        ! (gN/m2) litter lignin N
!   print *, clandunit(1),'clandunit'      ! index into landunit level quantities
!   print *, itypelun(1),'itypelun'      ! landunit type
!   print *, rootfr(1),'rootfr'      ! fraction of roots in each soil layer  (nlevgrnd)

!   print *,  fpi(1),'fpi'           ! fraction of potential immobilization (no units)
   print *,  cwdc_to_litr2c(1),'cwdc_to_litr2c'
   print *,  cwdc_to_litr3c(1),'cwdc_to_litr3c'
   print *,  litr1_hr(1),'litr1_hr'
   print *,  litr1c_to_soil1c(1),'litr1c_to_soil1c'
   print *,  litr2_hr(1),'litr2_hr'
   print *,  litr2c_to_soil2c(1),'litr2c_to_soil2c'
   print *,  litr3_hr(1),'litr3_hr'
   print *,  litr3c_to_soil3c(1),'litr3c_to_soil3c'
   print *,  soil1_hr(1),'soil1_hr'
   print *,  soil1c_to_soil2c(1),'soil1c_to_soil2c'
   print *,  soil2_hr(1),'soil2_hr'
   print *,  soil2c_to_soil3c(1),'soil2c_to_soil3c'
   print *,  soil3_hr(1),'soil3_hr'
   print *,  soil3c_to_soil4c(1),'soil3c_to_soil4c'
   print *,  soil4_hr(1),'soil4_hr'
!   print *,  cwdn_to_litr2n(1),'cwdn_to_litr2n'
!   print *,  cwdn_to_litr3n(1),'cwdn_to_litr3n'
!   print *,  potential_immob(1),'potential_immob'
!   print *,  litr1n_to_soil1n(1),'litr1n_to_soil1n'
!   print *,  sminn_to_soil1n_l1(1),'sminn_to_soil1n_l1'
!   print *,  litr2n_to_soil2n(1),'litr2n_to_soil2n'
!   print *,  sminn_to_soil2n_l2(1),'sminn_to_soil2n_l2'
!   print *,  litr3n_to_soil3n(1),'litr3n_to_soiln'
!   print *,  sminn_to_soil3n_l3(1),'sminn_to_soil3n_l3'
!   print *,  soil1n_to_soil2n(1),'soil1n_to_soil2n'
!   print *,  sminn_to_soil2n_s1(1),'sminn_to_soil2n_s1'
!   print *,  soil2n_to_soil3n(1),'soil2n_tosoil3n'
!   print *,  sminn_to_soil3n_s2(1),'sminn_tosoil3n_s2'
!   print *,  soil3n_to_soil4n(1),'soil3n_to_soil4n'
!   print *,  sminn_to_soil4n_s3(1),'sminn_to_soil4n_s3'
!   print *,  soil4n_to_sminn(1),'soil4n_to_sminn'
!   print *,  sminn_to_denit_l1s1(1),'sminn_to_denit_l1s1'
!   print *,  sminn_to_denit_l2s2(1),'sminn_to_denit_l2s2'
!   print *,  sminn_to_denit_l3s3(1),'sminn_to_denit_l3s3'
!   print *,  sminn_to_denit_s1s2(1),'sminn_to_denit_s1s2'
!   print *,  sminn_to_denit_s2s3(1),'sminn_to_denit_s2s3'
!   print *,  sminn_to_denit_s3s4(1),'sminn_to_denit_s3s4'
!   print *,  sminn_to_denit_s4(1),'sminn_to_denit_s4'
!   print *,  sminn_to_denit_excess(1),'sminn_to_denit_excess'
!   print *,  gross_nmin(1),'gross_nmin'
!   print *,  net_nmin(1),'net_nmin'

    write (55,163)dt,j,cwdc(1),litr1c(1),litr2c(1),litr3c(1),soil1c(1),soil2c(1),soil3c(1),soil4c(1),&
                     &cwdc_to_litr2c(1),cwdc_to_litr3c(1),litr1_hr(1),litr1c_to_soil1c(1),&
                     &litr2_hr(1),litr2c_to_soil2c(1),litr3_hr(1),litr3c_to_soil3c(1),&
                     &soil1_hr(1),soil1c_to_soil2c(1),soil2_hr(1),soil2c_to_soil3c(1),&
                     &soil3_hr(1),soil3c_to_soil4c(1),soil4_hr(1)
    163           format(I6,",",120(f12.4,","))
    end do
        


    end program DEMO_for_matrix_DecompMod

