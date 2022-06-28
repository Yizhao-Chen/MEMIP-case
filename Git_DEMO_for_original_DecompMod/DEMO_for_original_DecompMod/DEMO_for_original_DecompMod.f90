!  DEMO_for_original_DecompMod.f90 
!
!  FUNCTIONS:
!  DEMO_for_original_DecompMod - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: DEMO_for_original_DecompMod
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program DEMO_for_original_DecompMod
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
   
   
   !chen:数据输出
   open(55,file="C_outputs.csv")
   write(55,'(a500)')",runtime,cwdc,litr1c, litr2c,litr3c,soil1c,soil2c,soil3c,soil4c,&
                     &cwdc_to_litr2c,cwdc_to_litr3c,litr1_hr,litr1c_to_soil1c,&
                     &litr2_hr,litr2c_to_soil2c,litr3_hr,litr3c_to_soil3c,&
                     &soil1_hr,soil1c_to_soil2c,soil2_hr,soil2c_to_soil3c,&
                     &soil3_hr,soil3c_to_soil4c,soil4_hr"
   !end
   
    !&    N_up,N_fix,N_depo,N_leach,N_vol,N_def, Q_leaf,Q_wood,Q_root
   ! &   Q_fine,Q_coarse,Q_Micr,Q_Slow,Q_Pass,
   ! &   N_leaf,N_wood,N_root,N_fine,
   ! &   N_coarse,N_micr,N_slowC,N_pass,QNminer"

    !chen:定义初始状态变量 
    !chen: 碳库初始值
    cwdc(1) = 0.          ! (gC/m2) coarse woody debris C
    litr1c(1) = 50.        ! (gC/m2) litter labile C
    litr2c(1) = 30.        ! (gC/m2) litter cellulose C
    litr3c(1) = 20.       ! (gC/m2) litter lignin C
    soil1c(1) = 500.       ! (gC/m2) soil organic matter C (fast pool)
    soil2c(1) = 1000.       ! (gC/m2) soil organic matter C (medium pool)
    soil3c(1) = 2000.       ! (gC/m2) soil organic matter C (slow pool)
    soil4c(1) = 4000.       ! (gC/m2) soil organic matter C (slowest pool)
    !chen: 氮库初始值 DEMO不考虑氮循环
    cwdn(1)   = 30.       ! (gN/m2) coarse woody debris N
    litr1n(1) = 2.       ! (gN/m2) litter labile N
    litr2n(1) = 2.       ! (gN/m2) litter cellulose N
    litr3n(1) = 2.       ! (gN/m2) litter lignin N
    
    
    
    !chen:设置测试时间步长
    dt = 1.0
   dtd = 1.0    

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
    
      ! set initial values for potential C and N fluxes                            !chen:设置潜在碳氮通量的初始0值
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
 
   c=1  !设置column变量为定值，不进行column level loop
      
      

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
                
!==========================================================chen: revission 1, CLM4.0原来的C N分解过程=====================================================          
      ! calculate the non-nitrogen-limited fluxes                                     !chen:这里开始正式运行 第一部分（call之前）主要输出N固定的部分给CN_allocation,第二部分更新碳库
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.
!=============chen: 潜在分解速率计算,即不考虑N限制的分解速率，我将相关系数设为定值，这部分暂时不必考虑修改================================================
    do j =1,100
      ! CWD fragmentation -> litter pools
      cwdc_loss = cwdc(c) * ck_frag / dt
      cwdc_to_litr2c(c) = cwdc_loss * cwd_fcel
      cwdc_to_litr3c(c) = cwdc_loss * cwd_flig
      
      
      cwdn_loss = cwdn(c) * ck_frag / dt
      cwdn_to_litr2n(c) = cwdn_loss * cwd_fcel
      cwdn_to_litr3n(c) = cwdn_loss * cwd_flig

      ! litter 1 -> SOM 1
      if (litr1c(c) > 0.) then
         plitr1c_loss(c) = litr1c(c) * ck_l1 / dt
         ratio = 0.
         if (litr1n(c) > 0.) ratio = cn_s1/cn_l1(c)
         pmnf_l1s1(c) = (plitr1c_loss(c) * (1.0 - rf_l1s1 - ratio))/cn_s1
         
      end if  

      ! litter 2 -> SOM 2
      if (litr2c(c) > 0.) then
         plitr2c_loss(c) = litr2c(c) * ck_l2 / dt
         ratio = 0.
         if (litr2n(c) > 0.) ratio = cn_s2/cn_l2(c)
         pmnf_l2s2(c) = (plitr2c_loss(c) * (1.0 - rf_l2s2 - ratio))/cn_s2
      end if

      ! litter 3 -> SOM 3
      if (litr3c(c) > 0.) then
         plitr3c_loss(c) = litr3c(c) * ck_l3 / dt
         ratio = 0.
         if (litr3n(c) > 0.) ratio = cn_s3/cn_l3(c)
         pmnf_l3s3(c) = (plitr3c_loss(c) * (1.0 - rf_l3s3 - ratio))/cn_s3
      end if

      ! SOM 1 -> SOM 2
      if (soil1c(c) > 0.) then
         psoil1c_loss(c) = soil1c(c) * ck_s1 / dt
         pmnf_s1s2(c) = (psoil1c_loss(c) * (1.0 - rf_s1s2 - (cn_s2/cn_s1)))/cn_s2
      end if

      ! SOM 2 -> SOM 3
      if (soil2c(c) > 0.) then
         psoil2c_loss(c) = soil2c(c) * ck_s2 / dt
         pmnf_s2s3(c) = (psoil2c_loss(c) * (1.0 - rf_s2s3 - (cn_s3/cn_s2)))/cn_s3
      end if

      ! SOM 3 -> SOM 4
      if (soil3c(c) > 0.) then
         psoil3c_loss(c) = soil3c(c) * ck_s3 / dt
         pmnf_s3s4(c) = (psoil3c_loss(c) * (1.0 - rf_s3s4 - (cn_s4/cn_s3)))/cn_s4
      end if

      ! Loss from SOM 4 is entirely respiration (no downstream pool)
      if (soil4c(c) > 0.) then
         psoil4c_loss(c) = soil4c(c) * ck_s4 / dt
         pmnf_s4(c) = -psoil4c_loss(c)/cn_s4
      end if

      !print *,pmnf_l1s1(c),pmnf_l2s2(c),pmnf_l3s3(c),&
      !        pmnf_s1s2(c),pmnf_s2s3(c),pmnf_s3s4(c),pmnf_s4(c)
      ! Sum up all the potential immobilization fluxes (positive pmnf flux)
      ! and all the mineralization fluxes (negative pmnf flux)

      immob(c) = 0.
      ! litter 1 -> SOM 1
      if (pmnf_l1s1(c) > 0.) then
         immob(c) = immob(c) + pmnf_l1s1(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 -> SOM 2
      if (pmnf_l2s2(c) > 0.) then
         immob(c) = immob(c) + pmnf_l2s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 -> SOM 3
      if (pmnf_l3s3(c) > 0.) then
         immob(c) = immob(c) + pmnf_l3s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 -> SOM 2
      if (pmnf_s1s2(c) > 0.) then
         immob(c) = immob(c) + pmnf_s1s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 -> SOM 3
      if (pmnf_s2s3(c) > 0.) then
         immob(c) = immob(c) + pmnf_s2s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 -> SOM 4
      if (pmnf_s3s4(c) > 0.) then
         immob(c) = immob(c) + pmnf_s3s4(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4
      gross_nmin(c) = gross_nmin(c) - pmnf_s4(c)

      potential_immob(c) = immob(c)

 ! end column loop

   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.

!!!=========================================================chen: 潜在分解速率计算完成==================================================================

   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N

   dnp = 0.01
   fpi(1) = 0.5      !chen: 将N限制因子设为定值    
   c = 1

   

      ! upon return from CNAllocation, the fraction of potential immobilization
      ! has been set (cps%fpi). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.

      ! litter 1 fluxes (labile pool)
      if (litr1c(c) > 0.) then
          !print *, pmnf_l1s1(1)
         if (pmnf_l1s1(c) > 0.) then
            plitr1c_loss(c) = plitr1c_loss(c) * fpi(c)
            pmnf_l1s1(c) = pmnf_l1s1(c) * fpi(c)
            sminn_to_denit_l1s1(c) = 0.
         else
            sminn_to_denit_l1s1(c) = -dnp * pmnf_l1s1(c)
         end if
         
         litr1_hr(c) = rf_l1s1 * plitr1c_loss(c)                       !chen：litr1c异养呼吸输出变量
         litr1c_to_soil1c(c) = (1. - rf_l1s1) * plitr1c_loss(c)        !chen: litr1c至soil1c传输变量
         if (litr1n(c) > 0.) litr1n_to_soil1n(c) = plitr1c_loss(c) / cn_l1(c)
         sminn_to_soil1n_l1(c) = pmnf_l1s1(c)
         net_nmin(c) = net_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 fluxes (cellulose pool)
      if (litr2c(c) > 0.) then
         if (pmnf_l2s2(c) > 0.) then
            plitr2c_loss(c) = plitr2c_loss(c) * fpi(c)
            pmnf_l2s2(c) = pmnf_l2s2(c) * fpi(c)
            sminn_to_denit_l2s2(c) = 0.
         else
            sminn_to_denit_l2s2(c) = -dnp * pmnf_l2s2(c)
         end if
         litr2_hr(c) = rf_l2s2 * plitr2c_loss(c)                       !chen：litr2c异养呼吸输出变量
         litr2c_to_soil2c(c) = (1. - rf_l2s2) * plitr2c_loss(c)        !chen: litr2c至soil2c传输变量
         if (litr2n(c) > 0.) litr2n_to_soil2n(c) = plitr2c_loss(c) / cn_l2(c)
         sminn_to_soil2n_l2(c) = pmnf_l2s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 fluxes (lignin pool)
      if (litr3c(c) > 0.) then
         if (pmnf_l3s3(c) > 0.) then
            plitr3c_loss(c) = plitr3c_loss(c) * fpi(c)
            pmnf_l3s3(c) = pmnf_l3s3(c) * fpi(c)
            sminn_to_denit_l3s3(c) = 0.
         else
            sminn_to_denit_l3s3(c) = -dnp * pmnf_l3s3(c)
         end if
         litr3_hr(c) = rf_l3s3 * plitr3c_loss(c)                      !chen：litr3c异养呼吸输出变量
         litr3c_to_soil3c(c) = (1. - rf_l3s3) * plitr3c_loss(c)       !chen：litr3c至soil3c传输变量
         if (litr3n(c) > 0.) litr3n_to_soil3n(c) = plitr3c_loss(c) / cn_l3(c)
         sminn_to_soil3n_l3(c) = pmnf_l3s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 fluxes (fast rate soil organic matter pool)
      if (soil1c(c) > 0.) then
         if (pmnf_s1s2(c) > 0.) then
            psoil1c_loss(c) = psoil1c_loss(c) * fpi(c)
            pmnf_s1s2(c) = pmnf_s1s2(c) * fpi(c)
            sminn_to_denit_s1s2(c) = 0.
         else
            sminn_to_denit_s1s2(c) = -dnp * pmnf_s1s2(c)
         end if
         soil1_hr(c) = rf_s1s2 * psoil1c_loss(c)                     !chen：soil1c异养呼吸输出变量
         soil1c_to_soil2c(c) = (1. - rf_s1s2) * psoil1c_loss(c)      !chen：soil1c至soil2c传输变量
         soil1n_to_soil2n(c) = psoil1c_loss(c) / cn_s1
         sminn_to_soil2n_s1(c) = pmnf_s1s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 fluxes (medium rate soil organic matter pool)
      if (soil2c(c) > 0.) then
         if (pmnf_s2s3(c) > 0.) then
            psoil2c_loss(c) = psoil2c_loss(c) * fpi(c)
            pmnf_s2s3(c) = pmnf_s2s3(c) * fpi(c)
            sminn_to_denit_s2s3(c) = 0.
         else
            sminn_to_denit_s2s3(c) = -dnp * pmnf_s2s3(c)
         end if
         soil2_hr(c) = rf_s2s3 * psoil2c_loss(c)                     !chen：soil2c异养呼吸输出变量
         soil2c_to_soil3c(c) = (1. - rf_s2s3) * psoil2c_loss(c)      !chen：soil2c至soil3c传输变量
         soil2n_to_soil3n(c) = psoil2c_loss(c) / cn_s2
         sminn_to_soil3n_s2(c) = pmnf_s2s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 fluxes (slow rate soil organic matter pool)
      if (soil3c(c) > 0.) then
         if (pmnf_s3s4(c) > 0.) then
            psoil3c_loss(c) = psoil3c_loss(c) * fpi(c)
            pmnf_s3s4(c) = pmnf_s3s4(c) * fpi(c)
            sminn_to_denit_s3s4(c) = 0.
         else
            sminn_to_denit_s3s4(c) = -dnp * pmnf_s3s4(c)
         end if
         soil3_hr(c) = rf_s3s4 * psoil3c_loss(c)                     !chen：soil3c异养呼吸输出变量
         soil3c_to_soil4c(c) = (1. - rf_s3s4) * psoil3c_loss(c)      !chen：soil3c至soil4c传输变量
         soil3n_to_soil4n(c) = psoil3c_loss(c) / cn_s3
         sminn_to_soil4n_s3(c) = pmnf_s3s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4 fluxes (slowest rate soil organic matter pool)
      if (soil4c(c) > 0.) then
         soil4_hr(c) = psoil4c_loss(c)                               !chen：soil4c异养呼吸输出变量
         soil4n_to_sminn(c) = psoil4c_loss(c) / cn_s4
         sminn_to_denit_s4(c) = -dnp * pmnf_s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s4(c)
      end if
      
      
      !chen:各状态变量碳通量总结，屏幕输出和文件写入
      cwdc(1) = 0 + cwdc(1) - cwdc_loss
      litr1c(1) = 10.0 + litr1c(1) - plitr1c_loss(c)
      litr2c(1) = 5.0 + cwdc_to_litr2c(1) + litr2c(1) - plitr2c_loss(c)
      litr3c(1) = 5.0 + cwdc_to_litr3c(1) + litr3c(1) - plitr3c_loss(c)
      soil1c(1) = litr1c_to_soil1c(1) + soil1c(1) - psoil1c_loss(c)
      soil2c(1) = soil1c_to_soil2c(1) + litr2c_to_soil2c(1) + soil2c(1) - psoil2c_loss(c)
      soil3c(1) = soil2c_to_soil3c(1) + litr3c_to_soil3c(1) + soil3c(1) - psoil3c_loss(c)
      soil4c(1) = soil3c_to_soil4c(1) + soil4c(1) - psoil4c_loss(c)
        
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
    end program DEMO_for_original_DecompMod

