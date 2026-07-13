module SoilWaterMovementLateralMod

  use ColumnDataType     , only : col_es, col_ws, col_wf
  use VegetationDataType , only : veg_wf
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use elm_instMod        , only : waterflux_vars, waterstate_vars, temperature_vars
  use abortutils         , only : endrun
  use spmdMod            , only : masterproc, iam, npes, mpicom, comp_id


  implicit none
  save
  private

  public :: ComputeLateralUnsatFlux
  public :: SolveLateralSatFlow

contains
  !
  ! scatter_data_for_zengdecker
  ! ComputeLateralFlux
  !   ExtractUpOrDownHydroVarabiables
  !
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ComputeLateralUnsatFlux(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, num_h3dc, filter_h3dc, soilhydrology_vars, soilstate_vars, jwt, qflx_lateral_s)
    !
    ! !DESCRIPTION:
    ! Calculate watertable, considering aquifer recharge but no drainage.
    !
    ! !USES:
    use decompMod                 , only : bounds_type
    use elm_varctl                , only : use_var_soil_thick
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use shr_const_mod             , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use decompMod                 , only : bounds_type
    use elm_varcon                , only : wimp,grav,hfus,tfrz
    use elm_varcon                , only : e_ice,denh2o, denice
    use elm_varpar                , only : nlevsoi, max_patch_per_col, nlevgrnd
    use elm_time_manager          , only : get_step_size
    use column_varcon             , only : icol_roof, icol_road_imperv
    use TridiagonalMod            , only : Tridiagonal
    use SoilStateType             , only : soilstate_type
    use SoilHydrologyType         , only : soilhydrology_type
    use VegetationType            , only : veg_pp
    use ColumnType                , only : col_pp
    use TopounitType              , only : top_pp
    use spmdMod                   , only : masterproc, iam, npes, mpicom, comp_id
    use domainLateralMod          , only : ldomain_lateral
    use elm_varctl                , only : iulog
    use elm_varcon                , only : rpi
    use elm_varpar                , only : nh3dc_per_lunit
    use LandunitType              , only : lun_pp

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_h3dc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_h3dc(:) ! columnfilter for h3d points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    integer                  , intent(in)    :: jwt(bounds%begc:bounds%endc)
    real(r8)                 , intent(out)   :: qflx_lateral_s(bounds%begc:bounds%endc,1:nlevgrnd+1)
    !
    ! !LOCAL VARIABLES:
    integer                    :: fc,c,j, col_id_up, col_id_dn, l, k
    real(r8)                   :: h2osoi_vol_up, h2osoi_vol_dn
    real(r8)                   :: watsat_up, watsat_dn
    real(r8)                   :: bsw_up, bsw_dn
    real(r8)                   :: hksat_up, hksat_dn
    real(r8)                   :: fracice_up, fracice_dn
    real(r8)                   :: icefrac_up, icefrac_dn
    real(r8)                   :: smp_up, smp_dn
    real(r8)                   :: dzgmm, bswl, den, s1, s2, hkl, impedl, qflx_up_to_dn, theta
    real(r8)                   :: h2osoi_vol_up_1d(1:nlevgrnd), h2osoi_vol_dn_1d(1:nlevgrnd)
    real(r8)                   :: watsat_up_1d(1:nlevgrnd), watsat_dn_1d(1:nlevgrnd)
    real(r8)                   :: bsw_up_1d(1:nlevgrnd), bsw_dn_1d(1:nlevgrnd)
    real(r8)                   :: hksat_up_1d(1:nlevgrnd), hksat_dn_1d(1:nlevgrnd)
    real(r8)                   :: fracice_up_1d(1:nlevgrnd), fracice_dn_1d(1:nlevgrnd)
    real(r8)                   :: icefrac_up_1d(1:nlevgrnd), icefrac_dn_1d(1:nlevgrnd)
    real(r8)                   :: smp_up_1d(1:nlevgrnd), smp_dn_1d(1:nlevgrnd)
    integer                    :: c0, h3d_begc,h3d_endc
    logical                    :: up_soil_layer_saturated, dn_soil_layer_saturated
    logical                    :: up_local, dn_local

    !-----------------------------------------------------------------------

    associate(                                                &
         nbedrock       => col_pp%nlevbed                   , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)
         dz             => col_pp%dz                        , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         z              => col_pp%z                         , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zi             => col_pp%zi                        , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         h2osoi_liq     => col_ws%h2osoi_liq                , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice     => col_ws%h2osoi_ice                , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_vol     => col_ws%h2osoi_vol                , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         bsw            => soilstate_vars%bsw_col           , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         hksat          => soilstate_vars%hksat_col         , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         smp_l          => soilstate_vars%smp_l_col         , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
         watsat         => soilstate_vars%watsat_col        , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
         origflag       => soilhydrology_vars%origflag      , & ! Input:  constant
         fracice        => soilhydrology_vars%fracice_col   , & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)
         icefrac        => soilhydrology_vars%icefrac_col   , & ! Input:  [real(r8) (:,:) ]  fraction of ice
          hs_dx            =>    lun_pp%hs_dx                           ,      &
! Input   [real(r8) (:,:) ]  
          hs_dx_nod        =>    lun_pp%hs_dx_node                      ,      &
! Input   [real(r8) (:,:) ]  
          hs_slope         =>    col_pp%h3d_slope                       ,      &
! Input:  [real(r8) (:)   ] gridcell topographic slope (degree)
         zwt            => soilhydrology_vars%zwt_col        & ! Input: [real(r8) (:)   ]  water table depth (m)
         )


      if (use_var_soil_thick) then
         write(iulog,*)'use_var_soil_thick not supported for lateral flow'
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if
      qflx_lateral_s(bounds%begc:bounds%endc,1:nlevgrnd+1) = 0._r8
!print *,'chk---',num_h3dc,num_hydrologyc, nh3dc_per_lunit
   do fc = 1,num_h3dc,nh3dc_per_lunit  !loop for all soil columns that belong to same land unit
      c0 = filter_h3dc(fc)
      h3d_begc = c0
      h3d_endc = c0+nh3dc_per_lunit-1
      l  = col_pp%landunit(c0)

      do k = 1, nh3dc_per_lunit-1
         col_id_up = c0+(k+1)-1
         col_id_dn = c0+k-1
         den = hs_dx_nod(l,k)*1000._r8
         
         c = col_id_dn

         do j = 1, nlevgrnd


            up_soil_layer_saturated = (j > jwt(col_id_up) - 1)

            dn_soil_layer_saturated = (j > jwt(col_id_dn) - 1)

            if (up_soil_layer_saturated .or. dn_soil_layer_saturated) then
               ! Do not compute unsaturated lateral flux because WT in one of
               ! the column is shallower than 'j'-layer. Saturated lateral
               ! flux will be accounted for in SolveLateralSatFlow()
            else

               !hydraulic conductivity hkl(iconn,j) is
               !the lateral hydraulic conductivity is calculated using the geometric mean of the
               !neighbouring lateral cells and is approximated as 1000 times of the vertical hydraulic conductivity
               h2osoi_vol_up = h2osoi_vol(col_id_up,j)
               h2osoi_vol_dn = h2osoi_vol(col_id_dn,j)
               watsat_up     = watsat(col_id_up,j)
               watsat_dn     = watsat(col_id_dn,j)
               bsw_up        = bsw(col_id_up,j)
               bsw_dn        = bsw(col_id_dn,j)
               hksat_up      = hksat(col_id_up,j)
               hksat_dn      = hksat(col_id_dn,j)
               fracice_up    = fracice(col_id_up,j)
               fracice_dn    = fracice(col_id_dn,j)
               icefrac_up    = icefrac(col_id_up,j)
               icefrac_dn    = icefrac(col_id_dn,j)
               smp_up        = smp_l(col_id_up,j)
               smp_dn        = smp_l(col_id_dn,j)

               s1 = (h2osoi_vol_up + h2osoi_vol_dn)/(watsat_up + watsat_dn)
               s1 = min(1._r8, s1)

               bswl = (bsw_up + bsw_dn)/2._r8

               s2 = sqrt(hksat_up*hksat_dn)*s1**(2._r8*bswl + 3._r8)

               ! replace fracice with impedance factor, as in zhao 97,99
               if (origflag == 1) then
                  impedl = (1._r8-0.5_r8*(fracice_up + fracice_dn))
               else
                  impedl = 10._r8**(-e_ice*(0.5_r8*(icefrac_up + icefrac_dn)))
               endif

               hkl = impedl*s1*s2*10.0_r8
               theta = hs_slope(c)/180._r8*rpi
               qflx_up_to_dn = -(hkl*(smp_dn - smp_up)*cos(theta)/den + hkl*(sin(theta)))

               qflx_lateral_s(col_id_up,j) = qflx_lateral_s(col_id_up,j) - qflx_up_to_dn &
                          /hs_dx(l,k)

               qflx_lateral_s(col_id_dn,j) = qflx_lateral_s(col_id_dn,j)  + qflx_up_to_dn &
                          /hs_dx(l,k+1)
            endif
         enddo
      enddo
   enddo

    end associate

  end subroutine ComputeLateralUnsatFlux

  !-----------------------------------------------------------------------
  subroutine SolveLateralSatFlow(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, num_h3dc, filter_h3dc, soilhydrology_vars, soilstate_vars, jwt)
    !
    ! !DESCRIPTION:
    ! Calculate watertable while accounting for lateral flow
    !
    ! !USES:
    use decompMod                 , only : bounds_type
    use elm_varctl                , only : use_var_soil_thick
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use decompMod                 , only : bounds_type
    use elm_varpar                , only : nlevgrnd
    use elm_time_manager          , only : get_step_size
    use SoilStateType             , only : soilstate_type
    use SoilHydrologyType         , only : soilhydrology_type
    use ColumnType                , only : col_pp
    use domainLateralMod          , only : ldomain_lateral
    use elm_varcon                , only : rpi
    use LandunitType              , only : lun_pp
    use elm_varpar                , only : nh3dc_per_lunit
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(in)    :: num_h3dc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: filter_h3dc(:) ! columnfilter for h3d points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    integer                  , intent(inout) :: jwt(bounds%begc:bounds%endc)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,j, l, k, c0
    integer :: h3d_begc,h3d_endc
    integer :: col_id_up, col_id_dn
    integer :: step, nlevbed
    integer :: nstep
    logical :: up_local, dn_local
    real(r8) :: hksat_up, hksat_dn
    real(r8) :: den, qflx_up_to_dn
    real(r8) :: depth_up, depth_down
    real(r8) :: dtime, qlat_layer, qlat_tot, qlat_temp, s_y
    real(r8) :: rous, sy, trans, theta
    real(r8) :: qflx_lateral_s(bounds%begc:bounds%endc)

    !-----------------------------------------------------------------------

    associate(                                                &
         zi             => col_pp%zi                        , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         nlev2bed       => col_pp%nlevbed                   , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         h2osoi_liq     => col_ws%h2osoi_liq                , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)

         bsw            => soilstate_vars%bsw_col           , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         sucsat         => soilstate_vars%sucsat_col        , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         hksat          => soilstate_vars%hksat_col         , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         smp_l          => soilstate_vars%smp_l_col         , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
         watsat         => soilstate_vars%watsat_col        , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)

         origflag       => soilhydrology_vars%origflag      , & ! Input:  constant
         dz             => col_pp%dz                        , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zwt            => soilhydrology_vars%zwt_col       , & ! Input: [real(r8) (:)   ]  water table depth (m)
         hs_dx            =>    lun_pp%hs_dx                           ,      &
! Input   [real(r8) (:,:) ]  
         hs_dx_nod        =>    lun_pp%hs_dx_node                      ,      &
! Input   [real(r8) (:,:) ]  
         hs_slope         =>    col_pp%h3d_slope                       ,      &
! Input:  [real(r8) (:)   ] gridcell topographic slope (degree)
         wa             => soilhydrology_vars%wa_col         & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)

         )

      ! Water table changes due to qlateral in saturated GW
      nstep=1
      dtime = get_step_size()

      do fc = 1,num_h3dc,nh3dc_per_lunit  !loop for all soil columns that belong to same land unit
        c0 = filter_h3dc(fc)
        h3d_begc = c0
        h3d_endc = c0+nh3dc_per_lunit-1
        l  = col_pp%landunit(c0)

        qflx_lateral_s=0._r8

        do k = 1, nh3dc_per_lunit-1
           col_id_up = c0+(k+1)-1
           col_id_dn = c0+k-1
           c = col_id_dn

               ! local grid cell
            depth_up = zi(col_id_up,nlevgrnd) - zwt(col_id_up)  ! groundwater head(m)
            hksat_up = hksat(col_id_up,15)

            depth_down = zi(col_id_dn,nlevgrnd) - zwt(col_id_dn)
            hksat_dn = hksat(col_id_dn,15)

            den = hs_dx_nod(l,k)*1000._r8

            depth_up = max(depth_up, 0._r8)
            depth_down= max(depth_down, 0._r8)
            theta = hs_slope(c)/180._r8*rpi

            ! calculate transmissivity
            trans = 1.0_r8*sqrt(hksat_up*hksat_dn)*(depth_up+depth_down)/2._r8*1000._r8 ! (mm2/s)
            qflx_up_to_dn = -trans*((depth_down-depth_up)*1000._r8/den*cos(theta) + sin(theta))

            qflx_lateral_s(col_id_up) = qflx_lateral_s(col_id_up) - &
                    qflx_up_to_dn/1000._r8* dz(col_id_dn,j)/hs_dx(l,k+1)*cos(theta)

            qflx_lateral_s(col_id_dn) = qflx_lateral_s(col_id_dn) + &
                    qflx_up_to_dn/1000._r8* dz(col_id_dn,j)/hs_dx(l,k+1)*cos(theta)
         enddo
      enddo

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            nlevbed = nlev2bed(c)

            !scs: use analytical expression for aquifer specific yield
            rous = watsat(c,nlevbed) &
                 * ( 1._r8 - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
            rous=max(rous,0.02_r8)

            !--  water table is below the soil column  --------------------------------------

            qlat_temp = qflx_lateral_s(c)

            if(jwt(c) == nlevgrnd) then
               wa(c)  = wa(c) + qflx_lateral_s(c)  * dtime/nstep
               zwt(c) = zwt(c) - (qflx_lateral_s(c) * dtime)/nstep/1000._r8/rous
            else
               !-- water table within soil layers 1-15  -------------------------------------
               ! try to raise water table to account for qlat
               qlat_tot = qflx_lateral_s(c) * dtime/nstep
               if(qlat_tot > 0._r8) then !rising water table ! need to modify soil water content also, Han Qiu
                  do j = jwt(c)+1, 1,-1
                     !scs: use analytical expression for specific yield
                     s_y = watsat(c,j) &
                          * ( 1._r8 -  (1._r8 + 1.e3*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
                     s_y=max(s_y,0.02_r8)
                     qlat_layer=min(qlat_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                     qlat_layer=max(qlat_layer,0._r8)

                     h2osoi_liq(c,j) = h2osoi_liq(c,j) + qlat_layer

                     if(s_y > 0._r8) zwt(c) = zwt(c) - qlat_layer/s_y/1000._r8

                     qlat_tot = qlat_tot - qlat_layer
                     if (qlat_tot <= 0._r8) exit
                  enddo
               else ! deepening water table (negative qlat)
                  do j = jwt(c)+1, 15
                     !scs: use analytical expression for specific yield
                     s_y = watsat(c,j) &
                          * ( 1._r8 -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
                     s_y=max(s_y,0.02_r8)

                     qlat_layer=max(qlat_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                     qlat_layer=min(qlat_layer,0._r8)
                     h2osoi_liq(c,j) = h2osoi_liq(c,j) + qlat_layer
                     qlat_tot = qlat_tot - qlat_layer

                     if (qlat_tot >= 0._r8) then
                        zwt(c) = zwt(c) - qlat_layer/s_y/1000._r8
                        exit
                     else
                        zwt(c) = zi(c,j)
                     endif

                  enddo
                  if (qlat_tot > 0._r8) zwt(c) = zwt(c) - qlat_tot/1000._r8/rous
                  !rsub_top(c) = rsub_top(c) + qlat_tot/dtime

               endif

               !-- recompute jwt for following calculations  ---------------------------------
               ! allow jwt to equal zero when zwt is in top layer
               !jwt(c) = nlevbed
               jwt(c) = nlevgrnd
               do j = 1,nlevgrnd
                  if(zwt(c) <= zi(c,j)) then
                     jwt(c) = j-1
                     exit
                  end if
               enddo
            endif
         enddo

         call ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
              soilhydrology_vars, soilstate_vars)


    end associate

  end subroutine SolveLateralSatFlow

  !-----------------------------------------------------------------------
  subroutine ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars) 
    !
    ! !DESCRIPTION:
    ! Calculate watertable, considering aquifer recharge but no drainage.
    !
    ! !USES:
    use decompMod                 , only : bounds_type
    use elm_varctl                , only : use_var_soil_thick
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use shr_const_mod             , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use decompMod                 , only : bounds_type
    use elm_varcon                , only : wimp,grav,hfus,tfrz
    use elm_varcon                , only : e_ice,denh2o, denice
    use elm_varpar                , only : nlevsoi, max_patch_per_col, nlevgrnd
    use elm_time_manager          , only : get_step_size
    use column_varcon             , only : icol_roof, icol_road_imperv
    use TridiagonalMod            , only : Tridiagonal
    use SoilStateType             , only : soilstate_type
    use SoilHydrologyType         , only : soilhydrology_type
    use VegetationType            , only : veg_pp
    use ColumnType                , only : col_pp
    use TopounitType              , only : top_pp

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,i                                ! indices
    integer  :: k,k_zwt
    real(r8) :: sat_lev
    real(r8) :: s1,s2,m,b   ! temporary variables used to interpolate theta
    integer  :: sat_flag

    !-----------------------------------------------------------------------

    associate(                                                & 
         nbedrock           =>    col_pp%nlevbed            , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)           
         dz                 =>    col_pp%dz                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         z                  =>    col_pp%z                  , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                 =>    col_pp%zi                 , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         h2osoi_liq         =>    col_ws%h2osoi_liq         , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_ice         =>    col_ws%h2osoi_ice         , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_vol         =>    col_ws%h2osoi_vol         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         watsat             =>    soilstate_vars%watsat_col , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
         zwt                =>    soilhydrology_vars%zwt_col  & ! Output: [real(r8) (:)   ]  water table depth (m)                             
         )

      ! calculate water table based on soil moisture state
      ! this is a simple search for 1st layer with soil moisture 
      ! less than specified threshold (sat_lev)

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         ! initialize to depth of bottom of lowest layer
         zwt(c)=zi(c,nlevgrnd)

         ! locate water table from bottom up starting at bottom of soil column
         ! sat_lev is an arbitrary saturation level used to determine water table
         sat_lev=0.96

         k_zwt=nlevgrnd
         sat_flag=1 !will remain unchanged if all layers at saturation
         do k=nlevgrnd,1,-1
            h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                 + h2osoi_ice(c,k)/(dz(c,k)*denice)

            if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
               k_zwt=k
               sat_flag=0
               exit
            endif
         enddo
         if (sat_flag == 1) k_zwt=1

         ! if soil column above sat_lev, set water table to lower 
         ! interface of first layer
         if (k_zwt == 1) then
            zwt(c)=zi(c,1)
         else if (k_zwt < nlevgrnd) then
            ! interpolate between k_zwt and k_zwt+1 to find water table height
            s1 = (h2osoi_liq(c,k_zwt)/(dz(c,k_zwt)*denh2o) &
                 + h2osoi_ice(c,k_zwt)/(dz(c,k_zwt)*denice))/watsat(c,k_zwt)
            s2 = (h2osoi_liq(c,k_zwt+1)/(dz(c,k_zwt+1)*denh2o) &
                 + h2osoi_ice(c,k_zwt+1)/(dz(c,k_zwt+1)*denice))/watsat(c,k_zwt+1)
            m=(z(c,k_zwt+1)-z(c,k_zwt))/(s2-s1)*1.0_r8
            b=z(c,k_zwt+1)-m*s2
            zwt(c)=max(0._r8,m*sat_lev+b)

         else
            zwt(c)=zi(c,nlevgrnd)
         endif
      end do

    end associate

  end subroutine ThetaBasedWaterTable

end module SoilWaterMovementLateralMod
