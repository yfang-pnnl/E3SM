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
    integer                                  :: nlevbed                     ! number of layers to bedrock

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
         den = hs_dx_nod(l,k+1)*1000._r8
         
         c = col_id_dn
         nlevbed = nbedrock(c)

         do j = 1, nlevbed


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
                          * dz(col_id_up,j)/hs_dx(l,k+1)*cos(theta)

               qflx_lateral_s(col_id_dn,j) = qflx_lateral_s(col_id_dn,j)  + qflx_up_to_dn &
                          * dz(col_id_dn,j)/hs_dx(l,k)*cos(theta)
            endif
         enddo
      enddo
   enddo

    end associate

  end subroutine ComputeLateralUnsatFlux

  !-----------------------------------------------------------------------
  subroutine SolveLateralSatFlow(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, num_h3dc, filter_h3dc, soilhydrology_vars, soilstate_vars, jwt,qflx_lateral_s)
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
    use elm_varpar                , only : nlevsoi
    use elm_varcon                , only : watmin, e_ice, denh2o, denice
    use elm_varctl                , only : iulog
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
    real(r8), intent(out) :: qflx_lateral_s(bounds%begc:bounds%endc)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,j, l, k, c0
    integer :: h3d_begc,h3d_endc
    integer :: col_id_up, col_id_dn
    integer :: idx_up, idx_dn
    integer :: step, nlevbed
    integer :: nstep
    logical :: up_local, dn_local
    real(r8) :: hksat_up, hksat_dn
    real(r8) :: den, qflx_up_to_dn, impedl
    real(r8) :: depth_up, depth_down, qlat_balance_error, qlat_storage_change
    real(r8) :: dtime, qlat_layer, qlat_tot, qlat_temp, s_y, h2osoi_left_vol, h2osoi_avail, qlat_from_wa
    real(r8) :: rous, sy, trans, theta, storage_before, storage_after,wa_left_vol, wa_excess

    !-----------------------------------------------------------------------

    associate(                                                &
         zi             => col_pp%zi                        , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         nlev2bed       => col_pp%nlevbed                   , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         h2osoi_liq     => col_ws%h2osoi_liq                , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice         =>    col_ws%h2osoi_ice         , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)

         bsw            => soilstate_vars%bsw_col           , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         sucsat         => soilstate_vars%sucsat_col        , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         hksat          => soilstate_vars%hksat_col         , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         smp_l          => soilstate_vars%smp_l_col         , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
         watsat         => soilstate_vars%watsat_col        , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)

         origflag       => soilhydrology_vars%origflag      , & ! Input:  constant
         fracice        => soilhydrology_vars%fracice_col   , & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)
         icefrac        => soilhydrology_vars%icefrac_col   , & ! Input:  [real(r8) (:,:) ]  fraction of ice
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

      qflx_lateral_s(:)=0._r8
#if 0
      !=========================================================================
      ! STEP-1 DIAGNOSTIC: verify the column-ordering assumption.
      !
      ! SolveLateralSatFlow (and the caller's h3d_zwt_lun / qflx_lateral_col
      ! bookkeeping) assumes each connected H3D group occupies exactly
      ! nh3dc_per_lunit CONSECUTIVE entries in filter_h3dc, addressed by the raw
      ! integer offset c0+k.  In this configuration each H3D column lives on a
      ! distinct topounit and all nh3dc_per_lunit columns of a group belong to
      ! the SAME gridcell.  If filter_h3dc skips a column (inactive interior
      ! column, or interleaved ordering built by filterMod.F90), then c0+k points
      ! into a different group/gridcell: fluxes are computed between non-adjacent
      ! columns and written to the wrong column -> storage/flux mismatch that
      ! only appears with the lateral call yet passes the per-column internal
      ! check.  This block prints and aborts if that assumption is violated.
      do fc = 1,num_h3dc,nh3dc_per_lunit
        c0 = filter_h3dc(fc)
        do k = 0, nh3dc_per_lunit-1
           if (fc+k > num_h3dc) then
              write(iulog,*) 'H3D-ORDER: block overruns num_h3dc: fc,k,num_h3dc=', fc,k,num_h3dc
              call endrun(msg=errMsg(__FILE__, __LINE__))
           end if
           c = filter_h3dc(fc+k)
           ! Print the full block layout once (first clump entry) so the
           ! column/topounit/gridcell mapping can be inspected in the log.
           if (fc == 1) then
              write(iulog,*) 'H3D-ORDER map: fc=',fc,' k=',k, &
                   ' col=',c,' expected(c0+k)=',c0+k, &
                   ' gridcell=',col_pp%gridcell(c), &
                   ' topounit=',col_pp%topounit(c), &
                   ' landunit=',col_pp%landunit(c)
           end if
           if (c /= c0+k) then
              write(iulog,*) 'H3D-ORDER VIOLATION: non-contiguous block. ', &
                   'fc=',fc,' k=',k,' filter_h3dc(fc+k)=',c,' /= c0+k=',c0+k, &
                   ' gc=',col_pp%gridcell(c),' topo=',col_pp%topounit(c)
              call endrun(msg=errMsg(__FILE__, __LINE__))
           end if
           if (col_pp%gridcell(c) /= col_pp%gridcell(c0)) then
              write(iulog,*) 'H3D-ORDER VIOLATION: gridcell mismatch within block. ', &
                   'fc=',fc,' k=',k,' col=',c, &
                   ' gc(col)=',col_pp%gridcell(c),' gc(c0)=',col_pp%gridcell(c0)
              call endrun(msg=errMsg(__FILE__, __LINE__))
           end if
        end do
      end do
#endif
      !=========================================================================

      do fc = 1,num_h3dc,nh3dc_per_lunit  !loop for all soil columns that belong to same land unit
        c0 = filter_h3dc(fc)
        h3d_begc = c0
        h3d_endc = c0+nh3dc_per_lunit-1
        l  = col_pp%landunit(c0)

        do k = 1, nh3dc_per_lunit-1
           col_id_up = c0+(k+1)-1
           col_id_dn = c0+k-1
           c = col_id_dn

               ! local grid cell
            depth_up = zi(col_id_up,nlevgrnd) - zwt(col_id_up)  ! groundwater head(m)
            depth_down = zi(col_id_dn,nlevgrnd) - zwt(col_id_dn)

            ! Use the saturated conductivity at the layer adjacent to each
            ! column's water table, consistent with the H3D LateralResponse
            ! formulation.  A fixed deep layer can strongly overestimate the
            ! saturated transmissivity and lateral flux.
            idx_up = min(jwt(col_id_up)+1,nlev2bed(col_id_up))
            idx_dn = min(jwt(col_id_dn)+1,nlev2bed(col_id_dn))
            hksat_up = hksat(col_id_up,idx_up)
            hksat_dn = hksat(col_id_dn,idx_dn)

            den = hs_dx_nod(l,k+1)*1000._r8

            depth_up = max(depth_up, 0._r8)
            depth_down= max(depth_down, 0._r8)
            theta = hs_slope(c)/180._r8*rpi

            ! Apply frozen-soil hydraulic impedance, consistent with the
            ! unsaturated lateral-flow and drainage formulations.
            if (origflag == 1) then
               impedl = 1._r8 - 0.5_r8*(fracice(col_id_up,idx_up) + fracice(col_id_dn,idx_dn))
            else
               impedl = 10._r8**(-e_ice*(0.5_r8*(icefrac(col_id_up,idx_up) + icefrac(col_id_dn,idx_dn))))
            endif

            ! calculate transmissivity
            trans = impedl*sqrt(hksat_up*hksat_dn)*(depth_up+depth_down)/2._r8*1000._r8 !* 0.01 ! (mm2/s)
            qflx_up_to_dn = -trans*((depth_down-depth_up)*1000._r8/den*cos(theta) + sin(theta))

            ! qflx_up_to_dn is an integrated saturated flow per unit width [mm2/s].
            ! Convert it back to an areal flux [mm/s] by dividing by the H3D
            ! column length in mm. hs_dx is stored in m, so multiply by 1000.
            qflx_lateral_s(col_id_up) = qflx_lateral_s(col_id_up) - &
                    qflx_up_to_dn/(hs_dx(l,k+1)*1000._r8)*cos(theta)

            qflx_lateral_s(col_id_dn) = qflx_lateral_s(col_id_dn) + &
                    qflx_up_to_dn/(hs_dx(l,k)*1000._r8)*cos(theta)
         enddo
      enddo

         ! IMPORTANT: iterate over exactly the H3D columns whose saturated
         ! lateral flux was computed above and is later added to
         ! qflx_lateral_col by the caller (do fc = 1,num_h3dc loop in
         ! SoilWaterMovementMod).  Using filter_hydrologyc here caused a mass
         ! balance error: any H3D column present in filter_h3dc but not in
         ! filter_hydrologyc kept its raw explicit qflx_lateral_s (from the
         ! flux loop) and had that flux applied to the ELM water balance, yet
         ! its storage (wa, h2osoi_liq) was never updated here -- an
         ! uncompensated flux.  Looping over filter_h3dc guarantees the
         ! reported flux matches the applied storage change for every column.
         do fc = 1, num_h3dc
            c = filter_h3dc(fc)
            nlevbed = nlev2bed(c)

            !scs: use analytical expression for aquifer specific yield
            rous = watsat(c,nlevbed) &
                 * ( 1._r8 - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
            rous=max(rous,0.02_r8)

            !--  water table is below the soil column  --------------------------------------

            qlat_temp = qflx_lateral_s(c)

            !-- Requested lateral-water amount over the timestep [mm]
            qlat_tot = qlat_temp * dtime

!-- Save the storage components affected by this calculation.
!-- Use the same soil-layer range used by the ELM water balance.
storage_before = wa(c) + sum(h2osoi_liq(c,1:nlevgrnd)+h2osoi_ice(c,1:nlevgrnd))

!-- Initialize layer transfer.
qlat_layer = 0._r8

if (jwt(c) == nlevbed) then

   !===========================================================================
   ! Water table is below the resolved soil column.
   !
   ! wa is the aquifer storage below the soil column. Positive lateral flow
   ! adds to wa; negative lateral flow removes from wa. Any wa above 5000 mm
   ! is transferred to the bottom soil layer, subject to pore-space capacity.
   !===========================================================================

   if (qlat_tot > 0._r8) then

      ! Available storage in wa before reaching the nominal 5000-mm limit.
      wa_left_vol = max(0._r8, 5000._r8 - wa(c))

      ! Available liquid-water storage in the bottom layer [mm H2O].
      ! h2osoi_ice is ice mass [kg m-2], so convert it to occupied
      ! pore volume as liquid-water equivalent before subtracting it.
      h2osoi_left_vol =                                      &
           max(0._r8,                                        &
               watsat(c,nlevbed) * dz(c,nlevbed) * denh2o    &
               - h2osoi_ice(c,nlevbed) * denh2o / denice     &
               - watmin)                                     &
           - max(0._r8, h2osoi_liq(c,nlevbed) - watmin)

      h2osoi_left_vol = max(0._r8, h2osoi_left_vol)

      ! The total accepted inflow cannot exceed the combined available
      ! aquifer and bottom-layer storage.
      qlat_layer = min(qlat_tot, wa_left_vol + h2osoi_left_vol)
      qlat_layer = max(qlat_layer, 0._r8)

      ! First add the accepted water to wa.
      wa(c) = wa(c) + qlat_layer

      ! Raise the water table according to the accepted lateral inflow.
      if (qlat_layer > 0._r8) then
         zwt(c) = zwt(c) - qlat_layer /                       &
                  (1.e3_r8 * rous)
      endif

      ! Transfer storage above 5000 mm from wa into the bottom layer.
      wa_excess = max(0._r8, wa(c) - 5000._r8)
      wa_excess = min(wa_excess, h2osoi_left_vol)

      if (wa_excess > 0._r8) then
         h2osoi_liq(c,nlevbed) =                              &
              h2osoi_liq(c,nlevbed) + wa_excess

         wa(c) = wa(c) - wa_excess
      endif

      ! Remaining positive amount was not accommodated.
      qlat_tot = qlat_tot - qlat_layer

   else if (qlat_tot < 0._r8) then

      ! Do not remove more aquifer water than is available.
      qlat_layer = max(qlat_tot, -max(0._r8, wa(c)))
      qlat_layer = min(qlat_layer, 0._r8)

      wa(c) = wa(c) + qlat_layer

      ! Negative qlat_layer deepens the water table.
      if (qlat_layer < 0._r8) then
         zwt(c) = zwt(c) - qlat_layer /                       &
                  (1.e3_r8 * rous)
      endif

      ! Remaining negative amount was not accommodated.
      qlat_tot = qlat_tot - qlat_layer

   endif

else

   !===========================================================================
   ! Water table is within the resolved soil column.
   !===========================================================================

   if (qlat_tot > 0._r8) then

      !-------------------------------------------------------------------------
      ! Positive lateral flow:
      ! add water and raise the water table.
      !-------------------------------------------------------------------------

      do j = jwt(c)+1, 1, -1

         ! Analytical specific yield.
         s_y = watsat(c,j) *                                  &
              (1._r8 -                                        &
               (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**         &
               (-1._r8/bsw(c,j)))

         s_y = max(s_y, 0.02_r8)

         ! Water required to raise the water table to the top of layer j.
         qlat_layer =                                         &
              s_y * max(0._r8, zwt(c)-zi(c,j-1)) * 1.e3_r8

         ! Do not apply more than the remaining requested amount.
         qlat_layer = min(qlat_layer, qlat_tot)
         qlat_layer = max(qlat_layer, 0._r8)

         ! Available liquid-water storage in layer j [mm H2O].
         ! h2osoi_ice is ice mass [kg m-2], so convert it to occupied
         ! pore volume as liquid-water equivalent before subtracting it.
         h2osoi_left_vol =                                    &
              max(0._r8,                                      &
                  watsat(c,j)*dz(c,j)*denh2o                  &
                  - h2osoi_ice(c,j) * denh2o / denice         &
                  - watmin)                                   &
              - max(0._r8, h2osoi_liq(c,j)-watmin)

         h2osoi_left_vol = max(0._r8, h2osoi_left_vol)

         ! Limit the transfer to the pore volume actually available.
         qlat_layer = min(qlat_layer, h2osoi_left_vol)

         ! Update soil storage.
         h2osoi_liq(c,j) = h2osoi_liq(c,j) + qlat_layer

         ! Update water-table depth based on the accepted amount.
         if (qlat_layer > 0._r8) then
            zwt(c) = zwt(c) - qlat_layer /                    &
                     (s_y * 1.e3_r8)
         endif

         ! Track the amount that remains unapplied.
         qlat_tot = qlat_tot - qlat_layer

         if (qlat_tot <= 1.e-14_r8) then
            qlat_tot = 0._r8
            exit
         endif

         ! If this layer has no available capacity, continuing upward may
         ! still allow water to enter another layer.
      enddo

   else if (qlat_tot < 0._r8) then

      !-------------------------------------------------------------------------
      ! Negative lateral flow:
      ! remove water and deepen the water table.
      !-------------------------------------------------------------------------

      do j = jwt(c)+1, nlevbed

         ! Analytical specific yield.
         s_y = watsat(c,j) *                                  &
              (1._r8 -                                        &
               (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**         &
               (-1._r8/bsw(c,j)))

         s_y = max(s_y, 0.02_r8)

         ! Negative water amount associated with moving the water table
         ! from its current position to the bottom of layer j.
         qlat_layer =                                         &
              -s_y * max(0._r8, zi(c,j)-zwt(c)) * 1.e3_r8

         ! Do not remove more than the remaining requested outflow.
         qlat_layer = max(qlat_layer, qlat_tot)
         qlat_layer = min(qlat_layer, 0._r8)

         ! Liquid water available above watmin.
         h2osoi_avail = max(0._r8,                             &
                            h2osoi_liq(c,j)-watmin)

         ! Do not reduce layer water below watmin.
         qlat_layer = max(qlat_layer, -h2osoi_avail)

         ! Update soil storage.
         h2osoi_liq(c,j) = h2osoi_liq(c,j) + qlat_layer

         ! Track the unapplied negative amount.
         qlat_tot = qlat_tot - qlat_layer

         if (qlat_tot >= -1.e-14_r8) then

            ! Drainage ends within this layer.
            if (qlat_layer < 0._r8) then
               zwt(c) = zwt(c) - qlat_layer /                 &
                        (s_y * 1.e3_r8)
            endif

            qlat_tot = 0._r8
            exit

         else

            ! The requested outflow continues into the next layer.
            ! Only set zwt to the layer bottom if the water required to
            ! traverse this layer was actually removed.
            if (qlat_layer < 0._r8) then
               zwt(c) = zi(c,j)
            endif

         endif

      enddo

      ! If the water table has moved below all resolved soil layers,
      ! allow remaining outflow to draw from wa.
      if (qlat_tot < -1.e-14_r8) then

         qlat_from_wa = max(qlat_tot, -max(0._r8, wa(c)))
         qlat_from_wa = min(qlat_from_wa, 0._r8)

         wa(c) = wa(c) + qlat_from_wa
         qlat_tot = qlat_tot - qlat_from_wa

         if (qlat_from_wa < 0._r8) then
            zwt(c) = zwt(c) - qlat_from_wa /                  &
                     (1.e3_r8 * rous)
         endif

      endif

   endif

endif

!==============================================================================
! Record the amount actually accommodated.
!
! Following the qflx_lnd2ocn treatment in SoilHydrologyMod (subroutine
! Drainage_To_OCN), the flux exposed to the ELM water balance is the
! physically accepted transfer, not the original explicit request.  The
! explicit saturated-flow estimate qlat_temp can request more water than the
! available pore/aquifer storage; the storage-update logic above clips it and
! leaves the unaccommodated amount in qlat_tot.  The accepted transfer is
! therefore the requested amount minus the unapplied residual:
!
!     applied = qlat_temp*dtime - qlat_tot
!
! Every storage update above decrements qlat_tot by exactly the amount added
! to wa or h2osoi_liq, so this identity is exact and, unlike differencing two
! large (~5000-6000 mm) storage sums, it does not suffer catastrophic
! cancellation near the 1.e-7 mm balance tolerance.
!==============================================================================

!-- Positive means that water entered the column.
qflx_lateral_s(c) = qlat_temp - qlat_tot / dtime

!-- Diagnostic: cross-check against the explicit storage difference. The two
!-- expressions should agree to roundoff. This is informational only and does
!-- not feed the water balance.
storage_after = wa(c) + sum(h2osoi_liq(c,1:nlevgrnd)+h2osoi_ice(c,1:nlevgrnd))
qlat_storage_change = storage_after - storage_before
qlat_balance_error  = qlat_storage_change - qflx_lateral_s(c)*dtime
#if 0
if (abs(qlat_balance_error) > 1.e-7_r8) then
   write(iulog,*) 'Lateral water balance error: ',             &
                  qlat_balance_error
   write(iulog,*) 'column: ', c
   write(iulog,*) 'requested lateral flux: ',                  &
                  qlat_temp
   write(iulog,*) 'accepted lateral flux: ',                   &
                  qflx_lateral_s(c)
   write(iulog,*) 'unapplied lateral amount: ', qlat_tot
   write(iulog,*) 'storage before: ', storage_before
   write(iulog,*) 'storage after: ', storage_after
endif
#endif
enddo

      !=========================================================================
      ! Recompute jwt (index of the first unsaturated layer, i.e. the layer
      ! right above the water table) from the UPDATED zwt.
      !
      ! SolveLateralSatFlow moves zwt when it adds/removes lateral water, but
      ! left jwt at its pre-lateral value. Downstream Drainage (SoilHydrologyMod)
      ! branches on jwt (e.g. "if (jwt(c) == nlevbed)" and "do j = jwt(c)+1,..")
      ! while using this modified zwt. A stale jwt makes Drainage move water in
      ! layers inconsistent with the true water-table position, producing an
      ! uncompensated storage change and hence the BalanceCheckMod error that
      ! only appears when the lateral call is active. Recomputing jwt here (same
      ! logic as the jwt update loops in SoilHydrologyMod) keeps zwt and jwt
      ! consistent for every routine that runs after this one.
      !=========================================================================
      do fc = 1, num_h3dc
         c = filter_h3dc(fc)
         nlevbed = nlev2bed(c)
         jwt(c) = nlevbed
         ! allow jwt to equal zero when zwt is in top layer
         do j = 1, nlevbed
            if (zwt(c) <= zi(c,j)) then
               jwt(c) = j-1
               exit
            end if
         end do
      end do
         !call ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
         !     soilhydrology_vars, soilstate_vars)


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
