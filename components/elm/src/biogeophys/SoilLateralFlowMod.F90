module SoilLateralFlowMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lateral subsurface flow between grid cells, following Qiu et al. (2024),
  ! "Development of inter-grid-cell lateral unsaturated and saturated flow
  ! model in the E3SM Land Model (v2.0)", GMD 17, 143-167.
  !
  ! Lateral fluxes are computed explicitly from the state at the start of the
  ! soil-water solve:
  !
  ! - Unsaturated flow (Eq. 15): for each soil layer that lies entirely above
  !   the water table in both columns, a Darcy flux driven by the difference
  !   in matric potential plus elevation between the two columns. It enters
  !   the vertical (Zeng and Decker, 2009) Richards equation as a layer
  !   source term.
  ! - Saturated flow (Eq. 17): a Darcy flux driven by the difference in water
  !   table head (surface elevation minus water table depth), using a
  !   transmissivity based on the mean saturated thickness of the two columns.
  !   Each column's net saturated inflow is applied after the Richards solve
  !   at the water table (ApplySaturatedLateralFlux), as ELM does for
  !   baseflow: outflow drains the layers at and below the water table in
  !   turn, limited by the specific yield, and inflow fills the layers at and
  !   above the water table in turn. The layers below the water table stay
  !   saturated; the water table moves instead. (Zeng and Decker (2009) has
  !   no positive pressure below the water table, so a sink applied there
  !   inside the solve would desaturate the deep layers.)
  !
  ! Fluxes are first computed per connection and then summed for each column
  ! in the order of the neighbors' natural ids, so results do not depend on
  ! the number of MPI ranks.
  !
  ! Connections join the naturally vegetated columns of neighboring grid cells
  ! (one per grid cell). Water is conserved for any subgrid structure:
  !
  ! - Saturated flow treats the soil and crop columns of a grid cell as one
  !   aquifer: the water table head and transmissivity of a cell are the
  !   area-weighted means over its soil and crop columns, and the cell's net
  !   saturated inflow is applied at the same rate [mm/s] to each of them,
  !   i.e. divided by (cell area * summed weight of those columns).
  ! - Unsaturated flow uses the state of the naturally vegetated column only and
  !   is applied to that column, divided by (cell area * its weight).
  ! A connection is inactive if either cell's aquifer (or, for unsaturated
  ! flow, naturally vegetated) weight is below min_lateral_wt, to avoid very
  ! large per-area rates in nearly empty columns. Cells without a naturally
  ! vegetated column have no connections.
  !
  ! Saturated transmissivity: lateral_trans_method = 'elm' (Qiu et al., 2024;
  ! bottom-layer hksat times the saturated thickness inside the column) or
  ! 'asap' (Fan et al., 2007; Miguez-Macho et al., 2007; ASAP LATERALFLOW4:
  ! T = K*(fdepth + 1.5 - zwt) for zwt <= 1.5 m and K*fdepth*exp((1.5-zwt)/fdepth)
  ! below, with K = lateral_hk_anisotropy * hksat_lat_ref; arithmetic mean of
  ! the two cells). ASAP has saturated lateral flow only, so
  ! lateral_unsat_flow = .false. is recommended with 'asap'.
  !
  ! ComputeRiverGWExchange (use_gw2river) ports ASAP GW2RIVER: groundwater
  ! drains to the river when the water table is above the river surface, and
  ! the river feeds the aquifer otherwise (only with river water from MOSART,
  ! use_lnd_rof_gw_exchange). The exchange is added to the net saturated
  ! source applied by ApplySaturatedLateralFlux.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use elm_varpar        , only : nlevgrnd
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_ws, col_wf
  use SoilStateType     , only : soilstate_type
  use SoilHydrologyType , only : soilhydrology_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ComputeLateralFlux         ! lateral fluxes into owned columns
  public :: ComputeRiverGWExchange     ! groundwater-river exchange (ASAP GW2RIVER)
  public :: ApplySaturatedLateralFlux  ! apply the saturated lateral flux at the water table
  public :: ThetaBasedWaterTable    ! diagnose the water table from soil moisture
  !
  ! !PRIVATE DATA:
  real(r8), parameter :: sat_lev = 0.96_r8   ! saturation level that defines the theta-based water table
  real(r8), parameter :: min_lateral_wt = 0.05_r8 ! minimum column weight (of the cell) for a lateral connection
  real(r8), allocatable :: riv_seep_cum_grc(:)   ! river-to-aquifer seepage since the last MOSART coupling [mm, cell mean]
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ComputeLateralFlux(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, soilstate_vars, qflx_lat_layer, qflx_sat, zwt_lat)
    !
    ! !DESCRIPTION:
    ! Computes, for the columns in filter_hydrologyc, the unsaturated lateral
    ! flux into each soil layer and the net saturated lateral flux into the
    ! column [mm H2O/s]. Sets col_wf%qflx_lateral (outflow positive, both
    ! parts), qflx_lateral_unsat, qflx_lateral_sat and qflx_lat_layer (the
    ! unsaturated part; ApplySaturatedLateralFlux adds the saturated part).
    !
    ! The fluxes use a water table diagnosed from the current soil moisture
    ! (zwt_lat), not zwt itself: after the previous time step's water-table
    ! diagnosis, Drainage moves zwt (by the specific yield of the drained
    ! water, and caps it at 80 m), so zwt can show a saturated zone that the
    ! soil moisture does not have (e.g. zwt = 80 m in an unsaturated 82 m deep
    ! column).
    !
    ! Must be called on all MPI ranks (it performs a halo exchange), with one
    ! clump per MPI rank: 'bounds' (clump bounds) must cover all owned columns.
    ! Ghost columns are addressed through the processor bounds.
    !
#ifdef MOAB_LATERAL
    use decompMod               , only : get_proc_bounds
    use elm_varctl              , only : lateral_unsat_flow, lateral_hk_anisotropy, lateral_theta_watertable
    use elm_varctl              , only : lateral_trans_method
    use elm_varcon              , only : denh2o, denice
    use elm_varcon              , only : e_ice
    use landunit_varcon         , only : istsoil
    use ColumnConnectionSetType , only : c2c_connections
    use domainLateralMod        , only : NatVegColumnRealDataHaloExchange
#endif
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    real(r8)                 , intent(out)   :: qflx_lat_layer(bounds%begc:,1:)  ! unsaturated lateral flux into each soil layer [mm H2O/s]
    real(r8)                 , intent(out)   :: qflx_sat(bounds%begc:)           ! net saturated lateral flux into each column [mm H2O/s]
    real(r8)                 , intent(out)   :: zwt_lat(bounds%begc:)            ! water table depth used for the lateral fluxes [m]
    !
#ifdef MOAB_LATERAL
    ! !LOCAL VARIABLES:
    type(bounds_type)     :: bp              ! processor bounds (include ghost columns)
    integer               :: fc, c, g, j, k, iconn, cu, cd, nbu, nbd, nb, nconn, ncomp, cn, ic
    real(r8)              :: s1, bswl, hkl, impedl, smp_grad, face_area, area
    real(r8)              :: trans, sgn, wt, kl, tcol
    logical               :: asap
    real(r8), allocatable :: data(:,:)
    integer , allocatable :: jwt(:)
    integer , allocatable :: natcol(:)       ! naturally vegetated column of each owned grid cell (-1 if none)
    real(r8), allocatable :: zw(:)           ! water table depth used for the fluxes, owned and ghost columns [m]
    ! grid cell (aquifer) values, stored at the cell's naturally vegetated column; owned and ghost
    real(r8), allocatable :: zwc(:)          ! area-weighted water table depth of the soil and crop columns [m]
    real(r8), allocatable :: hsatc(:)        ! area-weighted saturated thickness inside the column ('elm') [m]
    real(r8), allocatable :: hkbc(:)         ! area-weighted bottom-layer hksat ('elm') [mm/s]
    real(r8), allocatable :: tc(:)           ! area-weighted transmissivity ('asap') [mm/s * m]
    real(r8), allocatable :: wcell(:)        ! summed weight (of the cell) of the soil and crop columns [-]
    real(r8), allocatable :: wnat(:)         ! weight (of the cell) of the naturally vegetated column [-]
    real(r8), allocatable :: sat_cell(:)     ! net saturated inflow per unit aquifer area, owned cells [mm H2O/s]
    real(r8)              :: vol(nlevgrnd)   ! volumetric water content [m3/m3]
    real(r8), allocatable :: qv_unsat(:,:)   ! [nconn, nlevgrnd] up-to-down volumetric flux per layer [mm H2O/s * m^2]
    real(r8), allocatable :: qv_sat(:)       ! [nconn] up-to-down volumetric saturated flux [mm H2O/s * m^2]
    real(r8)              :: unsat_in(nlevgrnd), sat_in
    integer , parameter   :: nextra = 7      ! zw, zwc, hsatc, hkbc, tc, wcell, wnat
    !-----------------------------------------------------------------------

    associate(                                                  &
         zi           => col_pp%zi                            , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         dz           => col_pp%dz                            , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         nlevbed      => col_pp%nlevbed                       , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_vol   => col_ws%h2osoi_vol                    , & ! Input:  [real(r8) (:,:) ] volumetric soil water (m3/m3)
         smp_l        => soilstate_vars%smp_l_col             , & ! Input:  [real(r8) (:,:) ] soil matric potential (mm)
         watsat       => soilstate_vars%watsat_col            , & ! Input:  [real(r8) (:,:) ] porosity
         hksat        => soilstate_vars%hksat_col             , & ! Input:  [real(r8) (:,:) ] saturated hydraulic conductivity (mm/s)
         bsw          => soilstate_vars%bsw_col               , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         origflag     => soilhydrology_vars%origflag          , & ! Input:  constant
         fracice      => soilhydrology_vars%fracice_col       , & ! Input:  [real(r8) (:,:) ] fractional impermeability
         icefrac      => soilhydrology_vars%icefrac_col       , & ! Input:  [real(r8) (:,:) ] fraction of ice
         zwt          => soilhydrology_vars%zwt_col           , & ! Input:  [real(r8) (:)   ] water table depth (m)
         conn         => c2c_connections                        &
         )

      call get_proc_bounds(bp)
      if (bounds%begc /= bp%begc .or. bounds%endc /= bp%endc) then
         call endrun(msg='ComputeLateralFlux: lateral flow requires one clump per MPI rank'//errMsg(__FILE__, __LINE__))
      end if
      asap = (trim(lateral_trans_method) == 'asap')

      ! --- Water table of owned columns, diagnosed from the current soil moisture
      allocate(zw(bp%begc_all:bp%endc_all))
      zw(:) = 0._r8
      do c = bp%begc, bp%endc
         zw(c) = zwt(c)
      end do
      do fc = 1, num_hydrologyc
         c  = filter_hydrologyc(fc)
         nb = nlevbed(c)
         if (lateral_theta_watertable) then
            do j = 1, nb
               vol(j) = col_ws%h2osoi_liq(c,j)/(dz(c,j)*denh2o) + col_ws%h2osoi_ice(c,j)/(dz(c,j)*denice)
            end do
            zw(c) = ThetaWaterTableDepth(nb, vol(1:nb), watsat(c,1:nb), col_pp%z(c,1:nb), zi(c,0:nb))
         else
            zw(c) = min(zwt(c), zi(c,nb))
         end if
      end do

      ! --- Grid cell (aquifer) values: area-weighted means over the soil and
      !     crop columns of each owned cell, stored at its naturally vegetated column
      allocate(natcol(bp%begg:bp%endg)) ; natcol(:) = -1
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (col_pp%itype(c) == istsoil) natcol(col_pp%gridcell(c)) = c
      end do

      allocate(zwc  (bp%begc_all:bp%endc_all)) ; zwc(:)   = 0._r8
      allocate(hsatc(bp%begc_all:bp%endc_all)) ; hsatc(:) = 0._r8
      allocate(hkbc (bp%begc_all:bp%endc_all)) ; hkbc(:)  = 0._r8
      allocate(tc   (bp%begc_all:bp%endc_all)) ; tc(:)    = 0._r8
      allocate(wcell(bp%begc_all:bp%endc_all)) ; wcell(:) = 0._r8
      allocate(wnat (bp%begc_all:bp%endc_all)) ; wnat(:)  = 0._r8

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
         cn = natcol(col_pp%gridcell(c))
         if (cn < 0) cycle
         nb = nlevbed(c)
         wt = col_pp%wtgcell(c)
         if (asap) then
            ! ASAP LATERALFLOW4 kcell, with wtd = -zw
            kl = lateral_hk_anisotropy * soilstate_vars%hksat_lat_ref_col(c)
            if (zw(c) <= 1.5_r8) then
               tcol = kl * (soilstate_vars%fdepth_col(c) + 1.5_r8 - zw(c))
            else
               tcol = kl * soilstate_vars%fdepth_col(c) * exp((1.5_r8 - zw(c))/soilstate_vars%fdepth_col(c))
            end if
            tc(cn) = tc(cn) + wt*tcol
         end if
         wcell(cn) = wcell(cn) + wt
         zwc(cn)   = zwc(cn)   + wt*zw(c)
         hsatc(cn) = hsatc(cn) + wt*max(zi(c,nb) - zw(c), 0._r8)
         hkbc(cn)  = hkbc(cn)  + wt*hksat(c,nb)
      end do
      do g = bp%begg, bp%endg
         cn = natcol(g)
         if (cn < 0) cycle
         wnat(cn) = col_pp%wtgcell(cn)
         if (wcell(cn) > 0._r8) then
            zwc(cn)   = zwc(cn)  /wcell(cn)
            hsatc(cn) = hsatc(cn)/wcell(cn)
            hkbc(cn)  = hkbc(cn) /wcell(cn)
            tc(cn)    = tc(cn)   /wcell(cn)
         end if
      end do

      ! --- Copy the state of ghost columns from the owning ranks (one MPI round)
      ! Components: h2osoi_vol, smp_l, fracice, icefrac (each 1:nlevgrnd), water table,
      ! and the grid cell values zwc, hsatc, hkbc, tc, wcell, wnat
      ncomp = 4*nlevgrnd + nextra
      ic    = 4*nlevgrnd
      allocate(data(bp%begc_all:bp%endc_all, ncomp))
      data(:,:) = 0._r8
      do c = bp%begc, bp%endc
         do j = 1, nlevgrnd
            data(c, 0*nlevgrnd + j) = h2osoi_vol(c,j)
            data(c, 1*nlevgrnd + j) = smp_l(c,j)
            data(c, 2*nlevgrnd + j) = fracice(c,j)
            data(c, 3*nlevgrnd + j) = icefrac(c,j)
         end do
         data(c, ic+1) = zw(c)
         data(c, ic+2) = zwc(c)
         data(c, ic+3) = hsatc(c)
         data(c, ic+4) = hkbc(c)
         data(c, ic+5) = tc(c)
         data(c, ic+6) = wcell(c)
         data(c, ic+7) = wnat(c)
      end do

      call NatVegColumnRealDataHaloExchange(bp, 'lateral_flow_state', data)

      do c = bp%endc + 1, bp%endc_all
         do j = 1, nlevgrnd
            h2osoi_vol(c,j) = data(c, 0*nlevgrnd + j)
            smp_l(c,j)      = data(c, 1*nlevgrnd + j)
            fracice(c,j)    = data(c, 2*nlevgrnd + j)
            icefrac(c,j)    = data(c, 3*nlevgrnd + j)
         end do
         zw(c)    = data(c, ic+1)
         zwc(c)   = data(c, ic+2)
         hsatc(c) = data(c, ic+3)
         hkbc(c)  = data(c, ic+4)
         tc(c)    = data(c, ic+5)
         wcell(c) = data(c, ic+6)
         wnat(c)  = data(c, ic+7)
      end do
      deallocate(data)

      ! --- Index of the deepest layer that lies entirely above the water table
      allocate(jwt(bp%begc_all:bp%endc_all))
      do c = bp%begc_all, bp%endc_all
         jwt(c) = 0
         if (nlevbed(c) <= 0 .or. nlevbed(c) > nlevgrnd) cycle
         jwt(c) = nlevbed(c)
         do j = 1, nlevbed(c)
            if (zw(c) <= zi(c,j)) then
               jwt(c) = j - 1
               exit
            end if
         end do
      end do

      ! --- Fluxes per connection (positive from the up to the down column)
      nconn = conn%nconn
      allocate(qv_unsat(nconn, nlevgrnd)) ; qv_unsat(:,:) = 0._r8
      allocate(qv_sat(nconn))             ; qv_sat(:)     = 0._r8

      do iconn = 1, nconn
         cu  = conn%col_id_up(iconn)
         cd  = conn%col_id_dn(iconn)
         nbu = nlevbed(cu)
         nbd = nlevbed(cd)

         ! unsaturated flux (Eq. 15) in layers above the water table of both
         ! naturally vegetated columns
         if (lateral_unsat_flow .and. wnat(cu) >= min_lateral_wt .and. wnat(cd) >= min_lateral_wt) then
            do j = 1, min(jwt(cu), jwt(cd))
               s1 = (h2osoi_vol(cu,j) + h2osoi_vol(cd,j)) / (watsat(cu,j) + watsat(cd,j))
               s1 = min(1._r8, max(s1, 0._r8))

               bswl = 0.5_r8*(bsw(cu,j) + bsw(cd,j))

               if (origflag == 1) then
                  impedl = 1._r8 - 0.5_r8*(fracice(cu,j) + fracice(cd,j))
               else
                  impedl = 10._r8**(-e_ice*(0.5_r8*(icefrac(cu,j) + icefrac(cd,j))))
               end if

               ! horizontal conductivity: geometric mean of the two columns
               hkl = lateral_hk_anisotropy * impedl * sqrt(hksat(cu,j)*hksat(cd,j)) &
                    * s1**(2._r8*bswl + 3._r8)

               ! gradient of (matric potential + elevation) along the land surface [-]
               smp_grad = (smp_l(cd,j) - smp_l(cu,j) + conn%dzg(iconn)*1.e3_r8) / (conn%dist(iconn)*1.e3_r8)

               face_area = conn%face_length(iconn) * 0.5_r8*(dz(cu,j) + dz(cd,j))   ! [m^2]

               qv_unsat(iconn, j) = -hkl * smp_grad * face_area
            end do
         end if

         ! saturated flux (Eq. 17) driven by the difference in water table head
         ! of the two cells' aquifers
         if (wcell(cu) < min_lateral_wt .or. wcell(cd) < min_lateral_wt) cycle

         if (asap) then
            ! ASAP: (kcell_i + kcell_j)*0.5 (arithmetic mean)
            trans = 0.5_r8*(tc(cu) + tc(cd))                                     ! [mm/s * m]
         else
            ! saturated thickness inside the column times the bottom-layer hksat
            trans = lateral_hk_anisotropy * sqrt(hkbc(cu)*hkbc(cd)) &
                 * 0.5_r8*(hsatc(cu) + hsatc(cd))                                ! [mm/s * m]
         end if

         ! head difference [m] = (elev_dn - zwt_dn) - (elev_up - zwt_up)
         qv_sat(iconn) = -trans * (conn%dzg(iconn) - zwc(cd) + zwc(cu)) / conn%dist(iconn) &
              * conn%face_length(iconn)
      end do

      ! --- Sum the fluxes for each owned cell, in a decomposition-independent order.
      !     Unsaturated: into the naturally vegetated column, per unit area of that
      !     column. Saturated: per unit area of the cell's soil and crop columns.
      qflx_lat_layer(bounds%begc:bounds%endc, :) = 0._r8
      qflx_sat(bounds%begc:bounds%endc)          = 0._r8
      zwt_lat(bounds%begc:bounds%endc)           = zw(bounds%begc:bounds%endc)

      allocate(sat_cell(bp%begc:bp%endc)) ; sat_cell(:) = 0._r8

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         col_wf%qflx_lateral(c)       = 0._r8
         col_wf%qflx_lateral_unsat(c) = 0._r8
         col_wf%qflx_lateral_sat(c)   = 0._r8
         col_wf%qflx_lat_layer(c,:)   = 0._r8

         if (col_pp%itype(c) /= istsoil) cycle
         if (conn%col_nconn(c) == 0) cycle

         nb          = nlevbed(c)
         unsat_in(:) = 0._r8
         sat_in      = 0._r8

         do k = conn%col_conn_beg(c), conn%col_conn_beg(c) + conn%col_nconn(c) - 1
            iconn = conn%col_conn_id(k)
            sgn   = conn%col_conn_sign(k)    ! +1: c is the down column (up-to-down flux enters c)
            if (sgn > 0._r8) then
               area = conn%downarea(iconn)
            else
               area = conn%uparea(iconn)
            end if

            if (wnat(c) >= min_lateral_wt) then
               do j = 1, nb
                  unsat_in(j) = unsat_in(j) + sgn*qv_unsat(iconn, j)/(area*wnat(c))
               end do
            end if
            if (wcell(c) >= min_lateral_wt) then
               sat_in = sat_in + sgn*qv_sat(iconn)/(area*wcell(c))
            end if
         end do

         do j = 1, nb
            qflx_lat_layer(c,j)          = unsat_in(j)
            col_wf%qflx_lat_layer(c,j)   = unsat_in(j)
            col_wf%qflx_lateral_unsat(c) = col_wf%qflx_lateral_unsat(c) + unsat_in(j)
         end do
         sat_cell(c) = sat_in
      end do

      ! the cell's saturated inflow goes to each of its soil and crop columns at the same rate
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
         cn = natcol(col_pp%gridcell(c))
         if (cn < 0) cycle
         qflx_sat(c)                = sat_cell(cn)
         col_wf%qflx_lateral_sat(c) = sat_cell(cn)
         col_wf%qflx_lateral(c)     = -(col_wf%qflx_lateral_unsat(c) + sat_cell(cn))
      end do

      deallocate(jwt, zw, qv_unsat, qv_sat, natcol, zwc, hsatc, hkbc, tc, wcell, wnat, sat_cell)

    end associate
#else
    call endrun(msg='ComputeLateralFlux: requires ELM built with -DMOAB_LATERAL'//errMsg(__FILE__, __LINE__))
#endif

  end subroutine ComputeLateralFlux

  !-----------------------------------------------------------------------
  subroutine ComputeRiverGWExchange(bounds, num_hydrologyc, filter_hydrologyc, dtime, &
       soilhydrology_vars, soilstate_vars, zwt_lat, qflx_sat, volrmch_grc)
    !
    ! !DESCRIPTION:
    ! Groundwater-river exchange, a port of ASAP GW2RIVER (module_wtable.f90)
    ! with wtd = -zwt. For each soil and crop column of a grid cell with a river
    ! (RIV_WIDTH * RIV_LENGTH > 0):
    !
    !   d_rs = maxdepth - rd               depth of the river surface [m]
    !   cond = width*length*hydcon/gw2river_bed_thickness
    !   A: zwt <  d_rs           : q_out = min(cond*(d_rs-zwt)/A, gw2river_max_rate)
    !   B: d_rs <= zwt < maxdepth: q_in  = max(min(min(cond*(zwt-d_rs)/A, max_rate), rd/dt), 0)*fr
    !   C: zwt >= maxdepth       : q_in  = max(min(hksat_top_ref, rd/dt), 0)*fr
    !
    ! where A is the area of the cell's soil and crop columns (the river bed is
    ! shared among them in proportion to their area, so each column sees the
    ! same conductance per unit area), fr = min(width*length/A, 1) and rd is the
    ! river water depth (from the MOSART main-channel volume with
    ! use_lnd_rof_gw_exchange, otherwise 0 so that B and C give 0). Cells with
    ! d_rs <= 0 (flooded) are skipped, as in ASAP.
    !
    ! q_out is reduced by ice impedance (as SIMTOP) and limited so that the
    ! saturated zone can supply it together with the lateral outflow; q_in of a
    ! cell is limited by the river water reported at the last MOSART coupling
    ! minus the seepage since then. q_in - q_out is added to qflx_sat, which
    ! ApplySaturatedLateralFlux applies at the water table. Sets
    ! col_wf%qflx_gw2riv (added to qflx_drain in Drainage) and col_wf%qflx_riv2gw.
    !
    ! !USES:
    use elm_varctl       , only : gw2river_bed_thickness, gw2river_max_rate
    use elm_varctl       , only : use_lnd_rof_gw_exchange, lnd_rof_coupling_nstep
    use elm_varcon       , only : denice, watmin, e_ice
    use elm_time_manager , only : get_nstep
    use GridcellType     , only : grc_pp
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    real(r8)                 , intent(in)    :: dtime                ! time step [s]
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    real(r8)                 , intent(in)    :: zwt_lat(bounds%begc:)     ! water table depth [m]
    real(r8)                 , intent(inout) :: qflx_sat(bounds%begc:)    ! net saturated source [mm H2O/s]
    real(r8)                 , intent(in)    :: volrmch_grc(bounds%begg:) ! MOSART main-channel volume [m3]
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, jwt, nb, nstep
    real(r8) :: wl, rd, d_rs, maxd, zw, hydcon, cond, fr, aaq, q_in, q_out, cap
    real(r8) :: s_y, avail, dzsum, icesum, imped, vol, avail_mm, scale
    real(r8), allocatable :: wcell(:)    ! summed weight of the cell's soil and crop columns [-]
    real(r8), allocatable :: demand(:)   ! river-to-aquifer seepage requested this step [mm, cell mean]
    !-----------------------------------------------------------------------

    associate(                                                    &
         zi         => col_pp%zi                                , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         dz         => col_pp%dz                                , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         nlevbed    => col_pp%nlevbed                           , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_liq => col_ws%h2osoi_liq                        , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice                        , & ! Input:  [real(r8) (:,:) ] ice (kg/m2)
         watsat     => soilstate_vars%watsat_col                , & ! Input:  [real(r8) (:,:) ] porosity
         sucsat     => soilstate_vars%sucsat_col                , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
         bsw        => soilstate_vars%bsw_col                   , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         fdepth     => soilstate_vars%fdepth_col                , & ! Input:  [real(r8) (:)   ] e-folding depth (m)
         hk_top     => soilstate_vars%hksat_top_ref_col         , & ! Input:  [real(r8) (:)   ] top-layer hksat (mm/s)
         width      => soilstate_vars%riv_width_col             , & ! Input:  [real(r8) (:)   ] river width (m)
         length     => soilstate_vars%riv_length_col            , & ! Input:  [real(r8) (:)   ] river length (m)
         maxdepth   => soilstate_vars%riv_depth_col               & ! Input:  [real(r8) (:)   ] bankfull depth (m)
         )

      if (.not. allocated(riv_seep_cum_grc)) then
         allocate(riv_seep_cum_grc(bounds%begg:bounds%endg))
         riv_seep_cum_grc(:) = 0._r8
      end if

      ! the river volume is updated at each MOSART coupling (as for use_lnd_rof_two_way);
      ! coupling times coincide with restart times, so the counter needs no restart field
      nstep = get_nstep()
      if (use_lnd_rof_gw_exchange) then
         if (mod(nstep, lnd_rof_coupling_nstep) == 1 .or. nstep <= 1 .or. lnd_rof_coupling_nstep == 1) then
            riv_seep_cum_grc(bounds%begg:bounds%endg) = 0._r8
         end if
      end if

      allocate(wcell (bounds%begg:bounds%endg)) ; wcell(:)  = 0._r8
      allocate(demand(bounds%begg:bounds%endg)) ; demand(:) = 0._r8

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         col_wf%qflx_gw2riv(c) = 0._r8
         col_wf%qflx_riv2gw(c) = 0._r8
         if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
         g = col_pp%gridcell(c)
         wcell(g) = wcell(g) + col_pp%wtgcell(c)
      end do

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
         g  = col_pp%gridcell(c)
         wl = width(c)*length(c)
         if (wl <= 0._r8 .or. wcell(g) <= 0._r8) cycle

         aaq  = grc_pp%area(g)*1.e6_r8*wcell(g)          ! area of the cell's soil and crop columns [m2]
         maxd = maxdepth(c)

         ! river water depth [m]
         rd = 0._r8
         if (use_lnd_rof_gw_exchange) then
            vol = volrmch_grc(g)
            if (vol > 0._r8 .and. vol < 1.e30_r8) rd = vol/wl
         end if

         d_rs = maxd - rd
         if (d_rs <= 0._r8) cycle                          ! flooded (ASAP skips this case)

         zw     = zwt_lat(c)
         hydcon = hk_top(c)*max(min(exp((1.5_r8 - maxd)/fdepth(c)), 1._r8), 0.1_r8)   ! [mm/s]
         cond   = wl*hydcon/gw2river_bed_thickness                                      ! [mm/s * m]
         fr     = min(wl/aaq, 1._r8)
         q_in   = 0._r8
         q_out  = 0._r8

         if (zw < d_rs) then
            ! A: water table above the river surface, groundwater drains to the river
            q_out = min(cond*(d_rs - zw)/aaq, gw2river_max_rate)

            ! ice impedance (as SIMTOP) and the layers at and below the water table
            nb  = nlevbed(c)
            jwt = nb
            do j = 1, nb
               if (zw <= zi(c,j)) then
                  jwt = j - 1
                  exit
               end if
            end do
            dzsum  = 0._r8
            icesum = 0._r8
            avail  = 0._r8
            do j = jwt+1, nb
               dzsum  = dzsum  + dz(c,j)
               icesum = icesum + dz(c,j)*min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*denice)/watsat(c,j))
               s_y = watsat(c,j) * (1._r8 - (1._r8 + 1.e3_r8*zw/sucsat(c,j))**(-1._r8/bsw(c,j)))
               s_y = max(s_y, 0.02_r8)
               avail = avail + min(s_y*(zi(c,j) - max(zw, zi(c,j-1)))*1.e3_r8, &
                    max(h2osoi_liq(c,j) - watmin, 0._r8))
            end do
            if (dzsum > 0._r8) then
               imped = 10._r8**(-e_ice*(icesum/dzsum))
            else
               imped = 1._r8
            end if
            q_out = imped*q_out
            ! the saturated zone must supply q_out together with any lateral outflow
            q_out = min(q_out, max(0.99_r8*avail/dtime + qflx_sat(c), 0._r8))
         else if (zw < maxd) then
            ! B: connected, water table below the river surface; the river feeds the aquifer
            cap  = min(cond*(zw - d_rs)/aaq, gw2river_max_rate)
            q_in = max(min(cap, rd*1.e3_r8/dtime), 0._r8)*fr
         else
            ! C: water table below the river bed, disconnected; infiltration at Ksat
            q_in = max(min(hk_top(c), rd*1.e3_r8/dtime), 0._r8)*fr
         end if

         col_wf%qflx_gw2riv(c) = q_out
         col_wf%qflx_riv2gw(c) = q_in
         demand(g) = demand(g) + col_pp%wtgcell(c)*q_in*dtime
      end do

      ! limit the seepage of each cell by the river water available since the last coupling
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (col_wf%qflx_riv2gw(c) <= 0._r8) cycle
         g = col_pp%gridcell(c)
         vol = volrmch_grc(g)
         if (.not. (vol > 0._r8 .and. vol < 1.e30_r8)) vol = 0._r8
         avail_mm = max(vol/(grc_pp%area(g)*1.e6_r8)*1.e3_r8 - riv_seep_cum_grc(g), 0._r8)
         if (demand(g) > avail_mm) then
            scale = avail_mm/demand(g)
            col_wf%qflx_riv2gw(c) = col_wf%qflx_riv2gw(c)*scale
         end if
      end do

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (.not. (col_pp%is_soil(c) .or. col_pp%is_crop(c))) cycle
         g = col_pp%gridcell(c)
         riv_seep_cum_grc(g) = riv_seep_cum_grc(g) + col_pp%wtgcell(c)*col_wf%qflx_riv2gw(c)*dtime
         qflx_sat(c) = qflx_sat(c) + col_wf%qflx_riv2gw(c) - col_wf%qflx_gw2riv(c)
      end do

      deallocate(wcell, demand)

    end associate

  end subroutine ComputeRiverGWExchange

  !-----------------------------------------------------------------------
  subroutine ApplySaturatedLateralFlux(bounds, num_hydrologyc, filter_hydrologyc, dtime, &
       soilhydrology_vars, soilstate_vars, qflx_sat, zwt_lat)
    !
    ! !DESCRIPTION:
    ! Applies each column's net saturated lateral flux at the water table,
    ! after the Richards solve, following ELM's baseflow treatment (Drainage):
    !
    ! - Outflow is removed from the layer that contains the water table and
    !   then from the layers below it in turn; each layer gives at most its
    !   specific yield times its saturated thickness, and the water table
    !   drops accordingly. If the saturated zone cannot supply the outflow,
    !   the rest is taken from the unsaturated layers above, bottom up.
    ! - Inflow fills the layer that contains the water table and then the
    !   layers above it in turn, each up to its effective porosity, and the
    !   water table rises accordingly. Any inflow left when the column is full
    !   is added to the top layer, whose excess goes to surface runoff in
    !   Drainage.
    !
    ! The water table starts from zwt_lat, the one the fluxes were computed
    ! with. The applied flux is added to col_wf%qflx_lat_layer.
    !
    ! !USES:
    use elm_varcon , only : denice, watmin
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    real(r8)                 , intent(in)    :: dtime                ! time step [s]
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    real(r8)                 , intent(in)    :: qflx_sat(bounds%begc:) ! net saturated lateral flux into each column [mm H2O/s]
    real(r8)                 , intent(in)    :: zwt_lat(bounds%begc:)  ! water table depth used for the lateral fluxes [m]
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, j, jwt, nb
    real(r8) :: tot        ! flux still to apply [mm H2O] (positive into the column)
    real(r8) :: dw         ! water added to a layer [mm H2O]
    real(r8) :: s_y        ! specific yield [-]
    real(r8) :: avail      ! water a layer can give [mm H2O]
    real(r8) :: cap        ! water a layer can take [mm H2O]
    real(r8), parameter :: tol = 1.e-8_r8   ! [mm H2O]
    !-----------------------------------------------------------------------

    associate(                                         &
         zi         => col_pp%zi                     , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         dz         => col_pp%dz                     , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         nlevbed    => col_pp%nlevbed                , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_liq => col_ws%h2osoi_liq             , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ] ice (kg/m2)
         watsat     => soilstate_vars%watsat_col     , & ! Input:  [real(r8) (:,:) ] porosity
         sucsat     => soilstate_vars%sucsat_col     , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
         bsw        => soilstate_vars%bsw_col        , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         zwt        => soilhydrology_vars%zwt_col      & ! Output: [real(r8) (:)   ] water table depth (m)
         )

      do fc = 1, num_hydrologyc
         c   = filter_hydrologyc(fc)
         nb  = nlevbed(c)
         tot = qflx_sat(c)*dtime
         if (tot == 0._r8) cycle
         zwt(c) = zwt_lat(c)

         ! index of the deepest layer that lies entirely above the water table
         jwt = nb
         do j = 1, nb
            if (zwt(c) <= zi(c,j)) then
               jwt = j - 1
               exit
            end if
         end do

         if (tot < 0._r8) then
            ! outflow: drain the layers at and below the water table in turn
            do j = jwt+1, nb
               s_y = watsat(c,j) * (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
               s_y = max(s_y, 0.02_r8)
               avail = min(s_y*(zi(c,j) - zwt(c))*1.e3_r8, max(h2osoi_liq(c,j) - watmin, 0._r8))
               dw = min(max(tot, -avail), 0._r8)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
               if (tot >= 0._r8) then
                  zwt(c) = zwt(c) - dw/s_y/1.e3_r8
                  exit
               else
                  zwt(c) = zi(c,j)
               end if
            end do
            ! saturated zone exhausted: take the rest from the layers above
            do j = min(jwt, nb), 1, -1
               if (tot >= 0._r8) exit
               dw = min(max(tot, -max(h2osoi_liq(c,j) - watmin, 0._r8)), 0._r8)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
            end do
            if (tot < -tol) then
               call endrun(msg='ApplySaturatedLateralFlux: lateral outflow exceeds the water in the column'// &
                    errMsg(__FILE__, __LINE__))
            end if
         else
            ! inflow: fill the layers at and above the water table in turn
            do j = min(jwt+1, nb), 1, -1
               cap = max(watsat(c,j) - h2osoi_ice(c,j)/(dz(c,j)*denice), 0.01_r8)*dz(c,j)*1.e3_r8
               cap = max(cap - h2osoi_liq(c,j), 0._r8)
               dw  = min(tot, cap)
               h2osoi_liq(c,j) = h2osoi_liq(c,j) + dw
               col_wf%qflx_lat_layer(c,j) = col_wf%qflx_lat_layer(c,j) + dw/dtime
               tot = tot - dw
               s_y = watsat(c,j) * (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
               s_y = max(s_y, 0.02_r8)
               zwt(c) = max(zi(c,j-1), zwt(c) - dw/s_y/1.e3_r8)
               if (tot <= 0._r8) exit
            end do
            if (tot > 0._r8) then
               h2osoi_liq(c,1) = h2osoi_liq(c,1) + tot
               col_wf%qflx_lat_layer(c,1) = col_wf%qflx_lat_layer(c,1) + tot/dtime
            end if
         end if
      end do

    end associate

  end subroutine ApplySaturatedLateralFlux

  !-----------------------------------------------------------------------
  subroutine ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, soilhydrology_vars, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Diagnose the water table depth from soil moisture (as in the CLM5
    ! theta-based method): find the deepest layer, searching up from the
    ! bottom of the hydrologically active column, whose saturation is at or
    ! below sat_lev, and interpolate between it and the layer below.
    !
    ! !USES:
    use elm_varcon , only : denh2o, denice
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer             :: c, fc, k, nb
    !-----------------------------------------------------------------------

    associate(                                         &
         dz         => col_pp%dz                     , & ! Input:  [real(r8) (:,:) ] layer thickness (m)
         z          => col_pp%z                      , & ! Input:  [real(r8) (:,:) ] layer depth (m)
         zi         => col_pp%zi                     , & ! Input:  [real(r8) (:,:) ] interface depth (m)
         nlevbed    => col_pp%nlevbed                , & ! Input:  [integer  (:)   ] number of hydrologically active layers
         h2osoi_liq => col_ws%h2osoi_liq             , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ] ice (kg/m2)
         h2osoi_vol => col_ws%h2osoi_vol             , & ! Output: [real(r8) (:,:) ] volumetric soil water (m3/m3)
         watsat     => soilstate_vars%watsat_col     , & ! Input:  [real(r8) (:,:) ] porosity
         zwt        => soilhydrology_vars%zwt_col      & ! Output: [real(r8) (:)   ] water table depth (m)
         )

      do fc = 1, num_hydrologyc
         c  = filter_hydrologyc(fc)
         nb = nlevbed(c)

         ! update h2osoi_vol from the bottom up to the first unsaturated layer
         do k = nb, 1, -1
            h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) + h2osoi_ice(c,k)/(dz(c,k)*denice)
            if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) exit
         end do

         zwt(c) = ThetaWaterTableDepth(nb, h2osoi_vol(c,1:nb), watsat(c,1:nb), z(c,1:nb), zi(c,0:nb))
      end do

    end associate

  end subroutine ThetaBasedWaterTable

  !-----------------------------------------------------------------------
  pure function ThetaWaterTableDepth(nb, vol, watsat, z, zi) result(zwt)
    !
    ! !DESCRIPTION:
    ! Water table depth diagnosed from soil moisture (as in the CLM5
    ! theta-based method): find the deepest layer, searching up from the
    ! bottom of the hydrologically active column, whose saturation is at or
    ! below sat_lev, and interpolate between it and the layer below. If the
    ! bottom layer is unsaturated, the water table is at the bottom of the
    ! column (no saturated zone).
    !
    ! !ARGUMENTS:
    integer  , intent(in) :: nb          ! number of hydrologically active layers
    real(r8) , intent(in) :: vol(:)      ! volumetric water content (1:nb) [m3/m3]
    real(r8) , intent(in) :: watsat(:)   ! porosity (1:nb)
    real(r8) , intent(in) :: z(:)        ! layer depth (1:nb) [m]
    real(r8) , intent(in) :: zi(0:)      ! interface depth (0:nb) [m]
    real(r8)              :: zwt         ! water table depth [m]
    !
    ! !LOCAL VARIABLES:
    integer  :: k, k_zwt
    logical  :: all_saturated
    real(r8) :: s1, s2, m, b
    !-----------------------------------------------------------------------

    k_zwt         = nb
    all_saturated = .true.
    do k = nb, 1, -1
       if (vol(k)/watsat(k) <= sat_lev) then
          k_zwt         = k
          all_saturated = .false.
          exit
       end if
    end do
    if (all_saturated) k_zwt = 1

    if (k_zwt == 1) then
       ! all layers below the first are saturated: water table at the bottom of layer 1
       zwt = zi(1)
    else if (k_zwt < nb) then
       ! interpolate between k_zwt and k_zwt+1
       s1  = vol(k_zwt  )/watsat(k_zwt  )
       s2  = vol(k_zwt+1)/watsat(k_zwt+1)
       m   = (z(k_zwt+1) - z(k_zwt))/(s2 - s1)
       b   = z(k_zwt+1) - m*s2
       zwt = max(0._r8, m*sat_lev + b)
    else
       ! bottom layer is unsaturated: water table at the bottom of the column
       zwt = zi(nb)
    end if

  end function ThetaWaterTableDepth

end module SoilLateralFlowMod
