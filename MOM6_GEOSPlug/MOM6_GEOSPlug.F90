#include "MAPL_Generic.h"

! GEOS default real kind
#define GeosKind      4
#define REAL_       real(kind=GeosKind)

module MOM6_GEOSPlugMod

!BOP
! !MODULE: MOM6_GEOSPlugMod -- to couple with MOM6.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a coupler for MOM.
! It uses ESMF AND MAPL. It has heavy dependencies on FMS and MOM.
!
! This should be built like MOM, so that its default reals
! are the same as for MOM.
!
! It does not use the configuration.
! Its time step is the clocks time step.
! Each run invocation runs one time step.
!

!USES:
  use ESMF
  use MAPL
  use MAPL_ConstantsMod,        only: MAPL_TICE

! FMS dependencies
  use field_manager_mod,        only: field_manager_init, field_manager_end
  use mpp_mod,                  only: mpp_exit

! MOM dependencies
  use MOM_diag_manager_infra,   only: MOM_diag_manager_init, MOM_diag_manager_end
  use MOM_coms_infra,           only: MOM_infra_init
  use MOM_io_infra,             only: io_infra_end

  use MOM_domain_infra,         only: get_domain_extent, &
                                      AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR

  use MOM_time_manager,         only: set_calendar_type, time_type, &
                                      set_time, set_date, &
                                      JULIAN

  use MOM_grid,                 only: ocean_grid_type

  use ocean_model_mod,          only: get_ocean_grid,                   &
                                      ocean_model_init,                 &
                                      ocean_model_init_sfc,             &
                                      update_ocean_model,               &
                                      ocean_model_end,                  &
                                      ocean_model_restart,              &
                                      ocean_model_data_get,             &
                                      ocean_public_type,                &
                                      ocean_state_type,                 &
                                      ocean_model_get_UV_surf,          &
                                      ocean_model_get_thickness,        &
                                      ocean_model_get_prog_tracer,      &
                                      ocean_model_get_prog_tracer_index

  use MOM_surface_forcing_gfdl, only: ice_ocean_boundary_type

! Nothing on the MOM side is visible through this module.

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
!EOP

! These are the MOM-side bulletin boards, where things are in
! MOM precision and on its grid

  type MOM_MAPL_Type
   type(ocean_public_type),       pointer   :: Ocean
   type(ice_ocean_boundary_type), pointer   :: Ice_ocean_boundary
   type(ocean_state_type),        pointer   :: Ocean_state
  end type MOM_MAPL_Type

! A wrapper-derived data type to connect our internal state with MOM
  type MOM_MAPLWrap_Type
     type(MOM_MAPL_Type), pointer :: Ptr
  end type MOM_MAPLWrap_Type

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGCM GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining
!   state specs and couplings among its children.  In addition, it creates the
!   children GCs (AGCM and OGCM) and runs their respective SetServices.

!EOP
!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals
!=============================================================================

    __Iam__('SetServices')

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//'SetServices'

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_Method_Initialize,   Initialize, _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_Method_Run,          Run,        _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_Method_Finalize,     Finalize,   _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_Method_WriteRestart, Record,     _RC)

!BOS

! !Import state:
#include "MOM6_GEOSPlug_Import___.h"

! !Export state:
#include "MOM6_GEOSPlug_Export___.h"

!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" , _RC)
    call MAPL_TimerAdd(GC,   name="RUN"        , _RC)
    call MAPL_TimerAdd(GC,   name="FINALIZE"   , _RC)

! Generic SetServices
! -------------------

    call MAPL_GenericSetServices    ( GC, _RC )

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!=============================================================================

!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalOcean wrapper

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),     intent(INOUT) :: GC     ! Gridded component
    type(ESMF_State),        intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK  ! The clock
    integer, optional,       intent(  OUT) :: RC     ! Error code:

!EOP

! Locals

    character(len=ESMF_MAXSTR)             :: COMP_NAME

    integer                                :: counts(7)
    integer                                :: Comm
    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed
    integer                                :: IM, JM, LM
    integer                                :: g_isc,g_iec,g_jsc,g_jec
    integer                                :: g_isd,g_ied,g_jsd,g_jed
    integer                                :: YEAR,MONTH,DAY,HR,MN,SC

! Locals with MOM types

    type(time_type)                        :: Time
    type(time_type)                        :: DT

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type(MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                => null()
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap

    type(ocean_grid_type),         pointer :: Ocean_grid              => null()

    REAL_, pointer                         :: TW  (:,:)        => null()
    REAL_, pointer                         :: SW  (:,:)        => null()
    REAL_, pointer                         :: AREA(:,:)        => null()
    REAL_, pointer                         :: MOM_2D_MASK(:,:) => null()
    REAL_, pointer                         :: SLV(:,:)         => null()

    real, allocatable                      :: Tmp2(:,:), Tmp3(:,:,:)

    REAL_, pointer, dimension(:, :, :)     :: DH, TL, SL

    integer                                :: DT_OCEAN
    character(len=7)                       :: wind_stagger     ! 'AGRID' or 'BGRID' or 'CGRID'
    integer                                ::iwind_stagger     !  AGRID  or  BGRID  or  CGRID
    integer                                :: i

    __Iam__('Initialize')

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//'Initialize'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"     )
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get the grid, configuration
!----------------------------

    call ESMF_GridCompGet( GC, grid=Grid,  _RC)

! Get the layout from the grid
!-----------------------------

    call ESMF_VMGetCurrent(VM, _RC)

! Set the time for MOM
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  _RC)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                       _RC)

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, _RC)

! Allocate this instance of the internal state and wrap
! -----------------------------------------------------

    allocate ( MOM_MAPL_internal_state, __STAT__)

    wrap%ptr => MOM_MAPL_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'MOM_MAPL_state', WRAP, _RC)

    allocate ( Boundary, __STAT__)
    allocate ( Ocean,    __STAT__)

    MOM_MAPL_internal_state%Ice_ocean_boundary => Boundary
    MOM_MAPL_internal_state%Ocean              => Ocean

! FMS initialization using the communicator from the VM
!------------------------------------------------------

    call ESMF_VMGet(VM, mpiCommunicator=Comm, _RC)

    call MOM_infra_init(Comm)

! Init MOM stuff
!---------------

    call field_manager_init
    call set_calendar_type ( JULIAN)
    call MOM_diag_manager_init

    DT   = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)

! Check run time wind stagger option set in AGCM.rc (GEOS config)
! to make sure it matches what is expected here
!----------------------------------------------------------------

    call MAPL_GetResource( MAPL, wind_stagger, Label="ocean_wind_stagger:", DEFAULT="AGRID", _RC)

    if ( trim(wind_stagger) == "AGRID") then
      iwind_stagger = AGRID
      if (MAPL_AM_I_Root()) print *, ' Surface stress stagger for MOM6: AGRID. Its value= ', AGRID
    elseif ( ( trim(wind_stagger) == "BGRID") .or. ( trim(wind_stagger) == "CGRID")) then
      print *, ' Surface stress stagger for MOM6: BGRID_NE or CGRID_NE. These options are not supported. Exiting!'
      ASSERT_(.false.)
    else
      print *, ' Surface stress stagger for MOM6 is invalid, stopping.'
      ASSERT_(.false.)
    endif

! Initialize ocean model
!-----------------------

    Ocean%is_ocean_pe = .true.
    call ocean_model_init  (Ocean, Ocean_state, Time, Time, iwind_stagger)
 
    MOM_MAPL_internal_state%Ocean_State => Ocean_State

    call ocean_model_init_sfc(Ocean_state, Ocean)

! Get the ocean grid and sizes of global and computational domains
!-----------------------------------------------------------------

    call get_ocean_grid (Ocean_state, Ocean_grid)
    isc  = Ocean_grid%isc; iec  = Ocean_grid%iec
    isd  = Ocean_grid%isd; ied  = Ocean_grid%ied

    jsc  = Ocean_grid%jsc; jec  = Ocean_grid%jec
    jsd  = Ocean_grid%jsd; jed  = Ocean_grid%jed

!   Notes:
!   ------
!   - We need indices for boundary fluxes: Boundary%...
!   - These arrays have "halos on this processor, centered at h points," see MOM_domain_infra.F90
    call get_domain_extent(Ocean_grid%Domain%mpp_domain, g_isc, g_iec, g_jsc, g_jec, g_isd, g_ied, g_jsd, g_jed)

! Check local sizes of horizontal dimensions
!--------------------------------------------
    call MAPL_GridGet(GRID, localCellCountPerDim=counts, _RC)

    IM=iec-isc+1
    JM=jec-jsc+1

    ASSERT_(counts(1)==IM)
    ASSERT_(counts(2)==JM)

    call get_ocean_grid(Ocean_state, Ocean_grid)
    LM=Ocean_grid%ke

! Check run time surface current stagger option set in MOM_input 
! to make sure it matches what is expected here
!---------------------------------------------------------------

    if (MAPL_AM_I_Root()) then
     if ( (Ocean%stagger == AGRID) .or. (Ocean%stagger == BGRID_NE)) then
       print *, ' Surface velocity stagger set in ocean model: (MOM6) AGRID or BGRID_NE. These are not supported, try CGRID_NE. Exiting!'
       ASSERT_(.false.)
     elseif (Ocean%stagger == CGRID_NE) then
       print *, ' Surface velocity stagger set in ocean model: (MOM6) CGRID_NE.'
     else
       print *, ' Surface velocity stagger set in ocean model: (MOM6) is invalid, stopping.'
       ASSERT_(.false.)
     endif
    endif

! Allocate MOM flux bulletin board.
!----------------------------------

    allocate ( Boundary% u_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% v_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% t_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% q_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% salt_flux       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% lw_flux         (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_vis_dir (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_vis_dif (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_nir_dir (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_nir_dif (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% lprec           (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% fprec           (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% runoff          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% calving         (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% stress_mag      (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% ustar_berg      (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% area_berg       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% mass_berg       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% runoff_hflx     (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% calving_hflx    (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% p               (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% mi              (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% ice_rigidity    (g_isd:g_ied,g_jsd:g_jed), &
               __STAT__)

! Initialize fluxes
!------------------

    Boundary%u_flux          = 0.0
    Boundary%v_flux          = 0.0
    Boundary%t_flux          = 0.0
    Boundary%q_flux          = 0.0
    Boundary%salt_flux       = 0.0
    Boundary%lw_flux         = 0.0
    Boundary%sw_flux_vis_dir = 0.0
    Boundary%sw_flux_vis_dif = 0.0
    Boundary%sw_flux_nir_dir = 0.0
    Boundary%sw_flux_nir_dif = 0.0
    Boundary%lprec           = 0.0
    Boundary%fprec           = 0.0
    Boundary%runoff          = 0.0
    Boundary%calving         = 0.0
    Boundary%stress_mag      = 0.0
    Boundary%ustar_berg      = 0.0
    Boundary%area_berg       = 0.0
    Boundary%mass_berg       = 0.0
    Boundary%runoff_hflx     = 0.0
    Boundary%calving_hflx    = 0.0
    Boundary%p               = 0.0
    Boundary%mi              = 0.0
    Boundary% ice_rigidity   = 0.0

! Profilers
! ---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"     )

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, _RC )

! Make sure exports neede by the parent prior to our run call are initialized
!----------------------------------------------------------------------------

    call MAPL_GetPointer(EXPORT, MOM_2D_MASK, 'MOM_2D_MASK', alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, TW,          'TW'  ,        alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, SW,          'SW'  ,        alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, AREA,        'AREA',        alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, SLV,         'SLV',         alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, DH,          'DH'  ,        alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, TL,          'T',           alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, SL,          'S',           alloc=.true., _RC)

! Get the 3-D MOM data
!---------------------
    allocate(Tmp3(IM,JM,LM), __STAT__)

    call ocean_model_get_thickness(Ocean_State, tmp3, isc, jsc)
    DH = real(tmp3, kind=GeosKind)

    call ocean_model_get_prog_tracer_index(Ocean_State,i,'temp')
    call ocean_model_get_prog_tracer(Ocean_State,i, tmp3, isc, jsc)
    TL = real(tmp3,kind=GeosKind) + MAPL_TICE

    call ocean_model_get_prog_tracer_index(Ocean_State,i,'salt')
    call ocean_model_get_prog_tracer(Ocean_State, i, tmp3, isc, jsc)
    SL = real(tmp3,kind=GeosKind)

! Get the 2-D MOM data
!---------------------
    allocate(Tmp2(IM,JM), __STAT__)

    call ocean_model_data_get(Ocean_State, Ocean, 'mask', Tmp2, isc, jsc)
    MOM_2D_MASK = real(Tmp2, kind=GeosKind)

    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', Tmp2, g_isc, g_jsc) ! this comes to us in deg C
    where(MOM_2D_MASK(:,:) > 0.0)
       TW = real(Tmp2, kind=GeosKind) + MAPL_TICE                                 ! because C to K was subtracted in MOM
    elsewhere
       TW = MAPL_UNDEF
    end where

    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', Tmp2, g_isc, g_jsc) ! comes to us in PSU
    where(MOM_2D_MASK(:,:) > 0.0)
       SW = real(Tmp2, kind=GeosKind)
    elsewhere
       SW = MAPL_UNDEF
    end where

    if(associated(area)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'area', Tmp2, isc, jsc)
       AREA = real(Tmp2, kind=GeosKind)
    end if

    deallocate(Tmp2, __STAT__)
    deallocate(Tmp3, __STAT__)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!=================================================================================

!BOP

! !IROUTINE: Run  -- Run method for External Ocean Model

! !INTERFACE:

  subroutine Run  ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC       ! Gridded component
    type(ESMF_State),    intent(INOUT) :: IMPORT   ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT   ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK    ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC       ! Error code:
    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals with ESMF and MAPL types

    type(MAPL_MetaComp),       pointer :: MAPL               => null()
    type(ESMF_Time)                    :: MyTime
    type(ESMF_TimeInterval)            :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                 => null()
    type(ocean_public_type),       pointer :: Ocean                    => null()
    type(ocean_state_type),        pointer :: Ocean_State              => null()
    type(ocean_grid_type),         pointer :: Ocean_grid               => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state  => null()
    type(MOM_MAPLWrap_Type)                :: wrap

!#include "MOM6_GEOSPlug_DeclarePointer___.h" ! Because these are "real(kind=GeosKind)" not using ACG.
! Exports
    REAL_, pointer                     :: TW    (:,:)        => null()
    REAL_, pointer                     :: SW    (:,:)        => null()
    REAL_, pointer                     :: UW    (:,:)        => null()
    REAL_, pointer                     :: VW    (:,:)        => null()
    REAL_, pointer                     :: UWB   (:,:)        => null()
    REAL_, pointer                     :: VWB   (:,:)        => null()
    REAL_, pointer                     :: UWC   (:,:)        => null()
    REAL_, pointer                     :: VWC   (:,:)        => null()
    REAL_, pointer                     :: SLV   (:,:)        => null()
    REAL_, pointer                     :: FRAZIL(:,:)        => null()
    REAL_, pointer                     :: MELT_POT(:,:)      => null()
    REAL_, pointer                     :: FRZMLT(:,:)        => null()
    REAL_, pointer                     :: MOM_2D_MASK  (:,:) => null()
    REAL_, pointer                     :: AREA  (:,:)        => null()
    REAL_, pointer                     :: T_Freeze (:,:)     => null()
    REAL_, pointer                     :: DH  (:,:,:)        => null()
    REAL_, pointer                     :: TL  (:,:,:)        => null()
    REAL_, pointer                     :: SL  (:,:,:)        => null()

! Imports
    REAL_, pointer                     :: TAUX(:,:)          => null()
    REAL_, pointer                     :: TAUY(:,:)          => null()
    REAL_, pointer                     :: PS  (:,:)          => null()
    REAL_, pointer                     :: PICE(:,:)          => null()
    REAL_, pointer                     :: LWFLX(:,:)         => null()
    REAL_, pointer                     :: SHFLX(:,:)         => null()
    REAL_, pointer                     :: QFLUX(:,:)         => null()
    REAL_, pointer                     :: RAIN(:,:)          => null()
    REAL_, pointer                     :: SNOW(:,:)          => null()
    REAL_, pointer                     :: SFLX(:,:)          => null()
    REAL_, pointer                     :: PENUVR(:,:)        => null()
    REAL_, pointer                     :: PENPAR(:,:)        => null()
    REAL_, pointer                     :: PENUVF(:,:)        => null()
    REAL_, pointer                     :: PENPAF(:,:)        => null()
    REAL_, pointer                     :: DRNIR(:,:)         => null()
    REAL_, pointer                     :: DFNIR(:,:)         => null()
    REAL_, pointer                     :: DISCHARGE(:,:)     => null()
    REAL_, pointer                     :: AICE(:,:)          => null()
    REAL_, pointer                     :: TAUXBOT(:,:)       => null()
    REAL_, pointer                     :: TAUYBOT(:,:)       => null()

! Temporaries
    real, allocatable                  :: U(:,:),  V(:,:), H(:,:,:)
    real, allocatable                  :: cos_rot(:,:)
    real, allocatable                  :: sin_rot(:,:)

    integer                            :: IM, JM, LM

    integer                            :: steady_state_ocean = 0       ! SA: Per Atanas T, "name" of this var is misleading
                                                                       ! We run ocean model only when it = 0

    character(len=7)                   :: pres_loading                 ! yes or no

    integer                            :: isc,iec,jsc,jec

    integer                            :: YEAR,MONTH,DAY,HR,MN,SC
    type(time_type)                    :: Time
    type(time_type)                    :: DT

    real                               :: pice_scaling = 1.0
    real                               :: dTf_dS
    integer                            :: DT_OCEAN
    integer                            :: tracer_index

    REAL_, pointer, dimension(:,:)     :: LATS  => null()
    REAL_, pointer, dimension(:,:)     :: LONS  => null()

    __Iam__('Run')

! Begin
!------

! Get the component name and set-up traceback handle.
! -----------------------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC( GC, MAPL, _RC)

    call MAPL_Get(MAPL,                      &
         INTERNAL_ESMF_STATE = INTERNAL,     &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
         _RC)

! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )

! Get the private internal state
!-------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS); VERIFY_(STATUS)
    MOM_MAPL_internal_state => WRAP%PTR

! Aliases to MOM types
!---------------------

    Boundary    => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean       => MOM_MAPL_internal_state%Ocean
    Ocean_State => MOM_MAPL_internal_state%Ocean_State

! Get domain size
!----------------
    call get_domain_extent(Ocean%Domain, isc, iec, jsc, jec)

    IM=iec-isc+1
    JM=jec-jsc+1

    call get_ocean_grid (Ocean_state, Ocean_grid)
    LM=Ocean_grid%ke

! Temporaries with MOM default reals
!-----------------------------------

    allocate(U(IM,JM   ),    __STAT__)
    allocate(V(IM,JM   ),    __STAT__)
    allocate(H(IM,JM,LM),    __STAT__)
    allocate(cos_rot(IM,JM), __STAT__)
    allocate(sin_rot(IM,JM), __STAT__)

! Get pointers to imports and exports
!------------------------------------
!#include "MOM6_GEOSPlug_GetPointer___.h" ! Because connectivities across MOM5 and MOM6 are messy, do these manually :(

    call MAPL_GetPointer(IMPORT, TAUX,     'TAUX'  ,    _RC)
    call MAPL_GetPointer(IMPORT, TAUY,     'TAUY'  ,    _RC)
    call MAPL_GetPointer(IMPORT, PS,       'PS'    ,    _RC)
    call MAPL_GetPointer(IMPORT, PICE,     'PICE'  ,    _RC)
    call MAPL_GetPointer(IMPORT, LWFLX,    'LWFLX'  ,   _RC)
    call MAPL_GetPointer(IMPORT, SHFLX,    'SHFLX'  ,   _RC)
    call MAPL_GetPointer(IMPORT, QFLUX,    'QFLUX'  ,   _RC)
    call MAPL_GetPointer(IMPORT, RAIN,     'RAIN'  ,    _RC)
    call MAPL_GetPointer(IMPORT, SNOW,     'SNOW'  ,    _RC)
    call MAPL_GetPointer(IMPORT, SFLX,     'SFLX'  ,    _RC)
    call MAPL_GetPointer(IMPORT, PENUVR,   'PENUVR'  ,  _RC)
    call MAPL_GetPointer(IMPORT, PENPAR,   'PENPAR'  ,  _RC)
    call MAPL_GetPointer(IMPORT, PENUVF,   'PENUVF'  ,  _RC)
    call MAPL_GetPointer(IMPORT, PENPAF,   'PENPAF'  ,  _RC)
    call MAPL_GetPointer(IMPORT, DRNIR,    'DRNIR'  ,   _RC)
    call MAPL_GetPointer(IMPORT, DFNIR,    'DFNIR'  ,   _RC)
    call MAPL_GetPointer(IMPORT, DISCHARGE,'DISCHARGE', _RC)
    call MAPL_GetPointer(IMPORT, AICE,     'AICE',      _RC)
    call MAPL_GetPointer(IMPORT, TAUXBOT,  'TAUXBOT',   _RC)
    call MAPL_GetPointer(IMPORT, TAUYBOT,  'TAUYBOT',   _RC)

! Get EXPORT pointers
!--------------------

    call MAPL_GetPointer(EXPORT, DH,    'DH'  ,   _RC)
    call MAPL_GetPointer(EXPORT, TL,    'T'   ,   _RC)
    call MAPL_GetPointer(EXPORT, SL,    'S'   ,   _RC)

    call MAPL_GetPointer(EXPORT, UW,    'UW'  ,   _RC)
    call MAPL_GetPointer(EXPORT, VW,    'VW'  ,   _RC)
    call MAPL_GetPointer(EXPORT, UWB,   'UWB' ,   _RC)
    call MAPL_GetPointer(EXPORT, VWB,   'VWB' ,   _RC)
    call MAPL_GetPointer(EXPORT, UWC,   'UWC' ,   _RC)
    call MAPL_GetPointer(EXPORT, VWC,   'VWC' ,   _RC)
    call MAPL_GetPointer(EXPORT, TW,    'TW'  ,   _RC)
    call MAPL_GetPointer(EXPORT, SW,    'SW'  ,   _RC)
    call MAPL_GetPointer(EXPORT, SLV,   'SLV',    _RC)

    call MAPL_GetPointer(EXPORT, FRAZIL,  'FRAZIL',   alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, MELT_POT,'MELT_POT', alloc=.true., _RC)
    call MAPL_GetPointer(EXPORT, FRZMLT,  'FRZMLT',                 _RC)
    call MAPL_GetPointer(EXPORT, T_Freeze,'T_Freeze',               _RC)

    call MAPL_GetPointer(EXPORT, MOM_2D_MASK, 'MOM_2D_MASK', _RC)
    call MAPL_GetPointer(EXPORT, AREA, 'AREA',               _RC)

! Fill in ocean boundary fluxes/forces
!-------------------------------------

    call MAPL_GetResource(MAPL, pice_scaling, Label="MOM_PICE_SCALING:", DEFAULT= 1.0, _RC)
    call MAPL_GetResource(MAPL, pres_loading, Label="pres_loading:",     DEFAULT="NO", _RC)
    call MAPL_GetResource(MAPL, dTf_dS,       Label="DTFREEZE_DS:",      DEFAULT= -0.054, _RC) !The derivative of the freezing temperature with salinity in deg_C/PSU

    ! NOTE: PICE that is available here is all = 0. This should be made realistic, for now it is from MOM5 legacy
    !       Need to study with zero pressure loading (CTL: as now), exp1 ( with PS only), exp2 (with PS and PICE), exp3 (PICE only).
    if ( pres_loading == "YES") then
      Boundary%P              (isc:iec,jsc:jec)= real(PS,          kind=KIND(Boundary%p)) + & ! Pressure of overlying atmospheric
                                  pice_scaling * real(PICE,        kind=KIND(Boundary%p))     ! and ice
    else
      Boundary%P              (isc:iec,jsc:jec)= &
                                   pice_scaling* real(PICE,        kind=KIND(Boundary%p)) ! Pressure of overlying ice only
    end if

    Boundary%lw_flux        (isc:iec,jsc:jec)= real(LWFLX,         kind=KIND(Boundary%p)) ! Long wave flux: both positive down
    Boundary%t_flux         (isc:iec,jsc:jec)= real(SHFLX,         kind=KIND(Boundary%p)) ! Sensible heat flux: both positive up
    Boundary%q_flux         (isc:iec,jsc:jec)= real(QFLUX,         kind=KIND(Boundary%p)) ! specific humidity flux [kg m-2 s-1] ( OR evaporation flux ?)
    Boundary%lprec          (isc:iec,jsc:jec)= real(RAIN,          kind=KIND(Boundary%p)) ! Liquid precipitation: both positive down
    Boundary%fprec          (isc:iec,jsc:jec)= real(SNOW,          kind=KIND(Boundary%p)) ! Frozen precipitation: both positive down
    Boundary%salt_flux      (isc:iec,jsc:jec)=-real(SFLX,          kind=KIND(Boundary%p)) ! Salt flux: MOM positive up, GEOS positive down
    Boundary%runoff         (isc:iec,jsc:jec)= real(DISCHARGE,     kind=KIND(Boundary%p)) ! mass flux of liquid runoff [kg m-2 s-1]

! All shortwave components are positive down in MOM and in GEOS
!--------------------------------------------------------------
    Boundary%sw_flux_vis_dir(isc:iec,jsc:jec)= real(PENUVR+PENPAR, kind=KIND(Boundary%p)) ! direct visible sw radiation        [W m-2]
    Boundary%sw_flux_vis_dif(isc:iec,jsc:jec)= real(PENUVF+PENPAF, kind=KIND(Boundary%p)) ! diffuse visible sw radiation       [W m-2]
    Boundary%sw_flux_nir_dir(isc:iec,jsc:jec)= real(DRNIR,         kind=KIND(Boundary%p)) ! direct  Near InfraRed sw radiation [W m-2]
    Boundary%sw_flux_nir_dif(isc:iec,jsc:jec)= real(DFNIR,         kind=KIND(Boundary%p)) ! diffuse Near InfraRed sw radiation [W m-2]

! Convert input stresses over water to MOM wind stagger
!------------------------------------------------------
    U = 0.0; V = 0.0
!   Using A-grid ice stress, note ice (CICE) stress has opposite sign to atmosphere
    U = real( TAUX*(1.-AICE) - TAUXBOT*AICE, kind=kind(U))
    V = real( TAUY*(1.-AICE) - TAUYBOT*AICE, kind=kind(V))

! Grid rotation angles - these could be saved in the first instance, rather doing every time
!   cos_rot = 1. ! A-grid
    call ocean_model_data_get(Ocean_State, Ocean, 'cos_rot', cos_rot, isc, jsc)
!   sin_rot = 0. ! A-grid
    call ocean_model_data_get(Ocean_State, Ocean, 'sin_rot', sin_rot, isc, jsc)

! Rotate input stress over water along i,j of tripolar grid, and combine with stress under ice
!---------------------------------------------------------------------------------------------
    Boundary%U_flux = 0.0;  Boundary%V_flux = 0.0 ! Initialize stress

    Boundary%U_flux  (isc:iec,jsc:jec)= real( (U*cos_rot - V*sin_rot), kind=KIND(Boundary%p))
    Boundary%V_flux  (isc:iec,jsc:jec)= real( (U*sin_rot + V*cos_rot), kind=KIND(Boundary%p))

! Set the time for MOM
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  _RC)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                       _RC)

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, _RC)

    DT   = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)

! Run MOM for one time step
!--------------------------

    ! set following to non-zero only if the coupled model becomes unstable 
    ! (inconsistent atmosphere or bad restart! or some instabilities) - per Atanas T
    call MAPL_GetResource(MAPL, steady_state_ocean, Label = "steady_state_ocean:", default = 0, _RC)

    if(steady_state_ocean == 0) then
      call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT)
    else
      ! SA: steady_state_ocean: Call update_ocean_model with additional args now available with MOM6
      !     though this will cause MOM_FATAL error
      call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT, .false., .false.)
    endif

! Get required Exports at GEOS precision
!---------------------------------------

!MOM 3D fields
!   Thickness
    call ocean_model_get_thickness(Ocean_State, H, isc, jsc)
    DH = real(H, kind=GeosKind)

!   Temp
    call ocean_model_get_prog_tracer_index(Ocean_State,tracer_index,'temp')
    call ocean_model_get_prog_tracer(Ocean_State,tracer_index, H, isc, jsc)
    TL = real(H,kind=GeosKind) + MAPL_TICE

!   Salt
    call ocean_model_get_prog_tracer_index(Ocean_State,tracer_index,'salt')
    call ocean_model_get_prog_tracer(Ocean_State, tracer_index, H, isc, jsc)
    SL = real(H,kind=GeosKind)

!MOM 2D fields
!   2D mask
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 'mask', U, isc, jsc)
    MOM_2D_MASK = real(U, kind=GeosKind)

!   surface (potential) temperature (K)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', U, isc, jsc) ! this comes to us in deg C
    where(MOM_2D_MASK(:,:) > 0.0)
       TW = real(U, kind=GeosKind) + MAPL_TICE                           ! because C to K was subtracted in MOM
    elsewhere
       TW = MAPL_UNDEF
    end where

!   surface salinity (PSU)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', U, isc, jsc) ! this comes to us in PSU
    where(MOM_2D_MASK(:,:) > 0.0)
       SW = real(U, kind=GeosKind)
    elsewhere
       SW = MAPL_UNDEF
    end where

!   sea level (m)
    U = 0.0
    if(associated(SLV)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', U, isc, jsc) ! this comes to us in m
       where(MOM_2D_MASK(:,:)>0.0)
          SLV = real(U, kind = GeosKind)
       elsewhere
          SLV = 0.0
       end where
    end if

!   frazil (W/m^2)
    U = 0.0
    if(associated(FRAZIL)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'frazil',   U, isc, jsc)  ! this comes to us in J/m2
       where(MOM_2D_MASK(:,:)>0.0)
          FRAZIL =  real( (U)/dt_ocean, kind = GeosKind) ! relying on fortran to promote the int (dt_ocean) to real
       elsewhere
          FRAZIL =  0.0
       end where
    endif

!   melt potential (W/m^2)
    U = 0.0
    if(associated(MELT_POT)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'melt_pot', U, isc, jsc)  ! this comes to us in J/m2
       where(MOM_2D_MASK(:,:)>0.0)
          MELT_POT = -real( (U)/dt_ocean, kind = GeosKind) ! relying on fortran to promote the int (dt_ocean) to real
       elsewhere
          MELT_POT =  0.0
       end where
       MELT_POT = MIN ( MELT_POT, 0.0) ! make sure melt potential is <= 0
    endif

!   freezing melt potential (W/m^2)
    if(associated(FRZMLT)) then
       if ( (.not.associated(FRAZIL)) .or. (.not.associated(MELT_POT))) then
         print *, 'You are asking for freeze melt potential, without asking for frazil and melt potential. You must ask for all. Exiting!'
         ASSERT_(.false.)
       endif

       where(MOM_2D_MASK(:,:)>0.0)
          FRZMLT = MELT_POT
       elsewhere
          FRZMLT = 0.0
       end where

       where(MOM_2D_MASK(:,:)>0.0 .and. FRAZIL>0.0)
          FRZMLT = FRAZIL
       end where

    end if

!   freezing temperature (deg C)
    if(associated(T_Freeze)) then
       where(MOM_2D_MASK(:,:)>0.0)
          T_Freeze = 0.0 + dTf_dS * SW ! exactly as in SIS2_ice_thm.F90
       elsewhere
          T_Freeze = -1.8
       end where
    end if

! currents (m/s)
!---------------
!   A-grid currents (for the atmospheric model)
    U = 0.0; V = 0.0
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'ua', U, isc, jsc) ! this comes to us in m/s
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'va', V, isc, jsc) ! this comes to us in m/s

    if(associated(UW ) .and. associated(VW )) then
      where(MOM_2D_MASK(:,:) > 0.0)
        UW = real(U, kind=GeosKind) * cos_rot + real(V, kind=GeosKind) * sin_rot
        VW = real(V, kind=GeosKind) * cos_rot - real(U, kind=GeosKind) * sin_rot
      elsewhere
        UW=0.0
        VW=0.0
      end where
    else
       print *, 'Both UW and VW MUST be allocated.'
       ASSERT_(.false.)
    endif



!   B-grid currents (for CICE dynamics)
    U = 0.0; V = 0.0
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'ub', U, isc, jsc) ! this comes to us in m/s
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'vb', V, isc, jsc) ! this comes to us in m/s

    if(associated(UWB )) then
      where(MOM_2D_MASK(:,:) > 0.0)
        UWB = real(U, kind=GeosKind)
      elsewhere
        UWB =0.0
      end where
    endif

    if(associated(VWB )) then
      where(MOM_2D_MASK(:,:) > 0.0)
        VWB = real(V, kind=GeosKind)
      elsewhere
        VWB =0.0
      end where
    end if

!   C-grid currents (for CICE dynamics)
    U = 0.0; V = 0.0
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'uc', U, isc, jsc) ! this comes to us in m/s
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'vc', V, isc, jsc) ! this comes to us in m/s

    if(associated(UWC )) then
       UWC = real(U, kind=GeosKind)
    endif

    if(associated(VWC )) then
       VWC = real(V, kind=GeosKind)
    end if

! Optional Exports at GEOS precision
!-----------------------------------
! none
!   3d exports with MOM6, other than depth, T, S (needed by NOBM), such as U, V, etc
!   will not be exported. If needed, write them directly from MOM6

    deallocate(U, V,            __STAT__)
    deallocate(cos_rot,sin_rot, __STAT__)
    deallocate(H,               __STAT__)

    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!BOP

!====================================================================

! !IROUTINE: Finalize        -- Finalize method for GuestOcean wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),           pointer :: MAPL
    type(ESMF_Time)                        :: MyTime
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(ice_ocean_boundary_type), pointer :: Boundary                => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals with MOM types

    type(time_type)                  :: Time
    integer                          :: YEAR,MONTH,DAY,HR,MN,SC

    __Iam__('Finalize')

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

! Get the private internal state
!-------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS); VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean
    Ocean_State => MOM_MAPL_internal_state%Ocean_State

! Set the times for MOM
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, _RC)

    call ESMF_TimeGet (MyTime,      &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         _RC)

    Time = set_date(YEAR,MONTH,DAY,HR,MN,SC)

    call ocean_model_end (Ocean, Ocean_State, Time) ! SA: this also calls ocean_model_save_restart(...)

    call MOM_diag_manager_end(Time )
    call field_manager_end
    call io_infra_end

    deallocate ( Boundary% u_flux          , &
                 Boundary% v_flux          , &
                 Boundary% t_flux          , &
                 Boundary% q_flux          , &
                 Boundary% salt_flux       , &
                 Boundary% lw_flux         , &
                 Boundary% sw_flux_vis_dir , &
                 Boundary% sw_flux_vis_dif , &
                 Boundary% sw_flux_nir_dir , &
                 Boundary% sw_flux_nir_dif , &
                 Boundary% lprec           , &
                 Boundary% fprec           , &
                 Boundary% runoff          , &
                 Boundary% calving         , &
                 Boundary% stress_mag      , &
                 Boundary% ustar_berg      , &
                 Boundary% area_berg       , &
                 Boundary% mass_berg       , &
                 Boundary% runoff_hflx     , &
                 Boundary% calving_hflx    , &
                 Boundary% p               , &
                 Boundary% mi              , &
                 Boundary% ice_rigidity    , &
                 __STAT__)

    deallocate ( Ocean,                   __STAT__)
    deallocate ( Boundary,                __STAT__)
    deallocate ( MOM_MAPL_internal_state, __STAT__)
!

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------

    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, _RC)

    call mpp_exit()

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

!====================================================================

! !IROUTINE: Record -- Record method for GuestOcean wrapper (write intermediate restarts)

! !INTERFACE:

  subroutine Record ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),     pointer :: MAPL
    type(MOM_MAPL_Type),     pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)          :: wrap
    type(ocean_state_type),  pointer :: Ocean_State             => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals
    character(len=14)                :: timeStamp
    logical                          :: doRecord

    __Iam__('Record')

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL")

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, MODE=MAPL_Write2Disk, _RC)

    if (doRecord) then

! Get the private internal state
!--------------------------------

       CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS); VERIFY_(STATUS)

       MOM_MAPL_internal_state => WRAP%PTR
       Ocean_State             => MOM_MAPL_internal_state%Ocean_State

       call MAPL_DateStampGet(clock, timeStamp, _RC)

! Write a restart
!-----------------

       call ocean_model_restart (Ocean_State, timeStamp); VERIFY_(STATUS)

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

end module MOM6_GEOSPlugMod

subroutine SetServices(GC, RC)
   use ESMF
   use MOM6_GEOSPlugMod, only : mySetservices=>SetServices
   type(ESMF_GridComp)  :: GC
   integer, intent(out) :: RC
   call mySetServices(GC, rc=RC)
end subroutine
