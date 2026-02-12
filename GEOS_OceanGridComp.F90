#include "MAPL_Generic.h"

module GEOS_OceanGridCompMod

!BOP
! !MODULE: GEOS_OceanGridCompMod -- Implements ESMF wrapper to invoke the DATASEA/MIT/MOM ocean models.

! !USES:

  use ESMF
  use MAPL
#ifdef BUILD_MIT_OCEAN
  use MIT_GEOS5PlugMod, only: MITSetServices => SetServices  ! this sets IRF
#endif
  use GEOS_DataSeaGridCompMod, only: DataSeaSetServices  => SetServices

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

  character(len=ESMF_MAXSTR)  :: OCEAN_NAME
  integer                     :: DO_DATASEA
  real                        :: OrphanDepth

! !DESCRIPTION:
!
!   {\tt GEOS\_OceanGridComp} is a light-weight gridded component that serves an
!   interface to ocean/data\_ocean components.
!
!EOP

  type :: T_PrivateState
     type(ESMF_Clock)  :: CLOCK
  end type T_PrivateState

  type :: T_PrivateState_Wrap
     type(T_PrivateState), pointer :: ptr
  end type T_PrivateState_Wrap

  integer ::          OCN
  integer ::          OCNd
  logical ::      DUAL_OCEAN
  logical ::     DO_DATA_ATM4OCN

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for Ocean

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices,
!       which sets the Run, Initialize, and Finalize services,
!       as well as allocating our instance of a generic state and putting it in the
!   gridded component (GC). Here we override all three methods and declare
!       the specs for the Imports and Export States (no MAPL controlled Internal State).
!
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Local vars
    type  (MAPL_MetaComp), pointer     :: MAPL
    type  (ESMF_Config)                :: CF
    integer                            :: iDUAL_OCEAN
    character(len=ESMF_MAXSTR)         :: charbuf_
    character(len=ESMF_MAXSTR)         :: sharedObj

    __Iam__('SetServices')

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Get constants from CF
! ---------------------

    call MAPL_GetResource (MAPL, DO_DATASEA,  Label="USE_DATASEA:" ,    DEFAULT=1,    _RC)
    call MAPL_GetResource (MAPL, OrphanDepth, Label="OGCM_TOP_LAYER:" , DEFAULT=10.0, _RC)

    call MAPL_GetResource (MAPL, DO_DATA_ATM4OCN,  Label="USE_DATA_ATM4OCN:", DEFAULT=.FALSE., _RC)

    if(DO_DATASEA/=0) then
       OCEAN_NAME="DATASEA"
       OCN = MAPL_AddChild(GC, NAME=OCEAN_NAME, SS=DataSeaSetServices, _RC)
    else
       call MAPL_GetResource ( MAPL, OCEAN_NAME, Label="OCEAN_NAME:", DEFAULT="MOM6", _RC)
       select case (trim(OCEAN_NAME))
          case ("MOM6")
             call MAPL_GetResource ( MAPL, sharedObj,  Label="MOM6_GEOSPLUG:",    DEFAULT="libMOM6_GEOSPlug.so",    _RC)
             OCN = MAPL_AddChild(OCEAN_NAME,'setservices_', parentGC=GC, sharedObj=sharedObj,  _RC)
#ifdef BUILD_MIT_OCEAN
          case ("MIT")
             OCN = MAPL_AddChild(GC, NAME=OCEAN_NAME, SS=MITSetServices, __RC__)
#endif
          case default
             charbuf_ = "OCEAN_NAME: " // trim(OCEAN_NAME) // " is not implemented, ABORT!"
             call WRITE_PARALLEL(charbuf_)
             VERIFY_(999)
       end select
    endif

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, _RC)
    DUAL_OCEAN = iDUAL_OCEAN /= 0

    OCNd = 0
    if (dual_ocean) then
       OCNd = MAPL_AddChild(GC, NAME="DATASEA", SS=DataSeaSetServices, _RC)
    endif

! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:
#include "GEOS_Ocean_Import___.h"

! !Export state:
#include "GEOS_Ocean_Export___.h"

! !Internal state:
#include "GEOS_Ocean_Internal___.h"

!ALT Note the FRACICE from datasea is inhereted (in AMIP or dual_ocean)

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'MASKO',                             &
         LONG_NAME          = 'ocean_mask',                        &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

! Exports of (a) child

    call MAPL_AddExportSpec (GC, SHORT_NAME = 'TW', CHILD_ID = OCN,  _RC)
    call MAPL_AddExportSpec (GC, SHORT_NAME = 'SW', CHILD_ID = OCN,  _RC)
    call MAPL_AddExportSpec (GC, SHORT_NAME = 'UW', CHILD_ID = OCN,  _RC)
    call MAPL_AddExportSpec (GC, SHORT_NAME = 'VW', CHILD_ID = OCN,  _RC)

    if(DO_DATASEA==0) then
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'DH',  CHILD_ID = OCN, _RC)
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'UWB', CHILD_ID = OCN, _RC)
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'VWB', CHILD_ID = OCN, _RC)
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'SLV', CHILD_ID = OCN, _RC)
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'T',   CHILD_ID = OCN, _RC)
       call MAPL_AddExportSpec (GC, SHORT_NAME = 'S',   CHILD_ID = OCN, _RC)

       if (trim(OCEAN_NAME) == "MOM6") then
          call MAPL_AddExportSpec (GC, SHORT_NAME = 'UWC', CHILD_ID   = OCN, _RC)
          call MAPL_AddExportSpec (GC, SHORT_NAME = 'VWC', CHILD_ID   = OCN, _RC)
       endif
    end if

!EOS

    if(DO_DATASEA==0) then
       call MAPL_TerminateImport    ( GC, SHORT_NAME=               &
           [character(len=9) :: 'TAUX  ','TAUY  ',                   &
              'PENUVR','PENPAR','PENUVF','PENPAF', 'DRNIR', 'DFNIR',  &
              'DISCHARGE','CALVING','LWFLX', 'SHFLX', 'QFLUX', 'RAIN', 'SNOW', &
              'SFLX','SWHEAT'],                                       &  ! do not terminate import of PEN_OCN since it is not used in the `plug'
              CHILD=OCN, _RC)
    end if

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, _RC)
! phase 1
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,    Run, _RC)
! phase 2 - this is only used in the predictor part of the replay for dual ocean
    if (DUAL_OCEAN) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run, _RC)
    endif

! terminate child's import for a temperature correction - we will fill it here, if we run in dual_ocean mode, otherwise nobody needs this variable
    if(DUAL_OCEAN) then
       call MAPL_TerminateImport    ( GC,   &
            SHORT_NAME = (/'DEL_TEMP'/),    &
            CHILD      = OCN,               &
            _RC)
    endif

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="INITIALIZE" , _RC)
    call MAPL_TimerAdd(GC, name="RUN"        , _RC)
    call MAPL_TimerAdd(GC, name="--ModRun"   , _RC)

! Generic SetServices--This creates the generic state and calls SetServices for children
!---------------------------------------------------------------------------------------

    call MAPL_GenericSetServices (GC, _RC )

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

! -----------------------------------------------------------------

!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalOcean wrapper

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),      intent(INOUT) :: GC     ! Gridded component
    type(ESMF_State),         intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),         intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),         intent(INOUT) :: CLOCK  ! The clock
    integer, optional,        intent(  OUT) :: RC     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: State
    type (ESMF_Grid)                    :: Grid
    type (T_PrivateState),    pointer   :: PrivateSTATE
    type (T_PrivateState_Wrap)          :: WRAP
    integer                             :: IM, JM, LM
    real                                :: DT

    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)
    type (ESMF_TimeInterval)            :: timeStep
    type (ESMF_Time)                    :: currTime

    real, pointer ::   MASK(:,:)
    real, pointer ::  MASKO(:,:)
    real, pointer :: MASK3D(:,:,:)
    real, pointer :: DH(:,:,:)

!=============================================================================

    __Iam__('Initialize')

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, grid=GRID, _RC )
    Iam = trim(COMP_NAME)//'::'//'Initialize'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, State, _RC)

! Profilers
!----------

    call MAPL_TimerOn(STATE,"INITIALIZE")
    call MAPL_TimerOn(STATE,"TOTAL"     )

! Get info from the Generic state
!--------------------------------

    call MAPL_Get(STATE, GIM=GIM, GEX=GEX, _RC)

! Allocate the private state...
!------------------------------

    allocate( PrivateSTATE, __STAT__)
    wrap%ptr => PrivateState

! And put it in the GC
!---------------------

    CALL ESMF_UserCompSetInternalState( GC, TRIM(OCEAN_NAME)//'_internal_state', WRAP, STATUS )
    VERIFY_(status)

! Initialize the PrivateState. First the time...
!-----------------------------------------------
    call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    _RC)             ! Get AGCM Heartbeat
    call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, _RC) ! set Default OCEAN_DT to AGCM Heartbeat

    CALL ESMF_TimeIntervalSet(timeStep, S=NINT(DT), _RC)
    call ESMF_ClockGet(CLOCK, currTIME=currTime, _RC)

!ALT: check with Max about moving the clock 1 step forward
    PrivateState%clock = ESMF_ClockCreate(NAME = TRIM(OCEAN_NAME)//"Clock", &
                         timeStep=timeStep, startTime=currTime, _RC)

! Initialize the Ocean Model.
!
!  This verifies the Grid and the Time against the private restarts.
!
!  Verifying that the Ocean grid and decomposition match those
!  inherited by the component. This is simply asserting that the
!  local im, jm, and lm are the same and making sure its internal
!  communication is consistent with the VM in the host's grid, probably
!  by initializing the ocean's internal communication with the
!  communicator that comes from the VM in the Grid's Layout.
!
!  Note thet tha bulletin board states , GIM(:) and GEX(:), have been created and
!  populated with no data fields. The ESMF arrays will be filled later.

!-----------------------------------------------------------------------

! Get sizes from my internal state
!---------------------------------
    call MAPL_Get(STATE, IM=IM, JM=JM, LM=LM, _RC)

! Once we know we have a valid  ESMF grid, we can call MAPL_GenericInitialize.
! This will allow us to use the built-in checkpoint/restart for our states.
!----------------------------------------------------------------------------

    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, _RC )
    call MAPL_TimerOn(STATE,"TOTAL"     )

    if(DO_DATASEA==0) then
       call MAPL_GetPointer(EXPORT, MASKO, 'MASKO'  , alloc=.true., _RC)

       select case (trim(OCEAN_NAME))
          case ("MIT")
             call MAPL_GetPointer(GEX(OCN), MASK3D, 'MIT_3D_MASK', _RC)
             MASK => MASK3D(:,:,1)
          case ("MOM6")
             call MAPL_GetPointer(GEX(OCN), MASK, 'MOM_2D_MASK', _RC)
       end select
       if(associated(MASKO)) MASKO = MASK
    end if

    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_TimerOff(STATE,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

! ========================================================

!BOP

! !IROUTINE: Run        -- Run method for ExternalModel wrapper

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
    type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: STATE
    type (ESMF_Time)                    :: EndTime
    type (ESMF_Time)                    :: MyTime,ct
    type (T_PrivateState),    pointer   :: PrivateSTATE
    type (T_PrivateState_Wrap)          :: WRAP
    type (ESMF_GridComp    ), pointer   :: GCS(:)
    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)

! Pointers to Imports

    real, pointer :: FROCEAN(:,:)
    real, pointer :: TAUXi(:,:)
    real, pointer :: TAUYi(:,:)
    real, pointer :: PENUVRi(:,:)
    real, pointer :: PENPARi(:,:)
    real, pointer :: PENUVFi(:,:)
    real, pointer :: PENPAFi(:,:)
    real, pointer :: DRNIRi(:,:)
    real, pointer :: DFNIRi(:,:)
    real, pointer :: HEATi(:,:,:)
    real, pointer :: DISCHARGEi(:,:)
    real, pointer :: CALVINGi(:,:)
    real, pointer :: LWFLXi(:,:)
    real, pointer :: SHFLXi(:,:)
    real, pointer :: QFLUXi(:,:)
    real, pointer :: SNOWi(:,:)
    real, pointer :: RAINi(:,:)
    real, pointer :: FHOCN(:,:)
    real, pointer :: FRESH(:,:)
    real, pointer :: FSALT(:,:)
    real, pointer :: PEN_OCN(:,:)

! Pointers to Exports

    real, pointer :: TS_FOUND (:,:)
    real, pointer :: SS_FOUND (:,:)
    real, pointer :: FRZMLTe(:,:)
    real, pointer :: T_Freeze_e(:,:)

! Diagnostics exports

    real, pointer :: RFLUX (:,:)
    real, pointer :: TAUXe (:,:)
    real, pointer :: TAUYe (:,:)
    real, pointer :: HEATe (:,:,:)
    real, pointer :: FROCEANe (:,:)
    real, pointer :: DISCHARGEe(:,:)
    real, pointer :: CALVINGe(:,:)
    real, pointer :: LWFLXe(:,:)
    real, pointer :: SWFLXe(:,:)
    real, pointer :: SHFLXe(:,:)
    real, pointer :: QFLUXe(:,:)
    real, pointer :: RAINe(:,:)
    real, pointer :: SNOWe(:,:)
    real, pointer :: SFLXe(:,:)
    real, pointer :: PEN_OCNe(:,:)

! Pointers to imports of child

    real, pointer :: TAUX(:,:)
    real, pointer :: TAUY(:,:)
    real, pointer :: PENUVR(:,:)
    real, pointer :: PENPAR(:,:)
    real, pointer :: PENUVF(:,:)
    real, pointer :: PENPAF(:,:)
    real, pointer :: DRNIR(:,:)
    real, pointer :: DFNIR(:,:)
    real, pointer :: HEAT(:,:,:)
    real, pointer :: DISCHARGE(:,:)
    real, pointer :: CALVING(:,:)
    real, pointer :: LWFLX(:,:)
    real, pointer :: SHFLX(:,:)
    real, pointer :: QFLUX(:,:)
    real, pointer :: RAIN(:,:)
    real, pointer :: SNOW(:,:)
    real, pointer :: SFLX(:,:)
    real, pointer :: FI(:,:)
    real, pointer :: FId(:,:)

! Pointers to exports of child

    real, pointer :: TW  (:,:)
    real, pointer :: SW  (:,:)
    real, pointer ::   MASK(:,:)
    real, pointer :: MASK3D(:,:,:)
    real, pointer :: FRZMLT(:,:)
    real, pointer :: T_Freeze(:,:) ! in deg C

    real, pointer :: TWd  (:,:)
    real, pointer :: DEL_TEMP (:,:)

! Locals

    integer           :: I,J,L
    integer           :: IM
    integer           :: JM
    integer           :: LM
    integer           :: NUM
    real, allocatable :: WGHT(:,:)
    real, pointer     :: WGHTi(:,:)
    real              :: DT, TAU_SST
    real              :: TAU_SST_UNDER_ICE
    real, pointer     :: LONS  (:,:)
    real, pointer     :: LATS  (:,:)
    real, parameter   :: OrphanSalinity=34.0  ! SA: ought to revisit this in context of OGCM rewrite. In reality one should update S as well.
    real              :: Tfreeze

    integer :: ID
    integer :: PHASE
    integer, allocatable :: PREDICTOR_CHLD(:)
    real, pointer        :: FR(:,:)    => null()
    real, pointer        :: FRI(:,:,:) => null()
    character(len=ESMF_MAXSTR) :: replayMode

    __Iam__('Run')

! Get the component's name and set-up traceback handle.
! -----------------------------------------------------

    call ESMF_GridCompGet( gc, NAME=COMP_NAME, currentPhase=PHASE, _RC)
    Iam = trim(COMP_NAME)//'::'//'Run'

    if (PHASE >= 10) PHASE = PHASE - 10 ! to be replaced by MAPL get_phase
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, _RC)

! Profilers
!----------

    call MAPL_TimerOn (STATE,"RUN"  )
    call MAPL_TimerOn (STATE,"TOTAL")

! Get child's import ad export to use as a bulletin board
!--------------------------------------------------------
    call MAPL_Get(STATE,    &
         GCS       = GCS,   &
         GIM       = GIM,   &
         GEX       = GEX,   &
         LONS      = LONS,  &
         LATS      = LATS,  &
         IM        = IM,    &
         JM        = JM,    &
         LM        = LM,    &
         _RC)


! Check the clocks to set set-up the "run-to" time
!-------------------------------------------------

    call ESMF_ClockGet( CLOCK, currTime=endTime, _RC)

! Get ocean model's private internal state
!---------------------------------

    CALL ESMF_UserCompGetInternalState( GC, TRIM(OCEAN_NAME)//'_internal_state', WRAP, STATUS )
    VERIFY_(STATUS)

    PrivateSTATE => WRAP%PTR

    call ESMF_ClockGet( PrivateState%CLOCK, currTime=myTime, _RC)

    if (myTime > EndTime) then
       call ESMF_ClockSet(PrivateState%Clock,direction=ESMF_DIRECTION_REVERSE, _RC)
       do
         call ESMF_ClockAdvance(PrivateState%Clock, _RC)
         call ESMF_ClockGet(PrivateState%Clock,currTime=ct, _RC)
         if (ct==endTime) exit
       enddo
       call ESMF_ClockSet(PrivateState%Clock, direction=ESMF_DIRECTION_FORWARD, _RC)
       call ESMF_ClockGet(PrivateState%CLOCK, currTime=myTime, _RC)
    end if

    if( MyTime <= EndTime ) then ! Time to run

! We get the ocean-land mask (now computed in Initialize of Plug)
! ---------------------------------------------------------------
       if(DO_DATASEA==0) then
          select case(trim(OCEAN_NAME))
             case ("MIT")
                call MAPL_GetPointer(GEX(OCN), MASK3D, 'MIT_3D_MASK', _RC)
                MASK => MASK3D(:,:,1)
             case ("MOM6")
                call MAPL_GetPointer(GEX(OCN), MASK, 'MOM_2D_MASK', _RC)
             end select
       else
          allocate(MASK3D(IM,JM,LM), __STAT__)
          MASK3D=1.0
          allocate(MASK(IM,JM), __STAT__)
          MASK=1.0
       end if

! Get ocean time step and misc. parameters
!-----------------------------------------

       call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    _RC)             ! Get AGCM Heartbeat
       call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, _RC) ! set Default OCEAN_DT to AGCM Heartbeat

! Get pointers to imports
!--------------------------------------------------------------------------------
       call MAPL_GetPointer(IMPORT, FROCEAN,    'FROCEAN',   _RC)
       call MAPL_GetPointer(IMPORT, TAUXi,      'TAUX'   ,   _RC)
       call MAPL_GetPointer(IMPORT, TAUYi,      'TAUY'   ,   _RC)
       call MAPL_GetPointer(IMPORT, PENUVRi,    'PENUVR' ,   _RC)
       call MAPL_GetPointer(IMPORT, PENPARi,    'PENPAR' ,   _RC)
       call MAPL_GetPointer(IMPORT, PENUVFi,    'PENUVF' ,   _RC)
       call MAPL_GetPointer(IMPORT, PENPAFi,    'PENPAF' ,   _RC)
       call MAPL_GetPointer(IMPORT, DRNIRi,     'DRNIR'  ,   _RC)
       call MAPL_GetPointer(IMPORT, DFNIRi,     'DFNIR'  ,   _RC)
       call MAPL_GetPointer(IMPORT, HEATi,      'SWHEAT',    _RC)
       call MAPL_GetPointer(IMPORT, DISCHARGEi, 'DISCHARGE', _RC)
       call MAPL_GetPointer(IMPORT, LWFLXi,     'LWFLX'  ,   _RC)
       call MAPL_GetPointer(IMPORT, SHFLXi,     'SHFLX'  ,   _RC)
       call MAPL_GetPointer(IMPORT, QFLUXi,     'QFLUX'  ,   _RC)
       call MAPL_GetPointer(IMPORT, SNOWi,      'SNOW'   ,   _RC)
       call MAPL_GetPointer(IMPORT, RAINi,      'RAIN'   ,   _RC)
       call MAPL_GetPointer(IMPORT, FHOCN,      'FHOCN'  ,   _RC)
       call MAPL_GetPointer(IMPORT, FRESH,      'FRESH'  ,   _RC)
       call MAPL_GetPointer(IMPORT, FSALT,      'FSALT'  ,   _RC)

       if(DO_DATA_ATM4OCN) then
          call MAPL_GetPointer(IMPORT, CALVINGi,  'CALVING', _RC)
       endif

! Get pointers from ImExState
!----------------------------
       if(DO_DATASEA==0) then
          call MAPL_GetPointer(GIM(OCN), TAUX,      'TAUX'  ,    _RC)
          call MAPL_GetPointer(GIM(OCN), TAUY,      'TAUY'  ,    _RC)
          call MAPL_GetPointer(GIM(OCN), PENUVR,    'PENUVR'  ,  _RC)
          call MAPL_GetPointer(GIM(OCN), PENPAR,    'PENPAR'  ,  _RC)
          call MAPL_GetPointer(GIM(OCN), PENUVF,    'PENUVF'  ,  _RC)
          call MAPL_GetPointer(GIM(OCN), PENPAF,    'PENPAF'  ,  _RC)
          call MAPL_GetPointer(GIM(OCN), DRNIR,     'DRNIR'  ,   _RC)
          call MAPL_GetPointer(GIM(OCN), DFNIR,     'DFNIR'  ,   _RC)
          call MAPL_GetPointer(GIM(OCN), HEAT,      'SWHEAT',    _RC)
          call MAPL_GetPointer(GIM(OCN), DISCHARGE, 'DISCHARGE', _RC)
          call MAPL_GetPointer(GIM(OCN), CALVING,   'CALVING',   _RC)
          call MAPL_GetPointer(GIM(OCN), LWFLX,     'LWFLX'  ,   _RC)
          call MAPL_GetPointer(GIM(OCN), SHFLX,     'SHFLX'  ,   _RC)
          call MAPL_GetPointer(GIM(OCN), QFLUX,     'QFLUX'  ,   _RC)
          call MAPL_GetPointer(GIM(OCN), RAIN,      'RAIN'  ,    _RC)
          call MAPL_GetPointer(GIM(OCN), SNOW,      'SNOW'  ,    _RC)
          call MAPL_GetPointer(GIM(OCN), SFLX,      'SFLX'  ,    _RC) ! and do not add import of PEN_OCN here since it is not used in the `plug'
       end if

       call MAPL_GetPointer(IMPORT,   PEN_OCN, 'PEN_OCN', _RC)

       call MAPL_GetPointer(GEX(OCN), TW,   'TW'  , alloc=.true., _RC)
       call MAPL_GetPointer(GEX(OCN), SW,   'SW'  , alloc=.true., _RC)

       if (dual_ocean) then
          call MAPL_GetPointer(GEX(OCNd), TWd, 'TW'  ,    alloc=.true., _RC)
          call MAPL_GetPointer(IMPORT,    FId, 'FRACICEd' ,             _RC)
       endif

       if(DO_DATASEA==0) then
          call MAPL_GetPointer(GEX(OCN), FRZMLT,   'FRZMLT', alloc=.true., _RC)
          if (trim(OCEAN_NAME) == "MOM6") then
           call MAPL_GetPointer(GEX(OCN), T_Freeze, 'T_Freeze',alloc=.true.,_RC)
          else
           nullify(T_Freeze)
          end if
       end if

! Get pointers to exports
!--------------------------------------------------------

       call MAPL_GetPointer(EXPORT, TS_FOUND,'TS_FOUND', _RC)
       call MAPL_GetPointer(EXPORT, SS_FOUND,'SS_FOUND', _RC)
       call MAPL_GetPointer(EXPORT, FRZMLTe, 'FRZMLT',   _RC)
       ! T_Freeze is only MOM6 in the StateSpecs
       if (trim(OCEAN_NAME) == "MOM6") then
          call MAPL_GetPointer(EXPORT, T_Freeze_e, 'T_Freeze', _RC)
       else
          nullify(T_Freeze_e)
       end if

! Diagnostics exports
!---------------------------------------------------------
       call MAPL_GetPointer(EXPORT, RFLUX,      'RFLUX' ,    _RC)
       call MAPL_GetPointer(EXPORT, FROCEANe,   'FROCEAN',   _RC)
       call MAPL_GetPointer(EXPORT, TAUXe,      'TAUX'   ,   _RC)
       call MAPL_GetPointer(EXPORT, TAUYe,      'TAUY'   ,   _RC)
       call MAPL_GetPointer(EXPORT, HEATe,      'SWHEAT' ,   _RC)
       call MAPL_GetPointer(EXPORT, DISCHARGEe, 'DISCHARGE', _RC)
       call MAPL_GetPointer(EXPORT, CALVINGe,   'CALVING',   _RC)
       call MAPL_GetPointer(EXPORT, LWFLXe,     'LWFLX'  ,   _RC)
       call MAPL_GetPointer(EXPORT, SWFLXe,     'SWFLX'  ,   _RC)
       call MAPL_GetPointer(EXPORT, SHFLXe,     'SHFLX'  ,   _RC)
       call MAPL_GetPointer(EXPORT, QFLUXe,     'QFLUX'  ,   _RC)
       call MAPL_GetPointer(EXPORT, RAINe,      'RAIN'  ,    _RC)
       call MAPL_GetPointer(EXPORT, SNOWe,      'SNOW'  ,    _RC)
       call MAPL_GetPointer(EXPORT, SFLXe,      'SFLX'  ,    _RC)
       call MAPL_GetPointer(EXPORT, PEN_OCNe,   'PEN_OCN',   _RC)


       if(associated(FROCEANe)) FROCEANe = FROCEAN

! Allocate space for temporary arrays
!------------------------------------

       allocate(WGHT(IM,JM), __STAT__)

! Weight for ocean grid
!----------------------
       wght=0.0
       where(MASK>0 .and. FROCEAN /= MAPL_UNDEF)
          WGHT = FROCEAN/MASK
       elsewhere
          WGHT = 0.0
       end where

#ifdef BUILD_MIT_OCEAN
       call MAPL_GetPointer(GIM(OCN), WGHTi, 'WGHT', __RC__)
       WGHTi = WGHT
#endif

       if(DO_DATASEA==0) then
! Copy imports into ImEx variables
!---------------------------------
          PENUVR    = PENUVRi    * WGHT
          PENPAR    = PENPARi    * WGHT
          PENUVF    = PENUVFi    * WGHT
          PENPAF    = PENPAFi    * WGHT
          DRNIR     = DRNIRi     * WGHT
          DFNIR     = DFNIRi     * WGHT
          if(DO_DATA_ATM4OCN) then
             DISCHARGE = DISCHARGEi
             CALVING   = CALVINGi
             ! may not be needed
             where(DISCHARGE < 0.0)
                DISCHARGE = 0.0
             endwhere
             where(CALVING < 0.0)
                CALVING = 0.0
             endwhere
          else
             DISCHARGE = DISCHARGEi * WGHT
             CALVING   = 0.0
          endif
          LWFLX     = LWFLXi     * WGHT
          QFLUX     = QFLUXi     * WGHT
          SHFLX     = (SHFLXi-FHOCN) * WGHT
          RAIN      = (RAINi+FRESH)  * WGHT
          SNOW      = SNOWi      * WGHT
          SFLX      = FSALT      * WGHT

! This stress forces the ocean, combined with sea ice bottom stress later
!------------------------------------------------------------------------
          TAUX = TAUXi * WGHT
          TAUY = TAUYi * WGHT

! Prepare radiative heating for ocean
!------------------------------------

          if(associated(RFLUX )) RFLUX  = 0.0
          select case (trim(OCEAN_NAME))
             case ("MIT", "DATASEA")
                do L=1,LM
                   HEAT(:,:,L) = HEATi(:,:,L)*WGHT
                   if(associated(RFLUX)) then
                      RFLUX = RFLUX + (1.0-MASK3D(:,:,L))*HEAT(:,:,L)
                   end if
                end do
             case ("MOM6")
                ! No 3D Mask from MOM6. Do nothing for now!
          end select

          if (associated(HEATe)) HEATe = HEAT
          if (associated(TAUXe)) TAUXe = TAUX
          if (associated(TAUYe)) TAUYe = TAUY
          if (associated(DISCHARGEe)) DISCHARGEe = DISCHARGE
          if (associated(CALVINGe)) CALVINGe = CALVING
          if (associated(LWFLXe)) LWFLXe = LWFLX
          if (associated(SWFLXe)) SWFLXe = PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR
          if (associated(SHFLXe)) SHFLXe = SHFLX
          if (associated(QFLUXe)) QFLUXe = QFLUX
          if (associated(RAINe)) RAINe = RAIN
          if (associated(SNOWe)) SNOWe = SNOW
          if (associated(SFLXe)) SFLXe = SFLX
          if (associated(PEN_OCNe)) PEN_OCNe = PEN_OCN
       end if !DO_DATASEA

! Loop the ocean model
!---------------------

       NUM = 0
       do while ( MyTime <= endTime )

! Run ocean for one time step (DT)
!---------------------------------

          call MAPL_TimerOff(STATE,"TOTAL")
          call MAPL_TimerOn (STATE,"--ModRun")

          if (.not. DUAL_OCEAN) then
             call MAPL_GenericRunChildren(GC, IMPORT, EXPORT, PrivateState%CLOCK, _RC)
          else
             if (PHASE == 1) then
                ! corrector
                call ESMF_GridCompRun( GCS(OCNd), importState=GIM(OCNd), &
                                       exportState=GEX(OCNd), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCNd, CLOCK=CLOCK, _RC )
                call ESMF_GridCompRun( GCS(OCN), importState=GIM(OCN), &
                     exportState=GEX(OCN), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCN, CLOCK=CLOCK, _RC )
                VERIFY_(STATUS)
             else
                ! predictor
                call ESMF_GridCompRun( GCS(OCNd), importState=GIM(OCNd), &
                     exportState=GEX(OCNd), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCNd, CLOCK=CLOCK, _RC )
             end if
          end if

          if (DUAL_OCEAN .and. PHASE == 1) then
             ! calculate temperature correction to send back to MOM
             call MAPL_GetPointer(GIM(OCN), DEL_TEMP, 'DEL_TEMP', _RC)
             call MAPL_GetPointer(GIM(OCNd), FI ,     'FRACICE' , _RC)

             call MAPL_GetResource(STATE,TAU_SST,           Label="TAU_SST:",           default=432000.0, _RC)
             call MAPL_GetResource(STATE,TAU_SST_UNDER_ICE, Label="TAU_SST_UNDER_ICE:", default=86400.0 , _RC)

             ! we should have valid pointers to TW and TWd by now

             DEL_TEMP = 0.0 ! we do not want uninitiazed variables
             where(MASK > 0.0 .and. FI < 0.05)
                ! what about relaxation
                DEL_TEMP = (TWd - TW)*DT/(DT+TAU_SST)
             end where

             where(MASK > 0.0 .and. FI >= 0.05 .and. FId > FI)
                ! 0.054 (C/psu) is the ratio between the freezing temperature and salinity of brine.
                ! -0.054*SW gives salinity dependent freezing temperature
                ! ideally this const should be from the ocean model, but doing so is difficult here
                DEL_TEMP = ((-0.054*SW+MAPL_TICE) - TW)*DT/(DT+TAU_SST_UNDER_ICE)

             end where

             ! put it back to MOM
             call ESMF_GridCompRun( GCS(OCN), importState=GIM(OCN), &
               exportState=GEX(OCN), clock=CLOCK, phase=2, userRC=STATUS )
             VERIFY_(STATUS)
          end if

          call MAPL_TimerOff(STATE,"--ModRun")
          call MAPL_TimerOn (STATE,"TOTAL")

! Bump the time in the internal state
!------------------------------------

          call ESMF_ClockAdvance( PrivateState%clock,                    _RC)
          call ESMF_ClockGet    ( PrivateState%clock, currTime= myTime , _RC)

          NUM = NUM + 1
       end do

       if(associated(SS_FOUND)) then
          SS_FOUND = OrphanSalinity
          where(WGHT > 0.0)
             SS_FOUND = SW
          end where
       end if

       if(associated(FRZMLTe)) then
          if(DO_DATASEA == 0) then
! assume frzmlt filled inside and get from ocean export
             where(WGHT > 0.0 )
                FRZMLTe = FRZMLT
             end where
          else
             FRZMLTe = 0.0
          end if
       end if

       if(associated(T_Freeze_e)) then
          if(DO_DATASEA == 0) then
             where(WGHT > 0.0 )
                T_Freeze_e = T_Freeze
             end where
          else
             T_Freeze_e = -1.8
          end if
       end if

       if (DUAL_OCEAN) then
          !ALT we might not have FI yet, so let get it again
          call MAPL_GetPointer(GIM(OCNd), FI , 'FRACICE'  , _RC)
          where(WGHT > 0.0)
             where(FI < 0.05)
                TS_FOUND = TWd
             elsewhere
                TS_FOUND = TW
             end where
          end where
       else
          where(WGHT > 0.0)
             TS_FOUND = TW
          end where
       end if

! Update orphan points
       if(DO_DATASEA == 0) then
          WGHT=FROCEAN*(1.0-MASK)
          Tfreeze=MAPL_TICE-0.054*OrphanSalinity ! in K

          where(wght>0.0)
             TS_FOUND=TS_FOUND+ &
                      DT*(LWFLXi+(PENUVRi+PENPARi+PENUVFi+PENPAFi+DRNIRi+DFNIRi - PEN_OCN)-SHFLXi-QFLUXi*MAPL_ALHL-MAPL_ALHF*SNOWi+FHOCN)/(OrphanDepth*MAPL_RHO_SEAWATER*MAPL_CAPWTR) ! explicit update in time
             FRZMLTe = (Tfreeze - TS_FOUND) * (MAPL_RHO_SEAWATER*MAPL_CAPWTR*OrphanDepth)/DT
             TS_FOUND=max(TS_FOUND, Tfreeze)
          end where
       end if

       deallocate(WGHT, __STAT__)

       if(DO_DATASEA/=0) then
          deallocate(MASK3D, __STAT__)
          deallocate(MASK,   __STAT__)
       end if

    end if ! Time to run

! Profilers
!----------

    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOff(STATE,"RUN"  )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

end module GEOS_OceanGridCompMod
