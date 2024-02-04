#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_DataSea -- A fake ocean surface

! !INTERFACE:

module GEOS_DataSeaGridCompMod

! !USES: 

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  logical          :: ocean_extData
  logical          :: ocean_sssData

! !DESCRIPTION:
! 
!   {\tt GEOS\_DataSea} is a gridded component that reads the 
!   ocean\_bcs file 
!   This module interpolates the SST and SSS data from 
!   either daily or monthly values to the correct time of the simulation.
!   Data are read only if the simulation time is not in the save interval.
!   It also sets surface currents US and VS to constant value (=0.) in 
!   Data Mode.
!

!EOP

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!  !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.
!
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases
    type (MAPL_MetaComp    ), pointer   :: MAPL => null()

!=============================================================================

    __Iam__('SetServices')

!****************************************************************************
! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(comp_name)//'::'//'SetServices'

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run, _RC)

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_GetResource (MAPL,   ocean_extData, Label="OCEAN_EXT_DATA:",   DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.

    ! This new SSS data feature will be ExtData based ONLY; No SSS data would set sss=30, as it was done with binary SST data
    call MAPL_GetResource (MAPL,   ocean_sssData,  Label="OCEAN_SSS_DATA:",  DEFAULT=.FALSE.,   _RC ) ! .TRUE. or .FALSE.

! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:
#include "GEOS_DataSea_Import___.h"

!  !Export state:
#include "GEOS_DataSea_Export___.h"

!EOS

    call MAPL_TimerAdd(GC,    name="RUN"     ,_RC)
    call MAPL_TimerAdd(GC,    name="-UPDATE" ,_RC)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices ( GC, _RC)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run stage for the DataSea component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )


! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the SST, sea ice and optionally SSS from dataset(s).

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_Time)                    :: CurrentTime
  character(len=ESMF_MAXSTR)          :: DATASeaFILE
  integer                             :: IFCST
  logical                             :: FCST
  integer                             :: adjSST
  real, pointer, dimension(:,:)       :: SST   => null()
  real, pointer, dimension(:,:)       :: SSS   => null()
  integer                             :: IM
  integer                             :: JM
  real                                :: TICE
  real                                :: CTB  ! Ocean-ice turbulent mixing coefficient (m/sec)
  real                                :: DT
  real                                :: RUN_DT

  real, pointer, dimension(:,:)       :: TNEW   => null()
  real, pointer, dimension(:,:)       :: F1     => null()

! Pointers to imports and exports
#include "GEOS_DataSea_DeclarePointer___.h"

   __Iam__('Run')

!****************************************************************************
!  Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, _RC)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )

! Get pointers to imports and exports
!------------------------------------
#include "GEOS_DataSea_GetPointer___.h"

! Set current time and calendar
!------------------------------

   call ESMF_ClockGet(CLOCK, currTime=CurrentTime, _RC)

   if (.not. ocean_extData) then
     ! Get the SST bcs file name from the resource file
     ! -------------------------------------------------
     call MAPL_GetResource(MAPL,DATASeaFILE,LABEL="DATA_SST_FILE:", RC=STATUS)
     VERIFY_(STATUS)
   endif

! In atmospheric forecast mode we do not have future SST and SSS
!--------------------------------------------------------------

   call MAPL_GetResource(MAPL,IFCST, LABEL="IS_FCST:",        default=0,    _RC)
   call MAPL_GetResource(MAPL,adjSST,LABEL="SST_ADJ_UND_ICE:",default=0,    _RC)

   FCST = IFCST==1

   call MAPL_Get(MAPL, IM=IM, JM=JM, _RC)

! SST is usually Reynolds/OSTIA SST or bulk SST
!------------------------------------------------

   allocate(SST(IM,JM), __STAT__)

! SSS is usually bulk SSS
!-------------------------

  if (ocean_sssData) then
    allocate(SSS(IM, JM), __STAT__)
  endif

!  Update data
!-------------

   call MAPL_TimerOn(MAPL,"-UPDATE" )

!  Read bulk SST from retrieval
!------------------------------

   if (ocean_extData) then
     sst = data_sst ! netcdf variable

     if (ocean_sssData) then  ! and bulk SSS (from retrieval)
       SSS = data_sss ! netcdf variable
     endif

   else ! binary
     call MAPL_ReadForcing(MAPL,'SST',DATASeaFILE, CURRENTTIME, sst, INIT_ONLY=FCST, _RC)
   endif

   call MAPL_TimerOff(MAPL,"-UPDATE" )

!  Update the exports
!--------------------

   if(associated(UW)) UW = 0.0
   if(associated(VW)) VW = 0.0

   TICE   = MAPL_TICE-1.8

   if (adjSST == 1) then
      SST = max(SST, TICE)
      SST = (1.-FRACICE)*SST+FRACICE*TICE
   endif

   if (adjSST == 2) then

      call MAPL_GetResource(MAPL,CTB    , LABEL="CTB:"   , default=1.0e-4, _RC)
      call MAPL_GetResource(MAPL,RUN_DT , LABEL="RUN_DT:",                 _RC)
      call MAPL_GetResource(MAPL,DT     , LABEL="DT:"    , default=RUN_DT, _RC)

      allocate(TNEW(size(TW,1),size(TW,2)), __STAT__)
      allocate(F1  (size(TW,1),size(TW,2)), __STAT__)

      TNEW=0.0
      F1  =0.0

!     ! SST below freezing point is set to freezing temperature
      TNEW   = max( SST,TICE)

      where(FRACICE == 1.0)
!     ! if fraction of ice is 1, set SST to freezing temperature        
        TNEW   =  TICE
      elsewhere
        F1=FRACICE*CTB/(2.0*(1.0-FRACICE))
        TNEW=(TNEW+TICE*F1*DT)/(1.0+F1*DT)
      end where

      SST = TNEW

      deallocate( TNEW, __STAT__)
      deallocate( F1,   __STAT__)
   endif

   if(associated(TW)) then
        TW = SST        ! SA: SST is in deg Kelvin, hence no need for abs(SST)
   end if

   if(associated(SW)) then
     if (ocean_sssData) then
       SW = SSS        ! SA: every SST data point must have SSS (in PSU) as well
     else
       SW = 30.0       ! SA: for now
     end if
   end if

! Clean-up
!---------

   deallocate(SST,     __STAT__)
   if (ocean_sssData) then
     deallocate(SSS,   __STAT__)
   endif

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN"  )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)
end subroutine RUN

end module GEOS_DataSeaGridCompMod
