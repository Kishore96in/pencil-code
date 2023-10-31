! $Id$
!
!  Add Coriolis force in the beta plane approximation to the hydro equation.
!  This assumes the expansion is being done about the equator.
!  The x-axis is in the direction of increasing colatitude, while the
!  z-axis is the radial direction.
!
!  31-Oct-2023: Kishore G. Added.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! NOTE: Omega is already defined in cdata.f90
  real :: R = 1 !Radius of the sphere
!
! run parameters
  namelist /special_run_pars/ &
    Omega, R
!
  contains
!***********************************************************************
    subroutine pencil_criteria_special
!
      lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
!
      df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + 2*Omega*x(l1:l2)*p%uu(:,2)/R
!
      df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + 2*Omega*p%uu(:,3) &
                          - 2*Omega*x(l1:l2)*p%uu(:,1)/R
!
      df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) - 2*Omega*p%uu(:,2)
!
    endsubroutine special_calc_hydro
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
