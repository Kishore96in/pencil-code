! $Id$
!
!  Special module that damps the velocity field outside a user-specified cuboid.
!  Useful in cases where you do not want waves reflected from the boundary to
!  interfere with the interior of the domain.
!
!  22-Oct-2024: Kishore G. Added.
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
  use Sub, only: step
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error, not_implemented
!
  implicit none
!
  include '../special.h'
!
  real :: x_1=impossible, x_2=impossible !corners of the cuboid outside which things should be damped
  real :: y_1=impossible, y_2=impossible
  real :: z_1=impossible, z_2=impossible
  real :: tau=1 !timescale over which the velocity should be damped
  real :: w=0 !width of the step function
  logical :: ldamp_rho=.false. !whether to damp the density to its initial profile
  logical :: ldamp_ss=.false. !whether to damp the entropy to its initial profile
!
  real, dimension (mx,my,mz) :: tauinv_prof
  real, dimension (mz) :: rho_prof, ss_prof
!
! run parameters
  namelist /special_run_pars/ &
    x_1, x_2, y_1, y_2, z_1, z_2, tau, w, ldamp_rho, ldamp_ss
!
  contains
! !***********************************************************************
    subroutine initialize_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (x_1==impossible) x_1=xyz0(1)
      if (y_1==impossible) y_1=xyz0(2)
      if (z_1==impossible) z_1=xyz0(3)
      if (x_2==impossible) x_2=xyz1(1)
      if (y_2==impossible) y_2=xyz1(2)
      if (z_2==impossible) z_2=xyz1(3)
!
      tauinv_prof = spread(spread(step(x,x_1,-w)+step(x,x_2,w),2,my),3,mz) &
                   *spread(spread(step(y,y_1,-w)+step(y,y_2,w),1,mx),3,mz) &
                   *spread(spread(step(z,z_1,-w)+step(z,z_2,w),1,mx),2,my) &
                   /tau
!
      if (ldamp_rho) then
        if (ilnrho==0) call fatal_error('initialize_special', 'could not find density variable to be damped')
        if (ldensity_nolog) call not_implemented('initialize_special', 'damping rho with ldensity_nolog=T')
        if (lreference_state) call not_implemented('initialize_special', 'damping rho with lreference_state=T')
        rho_prof = exp(f(l1,m1,:,ilnrho))
      endif
!
      if (ldamp_ss) then
        if (iss==0) call fatal_error('initialize_special', 'could not find entropy variable to be damped')
        if (pretend_lnTT) call not_implemented('initialize_special', 'damping entropy with pretend_lnTT=T')
        if (lreference_state) call not_implemented('initialize_special', 'damping entropy with lreference_state=T')
        ss_prof = f(l1,m1,:,iss)
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
      lpenc_requested(i_uu)=.true.
      if (ldamp_rho) lpenc_requested(i_rho)=.true.
      if (ldamp_ss) lpenc_requested(i_ss)=.true.
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
      df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - tauinv_prof(l1:l2,m,n)*p%uu(:,1)
      df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - tauinv_prof(l1:l2,m,n)*p%uu(:,2)
      df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) - tauinv_prof(l1:l2,m,n)*p%uu(:,3)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
!
      if (ldamp_rho) then
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - tauinv_prof(l1:l2,m,n)*(p%rho/rho_prof(n) - 1)
      endif
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
!
      if (ldamp_ss) then
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - tauinv_prof(l1:l2,m,n)*(p%ss - ss_prof(n))
      endif
!
    endsubroutine special_calc_energy
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
