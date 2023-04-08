! $Id$
!
! Adds a Fourier-filtered Lorentz force to the hydro equation (meant to be used along with magnetic by setting lkinematic=T and ljxb_as_aux=T there)
!
! April 2023: Kishore Gopalakrishnan (created by modifying nospecial.f90)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 6
!
!! Auxiliary variables required: 3 for jxb_filt, and 3 for jxb
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
! Temporary variables used while filtering jxb
  real, dimension(:,:,:,:), allocatable :: jxb_re, jxb_im
  
! Run parameters
  real :: k_min=-1, k_max=10
  
  namelist /special_run_pars/ &
    k_min, k_max
!
  contains
!***********************************************************************
  subroutine register_special
!
!  Set up indices for variables in special modules.
!
    use FArrayManager, only: farray_register_auxiliary
    
    if (lroot) call svn_id( &
          "$Id$")
    
    call farray_register_auxiliary('jxb_filt',ispecialvar,vector=3)
!
  endsubroutine register_special
!***********************************************************************
  subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
    use FArrayManager, only: farray_index_by_name
    
    real, dimension (mx,my,mz,mfarray) :: f
    
    if (farray_index_by_name('jxb') == -1) then
      call fatal_error("initialize_special", "Could not find jxb in f-array. Please set ljxb_as_aux=T in magnetic_init_pars.")
    endif
    
    if (.not.allocated(jxb_re)) allocate(jxb_re(nx,ny,nz,3))
    if (.not.allocated(jxb_im)) allocate(jxb_im(nx,ny,nz,3))
    
    call keep_compiler_quiet(f)
!
  endsubroutine initialize_special
!***********************************************************************
  subroutine finalize_special(f)
!
!  Called right before exiting.
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    deallocate(jxb_re)
    deallocate(jxb_im)
    
    call keep_compiler_quiet(f)
!
  endsubroutine finalize_special
!***********************************************************************
  subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
    lpenc_requested(i_rho1)=.true.
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
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
    integer :: ivec
    
    do ivec=1,3
      df(l1:l2,m,n,iuu+ivec-1) = df(l1:l2,m,n,iuu+ivec-1) + f(l1:l2,m,n,ispecialvar+ivec-1)*p%rho1
    enddo

  endsubroutine special_calc_hydro
!***********************************************************************
  subroutine special_after_timestep(f,df,dt_,llast)
!
!  Modify the f and df after df is updated.
!
    use Fourier, only: fft_xyz_parallel, kx_fft2, ky_fft2, kz_fft2
    
    logical, intent(in) :: llast
    real, dimension(mx,my,mz,mfarray), intent(inout) :: f
    real, dimension(mx,my,mz,mvar), intent(inout) :: df
    real, intent(in) :: dt_
    
    integer :: ikx,iky,ikz
    real :: k
    
    jxb_re = f(l1:l2,m1:m2,n1:n2,ijxb:ijxb+2)
    jxb_im = 0
    call fft_xyz_parallel(jxb_re, jxb_im)
    
    !See comments in special/maxwell.f90/compute_bb_from_aak_and_eek
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k = sqrt(kx_fft2(ikx+ipx*nx) + ky_fft2(iky+ipy*ny) + kz_fft2(ikz+ipz*nz))
          if (k<k_min .or. k>=k_max) then
            jxb_re(ikx,iky,ikz,:) = 0
            jxb_im(ikx,iky,ikz,:) = 0
          endif
        enddo
      enddo
    enddo
    
    call fft_xyz_parallel(jxb_re, jxb_im, linv=.true., lneed_im=.false.)
    
    f(l1:l2,m1:m2,n1:n2,ispecialvar:ispecialvar+2) = jxb_re
    
    call keep_compiler_quiet(df)
    call keep_compiler_quiet(dt_)
    call keep_compiler_quiet(llast)
!
  endsubroutine  special_after_timestep
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
  include '../special_dummies.inc'
!********************************************************************
!
endmodule Special
