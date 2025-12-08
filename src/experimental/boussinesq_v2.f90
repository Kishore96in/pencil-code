! $Id$
!
! 23-mar-2012/dintrans: coded
!
! Solve the Poisson equation for pressure when using the Boussinesq
! approximation (the so-called ``projection method'').
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! CPARAM logical, parameter :: lanelastic = .false.
! CPARAM logical, parameter :: lboussinesq = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!
! PENCILS PROVIDED rho; lnrho; rho1; glnrho(3); del2rho; del2lnrho
! PENCILS PROVIDED hlnrho(3,3); grho(3); glnrho2
! PENCILS PROVIDED del6lnrho; uij5glnrho(3); uglnrho; ugrho; sglnrho(3)
! PENCILS PROVIDED ekin; transprho
!
! PENCILS PROVIDED glnrhos(3)
! PENCILS PROVIDED totenergy_rel
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  logical :: lcalc_lnrhomean=.false.,lcalc_glnrhomean=.false.,lupw_lnrho=.false., &
             lremove_mean_temperature=.false.
!
  logical :: lwrite_debug=.false.
!
  real, dimension (mz) :: lnrhomz
  real, dimension (nz) :: glnrhomz
  real :: dx_2, dz_2
  real, pointer :: chi
  real, dimension(3) :: beta_glnrho_global=0., beta_glnrho_scaled=0.
!
  include '../density.h'
!
  namelist /density_run_pars/ lwrite_debug, lremove_mean_temperature
!
  real, pointer :: Pr
  integer :: igdu=0 !index for grad(div(u)) in the f-array
!
  contains
!***********************************************************************
    subroutine register_density
!
      use FArrayManager, only: farray_register_auxiliary
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call farray_register_auxiliary('gdu',igdu,communicated=.true., vector=3)
      if (lsphere_in_a_box) lgravr=.true.

      call put_shared_variable('beta_glnrho_scaled',beta_glnrho_scaled,caller='register_density')
      if (.not.lentropy) call put_shared_variable('beta_glnrho_global',beta_glnrho_global)
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState, only: select_eos_variable
      use DensityMethods, only: initialize_density_methods
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Boussinesq not implemented for entropy
!
      if (lenergy) then
        if (lentropy) call not_implemented("initialize_density","Boussinesq for entropy")
!
!  Boussinesq only implemented for ltemperature_nolog
!
        if (.not.ltemperature_nolog) call not_implemented("initialize_density", &
                                          "Boussinesq for log ltemperature")
      endif
!
!  Tell the equation of state that we're here and we don't have a
!  variable => isochoric (constant density).
!
      call select_eos_variable('lnrho',-1)   ! the same for rho?
!
!  For the implicit solver
!
      dx_2=1./dx**2
      dz_2=1./dz**2
!
      if (.not.lviscosity) then
        call get_shared_variable('Pr',Pr,caller='initialize_density')
        if (lroot) print*, 'Boussinesq: Pr=', Pr
      endif
!
      call initialize_density_methods
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!       f(:,:,:,igdu:igdu+2)=0.
!
!  Test of the Poisson solver
!
!      do n=n1,n2; do m=m1,m2
!        f(l1:l2,m,n,ipp)=-2.*sin(x(l1:l2))*cos(z(n))
!      enddo; enddo
!      write(10) f(l1:l2,m1:m2,n1:n2,ipp)
!      if (iorder_z==2) then
!        call inverse_laplacian_z_2nd(f(l1:l2,m1:m2,n1:n2,ipp))
!      else
!        call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))
!      endif
!      write(11) f(l1:l2,m1:m2,n1:n2,ipp)
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine read_density_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_run_pars, IOSTAT=iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine read_density_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_run_pars)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine density_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine density_after_boundary
!***********************************************************************
    subroutine pencil_criteria_density
!
!  All pencils that the Density module depends on are specified here.
!
      lpenc_requested(i_graddivu) = .true.
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ekin)) lpencil_in(i_u2)=.true.
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: lnrho0, rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f, p
!
      real, dimension (nx,3) :: gdu
!
! rho
      if (lpencil(i_rho)) p%rho=rho0
! lnrho
      if (lpencil(i_lnrho)) p%lnrho=lnrho0
! rho1
      if (lpencil(i_rho1)) p%rho1=1/rho0
! glnrho
      if (lpencil(i_glnrho)) p%glnrho=0.0
! grho
      if (lpencil(i_grho)) p%grho=0.0
! del6lnrho
      if (lpencil(i_del6lnrho)) p%del6lnrho=0.0
! hlnrho
      if (lpencil(i_hlnrho)) p%hlnrho=0.0
! sglnrho
      if (lpencil(i_sglnrho)) p%sglnrho=0.0
! uglnrho
      if (lpencil(i_uglnrho)) p%uglnrho=0.0
! ugrho
      if (lpencil(i_ugrho)) p%ugrho=0.0
! uij5glnrho
      if (lpencil(i_uij5glnrho)) p%uij5glnrho=0.0
! ekin
      if (lpencil(i_ekin)) p%ekin=0.5*p%u2
!
!     Populate auxiliary variables
!
      f(l1:l2,m,n,igdu:igdu+2) = p%graddivu
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine density_before_boundary(f)
!
      use Sub, only: remove_mean
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!     
      if (lremove_mean_temperature) call remove_mean(f,iTT)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine split_update_density(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_density
!***********************************************************************
    subroutine impose_density_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_density
!***********************************************************************
    subroutine get_slices_pressure(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pressure
!***********************************************************************
    subroutine get_init_average_density(f,init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density
!
    real, dimension (mx,my,mz,mfarray):: f
    real:: init_average_density
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine density_after_mn(f, df, mass_per_proc)
!
!  14-dec-09/dintrans: coded
!
      use Poisson, only: inverse_laplacian
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension(1), intent(in) :: mass_per_proc
!
      real, dimension (nx,ny,nz,3) :: correction
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(mass_per_proc)
!
      correction = f(l1:l2,m1:m2,n1:n2,igdu:igdu+2)
      call inverse_laplacian(correction)
      df(l1:l2,m1:m2,n1:n2,iux:iuz) = df(l1:l2,m1:m2,n1:n2,iux:iuz) - correction
!
    endsubroutine density_after_mn
!***********************************************************************
    subroutine dynamical_diffusion(uc)
!   
!  dummy routine
!  
      real, intent(in) :: uc
!  
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine boussinesq(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      call keep_compiler_quiet(f)
!
    endsubroutine boussinesq
!***********************************************************************
    function mean_density(f)
!
!  Return mean density as rho0 from eos
!
!  1-mar-15/MR: derived from mean_density in density
!
      use EquationOfState, only: rho0
!
      real :: mean_density
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      mean_density=rho0
!
      call keep_compiler_quiet(f)
!
    endfunction mean_density
!***********************************************************************
    subroutine update_char_vel_density(f)
!
!  Updates characteristic veelocity for slope-limited diffusion.
!  Most likely not yet a good method.
!
!  21-oct-15/MR: coded
!
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      if (lslope_limit_diff) f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) &
                            =f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) + rho0**2
!
    endsubroutine update_char_vel_density
!***********************************************************************
    subroutine impose_density_ceiling(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      call keep_compiler_quiet(f)

    endsubroutine impose_density_ceiling
!***********************************************************************
    subroutine calc_diagnostics_density(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_density
!***********************************************************************
    subroutine density_before_boundary_diagnostics(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
! 
      call keep_compiler_quiet(f)

    endsubroutine density_before_boundary_diagnostics
!***********************************************************************s
    subroutine write_z_stratification(f)

      real, dimension (mx,my,mz,mfarray) :: f
      call keep_compiler_quiet(f)

    endsubroutine write_z_stratification
!***********************************************************************
    subroutine load_variables_to_gpu_hydro
    endsubroutine load_variables_to_gpu_hydro
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(:) :: p_par

    call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Density
