! $Id$
!
!  This module sets up nonlinear streaming instability with
!  multiple particle species.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub, only: dragforce_equi_multispecies 
!
  implicit none
!
  include '../initial_condition.h'
!
! Module Variables
!
  real, dimension(npar_species) :: taus = 0.0, eps0 = 0.0, vpx0 = 0.0, vpy0 = 0.0
  real :: eta_vK = 0.0, ux0 = 0.0, uy0 = 0.0
!
! Input Parameters
!
  logical :: ltaus_log_center = .true.
  real :: logtausmin = -4.0, logtausmax = -1.0
  real :: dlnndlntaus = -4.0
  real :: dlnrhodlnr = -0.1
  real :: si_kx = 0.0, si_kz = 0.0
  complex :: si_ev = (0.0, 0.0)
!
  namelist /initial_condition_pars/ &
    ltaus_log_center, logtausmin, logtausmax, dlnndlntaus, dlnrhodlnr, si_kx, si_kz, si_ev
!
  contains
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
! Initialize any module variables which are parameter dependent.
!
! 10-mar-20/ccyang: coded
!
      use EquationOfState, only: cs0
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      integer :: i
      real :: dlogtaus
!
      call keep_compiler_quiet(f)
!
! Assemble stopping times.
!
      dlogtaus = (logtausmax - logtausmin) / real(npar_species)
      gettaus: if (ltaus_log_center) then
        taus = logtausmin + real((/ (i, i = 1, npar_species) /) - 0.5) * dlogtaus
        if (lroot) print *, "initialize_initial_condition: log(taus) = ", taus
        taus = 10.0**taus
      else gettaus
        taus = 0.5 * 10.0**logtausmin * (10.0**(real((/ (i, i = 0, npar_species - 1) /)) * dlogtaus) + &
                                         10.0**(real((/ (i, i = 1, npar_species) /)) * dlogtaus))
        if (lroot) print *, "initialize_initial_condition: taus = ", taus
      endif gettaus
!
! Evaluate the radial pressure gradient support.
!
      eta_vK = -0.5 * dlnrhodlnr * cs0
      if (lroot) print*, 'initialize_initial_condition: eta * v_K = ', eta_vK
!
! Find the density ratio for each species.
!
      eps0 = taus**(4.0 + dlnndlntaus)
      eps0 = eps_dtog / sum(eps0) * eps0
!
! Find the equilibrium velocities.
!
      call dragforce_equi_multispecies(npar_species, taus, eps0, eta_vK, vpx0, vpy0, ux0, uy0)
      eqvel: if (lroot) then
        print *, "initialize_initial_condition: ux0, uy0 = ", ux0, uy0
        print *, "initialize_initial_condition: vpx0 = ", vpx0
        print *, "initialize_initial_condition: vpy0 = ", vpy0
      endif eqvel
!
! Save the equilibrium velocities to a file.
!
      open(10, file="data/multisp_drag_eq.dat", form="unformatted", action="write")
      write(10) ux0, uy0, vpx0, vpy0
      close(10)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
! Initialize the velocity field.
!
!  10-mar-20/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
! Assign the equilibrium velocities to the gas
!
      f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux) + ux0
      f(l1:l2,m1:m2,n1:n2,iuy) = f(l1:l2,m1:m2,n1:n2,iuy) + uy0
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_xxp(f, fp)
!
! Initialize particles' positions.
!
! 10-mar-20/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(inout) :: fp
!
      integer, dimension(npar_species) :: ip
      integer :: npps, npx, npz, ix, iz, is, k
      real :: dxp, dzp, xp, yp, zp, dxp1, dzp1
      real :: ar, ai, a1, a2, a3, c1, c2
      real :: argx, argz
      real :: sinp, sinm, cosp, cosm
      real :: sinp2, sinm2, cosp2, cosm2
      real :: sin2kz, cos2kx, sin2kx
!
      call keep_compiler_quiet(f)
!
! Sanity checks
!
      proc: if (mod(npar, ncpus) /= 0 .or. npar_loc /= npar / ncpus) then
        if (lroot) print *, "initialize_initial_condition: npar, ncpus, npar_loc = ", npar, ncpus, npar_loc
        call fatal_error("initialize_initial_condition", "particles are not evenly distributed among processors")
      endif proc
!
      species: if (mod(npar_loc, npar_species) /= 0) then
        if (lroot) print *, "initialize_initial_condition: npar_loc, npar_species = ", npar_loc, npar_species
        call fatal_error("initialize_initial_condition", "particle species are not evenly divided")
      endif species
!
! Find the initial ID for each species.
!
      npps = npar_loc / npar_species
      ip = (/ ((k - 1) * (npar / npar_species) + iproc * npps, k = 1, npar_species) /)
!
! Find the spacing between particles.
!
      getnp: if (nzgrid > 1) then
        npx = nint(sqrt(Lxyz_loc(1) * real(npps) / Lxyz_loc(3)))
        npz = npps / npx
      else getnp
        npx = npps
        npz = 1
      endif getnp
!
      grid: if (npx * npz /= npps) then
        if (lroot) print *, "initialize_initial_condition: Lx_loc, Lz_loc = ", Lxyz_loc(1), Lxyz_loc(3)
        if (lroot) print *, "initialize_initial_condition: npps, npx, npz = ", npps, npx, npz
        call fatal_error("initialize_initial_condition", "cannot find equal spacing between particles")
      endif grid
!
      dxp = Lxyz_loc(1) / real(npx)
      dzp = Lxyz_loc(3) / real(npz)
      if (lroot) print *, "initialize_initial_condition: npx, npz = ", npx, npz
      if (lroot) print *, "initialize_initial_condition: dxp, dzp = ", dxp, dzp
!
! Compute repeated constants.
!
      c1 = si_kx**2 + si_kz**2
      c2 = c1**2
      coeff: if (c1 > 0.0) then
        c1 = 0.5 / c1
        c2 = 1.0 / c2
      endif coeff
!
      ar = real(si_ev)
      ai = aimag(si_ev)
      a1 = 0.25 * (ar**2 - ai**2)
      a2 = 0.5 * ar * ai
      a3 = 0.25 * (ar**2 + ai**2)
!
! Assign particle positions and IDs.
!
      k = 0
      yp = xyz0(2) + 0.5 * Lxyz(2)
!
      zloop: do iz = 1, npz
        zp = xyz0_loc(3) + (real(iz) - 0.5) * dzp
        argz = si_kz * zp
        sin2kz = sin(2.0 * argz)
!
        xloop: do ix = 1, npx
          xp = xyz0_loc(1) + (real(ix) - 0.5) * dxp
          argx = si_kx * xp
          cos2kx = cos(2.0 * argx)
          sin2kx = sin(2.0 * argx)
!
          sinp = sin(argx + argz)
          sinm = sin(argx - argz)
          cosp = cos(argx + argz)
          cosm = cos(argx - argz)
!
          sinp2 = sin(2.0 * (argx + argz))
          sinm2 = sin(2.0 * (argx - argz))
          cosp2 = cos(2.0 * (argx + argz))
          cosm2 = cos(2.0 * (argx - argz))
!
          sloop: do is = 1, npar_species
            dxp1 = -c1 * si_kx * (ar * (sinp + sinm) + ai * (cosp + cosm) &
                                - a1 * (sinp2 + sinm2) - a2 * (cosp2 + cosm2)) &
                   + c2 * si_kx**3 * (a2 * cos2kx + a1 * sin2kx)
            dzp1 = -c1 * si_kz * (ar * (sinp - sinm) + ai * (cosp - cosm) &
                                - a1 * (sinp2 - sinm2) - a2 * (cosp2 - cosm2)) &
                   + c2 * si_kz**3 * a3 * sin2kz
!
            k = k + 1
            fp(k,ixp) = xp + dxp1
            fp(k,iyp) = yp
            fp(k,izp) = zp + dzp1
            ip(is) = ip(is) + 1
            ipar(k) = ip(is)
          enddo sloop
        enddo xloop
      enddo zloop
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f, fp)
!
!  Initialize particles' mass and velocity.
!
!  10-mar-20/ccyang: coded
!
      use EquationOfState, only: rho0
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(inout) :: fp
!
      real, dimension(:), pointer :: tausp_species, tausp1_species
!
      real, dimension(npar_species) :: rhopj, vpx, vpy
!
      integer :: k, p
      real :: ux, uy, hgas, zp
!
      call keep_compiler_quiet(f)
!
!  Find the mass of each particle.
!
      rhopj = rho0 / real(npar / (npar_species * nxgrid * nygrid * nzgrid)) * eps0
      if (lroot) print *, "initial_condition_vvp: rhopj = ", rhopj
!
!  Assign the mass and velocity of each particle.
!
      ploop: do k = 1, npar_loc
        p = npar_species * (ipar(k) - 1) / npar + 1
        fp(k,irhopswarm) = rhopj(p)
        fp(k,ivpx) = fp(k,ivpx) + vpx0(p) 
        fp(k,ivpy) = fp(k,ivpy) + vpy0(p)
      enddo ploop
!
!  Override the stopping times in particles_dust.
!
      multisp: if (npar_species > 1) then
        call get_shared_variable("tausp_species", tausp_species)
        call get_shared_variable("tausp1_species", tausp1_species)
        tausp_species = taus / omega
        tausp1_species = omega / taus
        if (lroot) print *, "initial_condition_vvp: override tausp_species = ", tausp_species
      endif multisp
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
