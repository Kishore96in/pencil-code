! $Id: nograv.f90,v 1.21 2003-04-04 05:46:55 brandenb Exp $

module Gravity

!
!  Dummy model: no gravity
!

  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
  endinterface

  real :: z1,z2,zref,zgrav,gravz,zinfty  !(used by Entropy and Density)
  real :: lnrho_bot,ss_bot
  character (len=labellen) :: grav_profile='const'  !(used by Density)
  logical :: lself_gravity=.false.

  integer :: dummy              ! We cannot define empty namelists
  namelist /grav_init_pars/ dummy
  namelist /grav_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  9-jan-02/wolf: coded
! 28-mar-02/axel: adapted from grav_z
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_grav called twice')
      first = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nograv.f90,v 1.21 2003-04-04 05:46:55 brandenb Exp $")
!
      lgrav = .false.
      lgravz = .false.
      lgravr = .false.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up some variables for gravity; do nothing in nograv
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed from setup_grav
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!   9-jan-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df)
!
!  add nothing to duu/dt
!
! 28-mar-02/axel: adapted from grav_z
!
      use Cdata
!      use Mpicomm
      use Sub
      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
!
      if(ip==0) print*,f,df  !(keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz,pot,pot0)
!
!  gravity potential
!  28-mar-02/axel: adapted from grav_z
!
      use Cdata, only: mx,my,mz,lroot
!
      real, dimension (mx,my,mz) :: xx,yy,zz,pot
      real, optional :: pot0
!
      if (lroot) print*,'potential: should not have been called'
      pot = 0.
      pot0 = 0.
!
      if(ip==0) print*,xx(1,1,1),yy(1,1,1),zz(1,1,1)
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  gravity potential
!  28-mar-02/axel: adapted from grav_z
!
      use Cdata, only: nx,lroot
!
      real, dimension (nx) :: xmn,pot
      real :: ymn,zmn
      real, optional :: pot0
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
!
      if (lroot) print*,'potential: should not have been called'
      pot = 0.
      pot0 = 0.
!
      if(ip==0) print*,xmn,ymn,zmn,rmn,grav
    endsubroutine potential_penc
!***********************************************************************
    subroutine self_gravity(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mvar) :: f
!
      if(ip==0) print*,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine self_gravity
!***********************************************************************

endmodule Gravity
