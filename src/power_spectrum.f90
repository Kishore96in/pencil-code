! $Id$
!
! MODULE_DOC: reads in full snapshot and calculates power spetrum of u
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpower_spectrum = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
!
!***************************************************************
!    3-sep-02/axel+nils: coded
!    5-sep-02/axel: loop first over all points, then distribute to k-shells
!   23-sep-02/nils: adapted from postproc/src/power_spectrum.f90
!   14-mar-06/axel: made kx,ky,kz going only in integers. Works only for cubes.
!   11-nov-10/MR: intro'd flags for shell integration and z integration,
!   for that, changed namelist run_pars and corresp. read and write subroutines;
!   corresp. changes at the moment only in effect in power_xy
!
module power_spectrum
!
  use Cdata
  use Messages, only: svn_id, warning, fatal_error
!
  implicit none
!
  include 'power_spectrum.h'
!
  real :: pdf_max=30., pdf_min=-30., pdf_max_logscale=3.0, pdf_min_logscale=-3.
  logical :: lintegrate_shell=.true., lintegrate_z=.true., lcomplex=.false.
  logical :: lhalf_factor_in_GW=.false., lcylindrical_spectra=.false.
  logical :: lcyl_polar_spectra=.false.
  integer :: legendre_lmax=1.
  integer :: firstout = 0
!
  character (LEN=linelen) :: ckxrange='', ckyrange='', czrange=''
  integer, dimension(3,nk_max) :: kxrange=0, kyrange=0
  integer, dimension(3,nz_max) :: zrange=0
  integer :: n_spectra=0
  integer :: inz=0, n_segment_x=1, ndelx
!
  namelist /power_spectrum_run_pars/ &
      lintegrate_shell, lintegrate_z, lcomplex, ckxrange, ckyrange, czrange, &
      lcylindrical_spectra, inz, n_segment_x, lhalf_factor_in_GW, &
      pdf_max, pdf_min, pdf_min_logscale, pdf_max_logscale, &
      lcyl_polar_spectra, legendre_lmax
!
  contains
!***********************************************************************
    subroutine initialize_power_spectrum
!
      !!! the following warnings should become fatal errors
      if (nxgrid > nx) call warning ('power_spectrum', &
          "Part of the high-frequency spectrum are lost because nxgrid/= nx.")
      if (((dx /= dy) .and. ((nxgrid-1)*(nxgrid-1) /= 0)) .or. &
          ((dx /= dz) .and. ((nxgrid-1)*(nzgrid-1) /= 0))) &
          call warning ('power_spectrum', &
          "Shell-integration will be wrong; set dx=dy=dz to fix this.")
!
    endsubroutine initialize_power_spectrum
!***********************************************************************
    subroutine read_power_spectrum_run_pars(iostat)
!
! 05-feb-14/MR: added ordering of z ranges
! 12-mar-14/MR: changed merge_ranges into function
!
      use File_io, only: parallel_unit
      use General, only : parser, read_range, merge_ranges, quick_sort
!
      integer, intent(out) :: iostat
!
      integer :: i, iend_zrange
      character (LEN=20), dimension(nz_max) :: czranges
      integer, dimension(3) :: range
      integer, dimension(nz_max) :: iperm
      logical :: ldum
!
      read(parallel_unit, NML=power_spectrum_run_pars, IOSTAT=iostat)
      if (iostat /= 0) return
!
      kxrange(:,1) = (/1,nxgrid,1/)
      kyrange(:,1) = (/1,nygrid,1/)
!
      if ( lintegrate_shell .or. lintegrate_z ) lcomplex = .false.
!
      if ( .not.lintegrate_shell ) then
        call get_kranges( ckxrange, kxrange, nxgrid )
        call get_kranges( ckyrange, kyrange, nygrid )
      endif
!
      if ( .not.lintegrate_z ) then
!
        iend_zrange=0
        do i=1,parser( czrange, czranges, ',' )
!
          if ( read_range( czranges(i), range, (/1,nzgrid,1/) ) ) &
            ldum =  merge_ranges( zrange, iend_zrange, range )
            !!print*, 'iend_zrange, zrange(:,1:iend_zrange)=', iend_zrange, zrange(:,1:iend_zrange)
!
        enddo
!
        if (iend_zrange>0) then
!
! try further merging (not yet implemented)
!
          do i=iend_zrange-1,1,-1
          !!  ldum = merge_ranges( zrange, i, zrange(:,i+1), istore=iend_zrange )
          enddo
!
! sort ranges by ascending start value
!
          call quick_sort(zrange(1,1:iend_zrange),iperm)
          zrange(2:3,1:iend_zrange) = zrange(2:3,iperm(1:iend_zrange))
        else
!
! if no ranges specified: range = whole zgrid
!
          zrange(:,1) = (/1,nzgrid,1/)
        endif
      endif
!
      n_spectra = parser( xy_spec, xy_specs, ',' )
!
      do i=1,n_xy_specs_max
        if ( xy_specs(i) == 'u' ) then
          uxy_spec=.false.
        else if ( xy_specs(i) == 'jxb' ) then
          jxbxy_spec=.false.
        else if ( xy_specs(i) == 'b' ) then
          bxy_spec=.false.
        endif
      enddo
!
      if (uxy_spec  ) n_spectra = n_spectra+1
      if (bxy_spec  ) n_spectra = n_spectra+1
      if (jxbxy_spec) n_spectra = n_spectra+1
!
      if (n_segment_x < 1) &
        call fatal_error('read_power_spectrum_run_pars', &
                         'n_segment_x < 1')
      ndelx=nxgrid/n_segment_x

    endsubroutine read_power_spectrum_run_pars
!***********************************************************************
    subroutine get_kranges( ckrange, kranges, ngrid )
!
      use General, only : parser, read_range, merge_ranges
!
      character (LEN=*)      , intent(in) :: ckrange
      integer, dimension(:,:), intent(out):: kranges
      integer                , intent(in) :: ngrid
!
      integer :: nr, nre, i !--, ie
      character (LEN=20), dimension(size(kranges,2)) :: ckranges
!
      ckranges=''
!
      nr = parser( ckrange, ckranges, ',' )
      nre = nr
!
      do i=1,nr
!
        if ( read_range( ckranges(i), kranges(:,i), (/-ngrid/2,ngrid/2-1,1/) ) ) then
!
          if ( kranges(1,i)>=0 ) then
            kranges(1:2,i) = kranges(1:2,i)+1
          else
!
            if ( kranges(2,i)>=0 ) then
!
              if ( nre<nk_max ) then
!
                nre = nre+1
                kranges(:,nre) = (/1,kranges(2,i)+1,kranges(3,i)/)
!
                !!!call merge_ranges( kranges, i-1, kranges(:,nre) )
                !!!ldum =  merge_ranges( kranges, i-1, kranges(:,nre) )
                !!!call merge_ranges( kranges, nre-1, kranges(:,nre), nr+1 )
!
              else
                print*, 'get_kranges: Warning - subinterval could not be created!'
              endif
!
              kranges(2,i) = -1
!
            endif
!
            kranges(1:2,i) = ngrid + kranges(1:2,i) + 1
!
          endif
!
          kranges(2,i) = min(kranges(2,i),ngrid)
!
          !!!call merge_ranges( kranges, i-1, kranges(:,i) )
          !!!call merge_ranges( kranges, ie, kranges(:,i), nr+1 )
!
        endif
!
      enddo
!
    endsubroutine get_kranges
!***********************************************************************
    subroutine write_power_spectrum_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=power_spectrum_run_pars)
!
    endsubroutine write_power_spectrum_run_pars
!***********************************************************************
    subroutine power(f,sp,iapn_index)
!
!  Calculate power spectra (on spherical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fft_xyz_parallel
      use Mpicomm, only: mpireduce_sum
      use General, only: itoa
      use Sub, only: curli
!
  integer, intent(in), optional :: iapn_index
! integer, pointer :: inp,irhop,iapn(:)
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bb,oo
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  do ivec=1,3
     !
     if (trim(sp)=='u') then
        if (iuu==0) call fatal_error('power','iuu=0')
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     elseif (trim(sp)=='ud') then
        if (iuud(iapn_index)==0) call fatal_error('power','iuud=0')
        a1=f(l1:l2,m1:m2,n1:n2,iuud(iapn_index)+ivec-1)
     elseif (trim(sp)=='r2u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/2.)
     elseif (trim(sp)=='r3u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/3.)
     elseif (trim(sp)=='o') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iuu,oo,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=oo
           enddo
        enddo
     elseif (trim(sp)=='b') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bb,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=bb
           enddo
        enddo
     elseif (trim(sp)=='a') then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
     else
        print*,'There is no such sp=',trim(sp)
     endif
     b1=0.
!
!  Doing the Fourier transform
!
     call fft_xyz_parallel(a1,b1)
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
     do ikz=1,nz
        do iky=1,ny
           do ikx=1,nx
              k=nint(sqrt(kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
              if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                   +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
           enddo
        enddo
     enddo
     !
  enddo !(loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
    spectrum_sum=.5*spectrum_sum
    if (sp=='ud') then
       open(1,file=trim(datadir)//'/power_'//trim(sp)//'-'//&
            trim(itoa(iapn_index))//'.dat',position='append')
    else
      open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
    endif
!
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
    close(1)
  endif
!
    endsubroutine power
!***********************************************************************
    subroutine power_2d(f,sp)
!
!  Calculate power spectra (on circles) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fourier_transform_xz
      use Mpicomm, only: mpireduce_sum
      use Sub, only: curli
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bb
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=1) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
  do ivec=1,3
     !
    if (sp=='u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     elseif (sp=='b') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bb,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=bb
           enddo
        enddo
     elseif (sp=='a') then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
     else
        print*,'There are no such sp=',sp
     endif
     b1=0
!
!  Doing the Fourier transform
!
     !print*, 'ivec1=', ivec
     call fourier_transform_xz(a1,b1)    !!!! MR: causes error - ivec is set back from 1 to 0
     !print*, 'ivec2=', ivec
!    to be replaced by comp_spectrum( f, sp, ivec, ar, ai, fourier_transform_xz )
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over circles...'
     do ikz=1,nz
       do iky=1,ny
         do ikx=1,nx
           k=nint(sqrt(kx(ikx)**2+kz(ikz+ipz*nz)**2))
           if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
!           if (iky==16 .and. ikx==16) &
!           print*, 'power_2d:', ikx,iky,ikz,k,nk,a1(ikx,iky,ikz),b1(ikx,iky,ikz),spectrum(k+1)
         enddo
       enddo
     enddo
     !
  enddo !(loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'_2d.dat'
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power'//trim(sp)//'_2d.dat',position='append')
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
    close(1)
  endif
  !
  endsubroutine power_2d
!***********************************************************************
  subroutine comp_spectrum_xy( f, sp, ar, ai, ivecp )
!
! generates xy-spectrum of the component ivecp of the vector field, selected by sp
!
! 18-Jan-11/MR: outsourced from power_xy
!
    use Sub,      only: curli
    use General,  only: ioptest
    use Fourier,  only: fourier_transform_xy
!
    implicit none
!
    real, dimension(mx,my,mz,mfarray) :: f
    character (LEN=*)                 :: sp
    integer, optional                 :: ivecp
    real, dimension(nx,ny,nz)         :: ar, ai
!
    intent(in)  :: sp, f, ivecp
    intent(out) :: ar
    intent(out) :: ai
!
    real, dimension(nx) :: bb
    integer :: m,n,ind,ivec,i,la,le,res
!
    ivec = ioptest(ivecp,1)
!
    if (sp=='u') then
       if (iuu==0) call fatal_error('get_comp_spectrum','variable "u" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
    elseif (sp=='rho') then
       if ( ldensity_nolog ) then
         if (irho==0) call fatal_error('get_comp_spectrum','variable "rho" not existent')
         ind = irho
       else
         if (ilnrho==0) call fatal_error('get_comp_spectrum','variable "lnrho" not existent')
         ind = ilnrho
       endif
       if (ivec>1) return
       ar=f(l1:l2,m1:m2,n1:n2,ind)
    elseif (sp=='s') then
       if (iss==0) call fatal_error('get_comp_spectrum','variable "s" not existent')
       if (ivec>1) return
       ar=f(l1:l2,m1:m2,n1:n2,iss)
    elseif (sp=='b') then
        if (iaa==0) call fatal_error('get_comp_spectrum','variable "b" not existent')
        do n=n1-nghost,n2-nghost
          do m=m1-nghost,m2-nghost
             call curli(f,iaa,bb,ivec)
             ar(:,m,n)=bb
          enddo
       enddo
    elseif (sp=='a') then
       if (iaa==0) call fatal_error('get_comp_spectrum','variable "a" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
    elseif (sp=='jxb') then
       if (ijxb==0) call fatal_error('get_comp_spectrum','variable "jxb" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,ijxbx+ivec-1)
    else
       print*,'comp_spectrum_xy: Warning - There is no such sp=',sp
       return
    endif
!
    ai=0.
!
!  Doing the Fourier transform
!
    res=mod(nxgrid,n_segment_x)
    la=0
    do i=1,n_segment_x
      le=la+1; la=la+ndelx
      if (res>0) then
        la=la+1
        res=res-1
      endif
      call fourier_transform_xy(ar(la:le,:,:),ai(la:le,:,:))
    enddo
!
   return
!
   endsubroutine comp_spectrum_xy
!***********************************************************************
   subroutine power_xy(f,sp,sp2)
!
!  Calculate power spectra (on circles) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   11-nov-10/MR: extended to arbitrary combinations of shell/2d and z dependent/integrated spectra
!                 additional information about kind of spectrum + wavenumber vectors in output file;
!                 extended shell-integrated spectra to anisotropic boxes, extended k range to k_x,max^2 + k_y,max^2
!   18-jan-11/MR: modified for calculation of power spectra of scalar products
!   10-may-11/MR: modified for use with ranges in kx, ky, z; for output of true
!                 (complex and componentwise) instead of power spectra
!    5-may-14/MR: modifications for request of individual components of a vector field
!    4-nov-16/MR: correction: no k_x, k_y output for shell-integrated spectra
!
   use Mpicomm, only: mpireduce_sum, mpigather_xy, mpigather_and_out_real, mpigather_and_out_cmplx, &
                      mpimerge_1d, ipz, mpibarrier, mpigather_z
   use General, only: itoa, write_full_columns, get_range_no, write_by_ranges
!
  implicit none
!
  real, dimension(mx,my,mz,mfarray), intent(in) :: f
  character (len=*),                 intent(in) :: sp
  character (len=*), optional,       intent(in) :: sp2
!
  !integer, parameter :: nk=nx/2                      ! actually nxgrid/2 *sqrt(2.)  !!!
!
  integer :: i,il,jl,k,ikx,iky,ikz,ivec,nk,ncomp,nkx,nky,npz,nkl,iveca,cpos
  real,    dimension(nx,ny,nz)            :: ar,ai
  real,    dimension(:,:,:), allocatable  :: br,bi
  real,    allocatable, dimension(:)      :: spectrum1,spectrum1_sum, kshell
  real,    allocatable, dimension(:,:)    :: spectrum2,spectrum2_sum,spectrum2_global
  real,    allocatable, dimension(:,:,:)  :: spectrum3
  complex, allocatable, dimension(:,:,:,:):: spectrum3_cmplx
!
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nx,ny)  :: prod
  real                    :: prods
!
  character (len=80)   :: title
  character (len=fnlen):: filename
  character (len=3)    :: sp_field
  logical              :: l2nd
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  l2nd = .false.
  if (present(sp2)) l2nd = sp/=sp2
!
  if (l2nd) lcomplex = .false.
!
  if (l2nd) allocate(br(nx,ny,nz),bi(nx,ny,nz))
!
  cpos=0
!
!  to add further fields, modify here!
!
  if ( sp(1:1)=='u' .or. sp(1:1)=='b' .or. sp(1:1)=='a' .or. sp(1:1)=='s' ) then
    sp_field=sp(1:1)
    cpos=2
  elseif (len(trim(sp))>=3) then
    if ( sp(1:3)=='rho' .or. sp(1:3)=='jxb' ) then
      sp_field=sp(1:3)
      cpos=4
    endif
  endif
  if (cpos==0) &
    call fatal_error('power_xy','no implementation for field '//trim(sp))

  if ( sp_field=='u' .or. sp_field=='b' .or.  &
       sp_field=='a' .or. sp_field=='jxb' ) then  ! for vector fields
    if (len(trim(sp))>=cpos) then                 ! component specification expected
      ncomp=1
      select case (sp(cpos:cpos))
      case ('x')  ; iveca=1
      case ('y')  ; iveca=2
      case ('z')  ; iveca=3
      case default; call fatal_error('power_xy','no components other than x,y,z may be selected')
      end select
    else                                        ! no component specified -> all three components
      ncomp=3; iveca=1
    endif
  else
    ncomp=1; iveca=1
  endif
!
  if (lintegrate_shell) then
!
    title = 'Shell-integrated'
    nk = nint( sqrt( ((nxgrid+1)/Lx)**2+((nygrid+1)/Ly)**2 )*Lx/2 )+1
    allocate( kshell(nk) )
!
! To initialize variables with NaN, please only use compiler flags. (Bourdin.KIS)
!
    kshell = -1.0
!
    if (lintegrate_z) then
!
      title = trim(title)//' and z-integrated power'
      allocate( spectrum1(nk), spectrum1_sum(nk) )
!
      spectrum1=0.
      spectrum1_sum=0.
!
    else
!
      title = trim(title)//' and z-dependent power'
      allocate( spectrum2(nk,nz), spectrum2_sum(nk,nz) )
!
      if (lroot) then
        allocate( spectrum2_global(nk,nzgrid) )
      else
        allocate( spectrum2_global(1,1) )                  ! only a dummy
      endif
!
      spectrum2=0.
      spectrum2_sum=0.
!
    endif
!
  else if (lintegrate_z) then
!
    title = 'z-integrated power'
    allocate( spectrum2(nx,ny), spectrum2_sum(nx,ny) )
!
    if (lroot) then
      allocate( spectrum2_global(nxgrid,nygrid) )
    else
      allocate( spectrum2_global(1,1) )                  ! only a dummy
    endif
!
    spectrum2=0.
    spectrum2_sum=0.
!
  else
!
    title = 'z-dependent'
!
    if ( lcomplex ) then
      if ( ncomp>1 ) then
        title = trim(title)//' complex componentwise ('//trim(itoa(ncomp))//') '
      else
        title = trim(title)//' complex '
      endif
      allocate( spectrum3_cmplx(nx,ny,nz,ncomp) )
      spectrum3_cmplx=0.
    else
      title = trim(title)//' power'
      allocate( spectrum3(nx,ny,nz) )
      spectrum3=0.
    endif
!
  endif
!
  title = trim(title)//' spectrum w.r.t. x and y'
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)       !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)       !*2*pi/Ly
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part ar1-ar3; and put imaginary part, ai1-ai3, to zero
  !
  do ivec=iveca,iveca+ncomp-1
!
    call comp_spectrum_xy( f, sp_field, ar, ai, ivec )
    if (l2nd) call comp_spectrum_xy( f, sp2, br, bi, ivec )
!
!  integration over shells
!  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
!
    if (lroot .AND. ip<10) print*,'fft done; now collect/integrate over circles...'
!
!  Summing up the results from the different processors
!  The result is available only on root  !!??
!
    do ikz=1,nz
      if (lintegrate_shell) then
!
        do iky=1,ny
          do ikx=1,nx
!
            !!k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2))
            k=nint( sqrt( (kx(ikx)/Lx)**2+(ky(iky+ipy*ny)/Ly)**2 )*Lx ) ! i.e. wavenumber index k
                                                                        ! is |\vec{k}|/(2*pi/Lx)
            if ( k>=0 .and. k<=nk-1 ) then
!
              kshell(k+1) = k*2*pi/Lx
!
              if (l2nd) then
                prods = 0.5*(ar(ikx,iky,ikz)*br(ikx,iky,ikz)+ai(ikx,iky,ikz)*bi(ikx,iky,ikz))
              else
                prods = 0.5*(ar(ikx,iky,ikz)**2+ai(ikx,iky,ikz)**2)
              endif
!
              if (lintegrate_z) then
                spectrum1(k+1) = spectrum1(k+1)+prods*dz              ! equidistant grid required
              else
                spectrum2(k+1,ikz) = spectrum2(k+1,ikz) + prods
              endif
            endif
          enddo
        enddo
!
      else
!
        if (l2nd) then
          prod = ar(:,:,ikz)*br(:,:,ikz)+ai(:,:,ikz)*bi(:,:,ikz)
        elseif ( .not. lcomplex ) then
          prod = ar(:,:,ikz)**2+ai(:,:,ikz)**2
        endif
!
        if (lintegrate_z) then
          spectrum2(:,:)=spectrum2(:,:)+(0.5*dz)*prod                 ! equidistant grid required
        elseif ( lcomplex ) then
          spectrum3_cmplx(:,:,ikz,ivec-iveca+1)=cmplx(ar(:,:,ikz),ai(:,:,ikz))
        else
          spectrum3(:,:,ikz)=spectrum3(:,:,ikz)+0.5*prod
        endif
!
      endif
!
    enddo
!
  enddo !(of loop over ivec)
!
  if (lintegrate_shell .and. firstout<n_spectra .and. ipz==0) &        ! filling of the shell-wavenumber vector
    call mpimerge_1d(kshell,nk,12)
!
  if (lroot) then
!
!  on root processor, append global result to diagnostics file "power<field>_xy.dat"
!
    if ( sp2=='' ) then
      filename=trim(datadir)//'/power'//trim(sp)//'_xy.dat'
    else
      filename=trim(datadir)//'/power'//trim(sp)//'.'//trim(sp2)//'_xy.dat'
    endif
!
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',filename
!
    open(1,file=filename,position='append')
!
    if (lintegrate_shell) then
      nkl = nk
      if ( kshell(nk) == -1 ) nkl = nk-1
    endif
!
    if ( firstout<n_spectra ) then
!
      write(1,'(a)') title
!
      if (lintegrate_shell) then
!
        write(1,'(a)') 'Shell-wavenumbers k ('//trim(itoa(nkl))//'):'
        write(1,'(1p,8e15.7)') kshell(1:nkl)
!
      else

        nkx = get_range_no( kxrange, nk_max )
        nky = get_range_no( kyrange, nk_max )
!
        write(1,'(a)') 'Wavenumbers k_x ('//trim(itoa(nkx))//') and k_y ('//trim(itoa(nky))//'):'
!
        call write_by_ranges( 1, kx*2*pi/Lx, kxrange )
        call write_by_ranges( 1, ky*2*pi/Ly, kyrange )
!
      endif
!
      if (  zrange(1,1)>0 .and. &
           (zrange(1,1)>1 .or. zrange(2,1)<nzgrid .or. zrange(3,1)>1) ) then
!
        npz = get_range_no( zrange, nz_max )
!
        write(1,'(a)') 'z-positions ('//trim(itoa(npz))//'):'
        call write_by_ranges( 1, zgrid, zrange )
!
      endif
!
    endif
!
    firstout = firstout+1
!
    write(1,*) t
!
  endif
!
  if (lintegrate_shell) then
!
    if (lintegrate_z) then
      call mpireduce_sum(spectrum1,spectrum1_sum,nk)
    else
      call mpireduce_sum(spectrum2,spectrum2_sum,(/nk,nz/),12)
      call mpigather_z(spectrum2_sum,spectrum2_global,nk)
    endif
!
  else if (lintegrate_z) then
         call mpireduce_sum(spectrum2,spectrum2_sum,(/nx,ny/),3)
         call mpigather_xy( spectrum2_sum, spectrum2_global, 0 )
!
!  transposing output, as in Fourier_transform_xy; an unreverted transposition is performed
!  but no transposition when nygrid=1 (e.g., in 2-D setup for 1-D spectrum)
!
       elseif (lcomplex) then
         call mpigather_and_out_cmplx(spectrum3_cmplx,1,.not.(nygrid==1),kxrange,kyrange,zrange)
       else
         call mpigather_and_out_real(spectrum3,1,.not.(nygrid==1),kxrange,kyrange,zrange)
       endif
!
  if (lroot) then
!
    if (lintegrate_shell) then
!
      if (lintegrate_z) then
        write(1,'(1p,8e15.7)') spectrum1_sum(1:nkl)
      else
        do i=1,nz_max
          if ( zrange(1,i) > 0 ) then
            do jl=zrange(1,i), zrange(2,i), zrange(3,i)
              write(1,'(1p,8e15.7)') (spectrum2_global(il,jl), il=1,nkl)
            enddo
          endif
        enddo
!        print*, 'nach write'
      endif
!
    else
!
      if (lintegrate_z) &
        call write_by_ranges( 1, spectrum2_global, kxrange, kyrange, .true. )
                                                                     ! transposing output, as in fourier_transform_xy
                                                                     ! an unreverted transposition is performed
    endif
    close(1)
!
  endif
!
  call mpibarrier          ! necessary ?
!  print*, 'nach barrier:', iproc, ipy, ipz
!
  if (lintegrate_shell) then
!
    deallocate(kshell)
    if (lintegrate_z) then
      deallocate(spectrum1,spectrum1_sum)
    else
      deallocate(spectrum2,spectrum2_sum,spectrum2_global)
    endif
!
  else if (lintegrate_z) then
    deallocate(spectrum2,spectrum2_sum,spectrum2_global)
  elseif ( lcomplex ) then
    deallocate(spectrum3_cmplx)
  else
    deallocate(spectrum3)
  endif
  !
  endsubroutine power_xy
!***********************************************************************
  subroutine powerhel(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use Chiral, only: iXX_chiral, iYY_chiral
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, jkz, im, in, ivec, ivec_jj
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi, jji, b2, j2
  real, dimension(nx,3) :: bb, bbEP, jj
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nk,nzgrid) :: cyl_spectrum, cyl_spectrum_sum
  real, dimension(nk,nzgrid) :: cyl_spectrumhel, cyl_spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true., lwrite_krms_GWs=.false.
!
!  passive scalar contributions (hardwired for now)
!
  real, dimension(nx,3) :: gtmp1,gtmp2
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !
  if (lcylindrical_spectra) then
    cyl_spectrum=0.
    cyl_spectrum_sum=0.
    cyl_spectrumhel=0.
    cyl_spectrumhel_sum=0.
  endif
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      if (iuu==0) call fatal_error('powerhel','iuu=0')
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corresponds to vorticity)
        enddo
      enddo
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  neutral velocity power spectra (spectra of |un|^2 and on.un)
!
    elseif (sp=='neu') then
      if (iuun==0) call fatal_error('powerhel','iuun=0')
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuun,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corresponds to vorticity)
        enddo
      enddo
      b_re=f(l1:l2,m1:m2,n1:n2,iuun+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  magnetic power spectra (spectra of |B|^2 and A.B)
!
    elseif (sp=='mag') then
      if (iaa==0) call fatal_error('powerhel','iaa=0')
      if (lmagnetic) then
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,bbi,ivec)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
          enddo
        enddo
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
        a_im=0.
        b_im=0.
      else
        if (headt) print*,'magnetic power spectra only work if lmagnetic=T'
      endif
!
!  magnetic power spectra (spectra of |J|^2 and J.B) !!! should be J.A
!
    elseif (sp=='j.a') then
      if (iaa==0) call fatal_error('powerhel','iaa=0')
      if (lmagnetic) then
        do n=n1,n2
          do m=m1,m2
          call del2vi_etc(f,iaa,ivec,curlcurl=jji)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=jji  !(this corresponds to the current density)
          enddo
        enddo
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
        a_im=0.
        b_im=0.
      else
        if (headt) print*,'magnetic power spectra only work if lmagnetic=T'
      endif
!
!  current helicity spectrum (J.B)
!
    elseif (sp=='j.b') then
      if (iaa==0) call fatal_error('powerhel','iaa=0')
      if (lmagnetic) then
        do n=n1,n2
          do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          call del2vi_etc(f,iaa,ivec,curlcurl=jji)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corresponds to the magnetic field)
          b_re(:,im,in)=jji  !(this corresponds to the current density)
          enddo
        enddo
        a_im=0.
        b_im=0.
      else
        if (headt) print*,'magnetic power spectra only work if lmagnetic=T'
      endif
!
!  Gravitational wave power spectra (breathing mode; diagonal components of gij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWd') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','ihij=0 or igij=0')
      a_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec-1)  !(corresponds to hii)
      b_re=f(l1:l2,m1:m2,n1:n2,igij+ivec-1)  !(corresponds to gii)
      a_im=0.
      b_im=0.
!
!  Gravitational wave power spectra (off-diagonal components of gij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWe') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','igij=0 or igij=0')
      a_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec+2)  !(corresponds to hij)
      b_re=f(l1:l2,m1:m2,n1:n2,igij+ivec+2)  !(corresponds to gij)
      a_im=0.
      b_im=0.
!
!  Gravitational wave power spectra (breathing mode; diagonal components of hij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWf') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','ihij=0 or igij=0')
      a_re=f(l1:l2,m1:m2,n1:n2,igij+ivec-1)  !(corresponds to gii)
      b_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec-1)  !(corresponds to hii)
      a_im=0.
      b_im=0.
!
!  Gravitational wave power spectra (off-diagonal components of hij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWg') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','igij=0 or igij=0')
      a_re=f(l1:l2,m1:m2,n1:n2,igij+ivec+2)  !(corresponds to gij)
      b_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec+2)  !(corresponds to hij)
      a_im=0.
      b_im=0.
!
!  spectrum of u.b
!
    elseif (sp=='u.b') then
      if (iuu==0.or.iaa==0) call fatal_error('powerhel','iuu or iaa=0')
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='mgz') then
      if (ivec==3) then
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,bbi,ivec)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
          enddo
        enddo
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
        a_im=0.
        b_im=0.
      else
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
      endif
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='bb2') then
      if (ivec==3) then
        do n=n1,n2
          do m=m1,m2
            call curl(f,iaa,bb)
            call dot2(bb,b2)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=b2
          enddo
        enddo
        if (ilnrho/=0) then
          a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
        else
          a_re=0.
        endif
        a_im=0.
        b_im=0.
      else
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
      endif
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='jj2') then
      if (ivec==3) then
        do n=n1,n2
          do m=m1,m2
            call del2v_etc(f,iaa,curlcurl=jj)
            call dot2(jj,j2)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=j2
          enddo
        enddo
        if (ilnrho/=0) then
          a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
        else
          a_re=0.
        endif
        a_im=0.
        b_im=0.
      else
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
      endif
!
!  spectrum of uzs and s^2
!
    elseif (sp=='uzs') then
      if (ivec==3) then
        a_re=f(l1:l2,m1:m2,n1:n2,iuz)  !(this corresponds to uz)
        b_re=f(l1:l2,m1:m2,n1:n2,iss)  !(this corresponds to ss)
        a_im=0.
        b_im=0.
      else
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
      endif
!
!  magnetic energy spectra based on fields with Euler potentials
!
    elseif (sp=='bEP') then
      if (iXX_chiral/=0.and.iYY_chiral/=0) then
        do n=n1,n2
          do m=m1,m2
            call grad(f,iXX_chiral,gtmp1)
            call grad(f,iYY_chiral,gtmp2)
            call cross(gtmp1,gtmp2,bbEP)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=bbEP(:,ivec)  !(this corresponds to magnetic field)
            a_re(:,im,in)=.5*(f(l1:l2,m,n,iXX_chiral)*gtmp2(:,ivec) &
                             -f(l1:l2,m,n,iYY_chiral)*gtmp1(:,ivec))
          enddo
        enddo
        a_im=0.
        b_im=0.
      endif
!
!  Spectrum of uxj
!
    elseif (sp=='uxj') then
      do n=n1,n2
        do m=m1,m2
          if (ivec==1) ivec_jj=2
          if (ivec==2) ivec_jj=1
          if (ivec/=3) call del2vi_etc(f,iaa,ivec_jj,curlcurl=jji)
          im=m-nghost
          in=n-nghost
          if (ivec==1) b_re(:,im,in)=+jji
          if (ivec==2) b_re(:,im,in)=-jji
          if (ivec==3) b_re(:,im,in)=+0.
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      a_im=0.
      b_im=0.
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +b_re(ikx,iky,ikz)**2 &
               +b_im(ikx,iky,ikz)**2
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
!
!  allow for possibility of cylindrical spectral
!
    if (lcylindrical_spectra) then
      if (lroot .AND. ip<10) print*,'fft done; now integrate over cylindrical shells...'
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
            k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
            jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
            k=nint(sqrt(k2))
            if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
              cyl_spectrum(k+1,jkz)=cyl_spectrum(k+1,jkz) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
              cyl_spectrumhel(k+1,jkz)=cyl_spectrumhel(k+1,jkz) &
                 +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  end of loop through all points
!
            endif
          enddo
        enddo
      enddo
    endif
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !
  if (lcylindrical_spectra) then
    call mpireduce_sum(cyl_spectrum,cyl_spectrum_sum,(/nk,nzgrid/))
    call mpireduce_sum(cyl_spectrumhel,cyl_spectrumhel_sum,(/nk,nzgrid/))
  endif
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lcylindrical_spectra) then
      if (ip<10) print*,'Writing cylindrical power spectrum ',sp &
           ,' to ',trim(datadir)//'/cyl_power_'//trim(sp)//'.dat'
    !
      cyl_spectrum_sum=.5*cyl_spectrum_sum
      open(1,file=trim(datadir)//'/cyl_power_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrum_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) t
        write(1,'(1p,8e10.2)') cyl_spectrum_sum
      endif
      close(1)
      !
      open(1,file=trim(datadir)//'/cyl_powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrumhel_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) t
        write(1,'(1p,8e10.2)') cyl_spectrumhel_sum
      endif
      close(1)
    endif
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerhel
!***********************************************************************
  subroutine powerLor(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, im, in, ivec, stat
  real :: k2
  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(mx,my,mz,3) :: Lor
  real, dimension(:,:,:,:), allocatable :: tmpv, scrv
  real, dimension(:,:,:), allocatable :: c_re, c_im
  real, dimension(nx,ny,nz) :: a_re, a_im, b_re, b_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum, spectrum_sum, spectrum2, spectrum2_sum
  real, dimension(nk) :: spectrumhel, spectrumhel_sum, spectrum2hel, spectrum2hel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  Note, if lhydro=F, then f(:,:,:,1:3) does no longer contain
  !  velocity. In that case, we want the magnetic field instead.
  !
  if (.not.lhydro) then
    allocate(tmpv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate memory for tmpv')
    allocate(scrv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate memory for scrv')
    allocate(c_re(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate memory for c_re')
    allocate(c_im(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate memory for c_im')
  endif
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  spectrum2=0.
  spectrum2_sum=0.
  spectrum2hel=0.
  spectrum2hel_sum=0.
  !
  !  compute Lorentz force
  !
  do m=m1,m2
  do n=n1,n2
     aa=f(l1:l2,m,n,iax:iaz)
     call gij(f,iaa,aij,1)
     call gij_etc(f,iaa,aa,aij,bij)
     call curl_mn(aij,bb,aa)
     call curl_mn(bij,jj,bb)
     call cross_mn(jj,bb,jxb)
     Lor(l1:l2,m,n,:)=jxb
     if (.not.lhydro) tmpv(l1:l2,m,n,:)=bb
     if (.not.lhydro) scrv(l1:l2,m,n,:)=jj
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Lorentz force spectra (spectra of L*L^*)
!
    if (sp=='Lor') then
      b_re=Lor(l1:l2,m1:m2,n1:n2,ivec)
      if (lhydro) then
        a_re=f(l1:l2,m1:m2,n1:n2,ivec)
      else
        a_re=tmpv(l1:l2,m1:m2,n1:n2,ivec)
        c_re=scrv(l1:l2,m1:m2,n1:n2,ivec)
        c_im=0.
      endif
      a_im=0.
      b_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    if (.not.lhydro) call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!  Remember: a=B, b=Lor, c=J, so for nonhydro, we want a.b and c.b
!
            if (lhydro) then
              spectrum(k+1)=spectrum(k+1) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
            else
              spectrum(k+1)=spectrum(k+1) &
                 +c_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
              spectrum2(k+1)=spectrum2(k+1) &
                 +c_re(ikx,iky,ikz)**2 &
                 +c_im(ikx,iky,ikz)**2
              spectrum2hel(k+1)=spectrum2hel(k+1) &
                 +c_re(ikx,iky,ikz)*a_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*a_im(ikx,iky,ikz)
            endif
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  call mpireduce_sum(spectrum2,spectrum2_sum,nk)
  call mpireduce_sum(spectrum2hel,spectrum2hel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !  normal 2 spectra
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    !  additional 2 spectra
    !
    spectrum2_sum=.5*spectrum2_sum
    open(1,file=trim(datadir)//'/power_2'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum2_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum2_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_2'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum2hel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum2hel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  if (allocated(tmpv)) deallocate(tmpv)

  endsubroutine powerLor
!***********************************************************************
  subroutine powerLor_OLD(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec, stat
  real :: k2
  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(mx,my,mz,3) :: Lor
  real, dimension(:,:,:,:), allocatable :: tmpv, scrv
  real, dimension(:,:,:), allocatable :: c_re, c_im
  real, dimension(nx,ny,nz) :: a_re, a_im, b_re, b_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  Note, if lhydro=F, then f(:,:,:,1:3) does no longer contain
  !  velocity. In that case, we want the magnetic field instead.
  !
  if (.not.lhydro) then
    allocate(tmpv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for tmpv')
    allocate(scrv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for scrv')
    allocate(c_re(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for c_re')
    allocate(c_im(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for c_im')
  endif
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !
  !  compute Lorentz force
  !
  do m=m1,m2
  do n=n1,n2
     aa=f(l1:l2,m,n,iax:iaz)
     call gij(f,iaa,aij,1)
     call gij_etc(f,iaa,aa,aij,bij)
     call curl_mn(aij,bb,aa)
     call curl_mn(bij,jj,bb)
     call cross_mn(jj,bb,jxb)
     Lor(l1:l2,m,n,:)=jxb
     if (.not.lhydro) tmpv(l1:l2,m,n,:)=bb
     if (.not.lhydro) scrv(l1:l2,m,n,:)=jj
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Lorentz force spectra (spectra of L*L^*)
!
    if (sp=='Lor') then
      b_re=Lor(l1:l2,m1:m2,n1:n2,ivec)
      if (lhydro) then
        a_re=f(l1:l2,m1:m2,n1:n2,ivec)
      else
        a_re=tmpv(l1:l2,m1:m2,n1:n2,ivec)
        c_re=scrv(l1:l2,m1:m2,n1:n2,ivec)
        c_im=0.
      endif
      a_im=0.
      b_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    if (.not.lhydro) call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!  Remember: a=B, b=Lor, c=J, so for nonhydro, we want a.b and c.b
!
            if (lhydro) then
              spectrum(k+1)=spectrum(k+1) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
            else
              spectrum(k+1)=spectrum(k+1) &
                 +c_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
            endif
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  if (allocated(tmpv)) deallocate(tmpv)

  endsubroutine powerLor_OLD
!***********************************************************************
  subroutine powerEMF(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,3) :: EMF,JJJ,EMB,BBB
  real, dimension(nx,ny,nz) :: a_re,a_im, b_re,b_im, c_re,c_im, d_re,d_im
  real, dimension(nx,3) :: uu,aa,bb,jj,uxb,uxj
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
 real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !
  !  compute EMFentz force
  !
  do m=m1,m2
  do n=n1,n2
     uu=f(l1:l2,m,n,iux:iuz)
     aa=f(l1:l2,m,n,iax:iaz)
     call gij(f,iaa,aij,1)
     call gij_etc(f,iaa,aa,aij,bij)
     call curl_mn(aij,bb,aa)
     call curl_mn(bij,jj,bb)
     call cross_mn(uu,bb,uxb)
     call cross_mn(uu,jj,uxj)
     EMF(l1:l2,m,n,:)=uxb
     EMB(l1:l2,m,n,:)=uxj
     JJJ(l1:l2,m,n,:)=jj
     BBB(l1:l2,m,n,:)=bb
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Electromotive force spectra (spectra of L*L^*)
!
    if (sp=='EMF') then
      a_re=EMF(l1:l2,m1:m2,n1:n2,ivec)
      b_re=JJJ(l1:l2,m1:m2,n1:n2,ivec)
      c_re=EMB(l1:l2,m1:m2,n1:n2,ivec)
      d_re=BBB(l1:l2,m1:m2,n1:n2,ivec)
      a_im=0.
      b_im=0.
      c_im=0.
      d_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    call fft_xyz_parallel(c_re,c_im)
    call fft_xyz_parallel(d_re,d_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +c_re(ikx,iky,ikz)*d_re(ikx,iky,ikz) &
               +c_im(ikx,iky,ikz)*d_im(ikx,iky,ikz)
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerEMF
!***********************************************************************
  subroutine powerTra(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn, div_mn, multsv_mn, &
        h_dot_grad_vec
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,3) :: Adv, Str, BBB
  real, dimension(nx,ny,nz) :: a_re,a_im, b_re,b_im, c_re,c_im
  real, dimension(nx,3) :: uu, aa, bb, divu, bbdivu, bgradu, ugradb
  real, dimension(nx,3,3) :: uij, aij, bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !
  !  compute EMF transfer terms. Following Rempel (2014), we split
  !  curl(uxB) = -[uj*dj(Bi)+.5*Bi(divu)] -[-Bj*dj(ui)+.5*Bi(divu)]
  !            =  ---- advection --------  ------ stretching ------
  !
  do m=m1,m2
  do n=n1,n2
     uu=f(l1:l2,m,n,iux:iuz)
     aa=f(l1:l2,m,n,iax:iaz)
     call gij(f,iuu,uij,1)
     call gij(f,iaa,aij,1)
     call gij_etc(f,iaa,aa,aij,bij)
     call div_mn(uij,divu,uu)
     call curl_mn(aij,bb,aa)
     call multsv_mn(divu,bb,bbdivu)
     call h_dot_grad_vec(uu,bij,bb,ugradb)
     call h_dot_grad_vec(bb,uij,uu,bgradu)
     Adv(l1:l2,m,n,:)=+ugradb+.5*bbdivu
     Str(l1:l2,m,n,:)=-bgradu+.5*bbdivu
     BBB(l1:l2,m,n,:)=bb
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Electromotive force transfer spectra
!
    if (sp=='Tra') then
      a_re=BBB(l1:l2,m1:m2,n1:n2,ivec)
      b_re=Adv(l1:l2,m1:m2,n1:n2,ivec)
      c_re=Str(l1:l2,m1:m2,n1:n2,ivec)
      a_im=0.
      b_im=0.
      c_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*c_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*c_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerTra
!***********************************************************************
  subroutine powerGWs(f,sp,lfirstcall)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel, fourier_transform
    use Mpicomm, only: mpireduce_sum, mpigather_and_out_cmplx
    use SharedVariables, only: get_shared_variable
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
    use Special, only: special_calc_spectra
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=3) :: sp
  logical :: lfirstcall

  integer, parameter :: nk=nxgrid/2

  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real :: k2
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, allocatable, dimension(:) :: spectrum,spectrumhel
  real, allocatable, dimension(:) :: spectrum_sum,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  logical, save :: lwrite_krms_GWs=.false.
  real :: sign_switch, kk1, kk2, kk3
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! Select cases where spectra are precomputed
!
  if (iggXim>0.or.iggTim>0) then
    if (sp=='StT'.or.sp=='StX') then
      allocate(spectrum(nxgrid),spectrumhel(nxgrid))
      allocate(spectrum_sum(nxgrid),spectrumhel_sum(nxgrid))
    else
      allocate(spectrum(nk),spectrumhel(nk))
      allocate(spectrum_sum(nk),spectrumhel_sum(nk))
    endif
    call special_calc_spectra(f,spectrum,spectrumhel,lfirstcall,sp)
  else
    allocate(spectrum(nk),spectrumhel(nk))
!
!  Initialize power spectrum to zero. The following lines only apply to
!  the case where special/gravitational_waves_hij6.f90 is used.
!
    k2m=0.
    nks=0.
    spectrum=0.
    spectrumhel=0.
!
!  Define wave vector, defined here for the *full* mesh.
!  Each processor will see only part of it.
!  Ignore *2*pi/Lx factor, because later we want k to be integers
!
    kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
    ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
    kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
!
!  Gravitational wave tensor (spectra of g*g^* for gT and gX, where g=hdot)
!
    if (sp=='GWs') then
      if (iggX>0.and.iggT>0.and.iggXim==0.and.iggTim==0) then
        a_re=f(l1:l2,m1:m2,n1:n2,iggX)
        b_re=f(l1:l2,m1:m2,n1:n2,iggT)
      else
        call fatal_error('powerGWs','must have lggTX_as_aux=T')
      endif
      a_im=0.
      b_im=0.
!
!  Gravitational wave tensor (spectra of h*h^* for hT and hX)
!
    elseif (sp=='GWh') then
      if (ihhX>0.and.ihhXim==0) then
        a_re=f(l1:l2,m1:m2,n1:n2,ihhX)
        b_re=f(l1:l2,m1:m2,n1:n2,ihhT)
      else
        call fatal_error('powerGWs','must have lhhTX_as_aux=T')
      endif
      a_im=0.
      b_im=0.
!
!  Gravitational wave stress tensor (only if lStress_as_aux is requested)
!  Note: for aux_stress='d2hdt2', the stress is replaced by GW_rhs.
!
    elseif (sp=='Str') then
      if (iStressX>0.and.iStressXim==0) then
        a_re=f(l1:l2,m1:m2,n1:n2,iStressX)
        b_re=f(l1:l2,m1:m2,n1:n2,iStressT)
      else
        call fatal_error('powerGWs','must have lStress_as_aux=T')
      endif
      a_im=0.
      b_im=0.
    else
      call fatal_error('powerGWs','no valid spectrum (=sp) chosen')
    endif
!
!  Doing the Fourier transform
!
    !call fft_xyz_parallel(a_re,a_im)
    !call fft_xyz_parallel(b_re,b_im)
    call fourier_transform(a_re,a_im)
    call fourier_transform(b_re,b_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
! do ikz=1,nz
!   do iky=1,ny
!     do ikx=1,nx
    do iky=1,nz
      do ikx=1,ny
        do ikz=1,nx
          k2=kx(ikx+ipy*ny)**2+ky(iky+ipz*nz)**2+kz(ikz+ipx*nx)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  Switch sign for the same k vectors for which we also
!  switched the sign of e_X. Define (kk1,kk2,kk3) as short-hand
!
            kk1=kx(ikx+ipy*ny)
            kk2=ky(iky+ipz*nz)
            kk3=kz(ikz+ipx*nx)
!
            !kk1=kx(ikx+ipx*nx)
            !kk2=ky(iky+ipy*ny)
            !kk3=kz(ikz+ipz*nz)
!
!  possibility of swapping the sign
!
             sign_switch=1.
             if (kk3<0.) then
               sign_switch=-1.
             elseif (kk3==0.) then
               if (kk2<0.) then
                 sign_switch=-1.
               elseif (kk2==0.) then
                 if (kk1<0.) then
                   sign_switch=-1.
                 endif
               endif
             endif
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +a_re(ikz,ikx,iky)**2 &
               +a_im(ikz,ikx,iky)**2 &
               +b_re(ikz,ikx,iky)**2 &
               +b_im(ikz,ikx,iky)**2
            spectrumhel(k+1)=spectrumhel(k+1)+2*sign_switch*( &
               +a_im(ikz,ikx,iky)*b_re(ikz,ikx,iky) &
               -a_re(ikz,ikx,iky)*b_im(ikz,ikx,iky))
!
!  compute krms only once
!
            if (lwrite_krms_GWs) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
!
!  end from communicated versus computed spectra
!
  endif
!
!  open
!
  open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
!
!  XX
!
  if (sp=='StT'.or.sp=='StX') then
!
!  transposing output, as in Fourier_transform_xy; an unreverted transposition is performed
!  but no transposition when nygrid=1 (e.g., in 2-D setup for 1-D spectrum)
!
    call mpireduce_sum(spectrumhel,spectrumhel_sum,nxgrid)
    call mpireduce_sum(spectrum,spectrum_sum,nxgrid)
  else
    !
    !  Summing up the results from the different processors
    !  The result is available only on root
    !
    call mpireduce_sum(spectrum   ,spectrum_sum   ,nk)
    call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  endif
!
!  compute krms only once
!
  if (lwrite_krms_GWs) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms_GWs=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !  half factor or not?
    !  By default (lhalf_factor_in_GW=F), we have total(S) = gg2m.
    !  Otherwise we have total(S) = (1/2) * gg2m.
    !
    if (lhalf_factor_in_GW) then
      spectrum_sum=.5*spectrum_sum
      spectrumhel=.5*spectrumhel
    endif
    !
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    if ( all(sp.ne.(/'SCL','VCT','Tpq'/)) ) then
      open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do k = 1, nk
          write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
        enddo
      else
        write(1,*) t
        write(1,'(1p,8e10.2)') spectrumhel_sum
      endif
      close(1)
    endif
    !
    if (lwrite_krms_GWs) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms_GWs.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms_GWs=.false.
    endif
  endif
  !
  endsubroutine powerGWs
!***********************************************************************
  subroutine powerscl(f,sp,iapn_index,lsqrt)
!
!  Calculate power spectrum of scalar quantity (on spherical shells) of the
!  variable specified by `sp', e.g. spectra of cc, rho, etc.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
    use Fourier, only: fft_xyz_parallel
    use General, only: itoa
    use Mpicomm, only: mpireduce_sum
    use Sub, only: curli, grad
    use SharedVariables, only: get_shared_variable
!
  logical, intent(in), optional :: lsqrt
  integer, intent(in), optional :: iapn_index
  integer, pointer :: inp,irhop,iapn(:)
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz, ivec, im, in
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  real, dimension(nx) :: bbi
  real, dimension(nx,3) :: gLam
  character (len=2) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  For "kin", calculate spectra of <uk^2> and <ok.uk>
  !  For "mag", calculate spectra of <bk^2> and <ak.bk>
  !
  if (sp=='ro') then
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
  !
  !  spectrum of lnrho (or normalized enthalpy).
  !  Need to take log if we work with linear density.
  !
  elseif (sp=='lr') then
    if (ldensity_nolog) then
      a_re=alog(f(l1:l2,m1:m2,n1:n2,irho))
    else
      a_re=f(l1:l2,m1:m2,n1:n2,ilnrho)
    endif
  elseif (sp=='nd') then
    a_re=f(l1:l2,m1:m2,n1:n2,ind(iapn_index))
  elseif (sp=='np') then
    call get_shared_variable('inp', inp, caller='powerscl')
    a_re=f(l1:l2,m1:m2,n1:n2,inp)
  elseif (sp=='na') then
    call get_shared_variable('iapn', iapn, caller='powerscl')
    a_re=f(l1:l2,m1:m2,n1:n2,iapn(iapn_index))
  elseif (sp=='rp') then
    call get_shared_variable('irhop', irhop, caller='powerscl')
    a_re=f(l1:l2,m1:m2,n1:n2,irhop)
  elseif (sp=='TT') then
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
  elseif (sp=='ss') then
    a_re=f(l1:l2,m1:m2,n1:n2,iss)
  elseif (sp=='cc') then
    if (icc==0) call fatal_error('powerscl','icc=0, which is not allowed')
    a_re=f(l1:l2,m1:m2,n1:n2,icc)
  elseif (sp=='cr') then
    a_re=f(l1:l2,m1:m2,n1:n2,iecr)
  elseif (sp=='sp') then
    a_re=f(l1:l2,m1:m2,n1:n2,ispecialvar)
  elseif (sp=='Ssp') then
    a_re=sqrt(abs(f(l1:l2,m1:m2,n1:n2,ispecialvar)))
  elseif (sp=='mu') then
    a_re=f(l1:l2,m1:m2,n1:n2,ispecialvar2)
  elseif (sp=='hr') then
    a_re=0.
    do m=m1,m2
      do n=n1,n2
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=a_re(:,im,in)+bbi*f(l1:l2,m,n,iaa-1+ivec)
        enddo
      enddo
    enddo
    a_im=0.
  elseif (sp=='ha') then
    a_re=0.
    do m=m1,m2
      do n=n1,n2
        call grad(f,ispecialvar,gLam)
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=a_re(:,im,in)+bbi*(f(l1:l2,m,n,iaa-1+ivec)+&
              gLam(:,ivec))
        enddo
      enddo
    enddo
    a_im=0.
  endif
  a_im=0.
!
!  Allow for talking the square root defined for pos/neg arguments.
!
  if (present(lsqrt)) then
    a_re=sqrt(abs(a_re))*sign(a_re,1.)
  endif
!
!  Doing the Fourier transform
!
  call fft_xyz_parallel(a_re,a_im)
!
!  integration over shells
!
  if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
  do ikz=1,nz
    do iky=1,ny
      do ikx=1,nx
        k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
        if (k>=0 .and. k<=(nk-1)) then
          spectrum(k+1)=spectrum(k+1) &
             +a_re(ikx,iky,ikz)**2 &
             +a_im(ikx,iky,ikz)**2
        endif
      enddo
    enddo
  enddo
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    if (sp=='na' .or. sp=='nd') then
       open(1,file=trim(datadir)//'/power_'//trim(sp)//'-'//&
            trim(itoa(iapn_index))//'.dat',position='append')
    else
       open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    endif
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
    close(1)
  endif
  !
  endsubroutine powerscl
!***********************************************************************
  subroutine power_1d(f,sp,ivec,ivar)
!
!  Calculate power spectra of the variable specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!    27-apr-14/nishant: added inz to compute power_x at a given z
!
    use Fourier, only: fourier_transform_x
    use Mpicomm, only: mpireduce_sum, stop_it, transp
    use Sub, only: curli
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=1) :: sp
    integer :: ivec
    integer, optional :: ivar
!
    integer, parameter :: nk=nx/2
    integer :: ix,iy,iz,im,in,ikx,iky,ikz,nc
    real, dimension(nx,ny,nz) :: a1,b1,a2
    real, dimension(nx) :: bb
    real, dimension(:,:), allocatable :: spectrumx,spectrumx_sum
    real, dimension(nk) :: spectrumy,spectrumy_sum
    real, dimension(nk) :: spectrumz,spectrumz_sum
    character (len=fnlen) :: suffix
!
!  identify version
!
    if (lroot .AND. ip<10) call svn_id( &
        "$Id$")
!
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!
    if (sp=='u') then
      if (lhydro) then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
      else
        if (lroot) &
            print*, 'power_1d: must have hydro module for velocity power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='b') then
      if (lmagnetic) then
        do n=n1,n2; do m=m1,m2
          call curli(f,iaa,bb,ivec)
          im=m-nghost
          in=n-nghost
          a1(:,im,in)=bb
        enddo; enddo
      else
        if (lroot) &
            print*, 'power_1d: must have magnetic module for magnetic power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='a') then
      if (lmagnetic) then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
      else
        if (lroot) &
            print*, 'power_1d: must have magnetic module for magnetic power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='p') then
      if (present(ivar)) then
        if (ivar>0) then
          a1=f(l1:l2,m1:m2,n1:n2,ivar)
        else
          if (lroot) &
              print*, 'power_1d: ivar must be >0, ivar=', ivar
          call fatal_error('power_1d','')
        endif
      else
        call fatal_error('power_1d','ivar not set')
      endif
    else
      if (lroot) print*,'There is no such spectra variable: sp=',sp
      call fatal_error('power_1d','')
    endif
    b1=0
    a2=a1
!
    !!print*,'checking lcomplex, oned',lcomplex,oned
!
    if (lcomplex) then
      nc=2
    else
      nc=1
    endif
    allocate(spectrumx(nc,nk), spectrumx_sum(nc,nk) )

   !! print*,'nc=',nc
!
! Need to initialize
!
    spectrumx=0.
    spectrumx_sum=0.
    spectrumy=0.
    spectrumy_sum=0.
    spectrumz=0.
    spectrumz_sum=0.
!
!  Do the Fourier transform
!
    call fourier_transform_x(a1,b1)
!
!  Stop the run if FFT=nofft
!
    if (.not.lfft) &
        call stop_it('Need FFT=fft in Makefile.local to get spectra!')
!
!  Spectra in x-direction
!
!NS: added
   if (.not.lintegrate_z) then
    !print*,'NISHANT inz=',inz
    do ikx=1,nk; do iy=1,ny
      if (lcomplex) then
        spectrumx(:,ikx) = spectrumx(:,ikx) + &
            (/a1(ikx,iy,inz), b1(ikx,iy,inz)/)
      else
        spectrumx(1,ikx) = spectrumx(1,ikx) + &
            sqrt(a1(ikx,iy,inz)**2 + b1(ikx,iy,inz)**2)
      endif
    enddo; enddo
   else
    do ikx=1,nk; do iy=1,ny; do iz=1,nz
      if (lcomplex) then
        spectrumx(:,ikx) = spectrumx(:,ikx) + &
            (/a1(ikx,iy,iz), b1(ikx,iy,iz)/)
      else
        spectrumx(1,ikx) = spectrumx(1,ikx) + &
            sqrt(a1(ikx,iy,iz)**2 + b1(ikx,iy,iz)**2)
      endif
    enddo; enddo; enddo
   endif
!
!  Multiply all modes, except the constant mode, by two.
!
    spectrumx(:,2:nk)=2*spectrumx(:,2:nk)
!
!  Doing fourier spectra in all directions if onedall=T
!
    if (onedall) then
!
!  Spectra in y-direction
!
      if (nygrid/=1) then
        a1=a2
        b1=0
        call transp(a1,'y')
        call fourier_transform_x(a1,b1)
        do iky=1,nk; do ix=1,nxgrid/nprocy; do iz=1,nz
          spectrumy(iky) = spectrumy(iky) + &
              sqrt(a1(iky,ix,iz)**2 + b1(iky,ix,iz)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        spectrumy(2:nk)=2*spectrumy(2:nk)
      endif
!
!  Spectra in z-direction
!
      if (nzgrid/=1) then
        a1=a2
        b1=0
        call transp(a1,'z')
        call fourier_transform_x(a1,b1)
        do ikz=1,nk; do ix=1,nxgrid/nprocz; do iy=1,ny
          spectrumz(ikz) = spectrumz(ikz) + &
              sqrt(a1(ikz,iy,ix)**2 + b1(ikz,iy,ix)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        spectrumz(2:nk)=2*spectrumz(2:nk)
      endif
    endif
!
!  Summing up the results from the different processors
!  The result is available only on root
!
    call mpireduce_sum(spectrumx,spectrumx_sum,(/nc,nk/))
    if (onedall.and.nygrid/=1) call mpireduce_sum(spectrumy,spectrumy_sum,nk)
    if (onedall.and.nzgrid/=1) call mpireduce_sum(spectrumz,spectrumz_sum,nk)
!
!  on root processor, write global result to file
!  don't need to multiply by 1/2 to get \int E(k) dk = (1/2) <u^2>
!  because we have only taken the data for positive values of kx.
!
    if (ivec==1) then
      suffix='x_x.dat'
    elseif (ivec==2) then
      suffix='y_x.dat'
    elseif (ivec==3) then
      suffix='z_x.dat'
    else
      suffix='_x.dat'
    endif
!
!  Append to diagnostics file
!
    if (lroot) then
      if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
          'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
      open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
          position='append')
      write(1,*) t
!
      if (lcomplex) then
        write(1,'(1p,8("(",e10.2,",",e10.2,")"))') spectrumx_sum/(nygrid*nzgrid)
      else
        write(1,'(1p,8e10.2)') spectrumx_sum/(nygrid*nzgrid)
      endif
!

!
      close(1)
    endif
!
!  Save data for y and z spectra if onedall=.true.
!
    if (onedall) then
!
!  Save y data
!
      if (lroot .and. nygrid/=1) then
        if (ivec==1) then
          suffix='x_y.dat'
        elseif (ivec==2) then
          suffix='y_y.dat'
        elseif (ivec==3) then
          suffix='z_y.dat'
        else
          suffix='_y.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
            position='append')
        write(1,*) t
        write(1,'(1p,8e10.2)') spectrumy_sum/(nxgrid*nzgrid)
        close(1)
      endif
!
!  Save z data
!
      if (lroot .and. nzgrid/=1) then
        if (ivec==1) then
          suffix='x_z.dat'
        elseif (ivec==2) then
          suffix='y_z.dat'
        elseif (ivec==3) then
          suffix='z_z.dat'
        else
          suffix='_z.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*,'Writing power spectra of variable', sp,  &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
            position='append')
        write(1,*) t
        write(1,'(1p,8e10.2)') spectrumz_sum/(nxgrid*nygrid)
        close(1)
      endif
    endif
!
  endsubroutine power_1d
!***********************************************************************
    subroutine pdf(f,variabl,pdf_mean,pdf_rms)
!
!  Calculated pdf of scalar field.
!  This routine is in this module, because it is always called at the
!  same time when spectra are invoked (called in wsnaps).
!
!    2-dec-03/axel: coded
!
      use Sub, only: grad, dot2_mn
      use Mpicomm, only: mpireduce_sum_int
      use SharedVariables, only: get_shared_variable
!
  integer :: l,i_pdf
  integer, parameter :: n_pdf=3001
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (nx,3) :: gcc
  real, dimension (nx) :: pdf_var,gcc2
  integer, dimension (n_pdf) :: pdf_yy, pdf_yy_sum
  real :: pdf_mean, pdf_rms, pdf_dx, pdf_dx1, pdf_scl
  character (len=120) :: pdf_file=''
  character (len=*) :: variabl
  logical :: logscale=.false.
  integer, pointer :: ispecial
!
!  initialize counter and set scaling factor
!
   pdf_yy=0
   pdf_yy_sum=0
   pdf_scl=1./pdf_rms
!
!  m-n loop
!
   do n=n1,n2
   do m=m1,m2
!
!  select the right variable
!
     if (variabl=='rhocc') then
       pdf_var=exp(f(l1:l2,m,n,ilnrho))*f(l1:l2,m,n,ilncc)-pdf_mean
       logscale=.false.
     elseif (variabl=='cc') then
       pdf_var=f(l1:l2,m,n,ilncc)-pdf_mean
       logscale=.false.
     elseif (variabl=='lncc') then
       pdf_var=abs(f(l1:l2,m,n,ilncc)-pdf_mean)
       logscale=.true.
     elseif (variabl=='gcc') then
       call grad(f,ilncc,gcc)
       call dot2_mn(gcc,gcc2)
       pdf_var=sqrt(gcc2)
       logscale=.false.
     elseif (variabl=='lngcc') then
       call grad(f,ilncc,gcc)
       call dot2_mn(gcc,gcc2)
       pdf_var=sqrt(gcc2)
       logscale=.true.
     elseif (variabl=='special') then
       call get_shared_variable('ispecial', ispecial, caller='pdf')
       pdf_var=f(l1:l2,m,n,ispecial)
       logscale=.false.
     elseif (variabl=='lnspecial') then
       call get_shared_variable('ispecial', ispecial, caller='pdf')
       pdf_var=alog(f(l1:l2,m,n,ispecial))
       logscale=.false.
     endif
!
!  put in the right pdf slot
!
     if (logscale) then
       pdf_dx=(pdf_max_logscale-pdf_min_logscale)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*log10(pdf_scl*pdf_var(l))-pdf_min_logscale)
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure its inside array boundries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     else
       pdf_dx=(pdf_max-pdf_min)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*(pdf_scl*pdf_var(l)-pdf_min))
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure its inside array boundries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     endif
   enddo
   enddo
!
!  Communicate and append from root processor.
!
  call mpireduce_sum_int(pdf_yy,pdf_yy_sum,n_pdf)
  if (lroot) then
     pdf_file=trim(datadir)//'/pdf_'//trim(variabl)//'.dat'
     open(1,file=trim(pdf_file),position='append')
     if (logscale) then
       write(1,10) t, n_pdf, pdf_dx, pdf_max_logscale, pdf_min_logscale, pdf_mean, pdf_rms
     else
       write(1,10) t, n_pdf, pdf_dx, pdf_max, pdf_min, pdf_mean, pdf_rms
     endif
     write(1,11) pdf_yy_sum
     close(1)
  endif
!
10 format(1p,e12.5,0p,i6,1p,5e12.4)
11 format(8i10)
endsubroutine pdf
!***********************************************************************
    subroutine power_phi(f,sp)
!
! Power spectra in phi direction in spherical coordinates:
! I define power_phi of a variable 'u' in the following way:
! {\hat u}(r,\theta,k) \equiv FFT (u(r,\theta,k))
! power_phi(u) \equiv
!         \sum_{r,\theta} dr d\theta
!             {\hat u}(r,\theta,k)*{\hat u}(r,\theta,-k) r^2 sin(\theta)
! ---------------------------------------------------------------------
! As this subroutine is called at the end of a time-step df can be
! used for storing temporary data.
! The \phi direction is the z direction.
! ----------------------------------------------------------------------
!
      use Sub, only: curli
      use Mpicomm, only: stop_it, y2x, z2x
      use Fourier, only: fourier_transform_real_1
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1
  real, dimension(nx) :: bb
  real, dimension(nygrid/2) :: spectrumy,spectrumy_sum
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
  real, dimension(nygrid) :: aatempy
  real, dimension(nzgrid) :: aatemp
  real, dimension(2*nzgrid+15) :: fftpack_temp
  real :: nVol2d,spec_real,spec_imag
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!--------------Makes sense only in spherical coordinate system -----------
  if (.not.(lspherical_coords.or.lcylindrical_coords)) &
      call stop_it("power_phi works only in spherical or cylindrical coords")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  !
  nVol2d=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumy_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  do ivec=1,3
    !
    if (trim(sp)=='u') then
      a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
    elseif (trim(sp)=='b') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bb,ivec)
          im=m-nghost
          in=n-nghost
          a1(:,im,in)=bb
        enddo
      enddo
    elseif (trim(sp)=='a') then
      a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
    else
      print*,'There are no such sp=',trim(sp)
    endif
!
    ifirst_fft=1
    do l=1,nx
      if (lspherical_coords) then
        do m=1,ny
          do j=1,nprocy
            call z2x(a1,l,m,j,aatemp)
!
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
!AB: is nVol2d correctly initialized? Did this now above. OK?
!
            if (lroot) then
!             write(*,*)l,m,j,'got data shall fft'
              call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
              ifirst_fft = ifirst_fft+1
              spectrum(1)=(aatemp(1)**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              do ispec=2,nzgrid/2
                spec_real=aatemp(2*ispec-2)
                spec_imag=aatemp(2*ispec-1)
                spectrum(ispec)= 2.*(spec_real**2+spec_imag**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              enddo
              spectrum(nzgrid/2)=(aatemp(nzgrid)**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              spectrum_sum=spectrum_sum+spectrum
              nVol2d = nVol2d+r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
            else
              nVol2d=1.
            endif
          enddo ! loop over yproc
        enddo   ! loop over ny
      elseif (lcylindrical_coords) then
        do n=1,nz
          do j=1,nprocz
            call y2x(a1,l,n,j,aatempy)
!
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
!
            if (lroot) then
!             write(*,*)l,n,j,'got data shall fft'
              call fourier_transform_real_1(aatempy,nygrid,ifirst_fft,fftpack_temp)
              ifirst_fft = ifirst_fft+1
              spectrumy(1)=(aatempy(1)**2)&
                     *rcyl_weight(l)
              do ispec=2,nygrid/2
                spec_real=aatempy(2*ispec-2)
                spec_imag=aatempy(2*ispec-1)
                spectrumy(ispec)= 2.*(spec_real**2+spec_imag**2)&
                     *rcyl_weight(l)
              enddo
              spectrumy(nygrid/2)=(aatempy(nygrid)**2)&
                     *rcyl_weight(l)
              spectrumy_sum=spectrumy_sum+spectrumy
              nVol2d = nVol2d+rcyl_weight(l)
            else
              nVol2d=1.
            endif
          enddo ! loop over zproc
        enddo   ! loop over nz
      else
        call fatal_error('power_phi','neither spherical nor cylindrical')
      endif
    enddo     ! loop over nx
!
  enddo !(from loop over ivec)
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power_phi'//trim(sp)//'.dat'
    open(1,file=trim(datadir)//'/power_phi'//trim(sp)//'.dat',position='append')
    write(1,*) t
!
    if (lspherical_coords) then
      spectrum_sum=.5*spectrum_sum
      write(1,'(1p,8e10.2)') spectrum_sum/nVol2d
    elseif (lcylindrical_coords) then
      spectrumy_sum=.5*spectrumy_sum
      write(1,'(1p,8e10.2)') spectrumy_sum/nVol2d
    endif
    close(1)
  endif
  !
  endsubroutine power_phi
!***********************************************************************
  subroutine powerhel_phi(f,sp)
!
! Power spectra in phi direction in spherical coordinates:
! I define power_phi of a variable 'u' in the following way:
! {\hat u}(r,\theta,k) \equiv FFT (u(r,\theta,k))
! power_phi(u) \equiv
!         \sum_{r,\theta} dr d\theta
!             {\hat u}(r,\theta,k)*{\hat u}(r,\theta,-k) r^2 sin(\theta)
! ---------------------------------------------------------------------
! As this subroutine is called at the end of a time-step df can be
! used for storing temporary data.
! The \phi direction is the z direction.
! ----------------------------------------------------------------------
!
    use Fourier, only: fourier_transform_real_1
    use Mpicomm, only: z2x, stop_it
    use Sub, only: curli
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bbi
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
  real, dimension(nzgrid/2) :: spectrumhel,spectrumhel_sum
  real, dimension(nzgrid) :: aatemp,bbtemp
  real, dimension(2*nzgrid+15) :: fftpack_temp
  real :: nVol2d,spec_reala,spec_imaga,spec_realb,spec_imagb
  character (len=*) :: sp
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!--------------Makes sense only in spherical coordinate system -----------
  if (.not.lspherical_coords) call stop_it("powerhel_phi works only in spherical coordinates")
!
!  Define wave vector, defined here for the *full* mesh.
!  Each processor will see only part of it.
!  Ignore *2*pi/Lx factor, because later we want k to be integers
!
!
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
!
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
!
  do ivec=1,3
     !
     if (trim(sp)=='kin') then
       do n=n1,n2
         do m=m1,m2
           call curli(f,iuu,bbi,ivec)
           im=m-nghost
           in=n-nghost
           a1(:,im,in)=bbi  !(this corresponds to vorticity)
         enddo
       enddo
       b1=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1) !(this corresponds to velocity)
     elseif (trim(sp)=='mag') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bbi,ivec)
              im=m-nghost
              in=n-nghost
              b1(:,im,in)=bbi !(this corresponds to magnetic field)
           enddo
        enddo
        a1=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1) !(this corresponds to vector potential)
     else
        print*,'There are no such sp=',trim(sp)
     endif
!
     ifirst_fft=1
     do l=1,nx
       do m=1,ny
         do j=1,nprocy
           call z2x(a1,l,m,j,aatemp)
           call z2x(b1,l,m,j,bbtemp)
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
           if (lroot) then
!             write(*,*)l,m,j,'got data shall fft'
             call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
             call fourier_transform_real_1(bbtemp,nzgrid,ifirst_fft,fftpack_temp)
             ifirst_fft = ifirst_fft+1
             spectrum(1)=(bbtemp(1)*bbtemp(1))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrumhel(1)=(aatemp(1)*bbtemp(1))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             do ispec=2,nzgrid/2
               spec_reala=aatemp(2*ispec-2)
               spec_imaga=aatemp(2*ispec-1)
               spec_realb=bbtemp(2*ispec-2)
               spec_imagb=bbtemp(2*ispec-1)
               spectrum(ispec)= 2.*(spec_realb*spec_realb+spec_imagb*spec_imagb)&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
               spectrumhel(ispec)= 2.*(spec_reala*spec_realb+spec_imaga*spec_imagb)&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             enddo
             spectrumhel(nzgrid/2)=(aatemp(nzgrid)*bbtemp(nzgrid))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum(nzgrid/2)=(bbtemp(nzgrid)*bbtemp(nzgrid))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum_sum=spectrum_sum+spectrum
             spectrumhel_sum=spectrumhel_sum+spectrumhel
             nVol2d = nVol2d+r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
           else
             nVol2d=1.
           endif
         enddo ! loop over yproc
       enddo   ! loop over ny
     enddo     ! loop over nx
!
   enddo !(from loop over ivec)
!
!  append to diagnostics file
!
   if (lroot) then
     if (ip<10) print*,'Writing power spectrum ',sp &
       ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
!
     spectrum_sum=.5*spectrum_sum
     spectrumhel_sum=0.5*spectrumhel_sum
     open(1,file=trim(datadir)//'/power_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum
     close(1)
!
     open(1,file=trim(datadir)//'/powerhel_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrumhel_sum
     close(1)
   endif
  !
 endsubroutine powerhel_phi
!***********************************************************************
    subroutine power_vec(f,sp)
!
!  Calculate power spectra (on shperical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Sub, only: del2v_etc
      use Mpicomm, only: mpireduce_sum
      use Fourier, only: fourier_transform
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz,3) :: a1,b1
  real, dimension(nx,3) :: tmp_a1
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  if (trim(sp)=='j') then
     ! compute j = curl(curl(x))
     do n=n1,n2
       do m=m1,m2
         call del2v_etc(f,iaa,curlcurl=tmp_a1)
         a1(:,m-nghost,n-nghost,:) = tmp_a1
       enddo
     enddo
  else
     print*,'There are no such sp=',trim(sp)
  endif
  b1=0
!
!  Doing the Fourier transform
!
  do ivec=1,3
     call fourier_transform(a1(:,:,:,ivec),b1(:,:,:,ivec))
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
     do ikz=1,nz
        do iky=1,ny
           do ikx=1,nx
              k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
              if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                   +a1(ikx,iky,ikz,ivec)**2+b1(ikx,iky,ikz,ivec)**2
           enddo
        enddo
     enddo
     !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (lroot) then
     if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
          ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
     spectrum_sum=.5*spectrum_sum
     open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum
     close(1)
  endif
  !
  endsubroutine power_vec
!***********************************************************************
  subroutine anisoq_diag(f)
!
!  Anisotropic alpha2 dynamo diagostics.
!  Calculate azimuthally averaged spectra of
!  fluid (v.u), kinetic (u.o), magnetic (a.b), and current (b.j) helicities.
!  Here v is the vector potential of u if u is incompressible.
!
!  16-sep-20/hongzhe: added this subroutine, modified from powerhel subroutine
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use Chiral, only: iXX_chiral, iYY_chiral
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, jkz, im, in, ivec, ivec_jj
  real :: k2, ktot2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: u_re, u_im, o_re, o_im
  real, dimension(nx,ny,nz) :: a_re, a_im, b_re, b_im
  real, dimension(nx) :: bbi, jji, b2, j2
  real, dimension(nx,3) :: bb, bbEP, jj
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk,nzgrid) :: vu_spec, vu_spec_sum
  real, dimension(nk,nzgrid) :: uo_spec, uo_spec_sum
  real, dimension(nk,nzgrid) :: ab_spec, ab_spec_sum
  real, dimension(nk,nzgrid) :: bj_spec, bj_spec_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true., lwrite_krms_GWs=.false.
!
!  passive scalar contributions (hardwired for now)
!
  real, dimension(nx,3) :: gtmp1,gtmp2
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  vu_spec=0.
  vu_spec_sum=0.
  uo_spec=0.
  uo_spec_sum=0.
  ab_spec=0.
  ab_spec_sum=0.
  bj_spec=0.
  bj_spec_sum=0.
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    ! calculate u and o
    !
    if (iuu==0) call fatal_error('powerhel','iuu=0')
    do n=n1,n2
      do m=m1,m2
        call curli(f,iuu,bbi,ivec)
        im=m-nghost
        in=n-nghost
        o_re(:,im,in)=bbi  !(this corresponds to vorticity)
      enddo
    enddo
    u_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
    o_im=0.
    u_im=0.
    !
    ! calculate a, b, and j
    !
    if (iaa==0) call fatal_error('powerhel','iaa=0')
    if (lmagnetic) then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
      a_im=0.
      b_im=0.
    else
      if (headt) print*,'magnetic power spectra only work if lmagnetic=T'
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(u_re,u_im)
    call fft_xyz_parallel(o_re,o_im)
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  calculate cylindrical helicity spectra
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over cylindrical shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
          ktot2=k2+kz(ikz+ipz*nz)**2
          jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1) .and. ktot2>0.) then
            uo_spec(k+1,jkz)=uo_spec(k+1,jkz) &
              +u_re(ikx,iky,ikz)*o_re(ikx,iky,ikz) &
              +u_im(ikx,iky,ikz)*o_im(ikx,iky,ikz)
            ab_spec(k+1,jkz)=ab_spec(k+1,jkz) &
              +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
              +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
            vu_spec(k+1,jkz)=uo_spec(k+1,jkz)/ktot2    ! fluid helicity is calculated by u.o/k^2, only work for low Mach numbers
            bj_spec(k+1,jkz)=ab_spec(k+1,jkz)*ktot2    ! current helicity is a.b*k^2
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(vu_spec,vu_spec_sum,(/nk,nzgrid/))
  call mpireduce_sum(uo_spec,uo_spec_sum,(/nk,nzgrid/))
  call mpireduce_sum(ab_spec,ab_spec_sum,(/nk,nzgrid/))
  call mpireduce_sum(bj_spec,bj_spec_sum,(/nk,nzgrid/))
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  no multiplication by 1/2
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing cylindrical power spectrum to'  &
        ,trim(datadir)//'/powerhel_fg.dat files'
    !
    open(1,file=trim(datadir)//'/powerhel_vu.dat',position='append')
    if (lformat) then
      do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, vu_spec_sum(k,jkz)
        enddo
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') vu_spec_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_uo.dat',position='append')
    if (lformat) then
      do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, uo_spec_sum(k,jkz)
        enddo
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') uo_spec_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_ab.dat',position='append')
    if (lformat) then
      do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, ab_spec_sum(k,jkz)
        enddo
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') ab_spec_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_bj.dat',position='append')
    if (lformat) then
      do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, bj_spec_sum(k,jkz)
        enddo
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') bj_spec_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine anisoq_diag
!***********************************************************************
  subroutine corfunc_cyl(f,sp)
!
!  Calculate azimuthally averaged correlation functions v11,v33,v12.
!  Here v11=(w11+w22)/2, v33=w33, v12=w12/sqrt(-1),
!  where wij=<u_i^* u_j>, for i,j=1,2,3.
!  They are functions kr=norm(kx,ky,kz) and kz/kr.
!  For the moment u=velocity field only,
!  but possible to extend using 'sp'
!
!  2020-Oct-14/hongzhe:  added this subroutine
!  2020-Nov-23/hongzhe:  removed 3d outputs
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use General, only: plegendre
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, ikx, iky, ikz, ikr, ikmu
  integer, dimension(1) :: temploc
  integer, dimension(nk) :: nmu
  real, dimension (mx,my,mz,mfarray) :: f
  real, allocatable, dimension(:,:) :: kmu
  real, dimension(nx,ny,nz) :: ux_re, ux_im
  real, dimension(nx,ny,nz) :: uy_re, uy_im
  real, dimension(nx,ny,nz) :: uz_re, uz_im
  real, allocatable, dimension(:,:) :: vxx, vxy, vzz
  real, allocatable, dimension(:,:) :: vxx_sum, vxy_sum, vzz_sum
  real, allocatable, dimension(:,:) :: coeff_a, coeff_b, coeff_c
  real, allocatable, dimension(:,:) :: coeff_a_sum, coeff_b_sum, coeff_c_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al_a, legendre_al_a_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al_b, legendre_al_b_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al_c, legendre_al_c_sum
  real :: k2, mu, kmu2
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true., lwrite_krms_GWs=.false.
!
!  passive scalar contributions (hardwired for now)
!
  real, dimension(nx,3) :: gtmp1,gtmp2
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize
  !
  k2m=0.
  nks=0.
  do ikr=1,nk
    nmu(ikr)=2*(ikr-1)+1
  enddo
  allocate( kmu(nk,nmu(nk)) )
  kmu=0.
  do ikr=2,nk
    do ikmu=1, nmu(ikr)
      kmu(ikr,ikmu)=-1+2.*(ikmu-1)/(nmu(ikr)-1)
    enddo
  enddo
  allocate( vxx(nk,nmu(nk)) )
  allocate( vxx_sum(nk,nmu(nk)) )
  allocate( vzz(nk,nmu(nk)) )
  allocate( vzz_sum(nk,nmu(nk)) )
  allocate( vxy(nk,nmu(nk)) )
  allocate( vxy_sum(nk,nmu(nk)) )
  allocate( coeff_a(nk,nmu(nk)) )
  allocate( coeff_a_sum(nk,nmu(nk)) )
  allocate( coeff_b(nk,nmu(nk)) )
  allocate( coeff_b_sum(nk,nmu(nk)) )
  allocate( coeff_c(nk,nmu(nk)) )
  allocate( coeff_c_sum(nk,nmu(nk)) )
  vxx=0.
  vxx_sum=0.
  vxy=0.
  vxy_sum=0.
  vzz=0.
  vzz_sum=0.
  coeff_a=0.
  coeff_a_sum=0.
  coeff_b=0.
  coeff_b_sum=0.
  coeff_c=0.
  coeff_c_sum=0.
  legendre_al_a=0.
  legendre_al_a_sum=0.
  legendre_al_b=0.
  legendre_al_b_sum=0.
  legendre_al_c=0.
  legendre_al_c_sum=0.
  !
  !  calculate each components
  !
  if (sp=='kin') then
    if (iuu==0) call fatal_error('powerhel','iuu=0')
    ux_re=f(l1:l2,m1:m2,n1:n2,iuu+1-1)
    uy_re=f(l1:l2,m1:m2,n1:n2,iuu+2-1)
    uz_re=f(l1:l2,m1:m2,n1:n2,iuu+3-1)
    ux_im=0.
    uy_im=0.
    uz_im=0.
  endif
  !
  !  Doing the Fourier transform
  !
  call fft_xyz_parallel(ux_re,ux_im)
  call fft_xyz_parallel(uy_re,uy_im)
  call fft_xyz_parallel(uz_re,uz_im)
  !
  !  calculate correlation functions
  !
  if (lroot .AND. ip<10) print*,'fft done; now integrate over cylindrical shells...'
  do ikz=1,nz
    do iky=1,ny
      do ikx=1, nx
        k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
        ikr=nint(sqrt(k2))
        mu=kz(ikz+ipz*nz)/sqrt(k2)
        if (ikr>=0. .and. ikr<=(nk-1)) then
          temploc=minloc(abs(kmu(ikr+1,:)-mu))
          ikmu=temploc(1)
          vxx(ikr+1,ikmu)=vxx(ikr+1,ikmu)+0.5*( &
              ux_re(ikx,iky,ikz)**2+ux_im(ikx,iky,ikz)**2 &
              +uy_re(ikx,iky,ikz)**2+uy_im(ikx,iky,ikz)**2 )
          vzz(ikr+1,ikmu)=vzz(ikr+1,ikmu)+ &
              uz_re(ikx,iky,ikz)**2+uz_im(ikx,iky,ikz)**2
          vxy(ikr+1,ikmu)=vxy(ikr+1,ikmu)+ &
              ux_re(ikx,iky,ikz)*uy_im(ikx,iky,ikz) &
              -ux_im(ikx,iky,ikz)*uy_re(ikx,iky,ikz)
        endif
      enddo
    enddo
  enddo
  !
  ! compute legendre coefficients
  !
  do ikr=1,nk
    do ikmu=1,nmu(ikr)
      kmu2=kmu(ikr,ikmu)**2
      coeff_c(ikr,ikmu)=vxy(ikr,ikmu)/(2*pi*kmu(ikr,ikmu))
      if (kmu2==1.) then
        coeff_a(ikr,ikmu)=vxx(ikr,ikmu)/(pi*2.)
      else
        coeff_a(ikr,ikmu)=( 4.*(1-kmu2)*vxx(ikr,ikmu)-kmu2*vzz(ikr,ikmu) )/ &
            ( 2*pi*(1-kmu2)*(2+kmu2) )
        coeff_b(ikr,ikmu)=( -2.*(1-kmu2)*vxx(ikr,ikmu)+(1+kmu2)*vzz(ikr,ikmu) )/ &
            ( pi*(1-kmu2)**2*(2+kmu2) )
      endif
      do i=1,legendre_lmax+1
        legendre_al_a(i,ikr)=legendre_al_a(i,ikr)+ &
            1./nmu(ikr)*(2*i-1)/2*coeff_a(ikr,ikmu)* &
            sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
        legendre_al_b(i,ikr)=legendre_al_b(i,ikr)+ &
            1./nmu(ikr)*(2*i-1)/2*coeff_a(ikr,ikmu)* &
            sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
        legendre_al_c(i,ikr)=legendre_al_c(i,ikr)+ &
            1./nmu(ikr)*(2*i-1)/2*coeff_a(ikr,ikmu)* &
            sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
      enddo
    enddo
  enddo
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(vxx,vxx_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(vxy,vxy_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(vzz,vzz_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(coeff_a,coeff_a_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(coeff_b,coeff_b_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(coeff_c,coeff_c_sum,(/nk,nmu(nk)/))
  call mpireduce_sum(legendre_al_a,legendre_al_a_sum,(/legendre_lmax+1,nk/))
  call mpireduce_sum(legendre_al_b,legendre_al_b_sum,(/legendre_lmax+1,nk/))
  call mpireduce_sum(legendre_al_c,legendre_al_c_sum,(/legendre_lmax+1,nk/))
  !
  !  on root processor, write global result to file
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing two point correlations to'  &
        ,trim(datadir)//'/cor2_.dat files'
    open(1,file=trim(datadir)//'/cor2_la_a_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    do i=1,legendre_lmax+1; do ikr=1,nk
      write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_a_sum(i,ikr)
    enddo; enddo
    close(1)
    open(1,file=trim(datadir)//'/cor2_la_b_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    do i=1,legendre_lmax+1; do ikr=1,nk
      write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_b_sum(i,ikr)
    enddo; enddo
    close(1)
    open(1,file=trim(datadir)//'/cor2_la_c_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    do i=1,legendre_lmax+1; do ikr=1,nk
      write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_c_sum(i,ikr)
    enddo; enddo
    close(1)
  endif
    
  endsubroutine corfunc_cyl
!***********************************************************************
  subroutine k_omega_spectra(f,sp)
!
!  Calculate spectra with both (x,y,z) and t Fourier-transformed
!  For the moment only works with kinetic helicity
!  Use luut_as_aux=T and specify omega_fourier
!  in run.in under &hydro_run_pars
!
!  29-oct-20/hongzhe: added this subroutine
!  20-nov-20/hongzhe: can now also compute Legendre coefficients
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use General, only: plegendre

!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, jkz, im, in, ivec, ivec_jj
  integer :: ikr, ikmu
  integer, dimension(nk) :: nmu
  real, allocatable, dimension(:,:) :: kmu
  real :: k2,mu
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi, jji, b2, j2
  real, dimension(nx,3) :: bb, bbEP, jj
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk,nzgrid) :: cyl_spectrum, cyl_spectrum_sum
  real, dimension(nk,nzgrid) :: cyl_spectrumhel, cyl_spectrumhel_sum
  real, allocatable, dimension(:,:) :: cyl_polar_spec, cyl_polar_spec_sum
  real, allocatable, dimension(:,:) :: cyl_polar_spechel, cyl_polar_spechel_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al, legendre_al_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_alhel, legendre_alhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  integer, dimension(1) :: temp
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true., lwrite_krms_GWs=.false.
!
!  passive scalar contributions (hardwired for now)
!
  real, dimension(nx,3) :: gtmp1,gtmp2
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  !
  if (lcylindrical_spectra) then
    cyl_spectrum=0.
    cyl_spectrum_sum=0.
    cyl_spectrumhel=0.
    cyl_spectrumhel_sum=0.
!
! allow for polar representation
!
    if (lcyl_polar_spectra) then
      do ikr=1,nk
        nmu(ikr)=2*(ikr-1)+1
      enddo
      !
      allocate( kmu(nk,nmu(nk)) )
      allocate( cyl_polar_spec(nk,nmu(nk)) )
      allocate( cyl_polar_spec_sum(nk,nmu(nk)) )
      allocate( cyl_polar_spechel(nk,nmu(nk)) )
      allocate( cyl_polar_spechel_sum(nk,nmu(nk)) )
      !
      kmu=0.
      do ikr=2,nk
        do ikmu=1, nmu(ikr)
          kmu(ikr,ikmu) = -1+2.*(ikmu-1)/(nmu(ikr)-1)
        enddo
      enddo
      ! 
      cyl_polar_spec=0.
      cyl_polar_spec_sum=0.
      cyl_polar_spechel=0.
      cyl_polar_spechel_sum=0.
      legendre_al=0.
      legendre_al_sum=0.
      legendre_alhel=0.
      legendre_alhel_sum=0.
    endif
  endif
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      if (iuu==0) call fatal_error('k_omega_spectra','iuu=0')
      if (iuut==0) call fatal_error('k_omega_spectra','iuut=0')
      if (iuust==0) call fatal_error('k_omega_spectra','iuust=0')
      if (ioot==0) call fatal_error('k_omega_spectra','ioot=0')
      if (ioost==0) call fatal_error('k_omega_spectra','ioost=0')
      b_re=f(l1:l2,m1:m2,n1:n2,iuut+ivec-1)    ! the real part of u(\vec x,\omega)
      b_im=f(l1:l2,m1:m2,n1:n2,iuust+ivec-1)   ! the imaginary part of u(\vec x,\omega)
      a_re=f(l1:l2,m1:m2,n1:n2,ioot+ivec-1)    ! the real part of omega(\vec x,\omega)
      a_im=f(l1:l2,m1:m2,n1:n2,ioost+ivec-1)   ! the imaginary part of omega(\vec x,\omega)
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  Compute cylindrical spectrum
!
    if (lcylindrical_spectra) then
      if (lroot .AND. ip<10) print*,'fft done; now integrate over cylindrical shells...'
      do ikz=1,nz
      do iky=1,ny
      do ikx=1,nx
        k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
        jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
        k=nint(sqrt(k2))
        if (k>=0 .and. k<=(nk-1)) then
          cyl_spectrum(k+1,jkz)=cyl_spectrum(k+1,jkz) &
              +b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2
          cyl_spectrumhel(k+1,jkz)=cyl_spectrumhel(k+1,jkz) &
              +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
              +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
        endif
!  compute azimuthally averaged spectrum
!  but as a function of kr=norm(kx,ky,kz) and kz/kr
        if (lcyl_polar_spectra) then
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          ikr=nint(sqrt(k2))
          mu=kz(ikz+ipz*nz)/sqrt(k2)
          if (ikr>=0. .and. ikr<=(nk-1)) then
            temp=minloc(abs(kmu(ikr+1,:)-mu))
            ikmu=temp(1)
            cyl_polar_spec(ikr+1,ikmu)=cyl_polar_spec(ikr+1,ikmu) &
                +b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2
            cyl_polar_spechel(ikr+1,ikmu)=cyl_polar_spechel(ikr+1,ikmu) &
                +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
          endif
        endif
!
!  end of loop through all points
!
      enddo
      enddo
      enddo
!
! compute legendre coefficients
!
      if (lcyl_polar_spectra) then
        do ikr=1,nk; do ikmu=1,nmu(ikr)
          do i=1,legendre_lmax+1
            legendre_al(i,ikr)=legendre_al(i,ikr)+ &
                1./nmu(ikr)*(2*i-1)/2*cyl_polar_spec(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
            legendre_alhel(i,ikr)=legendre_alhel(i,ikr)+ &
                1./nmu(ikr)*(2*i-1)/2*cyl_polar_spechel(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
          enddo
        enddo; enddo
      endif
!
    endif
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  if (lcylindrical_spectra) then
    call mpireduce_sum(cyl_spectrum,cyl_spectrum_sum,(/nk,nzgrid/))
    call mpireduce_sum(cyl_spectrumhel,cyl_spectrumhel_sum,(/nk,nzgrid/))
    if (lcyl_polar_spectra) then
      call mpireduce_sum(cyl_polar_spec,cyl_polar_spec_sum,(/nk,nmu(nk)/))
      call mpireduce_sum(cyl_polar_spechel,cyl_polar_spechel_sum,(/nk,nmu(nk)/))
      call mpireduce_sum(legendre_al,legendre_al_sum,(/legendre_lmax+1,nk/))
      call mpireduce_sum(legendre_alhel,legendre_alhel_sum,(/legendre_lmax+1,nk/))
    endif
  endif
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (lcyl_polar_spectra) then
      !  outputs for debug
      !  legendre p
      open(1,file=trim(datadir)//'/plegendre.dat',position='append')
      write(1,*) t
      do i=1,legendre_lmax+1; do ikr=1,nk; do ikmu=1,nmu(ikr)
      write(1,'(2i4,1p,8e10.2,3p,8e10.2)') i-1,ikr-1,kmu(ikr,ikmu),plegendre(i-1,0,kmu(ikr,ikmu))
      enddo;enddo;enddo
      close(1)
      !  spectra for comparison
      open(1,file=trim(datadir)//'/cyl_omega_power_cyl_'//trim(sp)//'.dat',position='append')
      write(1,*) t
      do jkz=1, nzgrid; do k=1, nk
        write(1,'(2i4,3p,8e10.2)') k-1, jkz-1-nzgrid/2, cyl_spectrum_sum(k,jkz)
      enddo; enddo
      close(1)
      open(1,file=trim(datadir)//'/cyl_omega_power_polar_'//trim(sp)//'.dat',position='append')
      write(1,*) t
      do ikr=1,nk; do ikmu=1,nmu(ikr)
        write(1,'(i4,1p,8e10.2,3p,8e10.2)') ikr-1, kmu(ikr,ikmu), cyl_polar_spec_sum(ikr,ikmu)
      enddo; enddo
      close(1)
      !  legendre coefficients a_l, in the form (l,kr,a_l), l,kr=0,1,2,..., 
      open(1,file=trim(datadir)//'/legendre_al_'//trim(sp)//'.dat',position='append')
      write(1,*) t
      do i=1,legendre_lmax+1; do ikr=1,nk
      write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_sum(i,ikr)
      enddo;enddo
      close(1)
      open(1,file=trim(datadir)//'/legendre_alhel_'//trim(sp)//'.dat',position='append')
      write(1,*) t
      do i=1,legendre_lmax+1; do ikr=1,nk
      write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_alhel_sum(i,ikr)
      enddo;enddo
      close(1)
    endif
    if (lcylindrical_spectra) then
      if (ip<10) print*,'Writing cylindrical power spectrum ',sp &
           ,' to ',trim(datadir)//'/cyl_omega_power_'//trim(sp)//'.dat'
      open(1,file=trim(datadir)//'/cyl_omega_power_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrum_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) t
        write(1,'(1p,8e10.2)') cyl_spectrum_sum
      endif
      close(1)
      !
      open(1,file=trim(datadir)//'/cyl_omega_powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrumhel_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) t
        write(1,'(1p,8e10.2)') cyl_spectrumhel_sum
      endif
      close(1)
    endif
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine k_omega_spectra
!***********************************************************************
  subroutine power1d_plane(f,sp)
!
!  Calculate power and helicity spectra of planar-averaged 
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   9-nov-20/hongzhe: if this coincides with power_1d I will remove it
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use Chiral, only: iXX_chiral, iYY_chiral
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, ikx, iky, ikz, jkz, im, in, ivec
  integer :: k3, k
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
!
!  passive scalar contributions (hardwired for now)
!
  real, dimension(nx,3) :: gtmp1,gtmp2
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      if (iuu==0) call fatal_error('powerhel','iuu=0')
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corresponds to vorticity)
        enddo
      enddo
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  magnetic power spectra (spectra of |B|^2 and A.B)
!
    elseif (sp=='mag') then
      if (iaa==0) call fatal_error('powerhel','iaa=0')
      if (lmagnetic) then
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,bbi,ivec)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
          enddo
        enddo
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
        a_im=0.
        b_im=0.
      else
        if (headt) print*,'magnetic power spectra only work if lmagnetic=T'
      endif
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over the xy plane...'
    do ikz=1,nz
      k3=nint(kz(ikz+ipz*nz))
      if (k3>=0 .and. k3<=nk-1) then
       do iky=1,ny
          do ikx=1,nx
            spectrum(k3+1)=spectrum(k3+1) &
             +2*b_re(ikx,iky,ikz)**2 &
              +2*b_im(ikx,iky,ikz)**2
            spectrumhel(k3+1)=spectrumhel(k3+1) &
              +2*a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
              +2*a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
          enddo
        enddo
      endif
    enddo
    spectrum(1)=spectrum(1)/2
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/powerkz_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/powerkz_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhelkz_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
  endif
  !
  endsubroutine power1d_plane
endmodule power_spectrum
