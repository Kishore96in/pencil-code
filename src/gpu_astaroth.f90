! $Id$
!
! MODULE_DOC: This module contains GPU related types and functions to be used with the ASTAROTH nucleus.
!
! CPARAM logical, parameter :: lgpu = .true.
!
!**************************************************************************
!
module GPU
!
  use Cdata
  use General, only: keep_compiler_quiet, ioptest, loptest
  use Mpicomm, only: stop_it
  use iso_c_binding

  implicit none

  external initialize_gpu_c
  external finalize_gpu_c
  external rhs_gpu_c
  external copy_farray_c
  external pos_real_ptr_c

  include 'gpu.h'

  !integer(KIND=ikind8) :: pFarr_GPU_in, pFarr_GPU_out
  type(C_PTR) :: pFarr_GPU_in, pFarr_GPU_out
  
contains

!***********************************************************************
    subroutine initialize_GPU
!
      character(LEN=512) :: str
!
      str=''
      if (lanelastic) str=trim(str)//', '//'anelastic'
      if (lboussinesq) str=trim(str)//', '//'boussinesq'
      !if (lenergy) str=trim(str)//', '//'energy'
      if (ltemperature) str=trim(str)//', '//'temperature'
      if (lshock) str=trim(str)//', '//'shock'
      if (lgrav) str=trim(str)//', '//'gravity'
      if (lheatflux) str=trim(str)//', '//'heatflux'
      if (lhyperresistivity_strict) str=trim(str)//', '//'hyperresi_strict'
      if (lhyperviscosity_strict) str=trim(str)//', '//'hypervisc_strict'
      if (lADI) str=trim(str)//', '//'implicit_physics'
      if (llorenz_gauge) str=trim(str)//', '//'lorenz_gauge'
      if (ldustvelocity) str=trim(str)//', '//'dustvelocity'
      if (ldustdensity) str=trim(str)//', '//'dustdensity'
      if (ltestscalar) str=trim(str)//', '//'testscalar'
      if (ltestfield) str=trim(str)//', '//'testfield'
      if (ltestflow) str=trim(str)//', '//'testflow'
      if (linterstellar) str=trim(str)//', '//'interstellar'
      if (lcosmicray) str=trim(str)//', '//'cosmicray'
      if (lcosmicrayflux) str=trim(str)//', '//'cosmicrayflux'
      if (lshear) str=trim(str)//', '//'shear'
      if (lpscalar) str=trim(str)//', '//'pscalar'
      if (lascalar) str=trim(str)//', '//'ascalar'
      if (lradiation) str=trim(str)//', '//'radiation'
      if (lchemistry) str=trim(str)//', '//'chemistry'
      if (lchiral) str=trim(str)//', '//'chiral'
      if (ldetonate) str=trim(str)//', '//'detonate'
      if (lneutralvelocity) str=trim(str)//', '//'neutralvelocity'
      if (lneutraldensity) str=trim(str)//', '//'neutraldensity'
      if (lopacity) str=trim(str)//', '//'opacity'
      if (lpolymer) str=trim(str)//', '//'polymer'
      if (lpointmasses) str=trim(str)//', '//'pointmasses'
      if (lpoisson) str=trim(str)//', '//'poisson'
      if (lselfgravity) str=trim(str)//', '//'selfgravity'
      if (lsolid_cells) str=trim(str)//', '//'solid_cells'
      if (lspecial) str=trim(str)//', '//'special'
      if (lpower_spectrum) str=trim(str)//', '//'power_spectrum'
      if (lparticles) str=trim(str)//', '//'particles'

      if (str/='') call stop_it('No GPU implementation for module(s) "'//trim(str(3:))//'"')
!
      call initialize_gpu_c(pFarr_GPU_in,pFarr_GPU_out)
!print'(a,1x,Z0,1x,Z0)', 'pFarr_GPU_in,pFarr_GPU_out=', pFarr_GPU_in,pFarr_GPU_out
    endsubroutine initialize_GPU
!**************************************************************************
    subroutine gpu_init
!
      call init_gpu_c
!
    endsubroutine gpu_init
!**************************************************************************
    subroutine register_GPU(f)
!
      real, dimension(:,:,:,:), intent(IN) :: f

      call register_gpu_c(f)
!
    endsubroutine register_GPU
!**************************************************************************
    subroutine finalize_GPU
!
      call finalize_gpu_c
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,isubstep,early_finalize)
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
      logical,                            intent(IN)    :: early_finalize
!
      integer :: ll, mm, nn
      real :: val
      logical, save :: lvery_first=.true.

      goto 1
      val=1.
      do nn=1,mz
        do mm=1,my
          do ll=1,mx
            f(ll,mm,nn,iux)=val; val=val+1.
      enddo; enddo; enddo

      print*, 'vor integrate'
      do nn=1,3
        if (notanumber(f(:,:,nn,iux))) print*,'NaN in ux, lower z', nn
      enddo
      print*, '---------------'

1     continue
      call rhs_gpu_c(isubstep,lvery_first,early_finalize)
!
      lvery_first=.false.

      return
!
      if (.not.lroot) return
      do nn=1,mz   !  nghost+1,mz-nghost
        print*, 'nn=', nn
        do mm=1,my
          print'(22(1x,f7.0))',f(:,mm,nn,iux)
      enddo; enddo

      do nn=1,3
        if (notanumber(f(:,:,nn,iux))) print*,'NaN in ux, lower z', nn                
      enddo

    endsubroutine rhs_GPU
!**************************************************************************
    function get_ptr_GPU(ind1,ind2,lout) result(pFarr)

      use Cparam
      use iso_c_binding

      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:), pointer :: pFarr

      integer :: i2

      interface
        type(c_ptr) function pos_real_ptr_c(ptr,ind)
          import :: c_ptr, ikind8
          type(c_ptr) :: ptr
          integer :: ind
        endfunction
      endinterface

      i2 = ioptest(ind2,ind1)
      if (loptest(lout)) then
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_out,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1/))
      else
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_in,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1/))
      endif

    endfunction get_ptr_GPU
!**************************************************************************
    subroutine copy_farray_from_GPU(f)

      real, dimension (mx,my,mz,mfarray), intent(OUT) :: f

      call copy_farray_c(f(1,1,1,iux),f(1,1,1,iuy),f(1,1,1,iuz),f(1,1,1,ilnrho))

    endsubroutine copy_farray_from_GPU
!**************************************************************************
endmodule GPU
