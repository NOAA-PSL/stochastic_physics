!>@brief The module 'stochy_patterngenerator_mod' contains the derived type random_pattern
!! which controls the characteristics of the random pattern
!! This is a modified version of the original one stochy_patterngenerator.F90, where the
!! the random patterns are not properly normalized
module stochy_patterngenerator_mod

!> generate random patterns with specified temporal and spatial auto-correlation
!! in spherical harmonic space.
 use kinddef
 use spectral_transforms, only: len_trie_ls, len_trio_ls, ls_dim, ls_max_node
! use mersenne_twister_stochy, only: random_setseed,random_gauss,random_stat
 use mersenne_twister, only: random_setseed,random_gauss,random_stat
 ! DH* replacing this with mpi_wrapper changes results - look for value of iseed?
 use mpi_wrapper,only: is_rootpe, mp_bcst
 ! *DH
 implicit none
 private

 public :: computevarspec, setvarspect, &
  patterngenerator_init, patterngenerator_destroy, getnoise, getnoise_un, &
  patterngenerator_advance, patterngenerator_advance_jb, ndimspec, &
  chgres_pattern,computevarspec_r

! -----------------------------------------------
!>@brief Derived type 'random_pattern' contains the attrubites of a spefic random pattern
!>@details A seperate instance of this type is needed for each pattern
 type,public :: random_pattern                ! start type define
! -----------------------------------------------
    real(kind_dbl_prec), public :: lengthscale  ! length scale in m
    real(kind_dbl_prec), public :: tau
    real(kind_dbl_prec), public :: dt
    real(kind_dbl_prec), public :: phi
    real(kind_dbl_prec), public :: stdev
    real(kind_dbl_prec), allocatable, dimension(:), public :: varspectrum, varspectrum1d, lap
    integer, allocatable, dimension(:), public ::&
       degree,order,idx_e,idx_o
    integer, allocatable, dimension(:,:), public :: idx
    integer, public :: seed
    real(kind_dbl_prec), allocatable, dimension(:,:,:), public :: spec_e,spec_o
    type(random_stat), public :: rstate
! -----------------------------------------------------
 end type random_pattern             ! end type define
! -----------------------------------------------------

 integer :: nlons,nlats,ntrunc,ndimspec

 contains
!>@brief The subroutine 'patterngenerator_init' sets up the spherical harmonics
!
!>@details It populates array defining the zonal and total wavenumbers, amplitude,
!! temporaral and spatial correlations.
 subroutine patterngenerator_init(lscale, delt, tscale, stdev, iseed, rpattern,&
                                  nlon, nlat, jcap, ls_nodes, npatterns,&
                                  nlevs, varspect_opt,new_lscale,bn)
!\callgraph
   real(kind_dbl_prec), intent(in),dimension(npatterns) :: lscale,tscale,stdev
   real(kind=kind_dbl_prec), intent(in) :: delt
   integer, intent(in) :: nlon,nlat,jcap,npatterns,varspect_opt
   integer, intent(in) :: ls_nodes(ls_dim,3),nlevs
   logical, intent(in) :: new_lscale
   logical, optional, intent(in) :: bn ! berner normalization
   type(random_pattern), intent(out), dimension(npatterns) :: rpattern
   integer(8), intent(inout) :: iseed(npatterns)
   integer m,j,l,n,nm,nn,np,indev1,indev2,indod1,indod2
   integer(8) count, count_rate, count_max, count_trunc
   integer(8) :: iscale = 10000000000_8
   integer count4, ierr
   logical :: bn_local ! berner normalization
!   integer  member_id
   integer indlsod,indlsev,jbasev,jbasod
   include 'function_indlsod'
   include 'function_indlsev'
   bn_local=.false.
   if (present(bn)) bn_local=bn
   if (is_rootpe()) then
      if (present(bn)) then
         print*,'Berner norm=',bn
      else
         print*,'old norm'
      endif
   end if
   nlons = nlon
   nlats = nlat
   ntrunc = jcap
   ndimspec = (ntrunc+1)*(ntrunc+2)/2
!  propagate seed supplied from namelist to all patterns...
   if (iseed(1) .NE. 0) then
      do np=2,npatterns
         if (iseed(np).EQ.0) then
            iseed(np)=iseed(1)+np*100000000
         endif
      enddo
   endif

   do np=1,npatterns
      allocate(rpattern(np)%idx(0:ntrunc,0:ntrunc))
      allocate(rpattern(np)%idx_e(len_trie_ls))
      allocate(rpattern(np)%idx_o(len_trio_ls))
      allocate(rpattern(np)%spec_e(len_trie_ls,2,nlevs))
      allocate(rpattern(np)%spec_o(len_trio_ls,2,nlevs))
      rpattern(np)%idx_e = 0; rpattern(np)%idx_o = 0; rpattern(np)%idx = 0
      rpattern(np)%spec_e(:,:,:)=0.
      rpattern(np)%spec_o(:,:,:)=0.
      nm = 0
      do m=0,ntrunc
         do n=m,ntrunc
            nm = nm + 1
            rpattern(np)%idx(m,n) = nm
         enddo
      enddo
      do j = 1, ls_max_node
         l=ls_nodes(j,1) ! zonal wavenumber
         jbasev=ls_nodes(j,2)
         jbasod=ls_nodes(j,3)
         indev1 = indlsev(l,l)
         indod1 = indlsod(l+1,l)
         if (mod(l,2) .eq. mod(ntrunc+1,2)) then
            indev2 = indlsev(ntrunc+1,l)
            indod2 = indlsod(ntrunc  ,l)
         else
            indev2 = indlsev(ntrunc  ,l)
            indod2 = indlsod(ntrunc+1,l)
         endif
         n = l ! degree
         do nn=indev1,indev2
            if (n <= ntrunc .and. l <= ntrunc) then
              nm = rpattern(np)%idx(l,n)
              rpattern(np)%idx_e(nn) = nm
            endif
            n = n + 2
         enddo
         n = l+1
         do nn=indod1,indod2
            if (n <= ntrunc .and. l <= ntrunc) then
              nm = rpattern(np)%idx(l,n)
              rpattern(np)%idx_o(nn) = nm
            endif
            n = n + 2
         enddo
      enddo
      allocate(rpattern(np)%degree(ndimspec),rpattern(np)%order(ndimspec),rpattern(np)%lap(ndimspec))
      j = 0
      do m=0,ntrunc
        do n=m,ntrunc
          j = j + 1
          rpattern(np)%degree(j) = n
          rpattern(np)%order(j) = m
        end do
      end do
      rpattern(np)%lap = -rpattern(np)%degree*(rpattern(np)%degree+1.0)
      rpattern(np)%tau = tscale(np)
      rpattern(np)%lengthscale = lscale(np)
      rpattern(np)%dt = delt
      rpattern(np)%phi = exp(-delt/tscale(np))
      rpattern(np)%stdev = stdev(np)
      allocate(rpattern(np)%varspectrum(ndimspec))
      allocate(rpattern(np)%varspectrum1d(0:ntrunc))
      ! seed computed on root, then bcast to all tasks and set.
      if (is_rootpe()) then
!         read(ens_nam(2:3),'(i2)') member_id
!         print *,'ens_nam,member_id',trim(ens_nam),member_id
         if (iseed(np) == 0) then
           ! generate a random seed from system clock and ens member number
           call system_clock(count, count_rate, count_max)
           ! iseed is elapsed time since unix epoch began (secs)
           ! truncate to 4 byte integer
           count_trunc = iscale*(count/iscale)
           count4 = count - count_trunc !+ member_id
           print *,'using seed',count4
         else
           !count4 = iseed(np) + member_id
           ! don't rely on compiler to truncate integer(8) to integer(4) on
           ! overflow, do wrap around explicitly.
           !count4 = mod(iseed(np) + member_id + 2147483648_8, 4294967296_8) - 2147483648_8
           count4 = int(mod(int(iseed(np), 8) + 2147483648_8, 4294967296_8) - 2147483648_8)
           print *,'using seed',count4,iseed(np)!,member_id
         endif
      endif
      ! broadcast seed to all tasks.
      call mp_bcst(count4)
      rpattern(np)%seed = count4
      ! set seed (to be the same) on all tasks. Save random state.
      call random_setseed(rpattern(np)%seed,rpattern(np)%rstate)
      if (varspect_opt .lt. 0 .or. varspect_opt .gt. 1) then
         if (is_rootpe()) then
            print *,'WARNING: illegal value for varspect_opt (should be 0,1,or 2), using 0 (gaussian spectrum)...'
         endif
         call setvarspect(rpattern(np),0,new_lscale,bn_local)
      else
         call setvarspect(rpattern(np),varspect_opt,new_lscale,bn_local)
      endif
   enddo ! n=1,npatterns
 end subroutine patterngenerator_init


!>@brief The subroutine 'patterngenerator_destroy' dellaocate arrays
!>@details This is acutally never called
 subroutine patterngenerator_destroy(rpattern,npatterns)
!\callgraph
   type(random_pattern), intent(inout) :: rpattern(npatterns)
   integer, intent(in) :: npatterns
   integer n
   do n=1,npatterns
   deallocate(rpattern(n)%varspectrum,rpattern(n)%varspectrum1d)
   deallocate(rpattern(n)%degree,rpattern(n)%order,rpattern(n)%lap)
   deallocate(rpattern(n)%idx,rpattern(n)%idx_e,rpattern(n)%idx_o)
   enddo
 end subroutine patterngenerator_destroy

!>@brief The subroutine 'computevarspec' compute the globally integrated
!! variance from complex spectral coefficients
!>@details this is necessary to ensure the proper global variance
 subroutine computevarspec(rpattern,dataspec,var)
!\callgraph
    ! compute globally integrated variance from spectral coefficients
    complex(kind_dbl_prec), intent(in) :: dataspec(ndimspec)
    real(kind_dbl_prec), intent(out) ::  var
    type(random_pattern), intent(in) :: rpattern
    integer n
    var = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
           var = var + dataspec(n)*conjg(dataspec(n))
       else
           var = var + 0.5*dataspec(n)*conjg(dataspec(n))
       endif
    enddo
 end subroutine computevarspec

!>@brief The subroutine 'computevarspec_r' compute the globally integrated
!! variance from real spectral coefficients
!>@details this is necessary to ensure the proper global variance
 subroutine computevarspec_r(rpattern,dataspec,var)
!\callgraph
    ! compute globally integrated variance from spectral coefficients
    real(kind_dbl_prec), intent(in) :: dataspec(2*ndimspec)
    real(kind_dbl_prec), intent(out) ::  var
    type(random_pattern), intent(in) :: rpattern
    integer n
    var = 0.
    do n=1,ndimspec
       if (rpattern%order(n) .ne. 0) then
           var = var + dataspec(n)**2+dataspec(n+ndimspec)**2
       else
           var = var + 0.5*(dataspec(n)**2+dataspec(n+ndimspec)**2)
       endif
    enddo
 end subroutine computevarspec_r

!>@brief The subroutine 'getnoise' scales spectral cofficients with
!! white noise to the appropriate amplitude for speherical harmonincs
!! variance from real spectral c
 subroutine getnoise(rpattern,noise_e,noise_o)
!\callgraph
   real(kind_dbl_prec), intent(out) :: noise_e(len_trie_ls,2)
   real(kind_dbl_prec), intent(out) :: noise_o(len_trio_ls,2)
   ! generate white noise with unit variance in spectral space
   type(random_pattern), intent(inout) :: rpattern
   real(kind_dbl_prec) :: noise(2*ndimspec)
   integer nm,nn
   call random_gauss(noise,rpattern%rstate)
   noise(1) = 0.; noise(ndimspec+1) = 0.
   noise = noise*sqrt(1./ntrunc)
   noise_e = 0.; noise_o = 0.
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      noise_e(nn,1) = noise(nm)/sqrt(2.*rpattern%degree(nm)+1)
      noise_e(nn,2) = noise(ndimspec+nm)/sqrt(2.*rpattern%degree(nm)+1)
      if (rpattern%order(nm) .eq. 0) then
        noise_e(nn,1) = sqrt(2.)*noise_e(nn,1)
        noise_e(nn,2) = 0.
      endif
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      noise_o(nn,1) = noise(nm)/sqrt(2.*rpattern%degree(nm)+1)
      noise_o(nn,2) = noise(ndimspec+nm)/sqrt(2.*rpattern%degree(nm)+1)
      if (rpattern%order(nm) .eq. 0) then
        noise_o(nn,1) = sqrt(2.)*noise_o(nn,1)
        noise_o(nn,2) = 0.
      endif
   enddo
 end subroutine getnoise

!>@brief The subroutine 'getnoise' scales spectral cofficients with
!! white noise to the appropriate amplitude for speherical harmonincs
!! variance from real spectral c
 subroutine getnoise_un(rpattern,noise_e,noise_o)
!\callgraph
   real(kind_dbl_prec), intent(out) :: noise_e(len_trie_ls,2)
   real(kind_dbl_prec), intent(out) :: noise_o(len_trio_ls,2)
   ! generate white noise with unit variance in spectral space
   type(random_pattern), intent(inout) :: rpattern
   real(kind_dbl_prec) :: noise(2*ndimspec)
   integer nm,nn
   call random_gauss(noise,rpattern%rstate)
   noise(1) = 0.; noise(ndimspec+1) = 0.
   noise_e = 0.; noise_o = 0.
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      noise_e(nn,1) = noise(nm)
      noise_e(nn,2) = noise(ndimspec+nm)
      if (rpattern%order(nm) .eq. 0) then
        noise_e(nn,2) = 0.
      endif
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      noise_o(nn,1) = noise(nm)
      noise_o(nn,2) = noise(ndimspec+nm)
      if (rpattern%order(nm) .eq. 0) then
        noise_o(nn,2) = 0.
      endif
   enddo
 end subroutine getnoise_un

!>@brief The subroutine 'patterngenerator_advance' advance 1st-order autoregressive process
 subroutine patterngenerator_advance(rpattern,k,skeb_first_call)
!\callgraph

    ! advance 1st-order autoregressive process with
    ! specified autocorrelation (phi) and variance spectrum (spectrum)
    real(kind_dbl_prec) :: noise_e(len_trie_ls,2)
    real(kind_dbl_prec) :: noise_o(len_trio_ls,2)
    type(random_pattern), intent(inout) :: rpattern
    logical, intent(in) :: skeb_first_call
    integer j,l,n,nn,nm,k,k2
    call getnoise(rpattern,noise_e,noise_o)
    if (k.GT.1.AND.skeb_first_call) then
       k2=k-1
    else
       k2=k
    endif
    do nn=1,len_trie_ls
       nm = rpattern%idx_e(nn)
       if (nm == 0) cycle
       rpattern%spec_e(nn,1,k) =  rpattern%phi*rpattern%spec_e(nn,1,k2) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_e(nn,1)
       rpattern%spec_e(nn,2,k) =  rpattern%phi*rpattern%spec_e(nn,2,k2) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_e(nn,2)
    enddo
    do nn=1,len_trio_ls
       nm = rpattern%idx_o(nn)
       if (nm == 0) cycle
       rpattern%spec_o(nn,1,k) =  rpattern%phi*rpattern%spec_o(nn,1,k2) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_o(nn,1)
       rpattern%spec_o(nn,2,k) =  rpattern%phi*rpattern%spec_o(nn,2,k2) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*rpattern%varspectrum(nm)*noise_o(nn,2)
    enddo
 end subroutine patterngenerator_advance

!>@brief The subroutine 'patterngenerator_advance_jb' advance 1st-order autoregressive process
 subroutine patterngenerator_advance_jb(rpattern)
!\callgraph

    ! advance 1st-order autoregressive process with
    ! specified autocorrelation (phi) and variance spectrum (spectrum)
    real(kind_dbl_prec) :: noise_e(len_trie_ls,2)
    real(kind_dbl_prec) :: noise_o(len_trio_ls,2)
    type(random_pattern), intent(inout) :: rpattern
    integer j,l,n,nn,nm
    call getnoise_un(rpattern,noise_e,noise_o)
    !call getnoise(rpattern,noise_e,noise_o)
    !print*,'in advance_jb'
    !print*, 'random pattern stdev:', rpattern%stdev
    do nn=1,len_trie_ls
       nm = rpattern%idx_e(nn)
       if (nm == 0) cycle
       rpattern%spec_e(nn,1,1) = rpattern%phi*rpattern%spec_e(nn,1,1) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*noise_e(nn,1)*rpattern%varspectrum(nm)
       rpattern%spec_e(nn,2,1) = rpattern%phi*rpattern%spec_e(nn,2,1) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*noise_e(nn,2)*rpattern%varspectrum(nm)
    enddo
    do nn=1,len_trio_ls
       nm = rpattern%idx_o(nn)
       if (nm == 0) cycle
       rpattern%spec_o(nn,1,1) = rpattern%phi*rpattern%spec_o(nn,1,1) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*noise_o(nn,1)*rpattern%varspectrum(nm)
       rpattern%spec_o(nn,2,1) = rpattern%phi*rpattern%spec_o(nn,2,1) + &
       rpattern%stdev*sqrt(1.-rpattern%phi**2)*noise_o(nn,2)*rpattern%varspectrum(nm)
    enddo
 end subroutine patterngenerator_advance_jb

!>@brief The subroutine 'setvarspect' calculates the variance spectrum
! from a specified decorrelation lengthscale
 subroutine setvarspect(rpattern,varspect_opt,new_lscale,berner_normalize)
!\callgraph
 ! define variance spectrum (isotropic covariance)
 ! normalized to unit global variance
  type(random_pattern), intent(inout) :: rpattern
  integer, intent(in) :: varspect_opt
  logical, intent(in) :: new_lscale
  logical, intent(in) :: berner_normalize
  integer :: n, nm
  complex(kind_dbl_prec) noise(ndimspec)
  real(kind_dbl_prec) var,b_jb,rerth,inv_rerth_sq,pi,gamma_sum,deltaE
  pi    = 4.0d0*atan(1.0d0)
  rerth  =6.3712e+6      ! radius of earth (m)
  inv_rerth_sq=1.0/(rerth**2)
  !print*,'setvarspect,ntrunc',berner_normalize,ntrunc
  ! 1d variance spectrum (as a function of total wavenumber)
  if (varspect_opt == 0) then ! gaussian
     if (is_rootpe()) print*, 'Gaussian variance spectrum'
     ! rpattern%lengthscale is interpreted as an efolding length
     ! scale, in meters.
     ! scaling factors for spectral coeffs of white noise pattern with unit variance
     if (new_lscale) then
        !fix for proper lengthscale
        if (is_rootpe()) print*, 'Proper lengthscale condition'
        rpattern%varspectrum = exp((rpattern%lengthscale*0.25)**2*rpattern%lap*inv_rerth_sq)
        do n=0,ntrunc
           rpattern%varspectrum1d(n) = exp(-(rpattern%lengthscale*0.25)**2*float(n)*(float(n)+1.)*inv_rerth_sq)
        enddo
     else
        if (is_rootpe()) print*, 'Not proper lengthscale condition'
        rpattern%varspectrum = exp(rpattern%lengthscale**2*rpattern%lap/(4.*rerth**2))
        do n=0,ntrunc
           rpattern%varspectrum1d(n) = exp(-rpattern%lengthscale**2*(float(n)*(float(n)+1.))/(4.*rerth**2))
        enddo
        if (is_rootpe()) print*,'Finished'
     endif
  else if (varspect_opt == 1) then ! power law
     ! rpattern%lengthscale is interpreted as a power, not a length.
     !print*,'Power Law lengthscale',rpattern%lengthscale
     rpattern%varspectrum1d(0) = 0.0
     do n=1,ntrunc
        rpattern%varspectrum1d(n) = float(n)**(rpattern%lengthscale)
     enddo
     ! scaling factors for spectral coeffs of white noise pattern with unit variance
     rpattern%varspectrum(1) = 0
     rpattern%varspectrum(2:) = rpattern%degree(2:)**(rpattern%lengthscale)
  endif
  !print*, 'Berner Normalization:',berner_normalize

  if (.not. berner_normalize) then ! normalize the GEFS way
     noise = 0.
     do n=1,ndimspec
        if (rpattern%order(n) .ne. 0.) then
           noise(n) = cmplx(1.,1.)/sqrt(2.*rpattern%degree(n)+1)
        else
           noise(n) = sqrt(2.)/sqrt(2.*rpattern%degree(n)+1.)
        endif
     enddo
     noise(1) = 0 ! no global mean.
     ! make sure global mean variance is 1.
     noise = noise*sqrt(1./ntrunc)
     noise = rpattern%varspectrum*noise
     call computevarspec(rpattern,noise,var)
     rpattern%varspectrum = rpattern%varspectrum/sqrt(var)
     rpattern%varspectrum1d = rpattern%varspectrum1d/var

  else if (berner_normalize) then ! normalize by Berner et al. 2009
     do n=1,len_trie_ls
         nm = rpattern%idx_e(n)
         if (nm == 0) cycle
          if (rpattern%order(nm) .ne. 0.) then
            rpattern%varspectrum(nm) = rpattern%varspectrum(nm)/sqrt(ntrunc*(2.*rpattern%degree(nm)+1)*(1+rpattern%phi))
          else
            rpattern%varspectrum(nm) = rpattern%varspectrum(nm)*sqrt(2.)/sqrt(ntrunc*(2.*rpattern%degree(n)+1.)*(1+rpattern%phi))
          endif
     enddo
     do n=1,len_trio_ls
        nm = rpattern%idx_o(n)
        if (nm == 0) cycle
         if (rpattern%order(nm) .ne. 0.) then
           rpattern%varspectrum(nm) = rpattern%varspectrum(nm)/sqrt(ntrunc*(2.*rpattern%degree(nm)+1)*(1+rpattern%phi))
         else
           rpattern%varspectrum(nm) = rpattern%varspectrum(nm)*sqrt(2.)/sqrt(ntrunc*(2.*rpattern%degree(n)+1.)*(1+rpattern%phi))
         endif
     enddo
     do n=0,ntrunc
        rpattern%varspectrum1d(n) = rpattern%varspectrum1d(n)/sqrt(ntrunc*(2.*n+1.))
     enddo
     gamma_sum=0
     !print*,'summing gamma'
     do n=1,ntrunc
        !print*,n,gamma_sum,rpattern%varspectrum1d(n)**2
        !gamma_sum = gamma_sum+n*(n+1.)*(2*n+1)*(rpattern%varspectrum1d(n)**2)
        gamma_sum = gamma_sum+n*(n+1.)*(2.*n+2.)*(rpattern%varspectrum1d(n)**2)
     enddo
     !b_jb=sqrt(((1-rpattern%phi)*4*pi*(rerth**2))/(gamma_sum*sqrt(var))) !!!Changed phi to (1-phi)
     b_jb = sqrt(4*pi*(rerth**2)*(1-rpattern%phi)/gamma_sum)
     !print*,'normalizing by b_jb',gamma_sum,b_jb
     rpattern%varspectrum = rpattern%varspectrum*b_jb
     rpattern%varspectrum1d = rpattern%varspectrum1d*b_jb
  endif
 end subroutine setvarspect

!>@brief The subroutine 'chgres_pattern' truncates the spherical harmonics if
!! restarting from a higher-resolution pattern
 subroutine chgres_pattern(pattern2din,pattern2dout,ntruncin,ntruncout)
!\callgraph
   real(kind_dbl_prec), intent(in) :: pattern2din((ntruncin+1)*(ntruncin+2))
   real(kind_dbl_prec), intent(out) :: pattern2dout((ntruncout+1)*(ntruncout+2))
   integer, intent(in) :: ntruncin,ntruncout
   integer             :: m,n,nm,ndimsspecin,ndimsspecout
   integer,allocatable, dimension(:,:):: idxin
   allocate(idxin(0:ntruncin,0:ntruncin))
   ndimsspecin=(ntruncin+1)*(ntruncin+2)/2
   ndimsspecout=(ntruncout+1)*(ntruncout+2)/2
   nm = 0
   do m=0,ntruncin
      do n=m,ntruncin
         nm = nm + 1
         idxin(m,n) = nm
      enddo
   enddo
  ! chgres
  nm = 0
  do m=0,ntruncout
     do n=m,ntruncout
        nm = nm + 1
        if (m .le.  ntruncin .and.  n .le.  ntruncin) then
           pattern2dout(nm) = pattern2din(idxin(m,n))
           pattern2dout(ndimsspecout+nm) = pattern2din(ndimsspecin+idxin(m,n))
        endif
     enddo
  enddo
  deallocate(idxin)
end subroutine chgres_pattern

end module stochy_patterngenerator_mod
