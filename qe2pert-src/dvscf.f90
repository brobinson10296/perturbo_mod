!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jpark <jinsoop412@gmail.com>; jjzhou <jjchou.comphy@gmail.com>; 
! Comment:
!  adapted from QE/PHonon/PH/dfile_star
!  Formula for rotation of dvscf: Eq. C13 in Luis paper, PRB 97, 235146(2018)
!
!  N.B. sr in symm_base module are symmetry operation matrices.
!  sr is S operation in cartesian cooridnate, sr(i,j) = S_ij, 
!   and the same sr can be applied to both real and recip. space.
!
!  If if we transfer S from cartesian (c) coord. to crystal coord. (r).
!   then the matrix are different for real space and recip. space:
!  in  real sace,  x_c = A x_r, (A -> at in cell_base module) 
!  in recip space. k_c = B k_r  (B -> bg in cell_base module), and A^-1 = B^T
!  so,  
!  S in crystal coord. in real space is:  B^T * sr * A  (T means transpose)
!    S (A x_s) = A x_s'  -->  (A^-1 S A) x_s = x_s'
!   
!  S in crystal coord. in recip. space is:  A^T * sr * B
!    S (B k_s) = B k_s'  -->  (B^-1 S B) k_s = k_s'
!
!  the variable s in symm_base is: 
!    s(i,j) = (B^T * sr * A)_ji  (S in cryst. coord. in real space)
!  quivalently, 
!    s^T = (B^T * sr * A)  ==> s = A^T * sr^T * B. = A^T (sr^-1) *B
!  since sr is rotation operation (in cart.), sr^T = sr^-1
!  so, s(i,j) = (A^T (sr^-1) * B)_ij (S^-1 in cryst, coord. in recip. space)
!  in summary, s has different meaning when applied to real. and recip. space.
!
! Maintenance:
!===============================================================================

subroutine load_dvscf(iq, u, dvscf)
   use kinds,    only : dp
   use input_param, only: prefix, phdir
   use fft_base, only: dfftp
   use ions_base, only: nat
   use noncollin_module, only: nspin_mag
   implicit none
   integer,     intent(in) :: iq
   complex(dp), intent(in) :: u(3*nat, 3*nat)
   complex(dp), intent(out):: dvscf(dfftp%nnr, nspin_mag, 3*nat)
    
   logical :: exst
   integer :: rec_len, iunit, imod, nmodes, nrow
   character(len=256) :: fname
   complex(dp), allocatable :: dvscf_tmp(:,:,:)

   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
   complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
   ! functions:
   integer, external :: find_free_unit
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   fname = trim(phdir) // trim(prefix) // '.dvscf_q' // trim(int_to_char(iq))
 !  when the file is read, it must be locked - invoke mpi barrier later
   inquire(file=trim(fname), exist = exst)
   if(.not. exst) call errore('load_dvscf',trim(fname)//' does not exist!',1)
   !
   nmodes  = 3 * nat
   rec_len = 2 * dfftp%nnr * nspin_mag
   !
   ! Open the dvscf file for reading
   iunit = find_free_unit()
   call open_dvscf(iunit, fname, rec_len, nmodes, exst)
   !
   allocate( dvscf_tmp(dfftp%nnr, nspin_mag, nmodes) )
   dvscf_tmp = czero

   do imod = 1, nmodes
      !read in dvscf for all the irreps
      call davcio(dvscf_tmp(:,:,imod), rec_len, iunit, imod, -1)
   end do
   close(iunit)
   
   !Transform from the basis of the patterns to cartesian basis
   ! dvscf(:,:,i) = dvscf(:,:,i) + conjg(u(i,j))*dvscf_tmp(:,:,j) 
   !        for i,j loop over [1, nmodes]
   ! a more efficient way to do the transformation
   nrow = dfftp%nnr * nspin_mag
   call zgemm('N','C', nrow, nmodes, nmodes, cone, &
               dvscf_tmp, nrow, u, nmodes, czero, dvscf, nrow)

   deallocate( dvscf_tmp )
end subroutine load_dvscf


!-----------------------------------------------------------------------
!   rotate dvscf_xqirr to get dvscf_xq
!-----------------------------------------------------------------------
subroutine get_dvscf_star(xq, xq_irr, isym_inv, q2mq, dvscf_irr, dvscf_xq)
   use kinds,            only : dp
   use constants,        only : tpi
   use fft_base,         only : dfftp
   use cell_base,        only : at, bg
   use ions_base,        only : nat, tau
   use symm_base,        only : ft, irt, s, sr
   use noncollin_module, only : nspin_mag
   !
   implicit none
   ! input variables:
   ! xq     : q vector, in cartesian coordinate
   ! xq_irr : q_irr vector, in cartesian coordinate
   real(dp),intent(in) :: xq(3), xq_irr(3)
   ! index of symmetry operation that rotate q -> qirr, e.g. sr(:,:,isym_inv)
   !  is the symmetry operation S^-1 that S^-1 q = q_irr (so q = S q_irr)
   integer, intent(in) :: isym_inv
   ! if q2mq is ture, then S^-1 (-q) = q_irr, 
   ! q and q_irr is connect via S^-1 and time reversal symmetry.
   logical, intent(in) :: q2mq
   ! dvscf_irr : dvscf of the irredicible q (xq_irr) in cartesian coordinate
   complex(dp), intent(in) :: dvscf_irr(dfftp%nnr, nspin_mag, 3*nat)
   ! dvscf_irr : dvscf of xq in cartesian coordinate
   complex(dp), intent(out) :: dvscf_xq(dfftp%nnr, nspin_mag, 3*nat)
   !
   ! local variables
   integer :: na, i, j, nmodes, imode0, ijks, nrpt
   integer :: is, k, n, nn, ri, rj, rk, ipol, index0, nar, ftau(3)
   ! auxiliary xq\cdot\tau and \xq_s\cdot\tau
   real(dp) :: xq_tau, xqirr_tau, xq_tmp(3), xq_rot(3)
   character(len=6), external :: int_to_char
   !
   complex(dp) :: phase_xqirr
   complex(dp), allocatable :: dvscf_tmp(:,:,:), phase_xq(:)

   nmodes = 3*nat
   !
   if(isym_inv < 1 .or. isym_inv > 48) call errore('get_dvscf_star', &
      'illegal symmetry operation: '//int_to_char(isym_inv), 1)
   
   xq_rot = merge(-xq, xq, q2mq)
   !check symmetry operation: S^-1 q = q_irr; 
   ! both xq_rot, xq_irr and sr are in cartesian coordinate.
   xq_tmp = matmul(sr(:,:,isym_inv), xq_rot) - xq_irr

   if(sqrt(dot_product(xq_tmp, xq_tmp)) > 1.0E-5_dp) then
      write(*, '(9(f8.4,1x))') matmul(sr(:,:,isym_inv), xq_rot), xq_irr, xq_tmp
      call errore('dvscf','S^-1 xq is not equal to xq_irr',1)
   endif

   !NOTE: dfftp%nnr = dfftp%nr1x * dfftp%nr2x * dfftp%nr3x 
   ! since we enforce nproc = 1 for each pool
   allocate( dvscf_tmp(dfftp%nnr, nspin_mag, nmodes) )
   !
   ! Transform to crystalline coordinates (necessary in order to apply s)
   dvscf_tmp =cmplx(0._dp, 0._dp, kind=dp)
   do i = 1, nat
      na = (i-1) * 3
      do j = 1, 3
         dvscf_tmp(:,:,na+j) = dvscf_irr(:,:,na+1)*at(1,j) + &
                               dvscf_irr(:,:,na+2)*at(2,j) + &
                               dvscf_irr(:,:,na+3)*at(3,j)
      enddo
   enddo
   !
   ! take away the phase due to the q-point
   dvscf_xq = cmplx(0._dp, 0._dp, kind=dp)
   do i = 1, nat
     !
     xqirr_tau = tpi * dot_product(xq_irr(:), tau(:,i))
     phase_xqirr= cmplx( cos(xqirr_tau), sin(xqirr_tau), kind=dp)
     !
     do ipol=1,3
        imode0 = (i-1)*3 + ipol
        dvscf_tmp(:,:,imode0) = phase_xqirr * dvscf_tmp(:,:,imode0)
     enddo
   enddo
   !
   ! Now rotate the dvscf
   allocate( phase_xq(nat) )
   !
   do i = 1, nat  
      xq_tau = - tpi * dot_product(xq_rot(:), tau(:,i))
      phase_xq(i) = cmplx(cos(xq_tau), sin(xq_tau), kind=dp)
   enddo
   
   ftau(1) = nint( ft(1,isym_inv) * dfftp%nr1 )
   ftau(2) = nint( ft(2,isym_inv) * dfftp%nr2 )
   ftau(3) = nint( ft(3,isym_inv) * dfftp%nr3 )
   !
   !do is = 1, nspin_mag
   !  kloop : do k = 1, dfftp%nr3
   !    jloop : do j = 1, dfftp%nr2
   !      iloop : do i = 1, dfftp%nr1

   dvscf_xq = cmplx(0._dp, 0._dp, kind=dp)
   nrpt = dfftp%nr1 * dfftp%nr2 * dfftp%nr3

!$omp parallel do schedule(guided) default(shared) private(ijks, is, n, &
!$omp&  i, j, k, ri, rj, rk, nn, na, nar, index0, ipol, imode0)
   do ijks = 1, nspin_mag * nrpt
      is = (ijks-1) / nrpt + 1
      !
      n = mod( (ijks-1), nrpt ) + 1
      !n  = (i-1)  + (j-1)*dfftp%nr1  + (k-1)*dfftp%nr2*dfftp%nr1  + 1
      i = mod(n-1, dfftp%nr1) + 1
      j = mod( (n-1)/dfftp%nr1, dfftp%nr2) + 1
      k = (n-1) / (dfftp%nr1 * dfftp%nr2) + 1
           
      ! ruotaijk find the rotated of i,j,k with the inverse of S
      ! a.k.a ruotaijk perform [S^-1] (i, j, k) = (ri, rj, rk)
      call ruotaijk( s(:, :, isym_inv), ftau, i, j, k, &
                    dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
      !
      nn = (ri-1) + (rj-1)*dfftp%nr1 + (rk-1)*dfftp%nr2*dfftp%nr1 + 1
      !
      do na = 1, nat
        nar = irt(isym_inv, na)
        !irt is the list of rotated of each atom: nar = [S^-1] na
        index0 = (nar-1) * 3
        !
        do ipol = 1, 3
            imode0 = (na-1)*3 + ipol
            !N.B. s(a, b, isym_inv) ==> [S^-1]_ba in C13 of PRB97,235146
            dvscf_xq(n, is, imode0) = dvscf_xq(n, is, imode0) + &
                ( s(ipol, 1, isym_inv) * dvscf_tmp(nn, is, index0+1) + &
                  s(ipol, 2, isym_inv) * dvscf_tmp(nn, is, index0+2) + &
                  s(ipol, 3, isym_inv) * dvscf_tmp(nn, is, index0+3) )
        enddo
      enddo
   enddo
!$omp end parallel do
   !
   ! Add back the phase factor for the new q-point
   !
   dvscf_tmp=cmplx(0._dp, 0._dp, kind=dp)
   do na = 1, nat
     do ipol = 1, 3
       imode0 = (na-1)*3 + ipol
       dvscf_tmp(:,:,imode0 ) = dvscf_xq(:,:,imode0) * phase_xq(na)
     enddo
   enddo
   !
   ! Back to cartesian coordinates
   !
   dvscf_xq = cmplx(0._dp, 0._dp, kind=dp)
   do i = 1, nat
     imode0 = (i-1)*3
     do j = 1, 3
       dvscf_xq(:,:,imode0+j) = dvscf_tmp(:,:,imode0+1)*bg(j,1) + &
             dvscf_tmp(:,:,imode0+2)*bg(j,2) + dvscf_tmp(:,:,imode0+3)*bg(j,3)
     enddo
   enddo
   !
   !dvscf(q) = conjg(dvscf(-q))
   if(q2mq) dvscf_xq = conjg( dvscf_xq )
   !
   deallocate(dvscf_tmp, phase_xq)
   !
end subroutine get_dvscf_star


! Modified version of subroutine diropn in QE/Module/io_file.f90
subroutine open_dvscf(iunit, fname, reclen, tot_rec, exst)
   use kinds,            only: dp
   use io_global,        only: stdout
   implicit none
   integer, intent(in) :: iunit, reclen, tot_rec
   !input: unit of the file to open
   !input: length of the records
   !input: total number of records
   character(len=*), intent(in) :: fname
   logical, intent(out) :: exst
   ! output: if true the file exists

   real(dp):: dummy
   integer*8 :: unf_recl, file_size
   ! double precision to prevent integer overflow
   integer :: ios, direct_io_factor
   logical :: opnd
   !
   if (iunit < 0) call errore ('open_dvscf', 'wrong unit', 1)
   inquire( unit = iunit, opened = opnd )
   if (opnd) call errore ('open_dvscf', "can't open a connected unit", abs(iunit))

   inquire (file =trim(fname), exist = exst)
   if(.not. exst) call errore('open_dvscf',trim(fname)//' does not exist!',1)
   
   if (reclen == -1) call errore('open_dvscf','negative record length',1)

   ! the  record length in direct-access I/O is given by the number of
   ! real*8 words times direct_io_factor (may depend on the compiler)
   !
   INQUIRE (IOLENGTH=direct_io_factor) dummy
   unf_recl = direct_io_factor * int(reclen, kind=kind(unf_recl))
   if (unf_recl <= 0) call errore ('open_dvscf', 'wrong record length', 3)
   !
   inquire(file=trim(fname), size=file_size)
   if( file_size .ne. tot_rec*unf_recl ) then
      write(stdout,'(5x,a)') "Warning (open_dvscf): " //trim(fname)// &
         " file size does not match. (turn on -assume byterecl if use intel compiler)"
   endif

   ios = 0
   open (iunit, file=trim(adjustl(fname)), iostat=ios, form='unformatted', &
        status = 'unknown', access = 'direct', recl = unf_recl, action='read')
   if (ios /= 0) call errore ('open_dvscf', 'error opening '//trim(fname), iunit)
   return
end subroutine open_dvscf
