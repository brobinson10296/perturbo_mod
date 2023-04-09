!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  based on the workflow in Phonon/PH/elphon.f90
!  Calculation of the electron-phonon matrix elements el_ph_mat
!     <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
!
!   NOTE: Segment fault may occur if OpenMP private space runs out. 
!   We allocated memory for large arrays in each thread. however, 
!   the private stack space each thread has is limited (defined by OMP_STACKSIZE).
!   The default value is implementation dependent and is usually quite small. 
!   A segmentation fault may appear if running out of space. We can specify 
!   the size of the private space by 'export OMP_STACKSIZE=xxK (or xxM, xxG)' 
!   or use OpenMP runtime subroutines: kmp_{set, get}_stacksize_s().
!
! Maintenance:
!===============================================================================

module elph_matrix
   use kinds, only: dp
   use constants, only: tpi
   use gvect, only: ngm, g
   use fft_base, only: dfftp, dffts
   use fft_interfaces, only: fwfft, invfft
   use cell_base, only: tpiba2, omega, at, bg
   use atom,   ONLY : msh, rgrid
   use pwcom, only: npwx
   use gvecs, only: doublegrid
   use uspp_param, only: upf, nhm, nh
   use uspp,  only: nlcc_any, nkb, okvan
   use lsda_mod, only : nspin
   use spin_orb, only : lspinorb
   use noncollin_module, only: noncolin, npol, nspin_mag
   use ions_base, only: nat, ityp, tau, ntyp => nsp
   implicit none
   public

   complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)

   complex(dp), allocatable, save :: int1(:,:,:,:,:),    int2(:,:,:,:,:)
   complex(dp), allocatable, save :: int1_nc(:,:,:,:,:), int2_so(:,:,:,:,:,:)
   complex(dp), allocatable, save :: int3(:,:,:,:,:),    int3_nc(:,:,:,:,:)
   complex(dp), allocatable, save :: drc(:,:)
   
   integer, parameter, private :: npe = 1
   
   public :: compute_elph_matrix
contains
! 
subroutine compute_elph_matrix()
   use buffers, only: open_buffer, save_buffer, close_buffer
   use m_gth,   only : setlocq_gth
   USE Coul_cut_2D, ONLY: do_cutoff_2D  
   use Coul_cut_2D_ph, only: cutoff_fact_qg, cutoff_lr_Vlocq
   use fft_interfaces, only: fft_interpolate
   use io_global, only: stdout
   use symm_base, only: nsym, s, sr, invs, irt
   !
   use eqv, only: dmuxc
   ! xq: q in cartesian coordinate
   use qpoint, only: xq, eigqts
   !
   use input_param, only: lwannier, num_band
   use qe_mpi_mod,  only: mp_split_pools
   use electronic_data, only: nkstot, evc_sub, xk_tot, ktable, ngk_tot, &
      igk_k_tot, kdim, rot_wan, kpt2num, becp1_tot, nb_sub
   !
   use lattice_data, only: numq, xq_tot, pattern, xq_irr, q2irq_symop, &
      identity_op, iq2irrq, q2minus_q

   implicit none
   complex(dp), pointer :: dvscf(:,:,:), dvscfins(:,:)!, dvscf_tmp(:,:,:)
   complex(dp), allocatable :: uact(:), vkb_kq(:,:), ephmat(:,:,:), dvscf_irrq(:,:,:)
   complex(dp), allocatable :: dvpsi_noloc(:,:), dvpsi(:,:), aux1(:,:), aux2(:), aux(:)
   integer, external :: find_free_unit
   complex(dp), external :: zdotc
   
   logical :: exst
   integer, parameter :: current_spin = 1
   integer :: nmodes, iun_eph, eph_nword, i, iq, iq_st, iq_end, nq_loc, irec, nrow
   integer :: na, ib, jb, ik, ikq, ip, npwq, nt, current_irrq, irrq, is, isym, isym_inv
   real(dp) :: arg, xkq(3), xkqg(3), xq_irr_cart(3), prog_rate
   integer,  allocatable :: igk_kqg(:)
   real(dp), allocatable :: vlocq(:,:) ! vlocq(ngm, ntyp)
   real(dp), allocatable :: rtau(:,:,:)

   write(stdout,'(a)') "start elph_matrix:"

   !check if all xk+xq are on the grid
   do ik = 1, nkstot
   do iq = 1, numq
      xkq = xk_tot(:,ik) + xq_tot(:,iq)
      ikq = kpt2num(xkq, kdim)
      if( ikq < 0 ) call errore('elph_matrix',"k- and q-grid are incommensurate",1)
   enddo; enddo

   ! Set non linear core correction stuff
   nlcc_any = any( upf(1:ntyp)%nlcc )
   
   nmodes = 3*nat
   iun_eph = find_free_unit()
   eph_nword = nb_sub * nb_sub * nkstot
   !io_level = 1, save eph matrix to disk
   call open_buffer(iun_eph, 'ephmat', eph_nword, 1, exst)
   
   if(.not. allocated(eigqts)) allocate(eigqts(nat))
   if(.not. allocated(dmuxc) ) allocate(dmuxc(dfftp%nnr, nspin_mag, nspin_mag))
   
   call allocate_work_arrays()
   allocate( uact(nmodes), vlocq(ngm, ntyp) )
   allocate( ephmat(nb_sub, nb_sub, nkstot) )
   
   !work arrays
   allocate( dvscf_irrq(dfftp%nnr, nspin_mag, nmodes) )
   allocate( dvscf(dfftp%nnr, nspin_mag, nmodes) )
   !allocate( dvscf_tmp(dfftp%nnr, nspin_mag, nmodes) )
   if(doublegrid) allocate( dvscfins(dffts%nnr, nspin_mag) )
   ! pre-run fwfft to setup FFT tables for dffts, to avoid race condition 
   !   when calling fwfft and invffw in openMP parallel region.
   allocate( aux(dffts%nnr) );  aux = (0.0_dp,  0.0_dp)
   call fwfft('Wave', aux, dffts)
   deallocate( aux )

   !
   allocate( rtau(3,48,nat) )
   call sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)

   !setup, see PH/phq_setup.f90
   !3) Computes the derivative of the XC potential
   !explicit output to eqv: dmuxc
   call setup_dmuxc()
   ! setup all gradient correction stuff
   call setup_dgc()
   
   !mpi parallelization
   call mp_split_pools(numq, iq_st, iq_end, nq_loc)
   
   current_irrq = -1
   dvscf_irrq = cmplx(0._dp, 0._dp, kind=dp)

   do i = 1, nq_loc
      iq = iq_st + i - 1
      
      !Note: xq and eigqts are global variables defined in module qpoint
      ! cutoff_fact_qg() and cutoff_lr_Vlocq() requires xq in qpoint (in cartesian coord.)
      ! xq_tot(:,iq): the iq-th point in crystal coordinate
      xq = xq_tot(:,iq)
      !transfer to cartesian coordinate
      call cryst_to_cart(1, xq, bg, 1)

      !compute  exp(-i xq times tau)
      do na = 1, nat
         arg = dot_product(xq, tau(:,na)) * tpi
         eigqts(na) = cmplx(cos(arg), -sin(arg), kind=dp)
      enddo
      
      !initializatiion (q-dependent) (following QE/PHonon/PH/phq_init)
      ! a0) compute rhocore for each atomic-type if needed for nlcc
      if(nlcc_any)  call set_drhoc(xq, drc)
      
      ! b) the fourier components of the local potential at q+G
      vlocq(:,:) = 0.d0
      do nt = 1, ntyp
         if(upf(nt)%tcoulombp) then
            call setlocq_coul( xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
         else if (upf(nt)%is_gth) then
            call setlocq_gth( nt, xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
         else
            call setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
               upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
         endif
      enddo
      
      !(following QE/PHonon/PH/phq_init)
      ! for 2d calculations, we need to initialize the fact for the q+G 
      ! component of the cutoff of the COulomb interaction
      ! NOTE: implicit input via module qpoint: xq; output Coul_cut_2D_ph: cutoff_2D_qg
      if(do_cutoff_2D) call cutoff_fact_qg()
      !  in 2D calculations the long range part of vlocq(g) (erf/r part)
      ! was not re-added in g-space because everything is caclulated in
      ! radial coordinates, which is not compatible with 2D cutoff. 
      ! It will be re-added each time vlocq(g) is used in the code. 
      ! Here, this cutoff long-range part of vlocq(g) is computed only once
      ! by the routine below and stored
      ! NOTE: implicit input via module qpoint: xq; output Coul_cut_2D_ph: lr_Vlocq
      if(do_cutoff_2D) call cutoff_lr_Vlocq() 
      
      !in the case of us-pp
      ! compute int1, int2, int1_nc, int2_so for xq
      if(okvan) call ph_dvanqq(xq, eigqts, vlocq)
      !!end of initialization for the current q
   
      !self-consistent part:Eq. B30 of PRB 64 235118
      !
      irrq = iq2irrq(iq)
      if( irrq .ne. current_irrq ) then
         current_irrq = irrq
         !dvscf_irrq in cartesian coordinates
         call load_dvscf(irrq, pattern(:,:,irrq), dvscf_irrq)
         
         !write(stdout,'(a, 2i5, 3f12.8)') "Finish loading dvscf:", iq, irrq, xq
      endif
      
      ! xq_irr(:,irrq): the irrq-th irreducible q-point in crystal coordinate
      xq_irr_cart = xq_irr(:, current_irrq)
      ! transfer to cartesian coordinate
      call cryst_to_cart(1, xq_irr_cart, bg, 1)
      
      ! sr(:,:,isym) is the symmetry operation S that Sq = q_irr
      isym = q2irq_symop(iq)
      if( (isym .eq. identity_op) .and. (.not. q2minus_q(iq)) ) then
         write(stdout,'(a15, i3, 2x, 3f12.8)') "irreducible q: ", current_irrq, xq_irr_cart(:)
         dvscf(:,:,:) = dvscf_irrq(:,:,:)
         !
         !uact(:, :) = pattern(:, :, current_irrq)
      else
         !write(stdout,'(3f12.8,2x,3f12.8,2x,i2,1x,l)') xq, xq_irr_cart(:), isym, q2minus_q(iq)
         !dvscf in cartesian coordinates for xq, rotated from dvscf_irrq
         call get_dvscf_star(xq, xq_irr_cart, isym, q2minus_q(iq), dvscf_irrq, dvscf(1,1,1))

         !!rotate displacement pattern
         !isym_inv = invs(isym)
         !call rotate_mod( pattern(:,:,current_irrq), uact, sr(:,:,isym_inv), &
         !      irt, rtau, xq_irr_cart, nat, isym)
         !! u(q) = conjg( u(-q) )
         !if(q2minus_q(iq))  uact = conjg(uact)
      endif

      !!transfer to pattern coordinate
      !nrow = dfftp%nnr * nspin_mag
      !call zgemm('N','N', nrow, nmodes, nmodes, cone, &
      !         dvscf_tmp(1,1,1), nrow, uact, nmodes, czero, dvscf(1,1,1), nrow)
      
      do ip = 1, nmodes
         !preparing int3 for the 2nd term in Eq.B30 of PRB 64 235118
         ! compute int3 and int3_nc
         if(okvan) call ph_newdq(xq, eigqts, dvscf(1,1,ip), npe)
         !
         uact     = cmplx(0.0_dp, 0.0_dp, kind=dp)
         uact(ip) = cmplx(1.0_dp, 0.0_dp, kind=dp)
   
         !compute 1st term in Eq.B30 of PRB 64 235118
         IF (doublegrid) THEN
            DO is = 1, nspin_mag
               call fft_interpolate(dfftp, dvscf(:,is,ip), dffts, dvscfins(:,is))
            ENDDO
         ELSE
            dvscfins => dvscf(:,:,ip)
         ENDIF
         !
         ephmat = cmplx(0.E0_dp, 0.E0_dp, kind=dp)
         !
!$omp parallel default(shared) private(ik, xkqg, ikq, npwq, xkq, ib, jb, &
!$omp&   dvpsi, dvpsi_noloc, vkb_kq, aux1, aux2, igk_kqg)
         allocate( dvpsi_noloc(npwx*npol, num_band), dvpsi(npwx*npol, nb_sub) )
         allocate( aux1(dffts%nnr, npol), aux2(npwx*npol), vkb_kq(npwx, nkb) )
         allocate( igk_kqg(npwx) )
!$omp do schedule(guided)
         do ik = 1, nkstot
            !xk+xq could be xk_tot(:,ikq)+G, the G need to be taken into account
            xkqg = xk_tot(:,ik) + xq_tot(:,iq)
            ikq = ktable( kpt2num(xkqg, kdim) )
            !transfer from crystal to cartesian coordinate, xk+xq with possible G
            call cryst_to_cart(1, xkqg, bg, 1)
            npwq = ngk_tot(ikq)
            ! xk+xq removed possible G, xkq is on k-grid.
            xkq = xk_tot(:,ikq)
            ! to cartesian coordiante
            call cryst_to_cart(1, xkq, bg, 1)
            !generate igk_k for xk+xq
            igk_kqg(:) = 0
            call get_igk_kqg(npwq, igk_k_tot(:,ikq), xkq, xkqg, igk_kqg(1))
            !
            !compute dvscf_q*psi_{k,i}>
            !
            !non-self-consistent part: Eq. B29 of PRB 64 235118
            call ph_dvqpsi_us_local &
               (ik, xq, eigqts, npwq, igk_kqg(1), uact, dmuxc, vlocq, dvpsi, .false.)
            
            !compute vkb_kq(k+q+G)
            call pw_init_us_2(npwq, igk_kqg(1), xkqg, vkb_kq)
            
            !non-local part
            call ph_dvqpsi_us_noloc &
               (ik, xkqg, npwq, igk_kqg(1), uact, vkb_kq, dvpsi_noloc)
            
            !apply the rotation if compute in wannier gauge
            if(lwannier) then
               ! if lwannier = .true. then nb_sub is the same as num_wann
               ! dvpsi = dvpsi_noloc * rot_wan(1:num_band,1:num_wann,ik) + dvpsi
               call zgemm('N','N', npwx*npol, nb_sub, num_band, cone, dvpsi_noloc, &
                        npwx*npol, rot_wan(1,1,ik), num_band, cone, dvpsi, npwx*npol)
            else
               ! if lwannier = .false. then nb_sub is the same as num_band
               dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + dvpsi_noloc(:,1:nb_sub)
            endif
               
            if(okvan) then
               !add the 2nd term in Eq.B30 of PRB 64 235118
               call ph_adddvscf(npwq, npe, becp1_tot(ik), vkb_kq, dvpsi_noloc)
               !apply the rotation
               if(lwannier) then
                  ! dvpsi = dvpsi_noloc * rot_wan(1:num_band,1:num_wann,ik) + dvpsi
                  call zgemm('N','N', npwx*npol, nb_sub, num_band, cone, dvpsi_noloc, &
                        npwx*npol, rot_wan(1,1,ik), num_band, cone, dvpsi, npwx*npol)
               else
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + dvpsi_noloc(:,1:nb_sub)
               endif
            endif

            !NOTE: we enforce dffts%has_task_groups =.false. (only support 1 processor/pool )
            do ib = 1, nb_sub
               aux2(:) = czero
               !self-consistent part: first term of Eq. B30
               !call cft_wave (ik, evc_sub(1, ib, ik), aux1, +1)
               ! Inverse Fourier transform of a wavefunction: G-space to r-space
               call invfft_wave(ngk_tot(ik), igk_k_tot(1,ik), evc_sub(1,ib,ik), aux1)
               !
               !! NOTE: current_spin and isk only matter in lsda calc., 
               ! and the number of k-points is doubled:  the first half for spin up (current_spin=1),
               !  the second half for spin down (current_spin=2).
               ! for non-spin and non-colinear calc. current_spin (=1) does not matter.
               ! In this interface program, we drop support for lsda calc.
               ! for spin-polarized calculation, using non-colinear without soc instead
               call apply_dpot(dffts%nnr, aux1, dvscfins(1,1), current_spin)
               !
               ! Fourier transform of a wavefunction: r-space to G-space 
               !call cft_wave (ik, aux2, aux1, -1)
               call fwfft_wave(npwq, igk_kqg(1), aux2, aux1)
               !
               dvpsi(:,ib) = dvpsi(:,ib) + aux2(:)
            enddo
            
            ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
            do ib = 1, nb_sub
            do jb = 1, nb_sub
               ephmat(jb,ib,ik) = zdotc(npwq, evc_sub(1,jb,ikq), 1, dvpsi(1,ib), 1)
               if(noncolin) then
                  ephmat(jb,ib,ik) = ephmat(jb,ib,ik) + &
                    zdotc(npwq, evc_sub(npwx+1,jb, ikq), 1, dvpsi(npwx+1,ib), 1)
               endif
            enddo
            enddo
         enddo
!$omp end do
         deallocate(vkb_kq, dvpsi, dvpsi_noloc, aux1, aux2, igk_kqg)
!$omp end parallel
         
         !save result
         irec = (i - 1) * nmodes + ip
         call save_buffer(ephmat, eph_nword, iun_eph, irec)
         
         prog_rate = 1.0E2_dp * real(irec, dp) / real(nq_loc*nmodes, dp)
         write(stdout,'(a, 2i4, 8x, a, f7.2, a)') &
            "iq_loc, imode: ", i, ip, "progress: ", prog_rate, "%"
      enddo
   enddo
   call close_buffer(iun_eph, 'KEEP')
   
   !release memory
   call deallocate_work_arrays()
   if(doublegrid) deallocate(dvscfins)
   deallocate( dvscf_irrq, dvscf)
   deallocate(dmuxc, uact, eigqts, ephmat, vlocq, rtau)
end subroutine compute_elph_matrix

subroutine allocate_work_arrays()
   implicit none
   !work arrays
   !if .TRUE. at least one pseudo is Vanderbilt
   if (okvan) then
      allocate( int1(nhm, nhm, 3, nat, nspin_mag) )
      allocate( int2(nhm, nhm, 3, nat, nat) )
      allocate( int3(nhm, nhm, nat, nspin_mag, npe) )
      
      if( noncolin ) then
         allocate( int1_nc(nhm, nhm, 3, nat, nspin) )
         allocate( int3_nc(nhm, nhm, nat, nspin, npe) )
         if( lspinorb ) then
            allocate( int2_so(nhm, nhm, 3, nat, nat, nspin) )
         endif
      endif
   endif
   if(nlcc_any) allocate ( drc(ngm, ntyp) )
end subroutine allocate_work_arrays

subroutine deallocate_work_arrays()
   implicit none
   !work arrays
   !if .TRUE. at least one pseudo is Vanderbilt
   if (okvan) then
      deallocate( int1, int2, int3)
      if( noncolin ) then
         deallocate( int1_nc, int3_nc )
         if( lspinorb ) deallocate( int2_so )
      endif
   endif
   if(allocated(drc)) deallocate(drc)
end subroutine deallocate_work_arrays

end module elph_matrix
