!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!
! Maintenance:
!===============================================================================
module pert_param
   use pert_const, only: dp, ryd2mev, ryd2ev, kelvin2eV, bohr2ang, timeunit
   implicit none
   save
   character(len=80) :: prefix
   character(len=80) :: tmp_dir
   type, public :: temperature
      integer :: nt
      real(dp), pointer :: kt(:) => null()
      real(dp), pointer :: ef(:) => null()
   end type
     ! jjzhou
   character(len=80) :: calc_mode
   character(len=80) :: fklist
   character(len=80) :: fqlist
   character(len=80) :: ftemper !file store temperature, chemical potential, doping level
   !polar_split
   !'polar',   only polar part;
   !'rmpol', whole - polar (remaider part).
   !others,  the whole, no split;
   character(len=80) :: polar_split

!   logical :: lpolar
   logical :: find_efermi
   logical :: spinor ! if spinor is .true., then spin degeneracy will is added in DOS.
   logical :: hole   !.false.: electron, .true.: hole
   logical :: full_ite  ! if .true. solve BTE with both E- and T-field iteratively.
   !logical :: tdep   ! .true.: phonon from TDEP
   logical :: use_mem
   logical :: load_scatter_eph

!   logical :: boltz_restart ! .true. : enable restart in carrier dynamics simulations
   logical :: debug   ! .true. : debug mode
   integer :: boltz_kdim(3)
   integer :: boltz_qdim(3)  ! dimension for q-grid
   ! max number of steps in carrier dynamics or iterative steps in transport calc.
   ! if boltz_nstep = 0, RTA calculations, otherwise iterative solutions
   integer :: boltz_nstep
   ! band_min, band_max: define bands to be computed.
   integer :: band_min, band_max !, num_bands
   ! boltz_emin, boltz_emax, and boltz_kdim are energy window and kgrid for tranpsort calculation.
   ! boltz_de, in mev energy step in tetrahedron integration
   ! trans_thr,  threshold of iterative procedure.
   !    (converged if relative error of two consecutive steps is within trans_thr)
   real(dp) :: boltz_emin, boltz_emax, boltz_de, trans_thr
   ! parameter for the gaussian broadening of the delta function
   real(dp) :: delta_smear
   ! phonon energy cutoff, phonon with energy smaller than phfreq_cutoff will be excluded.
   real(dp) :: phfreq_cutoff

   integer :: ntemper
   real(dp), allocatable :: temper(:)
   real(dp), allocatable :: doping(:)
   real(dp), allocatable :: efermi(:)
   type(temperature) :: tempers

   !for carrier dynamics
   real(dp) :: boltz_efield(3)  ! External Electric Field in V/cm
   real(dp) :: time_step  !time step of carrier dynamcis,in femto-second
   integer  :: output_nstep  ! output results for every 'output_nstep'
   character(len=80) :: solver !options: 'euler', 'rk4'
   character(len=80) :: boltz_init_dist !options: 'restart', 'lorentz', 'fermi', 'gaussian'
   !used if boltz_init_dist is lorentz, fermi or gaussian:
   !lorentz : f(k) = 0.1 / (((Ek-e0)/smear)^2 + 1.0 )
   !gaussian: f(k) = gauss_height * exp( -((Ek-e0)/smear)^2 )
   !fermi   : f(k) = 1.0 / (exp((Ek-e0)/smear) + 1.0)
   real(dp) :: boltz_init_e0 !eV
   real(dp) :: boltz_init_smear !meV
   real(dp) :: E_fermi !eV BR ADDED
   real(dp) :: gauss_height !BR ADDED
   !character(len=80) :: boltz_init_fname !restart from the last step of 'fname'

   !for two phonon scattering rate calculation
   !real(dp) :: smear_2ph
   character(len=80) :: sampling !options: 'uniform', 'cauchy'
   real(dp) :: cauchy_scale ! scale parameter gamma for cauchy distribution
   integer :: nsamples ! number of samples to generate

contains
subroutine init_input_param()
   use io_files, only: check_tempdir
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: meta_ionode, meta_ionode_id, world_comm, mp_bcast, stdout
   implicit none
   logical :: ltmp1, ltmp2
   integer :: iunit, i, ios
   character(len=120) :: ctmp, msg
   !
   CHARACTER(LEN=256), EXTERNAL :: trimcheck

  ! namelist /perturbo/ prefix, calc_mode, fklist, fqlist, ftemper, hole, full_ite, debug, boltz_kdim, boltz_qdim, &
  !    band_min, band_max, boltz_emax, boltz_emin, tmp_dir, use_mem, boltz_de, trans_thr, &
  !    delta_smear, phfreq_cutoff, boltz_nstep, polar_split, load_scatter_eph, solver, &
  !    boltz_efield, time_step, boltz_init_dist, boltz_init_smear, boltz_init_e0, output_nstep, &
  !    sampling, cauchy_scale, nsamples

   namelist /perturbo/ prefix, calc_mode, fklist, fqlist, ftemper, hole, full_ite, debug, boltz_kdim, boltz_qdim, &
        band_min, band_max, boltz_emax, boltz_emin, tmp_dir, use_mem, boltz_de, trans_thr, &
        delta_smear, phfreq_cutoff, boltz_nstep, polar_split, load_scatter_eph, solver, &
        boltz_efield, time_step, boltz_init_dist, boltz_init_smear, boltz_init_e0, output_nstep, &
        sampling, cauchy_scale, nsamples, E_fermi, gauss_height !BR add
   !set default value
   prefix       = ''
   tmp_dir      = '.'
   calc_mode    = '' !default
   polar_split  = '' !default, on split
   fqlist       = ''
   fklist       = ''
   ftemper      = ''
   debug        = .false.
   !spinor      = .false. ! deprecated, will read from hdf5 file instead.
   !tdep         = .false.
   hole         = .false. ! default, electron carrier
   full_ite     = .false. !
   use_mem      = .true.  ! only used in boltz_scatter
   load_scatter_eph = .false.  ! read eph_g2 from files

   boltz_kdim   =  1
   boltz_qdim   =  0
   band_min     =  1    !the default will be 1
   band_max     =  9999999  !a very big number, will be set to nband
   boltz_emin   = -9999.0_dp !eV,
   boltz_emax   =  9999.0_dp !eV,
   boltz_nstep  =  0
   boltz_de     =  1.0_dp   !meV,
   trans_thr    =  2.E-3_dp !0.2% as default
   delta_smear  =  10.0_dp  !meV
   phfreq_cutoff= 1.0_dp    !meV
   !for dynamics
   boltz_efield = 0.0_dp !turn off E-field by default
   time_step    = 1.0_dp !fs
   solver       = 'rk4' !use rk4 by default
   boltz_init_dist = ''  !no default value, should be specified.
   boltz_init_e0   = -9999.0_dp !eV, should be specified
   boltz_init_smear  = 20.0_dp  !meV, 20 by default
   E_fermi = 0.0_dp !eV 0 by default should be specified BR ADDED
   gauss_height = 0.1 !by default BR ADDED
   output_nstep  = 1

   sampling = 'uniform' ! other options 'cauchy'
   cauchy_scale = 0.05
   nsamples = 100000 ! 0.1 million by default

   !readin parameters
   if(meta_ionode) then
      call input_from_file()

      read(5, perturbo, err=100, iostat=ios)
100   call errore('init_input_para','reading input namelist',abs(ios))

      tmp_dir = trimcheck(tmp_dir)
      !do some check
      !band_min and band_max should be no less than 1
      if(band_min > band_max .or. band_min < 1 .or. band_max < 1) then
         msg = "both band_min and band_max should > 0, and band_min > band_max"
         call errore('init_input_para', trim(msg), 1)
      endif

      fqlist     = adjustl(fqlist)
      fklist     = adjustl(fklist)
      ftemper   = adjustl(ftemper)
      calc_mode = adjustl(calc_mode)
      polar_split = adjustl(polar_split)

      if(any(boltz_kdim(1:3) < 1)) &
         call errore('init_input_para','illegal boltz_kdim',1)

      if( all(boltz_qdim(1:3) .eq. 0) ) then
         !by default, boltz_qdim = boltz_kdim
         boltz_qdim(:) = boltz_kdim(:)
      elseif( any(boltz_kdim(1:3) < 1) ) then
         call errore('init_input_para','boltz_qdim should all be positive!',1)
      elseif( any( mod(boltz_kdim(:), boltz_qdim(:)) .ne. 0 ) ) then
         call errore('init_input_para','boltz_qdim is incommensurate with boltz_kdim',1)
      endif

      if(boltz_emin>boltz_emax .or. boltz_de<1.0E-3_dp ) call errore &
         ('init_input_para','illegal boltz_emax, boltz_emin or boltz_de.',1)

      boltz_emin = boltz_emin/ryd2ev
      boltz_emax = boltz_emax/ryd2ev
      if(boltz_nstep < 0) boltz_nstep = 0
      ! from mev to Ryd
      boltz_de = boltz_de/ryd2mev
      phfreq_cutoff = phfreq_cutoff/ryd2mev
      delta_smear   = delta_smear/ryd2mev
      if(trans_thr < 1.0E-16) &
         call errore('init_input_param','trans_thr is too small or negative')

      !open temperature file
      iunit = find_free_unit()
      find_efermi  = .false.
      if(trim(prefix)  .eq. '') prefix = 'pert'
      if(trim(ftemper) .eq. '') then
         ntemper = 0
      else
         open(iunit,file=trim(ftemper),form='formatted',err=101,iostat=ios)
101      call errore('init_input_para','opening file '//trim(ftemper),ios)
         read(iunit,*,iostat=ios) ntemper, ctmp
         if(ios .ne. 0) call errore('init_input_param', &
            'reading ntemper in file '//trim(ftemper), 1)
         ctmp = trim( adjustl(ctmp) )
         call lower_case( ctmp )
         if( ctmp(1:1) .eq. 't' ) then
            find_efermi = .true.
         else if( ctmp(1:1) .eq. 'f' ) then
            find_efermi = .false.
         else
            find_efermi = .false.
            write(stdout,'(5x, a)') &
               "Warning: illegal mode in " // trim(ftemper) // ". Set to 'F' "
         endif
         !
         if(ntemper < 1) call errore('init_input_param', &
            '#. temperatures < 1 in '// trim(ftemper), 1)
         !
         allocate(temper(ntemper), doping(ntemper), efermi(ntemper))
         ! temper in K;  efermi in eV;   doping in cm^-3
         temper = 0.0_dp; efermi = 0.0_dp; doping = 0.0_dp
         do i = 1, ntemper
            read(iunit,*,iostat=ios) temper(i), efermi(i), doping(i)
            if(ios .ne. 0) call errore('init_input_para', &
               'reading temper in file '//trim(ftemper), 1)
         enddo
         close(iunit)
         temper = temper*kelvin2eV/ryd2ev
         efermi = efermi / ryd2ev
         ! do the conversion in a later stage since we don't know it's 2D or 3D system.
         !for 3D from #./cm^3 to #./bohr^3, for 2D from #./cm^2 to #./bohr^2
         !doping = doping*1.0E-24_dp*(bohr2ang)**3
      endif

      !for dynamics
      if(time_step < 0.0_dp) call errore('init_input_param','negative time step',1)
      !convert to Rydberg atomic unit
      boltz_init_e0 = boltz_init_e0 / ryd2ev
      boltz_init_smear = boltz_init_smear / ryd2mev
      E_fermi = E_fermi / ryd2ev !BR ADDED
      time_step = time_step / (timeunit*1.0E15_dp)
      !from e*V/cm to Rydberg atomic unit (e*E)
      ! convert to Rydberg atomic unit: eE, eV/cm -> E_Rydberg / Bohr
      boltz_efield(1:3)  = boltz_efield(1:3)*bohr2ang/ryd2ev*1.0E-8_dp
      boltz_init_dist = trim(adjustl(boltz_init_dist))
      solver = trim(adjustl(solver))
      if(trim(solver) .ne. 'euler' .and. solver .ne. 'rk4') &
         call errore('init_input_param',"solver should be 'euler' or 'rk4'.", 1)

      sampling = adjustl(sampling)
      if(trim(sampling) .eq. '')  sampling = 'uniform'
   endif

   !broadcast
   call mp_bcast(prefix, meta_ionode_id, world_comm)
   call mp_bcast(tmp_dir, meta_ionode_id, world_comm)
   call mp_bcast(calc_mode, meta_ionode_id, world_comm)
   call mp_bcast(fklist, meta_ionode_id, world_comm)
   call mp_bcast(fqlist, meta_ionode_id, world_comm)
   call mp_bcast(ftemper, meta_ionode_id, world_comm)
   call mp_bcast(polar_split, meta_ionode_id, world_comm)
   !
   call mp_bcast(find_efermi, meta_ionode_id, world_comm)
   call mp_bcast(hole, meta_ionode_id, world_comm)
   call mp_bcast(full_ite, meta_ionode_id, world_comm)
   call mp_bcast(debug, meta_ionode_id, world_comm)
   !call mp_bcast(spinor, meta_ionode_id, world_comm)
   !call mp_bcast(tdep, meta_ionode_id, world_comm)
   call mp_bcast(use_mem, meta_ionode_id, world_comm)
   call mp_bcast(load_scatter_eph, meta_ionode_id, world_comm)
   !
   call mp_bcast(boltz_kdim, meta_ionode_id, world_comm)
   call mp_bcast(boltz_qdim, meta_ionode_id, world_comm)
   call mp_bcast(boltz_nstep, meta_ionode_id, world_comm)
   call mp_bcast(band_min, meta_ionode_id, world_comm)
   call mp_bcast(band_max, meta_ionode_id, world_comm)
   !
   call mp_bcast(boltz_emax, meta_ionode_id, world_comm)
   call mp_bcast(boltz_emin, meta_ionode_id, world_comm)
   call mp_bcast(boltz_de, meta_ionode_id, world_comm)
   call mp_bcast(trans_thr, meta_ionode_id, world_comm)
   call mp_bcast(delta_smear, meta_ionode_id, world_comm)
   call mp_bcast(phfreq_cutoff, meta_ionode_id, world_comm)
   !
   call mp_bcast(ntemper, meta_ionode_id, world_comm)
   tempers%nt = ntemper
   if(ntemper > 0) then
      if(.not. allocated(temper)) allocate(temper(ntemper))
      if(.not. allocated(efermi)) allocate(efermi(ntemper))
      if(.not. allocated(doping)) allocate(doping(ntemper))
      call mp_bcast(temper, meta_ionode_id, world_comm)
      call mp_bcast(doping, meta_ionode_id, world_comm)
      call mp_bcast(efermi, meta_ionode_id, world_comm)

      if(associated(tempers%kt)) deallocate(tempers%kt)
      if(associated(tempers%ef)) deallocate(tempers%ef)
      allocate( tempers%kt(ntemper), tempers%ef(ntemper) )
      tempers%kt(:) = temper(:)
      tempers%ef(:) = efermi(:)
   endif
   call check_tempdir(tmp_dir, ltmp1, ltmp2)
   !
   !for dynamics
   call mp_bcast(boltz_efield, meta_ionode_id, world_comm)
   call mp_bcast(time_step, meta_ionode_id, world_comm)
   call mp_bcast(solver, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_dist, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_e0, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_smear, meta_ionode_id, world_comm)
   call mp_bcast(E_fermi, meta_ionode_id, world_comm) !BR add
   call mp_bcast(gauss_height, meta_ionode_id, world_comm) !BR add
   call mp_bcast(output_nstep, meta_ionode_id, world_comm)

   call mp_bcast(sampling,      meta_ionode_id, world_comm)
   call mp_bcast(cauchy_scale,  meta_ionode_id, world_comm)
   call mp_bcast(nsamples,      meta_ionode_id, world_comm)
end subroutine init_input_param


!convert string to lower case
subroutine lower_case(string)
   character(len=*) :: string
   integer  :: i, ic, nlen

   nlen = len(string)
   do i = 1, nlen
      ic = ichar( string(i:i) )
      if( ic >= 65 .and. ic < 90 ) string(i:i) = achar(ic+32)
   end do
end subroutine lower_case

end module pert_param
