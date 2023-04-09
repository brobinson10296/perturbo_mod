!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  save g(R_e, R_p), H(R_e), \Phi(R_p), and other basic info. to a hdf5 file.
!
! Maintenance:
!===============================================================================

subroutine save_perturbo_data()
   use kinds,     only: dp
   use ions_base, only: nat, tau
   use cell_base, only: at, bg
   use io_global, only: stdout, ionode
   use lattice_data, only: qdim, numq, xq_tot, dynmat
   use electronic_data, only: wannier_center_cryst, xk_tot, et_sub, et_corr, rot_wan, nkstot
   use input_param, only: num_wann, prefix, kdim, num_band, debug, eig_corr
   !
   use force_constant, only: lattice_ifc
   use electron_wannier, only: electron_wann
   use elph_matrix_wannier, only: elph_mat_wann
   !
   use epwan_hdf5_io, only: write_force_constant, write_electron_wannier
   use hdf5_utils
   implicit none
   integer(HID_T) :: file_id
   character(len=80) :: fname
   real(dp) :: ctau(3, nat)
   integer, external :: find_free_unit
   !
   type(lattice_ifc) :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph
   
   ctau(:,:) = tau(:,:)
   call cryst_to_cart(nat, ctau, bg, -1)
   !forceconstant
   call calc_lattice_ifc(phon, qdim, nat, at, ctau, numq, xq_tot, dynmat)
   !electron
   if( trim(eig_corr) .eq. '' ) then
      call calc_electron_wann(elec, kdim, num_wann, at, &
         wannier_center_cryst, num_band, nkstot, xk_tot, et_sub, rot_wan)
   else
      call calc_electron_wann(elec, kdim, num_wann, at, &
         wannier_center_cryst, num_band, nkstot, xk_tot, et_corr, rot_wan)
   endif
   !
   !write out basic info
   if(ionode) then
      write(stdout, '(5x,a)') "save electron and phonon data to file ..."
      !output to hdf5 file
      fname = trim(prefix)//"_epwan.h5"
      call hdf_open_file(file_id, trim(fname), status='NEW')
      !
      call write_basic_data_hdf5(file_id)
      !
      call write_force_constant(file_id, phon)
      !
      call write_electron_wannier(file_id, elec)
      !close file
      call hdf_close_file(file_id)
   endif
   !eph
   call calc_elph_mat_wann(elph)

end subroutine save_perturbo_data


subroutine write_basic_data_hdf5(file_id)
   use kinds, only: dp
   use ions_base, only: amass, ityp, nat, tau
   use cell_base, only: omega, at, bg, alat
   use lattice_data, only: zstar, epsil, qdim
   use electronic_data, only: wannier_center, wannier_center_cryst, noncolin
   use input_param, only: num_wann, kdim, polar_alpha, thickness_2d, system_2d
   use symm_base, only: s, nsym
   use hdf5_utils
   implicit none
   integer(HID_T), intent(in) :: file_id
   !local
   integer :: i, s2d
   integer(HID_T) :: group_id
   real(dp) :: tmp_mass(nat)

   do i = 1, nat
      tmp_mass(i) = amass(ityp(i))
   enddo

   call hdf_create_group(file_id, 'basic_data')
   call hdf_open_group(file_id, 'basic_data', group_id)
   !
   call hdf_write_dataset(group_id, "num_wann", num_wann)
   call hdf_write_dataset(group_id, "nat",       nat)
   call hdf_write_dataset(group_id, "kc_dim",  kdim)
   call hdf_write_dataset(group_id, "qc_dim",  qdim)
   call hdf_write_dataset(group_id, "volume",  omega)
   call hdf_write_dataset(group_id, "epsil",   epsil)
   call hdf_write_dataset(group_id, "alat",     alat)
   call hdf_write_dataset(group_id, "at",         at)
   call hdf_write_dataset(group_id, "bg",         bg)
   call hdf_write_dataset(group_id, "nsym",  nsym)
   call hdf_write_dataset(group_id, "symop", s(:,:,1:nsym))
   ! Ewald parameter for e-ph polar correction only!
   call hdf_write_dataset(group_id, "polar_alpha", polar_alpha)
   call hdf_write_dataset(group_id, 'spinor', merge(1, 0, noncolin))
   
   s2d = merge(1, 0, system_2d)
   call hdf_write_dataset(group_id, "system_2d", s2d)
   if(system_2d) call hdf_write_dataset(group_id, "thickness_2d", thickness_2d)

   call hdf_write_dataset(group_id, "mass",  tmp_mass)
   call hdf_write_dataset(group_id, "zstar", zstar(:,:,1:nat) )
   call hdf_write_dataset(group_id, "tau",   tau(:, 1:nat) )
   call hdf_write_dataset(group_id, "wannier_center", wannier_center(:,1:num_wann))
   call hdf_write_dataset(group_id, "wannier_center_cryst", wannier_center_cryst(:,1:num_wann))

   call hdf_close_group(group_id)
end subroutine write_basic_data_hdf5
