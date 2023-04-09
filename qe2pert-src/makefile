
include ../make.sys


PERTDIR=../pert-src

PERTOBJ = \
HDF5_utils.f90 \
qe_mpi_mod.f90 \
environment_setup.f90 \
wigner_seitz_cell.f90 \
polar_correction.f90 \
electron_wannier.f90 \
force_constant.f90 \
elph_matrix_wannier.f90 \
epwan_hdf5_io.f90 \


Q2PTMOD = \
input_param.f90 \
electronic_data.f90 \
lattice_data.f90 \
elph_matrix.f90 \


Q2PTSRC = \
get_igk_kqg.f90 \
pw_read_file.f90 \
pw_read_collected_to_evc.f90 \
pw_init_us_2.f90 \
ph_transform_int.f90 \
ph_transform_int3.f90 \
ph_dvqpsi_us_local.f90 \
ph_dvqpsi_us_noloc.f90 \
ph_dvanqq.f90 \
ph_adddvscf.f90 \
ph_newdq.f90 \
dvscf.f90 \
calc_electron_wann.f90 \
calc_lattice_ifc.f90 \
calc_elph_mat_wann.f90 \
save_perturbo_data.f90 \
collect_elph_matrix.f90 \
qe2pert.f90


OBJ = $(PERTOBJ:.f90=.o) $(Q2PTMOD:.f90=.o) $(Q2PTSRC:.f90=.o)

EXE = qe2pert.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

HDF5_utils.o: ${PERTDIR}/HDF5_utils.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

qe_mpi_mod.o: ${PERTDIR}/qe_mpi_mod.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

environment_setup.o: ${PERTDIR}/environment_setup.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

wigner_seitz_cell.o: ${PERTDIR}/wigner_seitz_cell.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

polar_correction.o: ${PERTDIR}/polar_correction.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

electron_wannier.o: ${PERTDIR}/electron_wannier.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

force_constant.o: ${PERTDIR}/force_constant.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

elph_matrix_wannier.o: ${PERTDIR}/elph_matrix_wannier.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

epwan_hdf5_io.o: ${PERTDIR}/epwan_hdf5_io.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

%.o: %.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

clean:
	rm -f *.o *~ *.d *.mod *.i

