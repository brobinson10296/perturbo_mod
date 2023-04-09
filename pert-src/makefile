
include ../make.sys


PERTMOD = \
HDF5_utils.f90 \
wigner_seitz_cell.f90 \
polar_correction.f90 \
electron_wannier.f90 \
force_constant.f90 \
elph_matrix_wannier.f90 \
epwan_hdf5_io.f90 \
pert_const.f90 \
qe_mpi_mod.f90 \
pert_utils.f90 \
pert_param.f90  \
pert_data.f90   \
pert_output.f90 \
band_structure.f90 \
phonon_dispersion.f90 \
elphon_coupling_matrix.f90 \
boltz_utils.f90 \
vector_list.f90 \
boltz_grid.f90  \
imaginary_selfenergy.f90 \
boltz_trans_output.f90 \
boltz_trans_mod.f90 \
boltz_scatter.f90 \
boltz_scatter_integral.f90 \
boltz_grid_neighbors.f90 \
boltz_external_field.f90 \
boltz_dynamics_solver.f90 \
boltz_dynamics_mod.f90


PERTSRC = \
environment_setup.f90 \
weight_dos.f90 \
load_data_pert.f90 \
calc_bands.f90 \
calc_ephmat.f90 \
calc_electron_mfp.f90 \
electron_imsigma.f90 \
boltz_setup.f90 \
transport.f90 \
trans_postproc.f90 \
carrier_dynamics_run.f90 \
carrier_dynamics_postproc.f90 \
perturbo.f90


OBJ = $(PERTMOD:.f90=.o) $(PERTSRC:.f90=.o)

EXE = perturbo_mod.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

%.o: %.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

boltz_grid.o: boltz_grid.f90 boltz_grid_generate.f90 boltz_grid_load.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

boltz_scatter.o: boltz_scatter.f90 setup_kq_pair.f90 compute_scatter_eph_g2.f90 load_scatter_eph_g2.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

clean:
	rm -f *.o *~ *.d *.mod *.i
