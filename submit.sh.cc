#!/bin/bash

#SBATCH --job-name="perturbo"
#SBATCH --mail-user=brianr5@illinois.edu
#SBATCH -N 4
#SBATCH --ntasks-per-node=2
#SBATCH -t 04:00:00
#SBATCH -p secondary
#SBATCH --exclude=ccc0093,ccc0094,ccc0095,ccc0049,ccc0275,golub041

module load mvapich2/2.3rc2-intel-18.0
LD_LIBRARY_PATH=/projects/sg/brianr5/libraries/hdf5/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR

echo "start: $(date)"

export OMP_NUM_THREADS=5
export MV2_USE_THREAD_WARNING=0

PREFIX='INPUT'
TASKS=4

mpirun -np $TASKS /projects/sg/brianr5/src/quantum_espresso/qe_6.5/perturbo_mod/bin/perturbo_mod.x -npools $TASKS -i $PREFIX.in > $PREFIX.out

echo "end: $(date)"

