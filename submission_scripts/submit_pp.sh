#!/bin/bash
#SBATCH --job-name="pert_run"
#SBATCH --mail-user=brianr5@illinois.edu
#SBATCH -N 2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH -t 04:00:00
#SBATCH --mem=128G
#SBATCH -p secondary
#SBATCH --exclude=golub342,ccc[0189,0218],ccc[0297-0298],ccc0286,ccc[0309-0310],ccc[0291,0293],ccc[0305,0337]

module load mvapich2/2.3rc2-intel-18.0
LD_LIBRARY_PATH=/projects/sg/brianr5/libraries/hdf5/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR

echo "start: $(date)"
echo $SLURM_JOB_NODELIST

export MV2_ENABLE_AFFINITY=0
PROCESSES=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MV2_USE_THREAD_WARNING=0

PREFIX='pert_pp'

mpirun -np $PROCESSES /projects/sg/brianr5/src/quantum_espresso/qe_6.5/perturbo_mod_beta/bin/perturbo_mod_beta.x -npools $PROCESSES -i $PREFIX.in > $PREFIX.out

echo "end: $(date)"
