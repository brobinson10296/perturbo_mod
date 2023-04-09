#!/bin/bash

echo "start: $(date)"

export OMP_NUM_THREADS=5
PREFIX='INPUT HERE'

#mpirun -np 1 ~/src/qe_6.5/perturbo_mod/bin/perturbo_mod.x -npools 1 -i $PREFIX.in > $PREFIX.out
mpirun -np 1 ~/src/qe_6.5/perturbo_mod/bin/perturbo_mod.x -npools 1 -i $PREFIX.in

echo "end: $(date)"
