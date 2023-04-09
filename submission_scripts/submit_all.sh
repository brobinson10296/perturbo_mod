#!/bin/bash

#submit pert_setup
setup_id=$(sbatch --parsable submit_setup.sh)

#submit pert_run
run_id=$(sbatch  --parsable --dependency=afterok:$setup_id  submit_run.sh)

#submit pert_pp
pp_id=$(sbatch  --parsable --dependency=afterok:$run_id  submit_pp.sh)
