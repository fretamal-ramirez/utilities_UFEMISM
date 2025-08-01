#!/bin/bash

#SBATCH -J test_ufe
#SBATCH -t 12:00:00
#SBATCH -n 16
#SBATCH -o out_run_UFEMISM
#SBATCH -e err_run_UFEMISM

mpiexec UFEMISM_program /home/x_frare/utilities_UFEMISM/config-files_tetralith/config_tetralith_PD_spinup.cfg
