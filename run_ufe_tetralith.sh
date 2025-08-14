#!/bin/bash

#SBATCH -J test_ufe
#SBATCH -t 12:00:00
#SBATCH -n 16
#SBATCH -o out_run_UFEMISM
#SBATCH -e err_run_UFEMISM

mpiexec UFEMISM_program /home/x_frare/utilities_UFEMISM/config-files_tetralith/config_tetralith_PMIP_AWIESM1.cfg

mpiexec UFEMISM_program /home/x_frare/utilities_UFEMISM/config-files_tetralith/config_tetralith_PMIP_INM-CM4-8.cfg

mpiexec UFEMISM_program /home/x_frare/utilities_UFEMISM/config-files_tetralith/config_tetralith_PMIP_MIROC-ES2L.cfg

mpiexec UFEMISM_program /home/x_frare/utilities_UFEMISM/config-files_tetralith/config_tetralith_PMIP_MPI-ESM1-2.cfg


