#!/bin/bash
#SBATCH --job-name=make_tests_job                       # Job name
#SBATCH --output=make_tests_output.txt     # File to save the output
#SBATCH --error=make_tests_error.txt       # File to save error messages
#SBATCH --ntasks=1                                      # Number of tasks
#SBATCH --cpus-per-task=46                               # Number of CPUs required (match with -j8)

# Run the make command
# make tests -j8 > files/output/make_tests_output.txt 2>&1
make tests -j8 > make_tests_output.txt 2>&1