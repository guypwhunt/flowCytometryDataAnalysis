#!/bin/sh
#SBATCH --time=05:05:00
#SBATCH --partition=cpu
#SBATCH --mem 500G

module load r/4.1.1-gcc-9.4.0-withx-rmath-standalone-python-3.8.12

echo "BYE"