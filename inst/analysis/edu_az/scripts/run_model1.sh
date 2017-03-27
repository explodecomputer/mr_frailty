#!/bin/bash

#PBS -N edu_az
#PBS -o job_reports/edu_az-output
#PBS -e job_reports/edu_az-error
#PBS -l walltime=4:00:00
#PBS -t 1-100
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
splits=100
outdir="${HOME}/repo/mr_frailty/inst/analysis/edu_az/scratch1"

R --no-save --args ${i} ${splits} ${outdir} < ${HOME}/repo/mr_frailty/inst/analysis/edu_az/scripts/model1.R

