#!/bin/bash

#PBS -N bmi_pd
#PBS -o job_reports/bmi_pd-output
#PBS -e job_reports/bmi_pd-error
#PBS -l walltime=4:00:00
#PBS -t 1-10
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
outdir="${HOME}/repo/mr_frailty/inst/analysis/bmi_pd/scratch4"

R --no-save --args ${i} ${splits} ${outdir} < ${HOME}/repo/mr_frailty/inst/analysis/bmi_pd/scripts/model4.R

