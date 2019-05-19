#!/bin/bash
#
## ucdbioinfo_supernova_pipeline run_supernova.slurm
## script_runs supernova based on parameters in a setup json file
## https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome
#
#SBATCH --job-name=run_supernova # Job name
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --time=30-0
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/run_supernova-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/run_supernova-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent


start=`date +%s`

export MRO_DISK_SPACE_CHECK=disable 

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=$( nproc ); else THREADS=${SLURM_NTASKS}; fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

#input
basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )
fastqs=$(eval jq .project.fastqs $1 | sed 's/^"\(.*\)"$/\1/')

#output
supernova_out="01-${id}-supernova_run"
read_count=$(eval jq .supernova.read_count $1)

#check if supernova has run before and a lock file exists
if [ -f $supernova_out/_lock ] ; then
    rm $supernova_out/_lock
fi
#software
module load supernova/2.1.1

#application call
call="MRO_DISK_SPACE_CHECK=disable supernova run \
      --id=${supernova_out} \
      --maxreads=${read_count} \
      --fastqs=${fastqs} \
      --localcores=${THREADS}"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
