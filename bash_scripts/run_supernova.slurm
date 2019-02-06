#!/bin/bash
#
## ucdbioinfo_supernova_pipeline run_supernova.slurm
## script_runs supernova
#
#SBATCH --job-name=run_supernova # Job name
#SBATCH --nodes=1
#SBATCH --time=12-0
#SBATCH --ntasks=48 # Number of cores
#SBATCH --mem=508000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/run_supernova-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/run_supernova-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=48; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=496; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )
fastqs=$(eval jq .project.fastqs $1 | sed 's/^"\(.*\)"$/\1/')

supernova_out="01-${id}-supernova_run"
read_count=$(eval jq .supernova.read_count $1)

module load supernova/2.1.1

call="supernova run \
      --id=${supernova_out} \
      --maxreads=${read_count} \
      --fastqs=${fastqs} \
      --localcores=${THREADS} \
      --localmem=${MEM}"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
