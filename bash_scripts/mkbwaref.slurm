#!/bin/bash
#
#SBATCH --job-name=mkbwaref # Job name
#SBATCH --nodes=1
#SBATCH --time=1-0
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=12000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/mkbwaref-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/mkbwaref-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=1; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=12; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )

#input
supernova_seqin="${basepath}/02-${id}-supernova_outs/${id}-supernova_outs-pseudohap2.1.fasta.gz"

#software
module load bwa/0.7.16a

#application call
call="bwa index ${supernova_seqin}"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
