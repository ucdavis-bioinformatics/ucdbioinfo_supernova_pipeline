#!/bin/bash
#
## ucdbioinfo_supernova_pipeline mkoutput_supernova.slurm
## produces fastq files from supernova run output
## https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome
## Assumes supernova run has completed sucessfully
#
#SBATCH --array=0-2
#SBATCH --job-name=mkoutput_supernova # Job name
#SBATCH --nodes=1
#SBATCH --time=1-0
#SBATCH --ntasks=8 # Number of cores
#SBATCH --mem=120000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/mkoutput_supernova-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/mkoutput_supernova-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=48; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=490; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )

## output style options are raw|megabubbles|pseudohap|pseudohap2
## minsize default is 500bp for raw, 1000bp for all others
## headers default is short
outstyle=(raw megabubbles pseudohap)
minsize=(500 1000 1000)
if [ -z ${SLURM_ARRAY_TASK_ID} ]; then echo SLURM_ARRAY_TASK_ID=0; fi

#input
amsdir="${basepath}/01-${id}-supernova_run/outs/assembly"

#output
supernova_seqout="${basepath}/02-${id}-supernova_outs"
mkdir -p ${supernova_seqout}
supernova_seqoutprefix="${supernova_seqout}/${id}-supernova_outs"

#software
module load supernova/2.1.1

#application call
cp ${amsdir}/report.txt ${supernova_seqout}/report.txt
cp ${amsdir}/summary.csv ${supernova_seqout}/summary.csv

call="supernova mkoutput \
--asmdir=${amsdir} \
--outprefix=${out_prefix}-${outstyle[SLURM_ARRAY_TASK_ID]} \
--style=${outstyle[SLURM_ARRAY_TASK_ID]} \
--minsize=${minsize[SLURM_ARRAY_TASK_ID]} \
--headers=full"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
