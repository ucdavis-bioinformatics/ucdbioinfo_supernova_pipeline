#!/bin/bash
#
## ucdbioinfo_supernova_pipeline proc10xG.slurm
## runs the process_10xReads.py script from the proc10xG repo
## https://github.com/ucdavis-bioinformatics/proc10xG
## Assumes only a single pair of fastq (R1/R2) files under the fastqs folder
#
#SBATCH --job-name=proc10xG # Job name
#SBATCH --nodes=1
#SBATCH --time=5-0
#SBATCH --ntasks=3 # Number of cores
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/proc10xG-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/proc10xG-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=3; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=2; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )
fastqs=$(eval jq .project.fastqs $1 | sed 's/^"\(.*\)"$/\1/')

#input
fastq1=$(eval ls -d -1 ${fastqs}/* | grep -x ".*R1_001.fastq.gz")
fastq2=$(eval ls -d -1 ${fastqs}/* | grep -x ".*R2_001.fastq.gz")

#output
proc10xG_out="${basepath}/01-${id}-proc10xG_reads"
mkdir -p ${proc10xG_out}
proc10xG_outprefix=${proc10xG_out}/${id}-proc10xG_reads
log_out=${proc10xG_outprefix}.log

#software
module load anaconda2
scriptpath="$( cd "$(dirname "$0")" ; pwd -P )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
scriptpath="/share/biocore/projects/Silva_Filho_M_USP/D_saccharalis_10X/ucdbioinfo_supernova_pipeline"
proc10xPath="${scriptpath}/proc10xG"

#application call
call="${proc10xPath}/process_10xReads.py \
 -1 ${fastq1} \
 -2 ${fastq2} \
 -o ${proc10xG_outprefix} -a 2> ${log_out}"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

exit
