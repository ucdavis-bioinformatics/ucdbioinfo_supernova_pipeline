#!/bin/bash
#
## ucdbioinfo_supernova_pipeline kat_reads.slurm
## runs kat on the proc10xG reads to profile kmer distribution of reads
## https://github.com/TGAC/KAT
## Assumes proc10xG.slurm has completed successfully
#
#SBATCH --array=0-2
#SBATCH --job-name=kat # Job name
#SBATCH --nodes=1
#SBATCH --time=10-9
#SBATCH --ntasks=48 # Number of cores
#SBATCH --mem=500000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/kat-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/kat-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=48; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=492; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )

kmers=(21 27 33)
if [ -z ${SLURM_ARRAY_TASK_ID} ]; then echo SLURM_ARRAY_TASK_ID=0; fi

#input
proc10xG_out="${basepath}/01-${id}-proc10xG_reads"
proc10xG_outprefix=${proc10xG_out}/${id}-proc10xG_reads

#output
kat_reads_out=${basepath}/02-${id}-kat_reads
mkdir -p ${kat_reads_out}
kat_reads_outprefix=${kat_reads_out}/${id}-kat_reads

#software
module load kat
module load anaconda2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="/share/biocore/projects/Tao_R_Kyoto_Japan/persimmon_10xGenome/ucdbioinfo_supernova_pipeline/bash_scripts"
proc10xPath="${DIR}/../proc10xG"

#application call
call="kat hist \
  -o ${kat_reads_outprefix}-k${kmers[SLURM_ARRAY_TASK_ID]} \
  -m ${kmers[SLURM_ARRAY_TASK_ID]} \
  -t $((THREADS - 1)) \
  <( ${proc10xPath}/filter_10xReads.py -1 ${proc10xG_outprefix}_R1_001.fastq.gz -2 ${proc10xG_outprefix}_R2_001.fastq.gz 2> ${kat_reads_outprefix}-filter-k${kmers[SLURM_ARRAY_TASK_ID]}.log)"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

exit
