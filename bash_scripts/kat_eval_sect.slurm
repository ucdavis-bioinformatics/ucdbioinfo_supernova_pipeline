#!/bin/bash
#
## ucdbioinfo_supernova_pipeline kat_reads.slurm
## runs kat on the proc10xG reads to profile kmer distribution of reads
## https://github.com/TGAC/KAT
## Assumes proc10xG.slurm has completed successfully
#
#SBATCH --job-name=kat-eval-sect # Job name
#SBATCH --nodes=1
#SBATCH --time=1-0
#SBATCH --ntasks=32 # Number of cores
#SBATCH --mem=200000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/kat-eval-sect-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/kat-eval-sect-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=32; else THREADS=${SLURM_NTASKS}; fi
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=200; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

pipeline=$(eval jq .pipeline.basepath $1 | sed 's/^"\(.*\)"$/\1/')
basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )

#input
proc10xG_out="${basepath}/01-${id}-proc10xG_reads"
mkdir -p ${proc10xG_out}
proc10xG_outprefix=${proc10xG_out}/${id}-proc10xG_reads

supernova_seqin1="${basepath}/02-${id}-supernova_outs/${id}-supernova_mkout-pseudohap2.1.fasta.gz"
supernova_seqin2="${basepath}/02-${id}-supernova_outs/${id}-supernova_mkout-pseudohap2.2.fasta.gz"

#output
assembly_eval_out=${basepath}/03-${id}-assembly_eval
mkdir -p ${assembly_eval_out}
assembly_eval_outprefix=${assembly_eval_out}/${id}-kat_eval

#software
module load anaconda3
module load kat
module load anaconda2
proc10xPath="${pipeline}/proc10xG"

call="kat sect \
  -t${THREADS} \
  -H10000000000 \
  -o ${assembly_eval_outprefix}-sect-h1_vs_pe \
  <( gunzip -c ${supernova_seqin1} ) \
  <( ${proc10xPath}/filter_10xReads.py -1 ${proc10xG_outprefix}_R1_001.fastq.gz -2 ${proc10xG_outprefix}_R2_001.fastq.gz )"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

exit
