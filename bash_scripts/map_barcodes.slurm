#!/bin/bash
#
#SBATCH --job-name=map_barcodes # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=20 # Number of cores
#SBATCH --mem=64000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=3-00
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurm_out/mapbarcodes-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/mapbarcodes-%N-%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=20; else THREADS=${SLURM_NTASKS}; fi
MAPTHREADS=$(expr ${THREADS} - 6)
SORTTHREADS=$(expr ${THREADS} - ${MAPTHREADS})
if [ -z ${SLURM_MEM_PER_NODE} ]; then MEM=64; else  MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024); fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

pipeline=$(eval jq .pipeline.basepath $1 | sed 's/^"\(.*\)"$/\1/')
basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )

#input
mapfasta="${basepath}/02-${id}-supernova_outs/${id}-supernova_mkout-pseudohap2.1.fasta.gz"
proc10xG_out="${basepath}/01-${id}-proc10xG_reads"
proc10xG_outprefix=${proc10xG_out}/${id}-proc10xG_reads

#output
assembly_eval_out=${basepath}/03-${id}-assembly_eval
mkdir -p ${assembly_eval_out}
assembly_eval_outprefix=${assembly_eval_out}/${id}-assembly_eval-bwa.bam

#software
module load anaconda2
module load bwa/0.7.16a
module load samtools/1.9

proc10xPath="${pipeline}/proc10xG"

#application call
call="bwa mem -t ${MAPTHREADS} -C \
 -R '@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA\tDS:Paired' \
 ${mapfasta} ${proc10xG_outprefix}_R1_001.fastq.gz ${proc10xG_outprefix}_R2_001.fastq.gz \
 | python ${proc10xPath}/samConcat2Tag.py \
 | samtools sort -m 768M --threads ${SORTTHREADS} \
 | samtools view -hb -o ${assembly_eval_outprefix} -"

echo $call
#eval $call

call="samtools index -@ ${THREADS} ${assembly_eval_outprefix}"
echo $call
#eval $call

call="samtools flagstat -@ ${THREADS} ${assembly_eval_outprefix} > ${assembly_eval_outprefix}.flagstat"
echo $call
#eval $call

#Only produce idxstats for primary mapped reads with q > 30  and in proper pair
call="samtools view -b -q 30 -f 0x2 -F 0x904 ${assembly_eval_outprefix} | samtools idxstats - > ${assembly_eval_outprefix}.idxstats"
echo $call
#eval $call

call="samtools stats -@ ${THREADS} ${assembly_eval_outprefix} > ${assembly_eval_outprefix}.stats"
echo $call
#eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
