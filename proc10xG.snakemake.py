#  HOW TO CALL THE SNAKEMAKE FILE:
#   snakemake -s proc10xG.snakemake.py -j 999 --cluster-config templates/cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}"


## ucdbioinfo_supernova_pipeline proc10xG.slurm
## runs the process_10xReads.py script from the proc10xG repo
## https://github.com/ucdavis-bioinformatics/proc10xG
## Assumes only a single pair of fastq (R1/R2) files under the fastqs folder
import os
import json

args = {}
sbatch_args = {}
configfile: "templates/keith.json"
#sbatchfile: "templates/cluster.json"

###########################################################################
# INPUT PARAMETERS
###########################################################################
args['kmers'] = config['kat_reads']['kmers']
args['pipeline'] = config['pipeline']['basepath']
args['basename'] = config['project']['basename']
args['id'] = config['project']['id']
args['fastqs'] = args['basename'] + '/' + config['project']['fastqs']


files = os.listdir(args['fastqs'])
for file in files:
    if "R1_001.fastq.gz" in file:
        args['fastq1'] = args['fastqs'] + '/' + file
    if "R2_001.fastq.gz" in file:
        args['fastq2'] = args['fastqs'] + '/' + file



###########################################################################
# OUTPUT PARAMETERS
###########################################################################
# TODO: SHOULD proc10xg be coming from the job name
# PROC10XG
args['proc10xg_out'] = args['basename'] + '/01-%s-%s_reads' % (args['id'], 'proc10xG')
args['proc10xg_outprefix'] = args['proc10xg_out'] + '/%s-%s' % (args['id'], 'proc10xG')
args['fastq1_proc10xg_out'] = args['proc10xg_outprefix'] + '_R1_001.fastq.gz'
args['fastq2_proc10xg_out'] = args['proc10xg_outprefix'] + '_R2_001.fastq.gz'


args['log_out'] = args['proc10xg_outprefix'] + '.log'
args['proc10xPath'] = args['pipeline'] + '/%s' % ('proc10xG')

# KAT READS
args['kat_reads_out'] = args['basename'] + '/02-%s-%s' % (args['id'], 'kat_reads')
args['kat_reads_outprefix'] = args['kat_reads_out'] + '/%s-%s' % (args['id'], 'kat_reads')



###########################################################################
# MODULE LOADS
###########################################################################
import socket
print (socket.gethostname())
shell.prefix("module load kat; module load anaconda2;")
print(json.dumps(args, indent=1))


###########################################################################
# KAT READS SBATCH
###########################################################################

args['cluster_time'] = config['kat_reads_sbatch']['main']['time']
args['cluster_account'] = config['kat_reads_sbatch']['main']['account']
args['cluster_partition'] = config['kat_reads_sbatch']['main']['partition']
args['cluster_nodes'] = config['kat_reads_sbatch']['main']['n']


print(json.dumps(args, indent=1))


rule kat_reads:
    input:
        proc10xg_out = args['log_out'],
        fastq1 = args['fastq1_proc10xg_out'],
        fastq2 = args['fastq2_proc10xg_out']

    params:
        proc10xg_outprefix = args['proc10xg_outprefix'],
        proc10xg_out = args['proc10xg_out'],
        proc10xg_path = args['proc10xPath'],
        kat_reads_out = args['kat_reads_out'],
        kat_reads_outprefix = args['kat_reads_outprefix'],
        log_out = args['log_out'],
        kmers = args['kmers'],
        outputs = expand(args['kat_reads_outprefix'] + '-' + '{kmer}', kmer = args['kmers'])

    run:
        import subprocess
        shell("module list")
        for kmer, output in zip(params.kmers, params.outputs):
            command = "sbatch -p %s -n %s -t %s ./kat_reads_call.sh %s %s %s %s/filter_10xReads.py %s %s" %(args['cluster_partition'], args['cluster_nodes'], args['cluster_time'], output, kmer, 48, params.proc10xg_path, args['fastq1_proc10xg_out'], args['fastq2_proc10xg_out'])
            print (command)
            os.system("mkdir %s" %(output))
            shell(command)


rule proc10xG:
    input:
        fastq1 = args['fastq1'],
        fastq2 = args['fastq2']
    params:
        proc10xg_outprefix = args['proc10xg_outprefix'],
        proc10xg_out = args['proc10xg_out'],
        proc10xg_path = args['proc10xPath'],
        log_out = args['log_out']
    output:
        log_out = args['log_out'],
        #fastq1_out = args['fastq1_proc10xg_out'],
        #fastq2_out = args['fastq2_proc10xg_out']
    shell:
        "python {params.proc10xg_path}/process_10xReads.py -1 {input.fastq1} -2 {input.fastq2} -o {params.proc10xg_outprefix} -a 2> {output}"