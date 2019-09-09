#  HOW TO CALL THE SNAKEMAKE FILE:
#  snakemake -s proc10xG.snakemake.py -j 999 --cluster-config templates/cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}"
#  TODO make a left side and right side cluster and json for proc10xG and run supernova respectively
#  this file, templates/cluster.json, templates/keith.json,
#  to specify only one rule just add it after -s "-s proc10xG.snakemake.py kat_reads"

## ucdbioinfo_supernova_pipeline proc10xG.slurm
## runs the process_10xReads.py script from the proc10xG repo
## https://github.com/ucdavis-bioinformatics/proc10xG
## Assumes only a single pair of fastq (R1/R2) files under the fastqs folder

import os
import json
import subprocess

args = {}
sbatch_args = {}
configfile: "templates/keith.json"

###########################################################################
# CHECK IF SRUN OR SBATCH
###########################################################################
if config["__default__"]["running_locally"]=="True":
    # print ("Running Locally")
    args["running_locally"] = True
    args['cluster_threads'] = 48
    args['cluster_memory'] = 492

else:
    print ("Running on cluster")
    #print ("My SLURM_JOB_ID: %s" %(os.environ['SLURM_JOB_ID']))
    #args['cluster_threads'] = os.environ['SLURM_NTASKS']
    #print ("My SLURM_JOB_ID: %s" %(os.environ['SLURM_JOB_ID']))

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

# PROC10XG
args['proc10xg_out'] = args['basename'] + '/01-%s-%s_reads' % (args['id'], 'proc10xG')
args['proc10xg_outprefix'] = args['proc10xg_out'] + '/%s-%s' % (args['id'], 'proc10xG_reads')
args['fastq1_proc10xg_out'] = args['proc10xg_outprefix'] + '_R1_001.fastq.gz'
args['fastq2_proc10xg_out'] = args['proc10xg_outprefix'] + '_R2_001.fastq.gz'
args['log_out'] = args['proc10xg_outprefix'] + '.log'
args['proc10xPath'] = args['pipeline'] + '/%s' % ('proc10xG')

# KAT READS
args['kat_reads_out'] = args['basename'] + '/02-%s-%s' % (args['id'], 'kat_reads')
args['kat_reads_outprefix'] = args['kat_reads_out'] + '/%s-%s' % (args['id'], 'kat_reads')

# RUN SUPERNOVA
args['supernova_out'] = args['basename'] + '/01-%s-%s' % (args['id'], 'supernova_run')
args['supernova_id'] = '01-%s-%s' % (args['id'], 'supernova_run')
args['supernova_read_count'] = config["supernova"]["read_count"]
args['supernova_out_dir'] = args['supernova_out'] + '/' + 'outs'
args['supernova_seqout'] = args['basename'] + '/02-%s-%s' %(args['id'], 'supernova_outs')
args['supernova_out_prefix'] = args['supernova_seqout'] + '/%s-%s' %(args['id'], 'supernova_mkout')


# MKBWA
args['supernova_seqin1'] = args['supernova_seqout'] + '/%s-supernova_mkout-pseudohap2.1.fasta.gz' % args['id']
args['supernova_seqin2'] = args['supernova_seqout'] + '/%s-supernova_mkout-pseudohap2.2.fasta.gz' % args['id']

# KAT COMP and SECT
args['kat_compsect_out'] = args['basename'] + '/03-%s-%s' % (args['id'], 'assembly_eval')
args['kat_comp1'] = args['kat_compsect_out'] + '/%s-kat_eval-h1_vs_pe' % args['id']
args['kat_comp2'] = args['kat_compsect_out'] + '/%s-kat_eval-all_vs_pe' % args['id']
args['kat_sect'] = args['kat_compsect_out'] + '/%s-kat_eval-sect-h1_vs_pe' % args['id']

# MAP BARCODES
args['assembly_eval_outprefix'] = args['kat_compsect_out'] + '/%s-assembly_eval-bwa.bam'
args['assembly_eval_flagstat'] = args['kat_compsect_out'] + '/%s-assembly_eval-bwa.bam.flagstat'
args['assembly_eval_idxstats'] = args['kat_compsect_out'] + '/%s-assembly_eval-bwa.bam.idxstats'
args['assembly_eval_stats'] = args['kat_compsect_out'] + '/%s-assembly_eval-bwa.bam.stats'


###########################################################################
# MODULE LOADS
###########################################################################

if args['running_locally']=="False":
    import socket
    print("Running Locally")
	print (socket.gethostname())

shell.prefix("set -o pipefail; ")
shell.prefix("module load kat; module load anaconda2; module load bwa/0.7.16a; module load samtools/1.9; module load supernova/2.1.1;")
shell("module list")


###########################################################################
# KAT READS SBATCH
###########################################################################
args['cluster_time'] = config['kat_reads_sbatch']['main']['time']
args['cluster_account'] = config['kat_reads_sbatch']['main']['account']
args['cluster_partition'] = config['kat_reads_sbatch']['main']['partition']
args['cluster_nodes'] = config['kat_reads_sbatch']['main']['n']
print(json.dumps(args, indent=1))

###########################################################################
# RULES
###########################################################################


rule kat_sect:
    input:
        seqin_1 = args['supernova_seqin1'],
        proc10xin_1 = args['fastq1_proc10xg_out'],
        proc10xin_2 = args['fastq2_proc10xg_out']
    output:
        kat_comp1_out = args['kat_sect']
    run:
        arg_list = args['cluster_threads'], args['kat_comp1'], args['proc10xPath'], \
                   args['fastq1_proc10xg_out'], args['fastq2_proc10xg_out'], args['supernova_seqin1']
        if args['running_locally']:
             command = "kat comp -t%s -I10000000000 -H10000000000 -o %s <( %s -1 %s -2 %s ) <( gunzip -c %s)" % arg_list
        else:
                #TODO get more parmeters from the .json
                command = "sbatch "
        print(command)
        shell(command)


rule kat_comp2:
    input:
        seqin_1 = args['supernova_seqin1'],
        proc10xin_1 = args['fastq1_proc10xg_out'],
        proc10xin_2 = args['fastq2_proc10xg_out']
    output:
        kat_comp1_out = args['kat_comp2']
    run:
        arg_list = args['cluster_threads'], args['kat_comp2'], args['proc10xPath'], \
                   args['fastq1_proc10xg_out'], args['fastq2_proc10xg_out'], args['supernova_seqin1'], args['supernova_seqin2']
        if args['running_locally']:
            command = "kat sect -t%s -I10000000000 -H10000000000 -o %s <( %s -1 %s -2 %s ) <( gunzip -c %s %s)" % arg_list
        else:
            #TODO get more parameters from the .json
		    command = "sbatch"
        print(command)
        shell(command)


rule kat_comp1:
    input:
        seqin_1 = args['supernova_seqin1'],
        proc10xin_1 = args['fastq1_proc10xg_out'],
        proc10xin_2 = args['fastq2_proc10xg_out']

    output:
        kat_comp1_out = args['kat_comp1']
    run:
        arg_list = args['cluster_threads'], args['kat_comp1'], args['proc10xPath'], \
                   args['fastq1_proc10xg_out'], args['fastq2_proc10xg_out'], args['supernova_seqin1']
        if args['running_locally']:
            command = "kat comp -t%s -I10000000000 -H10000000000 -o %s <( %s -1 %s -2 %s ) <( gunzip -c %s)" % arg_list
        else:
            #TODO get more parmeters from the .json
		    command = "sbatch "
        print(command)
        shell(command)


rule map_barcodes:
    input:
        seqin_1 = args['supernova_seqin1'],
        proc10xin_1 = args['fastq1_proc10xg_out'],
        proc10xin_2 = args['fastq2_proc10xg_out']
    output:
        bam_out = args['assembly_eval_outprefix'],
        bam_flagstat = args['assembly_eval_flagstat'],
        bam_idxstats = args['assembly_eval_idxstats'],
        bam_stats = args['assembly_eval_stats']
    run:
        #TODO mapthreads?
        arg_list = 48, args['id'], args['id'], args['supernova_seqin1'], args['fastq1_proc10xg_out'], \
                   args['fastq2_proc10xg_out'], str(args['proc10xPath']) + '/samConcat2Tag.py', args['assembly_eval_outprefix']
        command_bwa = "bwa mem -t %s -C -R '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tDS:Paired' %s %s %s | python %s | samtools sort -m 768M --threads %s | samtools view -hb -o %s -" % arg_list
        command_index = "samtools index -@ %s %s" %(str(48), args['assembly_eval_outprefix'])
        command_flagstat = "samtools flagstat -@ %s %s > %st" %(48, args['assembly_eval_outprefix'], args['assembly_eval_flagstat'])
        command_view = "samtools view -b -q 30 -f 0x2 -F 0x904 %s | samtools idxstats - > %s" %(args['assembly_eval_outprefix'], args['assembly_eval_idxstats'])
        command_stats = "samtools stats -@ %s %s > %s" %(48, args['assembly_eval_outprefix'], args['assembly_eval_stats'])

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
    output:
        kat_reads_out = args['kat_reads_out']

    run:
        for kmer, output in zip(params.kmers, params.outputs):
            arg_list = args['cluster_partition'], args['cluster_nodes'], args['cluster_time'], output, kmer, 48, params.proc10xg_path, args['fastq1_proc10xg_out'], args['fastq2_proc10xg_out']
            if args['running_locally']:
                command = "kat hist -o %s -m %s -t %s <(%s -1 %s -2 %s)" % arg_list[3:]
            else:
                command = "sbatch -p %s -n %s -t %s --wrap='kat hist -o %s -m %s -t %s <(%s -1 %s -2 %s)'" % arg_list
            print(command)
            os.system("mkdir %s" %(output))
            os.system(command)


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
        #log_out = args['log_out'],
        out_dir = args['proc10xg_out'],
        fastq1_out = args['fastq1_proc10xg_out'],
        fastq2_out = args['fastq2_proc10xg_out']
    run:
        arg_list = args['proc10xPath'], args['fastq1'], args['fastq2'], args['proc10xg_outprefix'], args['log_out']
        command = "`python %s/process_10xReads.py -1 %s -2 %s -o %s -a 2> %s`" % arg_list
        print(command)
        shell(command)


rule mkbwaref:
    input:
        bwa_seq = args['supernova_seqin1']
    output:
        bwa_out = str(args['supernova_seqin1']) + '.bwt'
    run:
        command = "bwa index %s" % args['supernova_seqin1']
        print(command)
        shell(command)


rule mkoutput_supernova:
    input:
       in_dir = args['supernova_out_dir']
    output:
       seqout = args['supernova_seqout'],
       bwa_seq = args['supernova_seqin1']
    run:
       for outstyle, minsize in zip(config['supernova']['outstyle'], config['supernova']['minsize']):
           arg_list = args['cluster_partition'], args['cluster_nodes'], args['cluster_time'], input.in_dir, args['supernova_out_prefix'] + '-' + outstyle, outstyle, minsize
           if args['running_locally']:
               command = "supernova mkoutput --asmdir=%s/assembly --outprefix=%s --style=%s --minsize=%s --headers=full" % arg_list[3:]
           else:
               command = "sbatch -p %s -n %s -t %s --wrap='supernova mkoutput --asmdir=%s/assembly --outprefix=%s --style=%s --minsize=%s --headers=full'" % arg_list
           print(command)
           shell(command)

rule run_supernova:
    input:
        fastq1 = args['fastq1'],
        fastq2 = args['fastq2']
    params:
        supernova_out = args['supernova_out'],
        read_count = args['supernova_read_count'],
        fastqs = args['fastqs']
        #TODO check logic of this add output
        #localcores = args['cluster_threads']
    output:
        out_dir = args['supernova_out_dir']
    run:
        arg_list = args['supernova_id'], args['supernova_read_count'], args['fastqs'], 48
        command = "supernova run --id=%s --maxreads=%s --fastqs=%s --localcores=%s" % arg_list
        print(command)
        shell(command)

rule Illumina_10x:
    output:
        fastq1 = args['fastq1'],
        fastq2 = args['fastq2']

# only useful if running locally?
rule all:
    input:
        rules.kat_reads.output,
        rules.proc10xG.output,
        rules.run_supernova.output,
        rules.mkoutput_supernova.output,
        rules.kat_comp1.output,
        rules.kat_comp2.output,
        rules.kat_sect.output,
        rules.mkbwaref.output,
        rules.map_barcodes.output,

# running_locally = "False", sbatch included in snakemake call and proc10xg runs on that job and call sbatch for kat reads (3x) (use leftside.json for proc10xG call)
rule left_side:
    input:
        rules.proc10xG.output,
        rules.kat_reads.output

# running_locally = "False", sbatch included in snakemake call and run supernova runs on that job and call sbatch for mkoutput (3x) (use rightside.json for run_supernova call)
rule right_side:
    input:
        rules.run_supernova.output,
        rules.mkoutput_supernova.output


# running_locally = "False", no sbatch included in call or included (doesnt matter) all portions submitted as an sbatch
rule bottom:
    input:
        rules.kat_comp1.output,
        rules.kat_comp2.output,
        rules.kat_sect.output,
        rules.mkbwaref.output,
        rules.map_barcodes.output

# or the user can just use the cluster.json and submit each job and the sbatch specific there individually