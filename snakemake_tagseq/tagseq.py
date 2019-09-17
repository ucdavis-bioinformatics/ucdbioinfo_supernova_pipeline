## ucdbioinfo_tagseq_pipeline
## Assumes only a single fastq (R1) for each sample directory in the fastqs folder

import os
import json
import glob

args = {}
sbatch_args = {}

# run module load snakemake
# source activate snakemake

# Check if the user is running all the samples as the "master command"
# snakemake -s tagseq.py master_rule runs all samples based on sbatch calls or locally in a sequence.

if 'master' in config.keys() and config['master'] == True:
    #TODO TO CALL FOR ALL SAMPLES FROM SOME FILE samples.txt
    SAMPLES = ['SampleAC3', 'SampleAC2', 'SampleAC1', 'SampleAC4', 'SampleAD3', 'SampleAD2', 'SampleAD1', 'SampleAD4']

# If not we are running one sample (likely from the rule master_rule) or sample=Sample1,Sample2,Sample3 as specified by the user
else:
    SAMPLES = str(config['sample']).split(',') # passed via the command --config sample="SampleAC3"

print ("SAMPLES:", SAMPLES)
#TODO take this in from the CLI argument
configfile: "templates/tagseq.json"

###########################################################################
# CHECK IF RUNNING LOCALLY OR ON SLURM CLUSTER
###########################################################################
if config["__default__"]["running_locally"]=="True":
    args["running_locally"] = True
    print ("Running Locally")

else:
    args["running_locally"] = False
    print ("Running on Cluster")

###########################################################################
# CORE SETUP
###########################################################################
args['basename'] = config['project']['basename']
args['id'] = config['project']['id']
args['fastqs'] = args['basename'] + '/' + config['project']['fastqs']
args['htsout'] = args['basename'] + '/' + '01-HTS_Preproc'
args['starout'] = args['basename'] + '/' + '02-STAR_alignment'

###########################################################################
# MODULES
###########################################################################
if args['running_locally']== False:
    import socket
    print ("SOCKET:", socket.gethostname())

shell.prefix("set -o pipefail; ")
shell.prefix("module load star/2.7.0e; module load htstream/1.1.0;")
shell("module list")

print(json.dumps(args, indent=1))
###########################################################################

rule master_rule:
    run:
        for sample in SAMPLES:
            if args['running_locally'] == False:
                command = f"sbatch \
                   --job-name={config['hts_star']['job-name'] + sample} \
                   --ntasks={config['hts_star']['ntasks']} \
                   --nodes={config['hts_star']['n']} \
                   --partition={config['hts_star']['partition']} \
                   --time={config['hts_star']['time']}  \
                   --mem={config['hts_star']['mem']} \
                   --output={config['hts_star']['output'].replace('.out', sample + '.out')} \
                   --error={config['hts_star']['error'].replace('.err', sample + '.err')} \
                   --mail-type={config['hts_star']['mail-type']} \
                   --mail-user={config['hts_star']['mail-user']} \
                   --wrap='snakemake -s tagseq.py all --config sample={sample}'"
            else:
                command = f"snakemake -s tagseq.py all --config sample={sample}"
            print(command)
            shell(command)



rule htspreproc:
    input:
        # glob.glob('%s/{sample}/{sample}*R1*' % args['fastqs'])
        '%s/{sample}/{sample}_L3_R1.fastq.gz' % args['fastqs']
    output:
        hts_log = '%s/{sample}/{sample}_htsStats.log' % args['htsout'],
        fastq = '%s/{sample}/{sample}_R1.fastq.gz' % args['htsout']
    params:
        # TODO take this from the JSON file
        ref = '/share/biocore/keith/workshop/rnaseq_examples/References/human_rrna.fasta'
    run:
        # shell("mkdir %s/{wildcards.sample}" % args['htsout'])
        stats1 = "hts_Stats -L {output.hts_log} -U {input}"
        seq_screen1 = "hts_SeqScreener -A -L {output.hts_log}"
        seq_screen2 = "hts_SeqScreener -s {params.ref} -r -A -L {output.hts_log}"
        adapter_trim = "hts_AdapterTrimmer -n -A -L {output.hts_log}"
        cut_trim1 = "hts_CutTrim -n -a 20 -A -L {output.hts_log}"
        qwindow_trim = "hts_QWindowTrim -n -A -L {output.hts_log}"
        ntrimmer = "hts_NTrimmer -n -A -L {output.hts_log}"
        cut_trim2 = "hts_CutTrim -n -m 50 -A -L {output.hts_log}"
        stats2 = "hts_Stats -A -L {output.hts_log} -f %s/{wildcards.sample}/{wildcards.sample}" % (args['htsout'])
        master_list = [stats1, seq_screen1, seq_screen2, adapter_trim, cut_trim1, qwindow_trim, ntrimmer, cut_trim2, stats2]
        command = ' | '.join(master_list)
        print(command)
        shell(command)


rule star:
    input:
        fastq = '%s/{sample}/{sample}_htsStats.log' % args['htsout']
    output:
        bam = '%s/{sample}/{sample}_Aligned.sortedByCoord.out.bam' % args['starout'],
        log = '%s/{sample}/{sample}_Log.out' % args['starout'],
        prog = '%s/{sample}/{sample}_Log.progress.out' % args['starout'],
        stderr = '%s/{sample}/{sample}-STAR.stderr' % args['starout'],
        stdout = '%s/{sample}/{sample}-STAR.stdout' % args['starout'],
        temp = '%s/{sample}/{sample}__STARtmp' % args['starout']
    params:
        ref = "/share/biocore/keith/workshop/rnaseq_examples/References/star.overlap100.gencode.v31",
        prefix = "%s/{wildcards.sample}/{wildcards.sample}_" % (args['starout'])
    run:
        command = "STAR --runThreadN 8 --genomeDir {params.ref} \
                         --outSAMtype BAM SortedByCoordinate \
                         --readFilesCommand zcat \
                         --readFilesIn {input.fastq} \
                         --quantMode GeneCounts \
                         --outFileNamePrefix {params.prefix} \
                         > {output.stdout} 2> {output.stderr}"
        print(command)
        shell(command)

rule all:
    input:
        expand('%s/{sample}/{sample}_Aligned.sortedByCoord.out.bam'  % args['starout'], sample=SAMPLES),
        expand('%s/{sample}/{sample}_Log.out'  % args['starout'], sample=SAMPLES)
