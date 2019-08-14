##!/bin/bash

## ucdbioinfo_supernova_pipeline proc10xG.slurm
## runs the process_10xReads.py script from the proc10xG repo
## https://github.com/ucdavis-bioinformatics/proc10xG
## Assumes only a single pair of fastq (R1/R2) files under the fastqs folder




# bashCommand = "cwm --rdf test.rdf --ntriples > test.nt"
# import subprocess
# process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)


# read config info into this namespace
configfile: "../templates/keith.json"
print (configfile['samples'])

rule proc10xG:
    shell:
        "./../bash_scripts/proc10xG.slurm {configfile}"
