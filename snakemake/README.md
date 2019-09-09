# Running UCD Bioinformatics Supernova Pipeline with Snakemake


<img src="dag.pdf" alt="supernova" width="800px"/>

##  HOW TO CALL THE SNAKEMAKE FILE (GENERAL):
1.  `snakemake -s snakemake.py -j 999 --cluster-config templates/cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}"`
2.  Files needed: `snakemake.py`, `templates/cluster.json`, `templates/keith.json`
3.  Be sure to specify `running_locally` in the `template/keith.json` file. 
    - If `running_locally` == 'True':
        + `templates/cluster.json` is not needed if running locally only.
        + don't need to worry about the sbatch parameters in the `templates/keith.json` file
        + Follow rules under **RUNNING LOCALLY** below.
    - Else if `running_locally` == 'False':
        + be sure to set all sbatch parameters in the `templates/keith.json` file
        + Follow rules under **RUNNING ON CLUSTER** below.
        
## RUNNING ON CLUSTER
***
1.  Left Side
    - running_locally = "False", sbatch included in snakemake call and proc10xg runs on that job and call sbatch for kat reads (3x) (use leftside.json for proc10xG call)
    - `snakemake -s snakemake.py -j 999 --cluster-config templates/cluster_left.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}" left_side`
2.  Right Side
    - running_locally = "False", sbatch included in snakemake call and run supernova runs on that job and call sbatch for mkoutput (3x) (use rightside.json for run_supernova call)
    - `snakemake -s snakemake.py -j 999 --cluster-config templates/cluster_left.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}" right_side`
3.  Bottom 
    - running_locally = "False", no sbatch included in call or included (doesnt matter) all portions submitted as an sbatch
    - `snakemake -s snakemake.py -j 999 --cluster-config templates/cluster_left.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}" bottom`
***


<img src="dag_circled.pdf" alt="supernova_circled" width="800px"/>


## RUNNING LOCALLY

***
1. All
    - running_locally = "True"
    - `snakemake -s snakemake.py -j 999 --cluster-config templates/cluster_left.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}" left_side`
***


3.  to specify only one rule just add it after -s "-s proc10xG.snakemake.py kat_reads"



##  TODO 
- make a left side and right side cluster and json for proc10xG and run supernova respectively
