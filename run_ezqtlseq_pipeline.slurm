#!/bin/bash
# run snakemake EzQTLseq workflow

#SBATCH --job-name=EzQTLseq
#SBATCH --time="02-20:00:00"
## SBATCH --partition=workq
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --output=EzQTLseq_run.%N.%j.out

########################## on genotoul ###############################
## uncomment the 2 next lines for genotoul
#module load system/Python-3.6.3
#module load system/singularity-3.5.3
######################################################################

########################## on ifb-core ###############################
## uncomment the next line
#module load snakemake
######################################################################

# Full cleanup
#rm -fr *.pdf Bulks *.log logs *.fasta.* Parents .snakemake *.out

rm -fr .snakemake

export SINGULARITY_BINDPATH="/work2/project/gafl"

CONFIG=config.yaml
RULES=Snakefile.singularity.smk
CLUSTER="sbatch -J {cluster.jobname} -t {cluster.time}  -c {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb} --output {cluster.output}"
CLUSTER_CONFIG=cluster.json

MAX_JOBS=500

# log directory is mandatory (see $CLUSTER) else slurm jobs failed but not the master job
mkdir -p logs

# generate the dag file
snakemake --configfile $CONFIG -s $RULES --dag | dot -Tpdf > dag.pdf
snakemake --configfile $CONFIG -s $RULES --rulegraph | dot -Tpdf > dag_rules.pdf

# run the dry run (simulation)
snakemake --configfile $CONFIG -s $RULES -np -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER" >pipeline_dryrun.out

# uncomment the following line for the real run
#snakemake --configfile $CONFIG -s $RULES -p -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER"




# If latency problem add
# --latency-wait 60

#generate a final report
#snakemake --configfile $CONFIG -s $RULES --report smk_report.html
#ERROR doesn't work

exit 0

