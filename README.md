 EzQTLseq workflow v1.0

GAFL 2020


1) create symbolic link in ref_genome
ln -s <path>/<Mygenome.fasta> ref_genome/

2) Edit config.yaml
     - adapt the parameters
     - add your samples

3) edit run_ezqtlseq_pipeline.slurm for Genotoul or ifb-core

4) run the pipeline with:
   sbatch run_ezqtlseq_pipeline.slurm




