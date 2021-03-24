from os.path import join, exists
import os, functools
from snakemake.utils import min_version

min_version("5.3.0")

WORKDIR = os.getcwd()


# this means that we will try to link from the ref directory rather than creating directly the indices
ruleorder: link_ref_files > fasta_to_sa > fasta_to_fai > fasta_to_dict

report: "report/workflow.rst"
reportA: "report/bam.rst"

###### SUMMARY RULES #############################
rule all:
    input:
        #rules.fastqc_reads_all.output,
        #rules.multiqc_reads_all.output,
        "Bulks/qtlsT.csv"

rule parents:
    input:
        "Parents/p2.on.p1.fb.vcf.gz.tbi",
        "reports/Parents/p1.on.ref.1.stats.html",
        "reports/Parents/p1.on.ref.stats.html",
        "Parents/p2.fasta",
        "Parents/p1.on.p1.fb.table",
        "Parents/p2.on.p1.fb.table",
        "reports/Parents/p1.on.ref.1.stats.html",
        "reports/Parents/p1.on.ref.stats.html"

###### END SUMMARY RULES ############################

###### FUNCTIONS FOR READS ##########################

# wildcards.X = R|S|p1|p2
def get_R2_reads_files(wildcards):
    if config["dedup_reads"]:
        return expand(config[wildcards.X + "_reads_dir"]+"/{read_file2}", read_file2 = config[wildcards.X + "_reads2"])
    return config[wildcards.X + "_reads_dir"]

## wildcards.X = R|S|p1|p2
#def get_reads_files(wildcards):
#    return expand(config[wildcards.X + "_reads_dir"]+"/{read_file}", read_file = config[wildcards.X + "_reads"])

# wildcards.X = R|S|p1|p2
def get_reads1_files(wildcards):
    return expand(config[wildcards.X + "_reads_dir"]+"/{read_file1}", read_file1 = config[wildcards.X + "_reads1"])
# wildcards.X = R|S|p1|p2
def get_reads2_files(wildcards):
    return expand(config[wildcards.X + "_reads_dir"]+"/{read_file2}", read_file2 = config[wildcards.X + "_reads2"])


###### END FUNCTIONS FOR READS #############################

###### REF SEQ RULES    #########################################

rule link_ref_files:
    input:
        conf_ref = os.path.join(WORKDIR, config["ref_seq_path"]+"{suf}"),
    output:
        "Parents/ref.fasta{suf}"
    shell:
        "ln -s {input.conf_ref} {output}"

# needed because {suf} cannot be empty in the previous rule
rule link_ref_fasta:
    input:
        conf_ref = os.path.join(WORKDIR, config["ref_seq_path"]),
    output:
        "Parents/ref.fasta"
    shell:
        "ln -s {input.conf_ref} {output}"
###### END REF SEQ RULES    #########################################

include: "./parents.singularity.smk"
include: "./bulks.singularity.smk"
include: "./spet.singularity.smk"

###### MAPPING RULES    #############################################

# map reads -> remove non primary reads
#-> remove reads that contain an S in the cigar -> sort the reads
rule map_with_bwa:
    input:
        ref = "Parents/{ref}.fasta",
        ref_amb = "Parents/{ref}.fasta.amb",
        ref_ann = "Parents/{ref}.fasta.ann",
        ref_bwt = "Parents/{ref}.fasta.bwt",
        ref_pac = "Parents/{ref}.fasta.pac",
        ref_sa = "Parents/{ref}.fasta.sa",
        reads1 = get_reads1_files
    output:
        "{dir}/{X}.on.{ref}.1.sam"
    log:
        "logs/{dir}/{X}.on.{ref}.1.sam.log"
    wildcard_constraints:
        X="p1|p2|R|S",
        ref="ref|p1",
        dir="Parents|Bulks"
    params:
        reads2 = get_reads2_files,
        rg = r"@RG\tID:{X}.on.{ref}\tSM:{X}.on.{ref}\tPL:illumina",
        threads = 10
    shell:
        """
        r2='{params.reads2}'
        if [ -z "$r2" ] || [[ $r2 =~ "None" ]]
        then
           /work2/project/gafl/tools/containers/bwa_0.7.17.sif mem -R '{params.rg}' -t {params.threads} -o {output} {input.ref} <(cat {input.reads1}) > {log} 2>&1
           exit 0
        fi
        /work2/project/gafl/tools/containers/bwa_0.7.17.sif mem -R '{params.rg}' -t {params.threads} -o {output} {input.ref} <(cat {input.reads1}) <(cat {params.reads2}) > {log} 2>&1
        """

rule filter1:
    input:
        "{dir}/{X}.on.{ref}.1.sam"
    output:
        "{dir}/{X}.on.{ref}.2.sam"
    log:
        "logs/{dir}/{X}.on.{ref}.2.sam.log"
    wildcard_constraints:
        X="p1|p2|R|S",
        ref="ref|p1",
        dir="Parents|Bulks"
    params:
        rg = r"@RG\tID:{X}.on.{ref}\tSM:{X}.on.{ref}\tPL:illumina",
        threads = 1
    shell:
        "/work2/project/gafl/tools/containers/samtools_v1.9.sif view -h -F 2048 {input} > {output} 2> {log}"

rule filter2:
    input:
        "{dir}/{X}.on.{ref}.2.sam"
    output:
        "{dir}/{X}.on.{ref}.3.bam"
    log:
        "logs/{dir}/{X}.on.{ref}.3.bam.log"
    wildcard_constraints:
        X="p1|p2|R|S",
        ref="ref|p1",
        dir="Parents|Bulks"
    params:
        rg = r"@RG\tID:{X}.on.{ref}\tSM:{X}.on.{ref}\tPL:illumina",
        threads = 1
    shell:
        "singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif awk '($0 ~ /^@/) || ($6 !~ /S/)' {input} | /work2/project/gafl/tools/containers/samtools_v1.9.sif sort -O bam -o {output} > {log} 2>&1"

# piping the result of cat into the command won't work for some reason.
#I have to create a temp file with the cat output and delete it when done.
rule concat_R2:
    input:
        n6_reads = get_R2_reads_files
    output:
        temp("{dir}/{X}.R2.fastq.gz")
    wildcard_constraints:
        X="p1|p2|R|S",
        dir="Parents|Bulks"
    params:
        dedup_reads = config["dedup_reads"],
    shell: """
        if [ {params.dedup_reads} == False ]; then
            touch {output}
        else
            cat {input.n6_reads} > {output}
        fi
        """

#remove duplicates using N6 information
rule filter3:
    input:
        bam = "{dir}/{X}.on.{ref}.3.bam",
        n6_reads = "{dir}/{X}.R2.fastq.gz"
    output:
        keep_bam = "{dir}/{X}.on.{ref}.bam",
        discard_bam = temp("{dir}/{X}.on.{ref}.sorted.markdup.bam"),
    log:
        "logs/{dir}/{X}.on.{ref}.bam.log"
    wildcard_constraints:
        X="p1|p2|R|S",
        ref="ref|p1",
        dir="Parents|Bulks"
    params:
        prefix = "{dir}/{X}.on.{ref}",
        dedup_reads = config["dedup_reads"],
    shell: """
        if [ {params.dedup_reads} == False ]; then
            cp {input.bam} {output.keep_bam} && touch {output.discard_bam} && echo 'no dedup, so just copied the bam' > {log}
        else
            /work2/project/gafl/tools/containers/nudup_v2.2.sif -f {input.n6_reads} -o {params.prefix} {input.bam} && mv {params.prefix}.sorted.dedup.bam {output.keep_bam}  && mv {params.prefix}_dup_log.txt {log}
        fi
        """
###### END MAPPING RULES    ##########################################

###### MAPPING STATS RULES  #############################################

rule bam_stat:
    input:
        "{dir}/{X}.on.{ref}{n}bam"
    output:
        stats = report("reports/{dir}/{X}.on.{ref}{n}stats.html", caption="report/bam.rst", category = "bam"),
        flagstats = report("reports/{dir}/{X}.on.{ref}{n}flagstats.txt", caption="report/bam.rst", category = "bam")
    log:
        "logs/{dir}/{X}.on.{ref}{n}stats.log"
    wildcard_constraints:
        dir="Parents|Bulks",
        X="p1|p2|R|S",
        ref="ref|p1",
        n = ".|.\d.", # just a . or .n.
    shell:
        """
        /work2/project/gafl/tools/containers/samtools_v1.9.sif stats {input} | singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif plot-bamstats -p reports/{wildcards.dir}/{wildcards.X}.on.{wildcards.ref}{wildcards.n}stats - > {log} 2>&1
        /work2/project/gafl/tools/containers/samtools_v1.9.sif flagstat {input} > {output.flagstats} 2> {log}
        """

rule sam_stat:
    input:
        "{dir}/{X}.on.{ref}{n}sam"
    output:
        stats = report("reports/{dir}/{X}.on.{ref}{n}stats.html", caption="report/bam.rst", category = "bam"),
        flagstats = report("reports/{dir}/{X}.on.{ref}{n}flagstats.txt", caption="report/bam.rst", category = "bam")
    log:
        "logs/{dir}/{X}.on.{ref}{n}stats.log"
    wildcard_constraints:
        dir="Parents|Bulks",
        X="p1|p2|R|S",
        ref="ref|p1",
        n = ".|.\d.", # just a . or .n.
    shell:
        """
        /work2/project/gafl/tools/containers/samtools_v1.9.sif stats {input} | singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif plot-bamstats -p reports/{wildcards.dir}/{wildcards.X}.on.{wildcards.ref}{wildcards.n}stats - > {log} 2>&1
        /work2/project/gafl/tools/containers/samtools_v1.9.sif flagstat {input} > {output.flagstats} 2> {log}
        """

rule bam_stats:
    input:
        "{path}/{pref}.{suff}"
    wildcard_constraints:
        suff="sam|bam",
    output:
        stats = report("{path}/reports/bam/{pref}.{suff}.stats.html", caption="report/bam.rst", category = "bam"),
        flagstats = report("{path}/reports/bam/{pref}.{suff}.flagstats.txt", caption="report/bam.rst", category = "bam")
    shell:
        """
        /work2/project/gafl/tools/containers/samtools_v1.9.sif stats {input} | singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif plot-bamstats -p {wildcards.path}/reports/bam/{wildcards.pref}.{wildcards.suff}.stats -
        /work2/project/gafl/tools/containers/samtools_v1.9.sif flagstat {input} > {output.flagstats}
        """

###### END MAPPING STATS RULES  #############################################

###### FILE FORMAT CONVERSION RULES  #############################################

rule vcf_to_table:
    input:
        ref = "Parents/p1.fasta",
        vcf = "{dir}/{file}.vcf"
    output:
        "{dir}/{file}.table"
    log:
        "logs/{dir}/{file}.table.log"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif VariantsToTable -V {input.vcf} -R {input.ref} -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -O {output} > {log} 2>&1"

rule bgzip:
    input:
        "{sthg}.vcf"
    output:
        "{sthg}.vcf.gz"
    shell:
        "singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif bgzip -c {input} > {output}"

##### Index generation rules #####
rule fasta_to_sa:
    input:
        "{sthg}.fasta"
    output:
        amb = "{sthg}.fasta.amb",
        ann= "{sthg}.fasta.ann",
        bwt = "{sthg}.fasta.bwt",
        pac = "{sthg}.fasta.pac",
        sa = "{sthg}.fasta.sa"
    shell:
        "/work2/project/gafl/tools/containers/bwa_0.7.17.sif index {input}"

rule bam_to_bai:
    input:
        "{sthg}.bam"
    output:
        "{sthg}.bam.bai"
    shell:
        "/work2/project/gafl/tools/containers/samtools_v1.9.sif index {input} {output}"

rule fasta_to_fai:
    input:
        "{sthg}.fasta"
    output:
        "{sthg}.fasta.fai"
    shell:
        "/work2/project/gafl/tools/containers/samtools_v1.9.sif faidx {input}"

rule gz_to_tbi:
    input:
        "{sthg}.vcf.gz"
    output:
        "{sthg}.vcf.gz.tbi"
    shell:
        "singularity exec /work2/project/gafl/tools/containers/samtools_v1.9.sif tabix {input}"

rule vcf_to_idx:
    input:
        "{sthg}.vcf",
    output:
        "{sthg}.vcf.idx"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif IndexFeatureFile -I {input}"

rule fasta_to_dict:
    input:
        ref = "{sthg}.fasta",
    output:
        "{sthg}.dict"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif CreateSequenceDictionary -R {input.ref} -O {output}"

###### END FILE FORMAT CONVERSION RULES  #############################################
