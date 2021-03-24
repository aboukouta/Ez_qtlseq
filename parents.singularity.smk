
# Sequential rules

rule call_fb_for_parent:
    input:
        bam = "Parents/p{i}.on.{ref}.bam",
        bai = "Parents/p{i}.on.{ref}.bam.bai",
        ref = "Parents/{ref}.fasta",
        ref_fai = "Parents/{ref}.fasta.fai",
    output:
        "Parents/p{i}.on.{ref}.fb.vcf"
    log:
        "logs/Parents/p{i}.on.{ref}.fb.vcf.log"
    params:
        minCoverage = config["parentsMinCoverage"],
        minObservationsToCall = config["parentsMinObservationsToCall"],
        hapLen = "--haplotype-length 0",
        cut_size = 1000000,
        threads = 32,
    wildcard_constraints:
        pi="p1|p2",
        ref="ref|p1"
    shell:
        "singularity exec /work2/project/gafl/tools/containers/freebayes_v1.3.1.sif freebayes-parallel <( singularity exec /work2/project/gafl/tools/containers/freebayes_v1.3.1.sif fasta_generate_regions.py {input.ref_fai} {params.cut_size}) {params.threads} --genotype-qualities {params.hapLen} -C {params.minObservationsToCall} --min-coverage {params.minCoverage} -f {input.ref} {input.bam} > {output} 2> {log}"

rule p1_fasta:
    input:
        vcf = "Parents/p1.on.ref.fb.vcf.gz",
        tbi = "Parents/p1.on.ref.fb.vcf.gz.tbi",
        ref = "Parents/ref.fasta",
    output:
        "Parents/p1.fasta"
    log:
        "logs/Parents/p1.fasta.log"
    shell:
        "/work2/project/gafl/tools/containers/bcftools_v1.10.2.sif consensus --fasta-ref {input.ref} {input.vcf} -H R -o {output} > {log} 2>&1"   #" a modifier bcftools_v1.11.sif"

# generate p2 sequence incorporating all homozyg snps of p2
rule p2_fasta:
    input:
        vcf = "Parents/p2.on.p1.fb.vcf.gz",
        tbi = "Parents/p2.on.p1.fb.vcf.gz.tbi",
        ref = "Parents/p1.fasta",
    output:
        fasta = "Parents/p2.fasta",
        chain = "Parents/p2.chain"
    log:
        "logs/Parents/p2.fasta.log"
    shell:
        "/work2/project/gafl/tools/containers/bcftools_v1.10.2.sif consensus --fasta-ref {input.ref} {input.vcf} -H R -o {output.fasta} -c {output.chain} > {log} 2>&1"

rule keep_alt_homozyg:
    input:
        ref = "Parents/p1.fasta",
        ref_fai = "Parents/p1.fasta.fai",
        ref_dict = "Parents/p1.dict",
        p2_vcf = "Parents/p2.on.p1.fb.vcf",
    output:
        "Parents/p2.on.p1.fb.hz.vcf"
    log:
        "logs/Parents/p2.on.p1.fb.hz.vcf.log"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif SelectVariants -V {input.p2_vcf} -select \'vc.getHomVarCount() >= 1\' -R {input.ref} -O {output} > {log} 2>&1 "

rule exclude_p1_var:
    input:
        ref = "Parents/p1.fasta",
        ref_dict = "Parents/p1.dict",
        p1_vcf = "Parents/p1.on.p1.fb.vcf",
        p1_ind = "Parents/p1.on.p1.fb.vcf.idx",
        p2_vcf = "Parents/p2.on.p1.fb.hz.vcf",
        p2_ind = "Parents/p2.on.p1.fb.hz.vcf.idx"
    output:
        "Parents/p2.on.p1.fb.final.vcf"
    log:
        "logs/Parents/p2.on.p1.fb.final.vcf.log"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif SelectVariants -V {input.p2_vcf} --discordance {input.p1_vcf} -R {input.ref} -O {output} > {log} 2>&1"

# this is not used for the moment
# bulks still rely on p2.on.p1.fb.final.vcf
rule call_fb_p1p2_on_p1:
    input:
        bam1 = "Parents/p1.on.p1.bam",
        bam2 = "Parents/p2.on.p1.bam",
        ref = "Parents/p1.fasta",
        ref_idx = "Parents/p1.fasta.fai",
    output:
        "Parents/p12.on.p1.fb.vcf"
    log:
        "logs/Parents/p12.on.p1.fb.vcf"
    params:
        minCoverage = config["parentsMinCoverage"],
        minObservationsToCall = config["parentsMinObservationsToCall"],
        hapLen = "--haplotype-length 0"
    shell:
        "/work2/project/gafl/tools/containers/freebayes_v1.3.1.sif --genotype-qualities {params.hapLen} -C {params.minObservationsToCall} --min-coverage {params.minCoverage} -f {input.ref} \
         -b {input.bam1} -b {input.bam2} > {output} 2> {log}"

