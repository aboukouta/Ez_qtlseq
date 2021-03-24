rule call_fb_for_bulks:
    input:
        r_bam = "Bulks/R.on.p1.bam",
        r_bai = "Bulks/R.on.p1.bam.bai",
        s_bam = "Bulks/S.on.p1.bam",
        s_bai = "Bulks/S.on.p1.bam.bai",
        ref = "Parents/p1.fasta",
        ref_fai = "Parents/p1.fasta.fai",
    output:
        "Bulks/RS.fb.vcf"
    log:
        "logs/Bulks/RS.fb.vcf.log"
    params:
        pooledC = "--pooled-continuous",
        hapLen = "--haplotype-length 0",
        cut_size = 1000000,
        threads = 10,
    shell:
        "singularity exec /work2/project/gafl/tools/containers/freebayes_v1.3.1.sif freebayes-parallel <( singularity exec /work2/project/gafl/tools/containers/freebayes_v1.3.1.sif fasta_generate_regions.py {input.ref_fai} {params.cut_size}) {params.threads} --genotype-qualities {params.hapLen} {params.pooledC} -f {input.ref} -b {input.r_bam} -b {input.s_bam} > {output} 2> {log}"

rule intersect_with_parents:
    input:
        ref = "Parents/p1.fasta",
        vcf = "Bulks/RS.fb.vcf",
        vcf_index = "Bulks/RS.fb.vcf.idx",
        parents_vcf = "Parents/p2.on.p1.fb.final.vcf",
        parents_vcf_index = "Parents/p2.on.p1.fb.final.vcf.idx"
    output:
        "Bulks/RSp.fb.vcf"
    log:
        "logs/Bulks/RSp.fb.vcf.log"
    shell:
        "/work2/project/gafl/tools/containers/gatk4_v4.1.4.1.sif SelectVariants -V {input.vcf} --concordance {input.parents_vcf} -R {input.ref} -O {output} > {log} 2>&1"

rule do_qtlseqR:
    input:
        snps = "Bulks/RSp.fb.table",
    output:
        jpg6 = "Bulks/6.Takagi.jpg",
        qtlsT = "Bulks/qtlsT.csv",
        qtlsG = "Bulks/qtlsG.csv",
        filtered = "Bulks/filtered.csv",
        stats = "Bulks/stats.txt"
    params:
        R_bulk_size      = config["R_bulk_size"],
        S_bulk_size      = config["S_bulk_size"],
        nb_takagi_reps      = config["nb_takagi_reps"],
        filter_threshold      = config["filter_threshold"],
        window_size      = config["window_size"],
        false_discovery_rate_G      = config["false_discovery_rate_G"],
        min_depth_in_bulk      = config["min_depth_in_bulk"],
        max_depth_in_bulk      = config["max_depth_in_bulk"]
    shell:
        """
        singularity exec /work2/project/gafl/tools/containers/R_QTLseqr_v0.7.5.2.sif Rscript scripts/do_qtl.R \
        {input.snps} \
        {output.stats} \
        {output.filtered} \
        {output.qtlsT} \
        {output.qtlsG} \
        {params.R_bulk_size} \
        {params.S_bulk_size} \
        {params.nb_takagi_reps} \
        {params.filter_threshold} \
        {params.window_size} \
        {params.false_discovery_rate_G} \
        {output.jpg6} \
        {params.min_depth_in_bulk} \
        {params.max_depth_in_bulk}
        """


