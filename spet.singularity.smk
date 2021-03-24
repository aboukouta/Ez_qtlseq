# convert bam to bed keeping cigar information
# *.bed: 1 = chrom, 2 = start, 3 = end, 4 = readname, 5 = quality, 6 = strand, 7 = cigar
rule bam_to_bed:
    input:
        "Parents/{X}.bam"
    output:
        temp("Parents/{X}.bed")
    wildcard_constraints:
        X="p2.on.p1"
    shell:
        "/work2/project/gafl/tools/containers/bedtools_v2.29.2.sif bamtobed -cigar -i {input} | singularity exec /work2/project/gafl/tools/containers/bedtools_v2.29.2.sif sort -k 1,1 -k2,2n -k3,3n  > {output} "

# Merge enforcing strandedness (-s)
# replace $4 which is the count of reads in a region by a name of type SLmic1.0_ch05:389473-389619.zoneid
# the count goes in $7
rule bed_merge_filter:
    input:
        "Parents/{X}.bed"
    output:
        "Parents/{X}.merged.bed"
    wildcard_constraints:
        X="p2.on.p1"
    params:
        minNbReads = 100
    shell:
        "/work2/project/gafl/tools/containers/bedtools_v2.29.2.sif merge -i {input} -s -c 4,5,6 -o count,max,first | singularity exec /work2/project/gafl/tools/containers/bedtools_v2.29.2.sif awk '{{ if ($4 > {params.minNbReads}) print }}' | \
         singularity exec /work2/project/gafl/tools/containers/bedtools_v2.29.2.sif awk '{{FS=OFS=\"\t\"}} {{$7=$4; $4=$1\":\"$2\"-\"$3\".\"FNR; print $0}}' > {output}"

rule keep_spets_with_snps:
    input:
        bed = "Parents/p2.on.p1.merged.bed",
        vcf = "Parents/p2.on.p1.fb.final.vcf.gz"
    output:
        "Parents/p1.regions.snps.bed"
    shell:
        "/work2/project/gafl/tools/containers/bedtools_v2.29.2.sif intersect -u -a {input.bed} -b {input.vcf} >  {output}"

# When there are indels, CrossMap will generate 2 bed entries on each side of the indel.
# With bedtools merge we can reunite such entries if they are less than mergeSplit apart
rule convert_spets_to_p2_coords:
    input:
        p1_bed = "Parents/p1.regions.snps.bed",
        chain = "Parents/p2.chain"
    output:
        "Parents/p2.regions.snps.split.bed",
    shell:
        "/work2/project/gafl/tools/containers/crossmap_v0.4.1.sif bed {input.chain} {input.p1_bed} {output}"

rule merge_p2_splits:
    input:
        "Parents/p2.regions.snps.split.bed",
    output:
        "Parents/p2.regions.snps.bed",
    params:
        mergeSplit = 50
    shell:
        "/work2/project/gafl/tools/containers/bedtools_v2.29.2.sif merge -i {input} -d {params.mergeSplit} -c 4,5,6,7 -o first,first,first,first > {output}"

rule get_regions_fasta:
    input:
        fasta = "Parents/{X}.fasta",
        fai = "Parents/{X}.fasta.fai",
        bed = "Parents/{X}.regions.snps.bed"
    output:
        "Parents/{X}.regions.snps.fasta"
    shell:
        "/work2/project/gafl/tools/containers/bedtools_v2.29.2.sif getfasta -fi {input.fasta} -bed {input.bed} -name  | singularity exec /work2/project/gafl/tools/containers/bedtools_v2.29.2.sif sed  's/>.*/&_{wildcards.X}/' > {output}"

rule interleave:
    input:
        p1 = "Parents/p1.regions.snps.fasta",
        p2 = "Parents/p2.regions.snps.fasta"
    output:
        "Parents/p1p2.regions.snps.fasta"
    shell:
        "cat {input.p1} {input.p2} | /work2/project/gafl/tools/containers/bioawk_v1.0.sif -c fastx '{{print}}' | singularity exec /work2/project/gafl/tools/containers/bioawk_v1.0.sif sort -V | /work2/project/gafl/tools/containers/bioawk_v1.0.sif '{{print \">\"$1;print $2}}' > {output}"



# groupby exact regions (1,2,3,6) (chrom,start,end,strand) and keep: (4) the first read name (5) the max quality (6) the first strand, (7) the count of reads
# *.grouped.bed: 1 = chrom, 2 = start, 3 = end, 4 = strand, 5 = readname, 6 = maxquality, 7 = strand, 8 = read count
# rule bed_group:
#     input:
#         "Parents/{X}.bed"
#     output:
#         "Parents/{X}.grouped.bed"
#     wildcard_constraints:
#         X="p2.on.p1"
#     shell:
#         "bedtools groupby -i {input} -g 1,2,3,6 -c 4,5,6,7 -ops first,max,first,count > {output}"

# remove column 4 which is the strand at the wrong place
# keep only regions for which the maxquality is at least 50 and for which there are at least 10 reads
# *.grouped.filtered.bed: 1 = chrom, 2 = start, 3 = end, 4 = readname, 5 = maxquality, 6 = strand, 7 = read count
# rule bed_filter:
#     input:
#         "Parents/{RorH}.grouped.bed"
#     wildcard_constraints:
#         RorH="RdB|H33"
#     output:
#         "Parents/{RorH}.grouped.filtered.bed"
#     params:
#         minNbReads = 10,
#         minMaxQuality = 49
#     shell:
#         "cut -f 1,2,3,5,6,7,8 {input} | awk '{{ if ($5 > {params.minMaxQuality} && $7 > {params.minNbReads}) print }}' > {output}"


#########################
# intersection
#########################
# Replace read names with the id of the region (the line number in this case)
# intersect.bed: 1 = chrom, 2 = start,  3 = end,  4 = id,  5 = maxquality,  6 = strand, 7 = read count
# intersect.bed: 8 = chrom, 9 = start, 10 = end, 11 = id, 12 = maxquality, 13 = strand, 14 = read count
# intersect.bed  15 = overlap, 16 = inter_start, 17 = inter_end
# in case where $9 or $10 = -1 (region found on RdB but not on H33), use the $2 and $3 column
# otherwise, $16 = max($2,$9) and $17 = min($3,$10) so we have the real intersection in cols 9,10
# rule intersect:
#     input:
#         a = "Parents/RdB.merged.bed",
#         b = "Parents/H33.merged.bed"
#     output:
#         "Parents/intersect.bed"
#     shell:
#         "bedtools intersect -a {input.a} -b {input.b} -wao | awk '{{FS=OFS=\"\t\"}} {{$4=$11=FNR; print $0, $9==-1?$2:$2<$9?$9:$2, $10==-1?$3:$3<$10?$3:$10}}' > {output}"
#
# # keep the strand in col 6 for the -s option in the next command
# rule get_raw_bed:
#     input:
#         "Parents/intersect.bed"
#     output:
#         "Parents/intersect.raw.bed"
#     shell:
#         "awk '{{FS=OFS=\"\t\"}} {{print $1,$16,$17,$4,$5,$6}}' {input} > {output}"


# rule blast_regions_on_primers:
#     input:
#         db = join(probes_dir,"all_probes.fasta"),
#         fasta = "Parents/regions.fasta"
#     output:
#         "Parents/regions.on.primers.xml"
#     shell:
#         "blastn -db {input.db} -query {input.fasta} -word_size 10 -outfmt 5 -dust no > {output}"
#
# # make correspondence between region and designed primers. In some (rare) cases, there is more than one match which
# # means that there are up to
# rule add_primer_info:
#     input:
#         bed = "Parents/intersect.bed",
#         blast = "Parents/regions.on.primers.xml"
#     output:
#         "Parents/intersect.probes.bed",
#     run:
#         primers.append_primer_info(input.bed,input.blast,output[0])
#
# # tag regions with "keep" or "reject" after column 19 (id column)
# # keep when
# # minNbReads <= #reads <= maxNbReads
# # minRdBtoH33ratio <= RdB/H33 <= maxRdBtoH33ratio
# # have been found on both parents (which is implicit as the ratio is -1 if no match on H33)
# # have been matched to at least one primer [there are 2 regions with no match at all and 10 with 2 matches]
# rule filter_regions:
#     input:
#         bed = "Parents/intersect.probes.bed",
#         rdb_bam = "Parents/RdB.bam",
#         rdb_bai = "Parents/RdB.bam.bai",
#         h33_bam = "Parents/H33.bam",
#         h33_bai = "Parents/H33.bam.bai"
#     output:
#         "Parents/intersect.probes.filtered.bed"
#     params:
#         minNbReads = 100,
#         maxNbReads = 5000,
#         minRdBtoH33ratio = 0.5,
#         maxRdBtoH33ratio = 4
#     run:
#         coverages.update_counts_and_filter(input.bed, output[0], input.rdb_bam, input.h33_bam, params.minNbReads,params.maxNbReads,params.minRdBtoH33ratio,params.maxRdBtoH33ratio)
#
# # consider only lines with "keep" tag
# # remove the first 40 or the last 40 bases (primer size) depending on the strand
# # keep the id of the region (col 4)
# rule get_mask:
#     input:
#         "Parents/intersect.probes.filtered.bed",
#     output:
#         "Parents/mask.bed",
#     shell:
#         "awk '{{OFS=\"\t\"}}  {{ if ($21 == \"keep\") {{ if ($6 == \"+\") $18=$18+40; else $19=$19-40; print $1,$18,$19,$4 }} }}' {input} > {output}"
#
# rule per_chrom_stat:
#     input:
#         "Parents/intersect.probes.filtered.bed",
#     output:
#         "Parents/per_chrom.txt",
#     shell:
#         "grep 'keep' {input} | cut -f 1 | uniq -c > {output}"
