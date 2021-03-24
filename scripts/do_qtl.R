library(QTLseqr)
library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

snpfile               = args[1]    #snakemake@input[["snps"]]
statsfile             = args[2]    #snakemake@output[["stats"]]
filteredfile          = args[3]    #snakemake@output[["filtered"]]
qtlTfile              = args[4]    #snakemake@output[["qtlsT"]]
qtlGfile              = args[5]    #snakemake@output[["qtlsG"]]
bulkSizeR             = as.numeric(args[6])    #snakemake@config[["R_bulk_size"]]
bulkSizeS             = as.numeric(args[7])    #snakemake@config[["S_bulk_size"]]
nbReps                = as.numeric(args[8])    #snakemake@config[["nb_takagi_reps"]]
filterThreshold       = as.numeric(args[9])    #snakemake@config[["filter_threshold"]]
windowSize            = as.numeric(args[10])   #snakemake@config[["window_size"]]
falseDiscoveryRateG   = as.numeric(args[11])   #snakemake@config[["false_discovery_rate_G"]]
graph6                = args[12]   #snakemake@output[["jpg6"]]
minTotalDepth         = 2*as.numeric(args[13]) #2*snakemake@config[["min_depth_in_bulk"]]
maxTotalDepth         = 2*as.numeric(args[14]) #2*snakemake@config[["max_depth_in_bulk"]]
minSampleDepth        =  as.numeric(args[13])  #snakemake@config[["min_depth_in_bulk"]]

pref    = dirname(graph6)
graph1  = file.path(pref,"1.DP.dist.jpg")
graph2  = file.path(pref,"2.REF_FRQ.dist.jpg")
graph3a = file.path(pref,"3.SNP_LOW.jpg")
graph3b = file.path(pref,"3.SNP_HIGH.jpg")
graph4a = file.path(pref,"4.Gp.dist.a.jpg")
graph4b = file.path(pref,"4.Gp.dist.b.jpg")
graph5  = file.path(pref,"5.SNP.density.jpg")
graph7  = file.path(pref,"7.Gp.jpg")

#### non-configurable parameters ####
# graphWidth = 1500
# graphHeight = 600
refAlleleFreq = -0.001
depthDifference = 10000
minGQ = 0

dataSNP <-
  importFromGATK(
    file = file.path(snpfile),
    highBulk = "R.on.p1",
    lowBulk = "S.on.p1"
  )


ggplot(data = dataSNP) +
 geom_histogram(aes(x = DP.HIGH + DP.LOW))
ggsave(graph1)

ggplot(data = dataSNP) +
  geom_histogram(aes(x = REF_FRQ))
ggsave(graph2)

ggplot(data = dataSNP) +
  geom_histogram(aes(x = SNPindex.LOW))
ggsave(graph3a)

ggplot(data = dataSNP) +
  geom_histogram(aes(x = SNPindex.HIGH))
ggsave(graph3b)

per_chrom_count <- dataSNP %>% count(CHROM, name="num_snps_before_filter")

df_filt <- filterSNPs(SNPset = dataSNP,refAlleleFreq = refAlleleFreq, minTotalDepth = minTotalDepth, maxTotalDepth = maxTotalDepth, 
                      minSampleDepth = minSampleDepth, minGQ = minGQ, verbose = TRUE)

per_chrom_count <- right_join(per_chrom_count, df_filt %>% count(CHROM, name = "num_snps_after_filter"))
per_chrom_count <- per_chrom_count %>% add_row(CHROM = "Total", 
                            num_snps_before_filter = sum(per_chrom_count$num_snps_before_filter), 
                            num_snps_after_filter = sum(per_chrom_count$num_snps_after_filter))

write_tsv(per_chrom_count, statsfile)

write_csv(select(df_filt %>% mutate(REF_FRQ = format(REF_FRQ, digits =1), deltaSNP = format(deltaSNP * 100, digits=0)),-c(8,9,10,14,15,16)),filteredfile)

df_filt <- runQTLseqAnalysis(df_filt, windowSize = windowSize, popStruc = "F2", bulkSize = c(bulkSizeR,bulkSizeS), replications = nbReps, intervals = c(95, 99))

df_filt <- runGprimeAnalysis(df_filt, windowSize = windowSize, outlierFilter = "deltaSNP",filterThreshold = filterThreshold)

plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
ggsave(graph4a)

plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP", filterThreshold = filterThreshold)
ggsave(graph4b)

plotQTLStats(SNPset = df_filt, var = "nSNPs")
ggsave(graph5, width = 50, height = 10, units = "cm")

plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
ggsave(graph6, width = 50, height = 10, units = "cm")

plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = falseDiscoveryRateG)
ggsave(graph7, width = 50, height = 10, units = "cm")

result = tryCatch({
  getQTLTable(
    SNPset = df_filt,
    method = "Gprime",
    export = TRUE,
    fileName = qtlGfile
  )
}, warning = function(warning_condition) {
  file.create(qtlGfile)
  message("warning Gprim")
}, error = function(error_condition) {
  message("error Gprime")
})

result = tryCatch({
  getQTLTable(
    SNPset = df_filt,
    method = "QTLseq",
    export = TRUE,
    fileName = qtlTfile
  )
}, warning = function(warning_condition) {
  file.create(qtlTfile)
  message("warning Takagi")
}, error = function(error_condition) {
  message("error Takagi")
})




#c1 = c("SLmic1.0_ch02", "SLmic1.0_ch04", "SLmic1.0_ch05","SLmic1.0_ch06")
# c1 = c("SLmic1.0_ch04")
# 
# plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = fdrThreshold, subset = c1)
# ggsave(graph8)
# 
# plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE, subset = c1)
# ggsave(graph9)
