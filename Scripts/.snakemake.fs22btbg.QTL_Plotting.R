
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('QTL_VCF_to_Table/QTL_Table.csv'),
    output = list('QTL_Plots/DP_Filtering Data.pdf', 'QTL_Plots/REF Frequency for Filtering Data.pdf', 'QTL_Plots/SNP Index for Filtering Data.pdf', 'QTL_Plots/GPrime Distribution with Hampel Outlier Filter.pdf', 'QTL_Plots/GPrime Distribution with deltaSNP Outlier Filter.pdf', 'QTL_Plots/SNP Density Plot.pdf', 'QTL_Plots/Delta SNP Index Plot with Intervals.pdf', 'QTL_Plots/GPrime Value Plot.pdf'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("vcffilter" = list("filters" = '-i "GT[0]=\'0/0\' && GT[1]=\'1/1\' || GT[0]=\'1/1\' && GT[1]=\'0/0\'"'), "vcf_file" = 'pools_filtered_removedjunk', "QTL_Parser" = list("VCF_File_Name_or_Path" = 'Homozygous_Filtered_VCF/Homozygous_pools_filtered_removedjunk.vcf.gz', "Number_of_Chromosomes" = 12, "HighBulk" = 'rice_pools_res', "LowBulk" = 'rice_pools_sus', "WindowSize" = 1000000, "High_Bulk_Size" = 178, "Low_Bulk_Size" = 280, "REF_Allel_Frequency" = 0.2, "Min_Total_Depth" = 20, "Max_Total_Depth" = 100, "Depth_Difference" = 80, "Min_Sample_Depth" = 20, "Min_GQ" = 80, "Filter_Threshold" = 0.1, "FDR_q" = 0.01), "AllelPlots" = list("vcffile" = 'freebayes_D2.filtered', "Samples" = '-s con-all,D2,D2_F2_tt,D2_F2_TT', "data" = '-f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%RO\\t%AO\\t[,%GT]\\t[,%GQ]\\t[,%RO]\\t[,%AO]\\n"')),
    rule = 'QTL_Plotting'
)
######## Original script #########
library(yaml)
library(data.table)
library("QTLseqr")
library("ggplot2")
config <- yaml.load_file("config.yaml")
QTL_Plotting <- function(HighBulk,LowBulk, No_of_Chromosomes){
  options(warn=-1)
  table_file <- 'QTL_VCF_to_Table/QTL_Table.csv'
  chrom_names<-as.list(unique(fread(table_file,sep = "\t", select = c("CHROM") ))) 
  #chrom_vcf_names<- as.character.Array(droplevels(scan_vcf$`*:*-*`$rowRanges@seqnames@values[1:No_of_Chromosomes]))
  df <- importFromTable(file = table_file,
                        highBulk = HighBulk,
                        lowBulk = LowBulk,
                        chromList = chrom_names$CHROM[1:No_of_Chromosomes],sep = '\t')
  
  ggplot(data = df)+geom_histogram(aes(x = DP.HIGH + DP.LOW))+xlim(0,1000)
  ggsave(filename="QTL_Plots/DP_Filtering Data",device= "pdf", width=20, height=10)
  ggplot(data = df)+geom_histogram(aes(x = REF_FRQ))
  ggsave(filename="QTL_Plots/REF Frequency for Filtering Data",device = "pdf", width=20, height=10)
  ggplot(data = df)+geom_histogram(aes(x = SNPindex.HIGH))
  ggsave(filename = "QTL_Plots/SNP Index for Filtering Data",device="pdf", width=20, height=10)
  
  
  df_filt <-
    filterSNPs(
      SNPset = df,
      refAlleleFreq = config$QTL_Parser$REF_Allel_Frequency,
      minTotalDepth = config$QTL_Parser$Min_Total_Depth,
      maxTotalDepth = config$QTL_Parser$Max_Total_Depth,
      depthDifference = config$QTL_Parser$Depth_Difference,
      minSampleDepth = config$QTL_Parser$Min_Sample_Depth,
      minGQ = config$QTL_Parser$Min_GQ,
      verbose = TRUE
    )
  
  df_filt <- runQTLseqAnalysis(df_filt,
                               windowSize = config$QTL_Parser$WindowSize,
                               popStruc = "F2",
                               bulkSize = c(config$QTL_Parser$High_Bulk_Size, config$QTL_Parser$Low_Bulk_Size),
                               replications = 10000,
                               intervals = c(95, 99))
                               
  df_filt <- runGprimeAnalysis(df_filt,windowSize = config$QTL_Parser$WindowSize,outlierFilter = "deltaSNP",filterThreshold = config$QTL_Parser$Filter_Threshold)
  
  plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
  ggsave(filename = "QTL_Plots/GPrime Distribution with Hampel Outlier Filter",device = "pdf", width=20, height=10)
  
  plotGprimeDist(SNPset =df_filt, outlierFilter = "deltaSNP", filterThreshold = config$QTL_Parser$Filter_Threshold)
  ggsave(filename = "QTL_Plots/GPrime Distribution with deltaSNP Outlier Filter",device = "pdf", width=20, height=10)
  
  p1 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")
  ggsave(filename = "QTL_Plots/SNP Density Plot",plot=p1,device = "pdf", width=20, height=10)
  
  p2 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
  ggsave(filename = "QTL_Plots/Delta SNP Index Plot with Intervals",plot=p2,device = "pdf", width=20, height=10)
  p3 <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = config$QTL_Parser$FDR_q)
  ggsave(filename = "QTL_Plots/GPrime Value Plot",plot=p3,device = "pdf", width=20, height=10)
  
}

QTL_Plotting(config$QTL_Parser$HighBulk,config$QTL_Parser$LowBulk,config$QTL_Parser$Number_of_Chromosomes)
