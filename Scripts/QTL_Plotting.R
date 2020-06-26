library(yaml)
library(data.table)
library(devtools)

install_github("bmansfeld/QTLseqr",upgrade_dependencies=  TRUE)
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
