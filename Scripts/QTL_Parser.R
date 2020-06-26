install.packages('vcfR', repos='http://cran.us.r-project.org')
library(VariantAnnotation)
#library(vcfR)
library(yaml)
config <- yaml.load_file("config.yaml")
QTL_VCF_to_Table_Parser<-function(vcf_file_path,No_of_Chromosomes,HighBulk,LowBulk){
  
rice_vcf <- readVcf(vcf_file_path)
scan_vcf<- scanVcf(vcf_file_path)
ADREF_HighBulk <- paste0("AD_REF.",HighBulk)
ADALT_HighBulk <- paste0("AD_ALT.",HighBulk)
ADREF_LowBulk <- paste0("AD_REF.",LowBulk)
ADALT_LowBulk <- paste0("AD_ALT.",LowBulk)
GQ_HighBulk <- paste0("GQ.",HighBulk)
GQ_LowBulk <- paste0("GQ.",LowBulk)
table_colnames<- c("CHROM","POS","REF","ALT",ADREF_HighBulk,ADALT_HighBulk,GQ_HighBulk,ADREF_LowBulk,ADALT_LowBulk,GQ_LowBulk)
  
  
vcf_len<- scan_vcf$`*:*-*`$rowRanges@elementMetadata@nrows
chrom_vcf_names<- as.character.Array(droplevels(scan_vcf$`*:*-*`$rowRanges@seqnames@values[1:No_of_Chromosomes]))
QTL_Rice <- as.data.frame(matrix(ncol=10,nrow = vcf_len))
QTL_Rice$V2<-as.array(rice_vcf@rowRanges@ranges@start)
QTL_Rice$V1<-unlist(lapply(strsplit(as.array(rice_vcf@rowRanges@ranges@NAMES),':'),'[[',1))
QTL_Rice$V3<-scan_vcf$`*:*-*`$REF
QTL_Rice$V4<-scan_vcf$`*:*-*`$ALT
QTL_Rice$V7<- scan_vcf$`*:*-*`$GENO$GQ[,HighBulk]
QTL_Rice$V10<- scan_vcf$`*:*-*`$GENO$GQ[,LowBulk]
QTL_Rice$V5<-lapply(scan_vcf$`*:*-*`$GENO$AD[,HighBulk],'[[',1)
QTL_Rice$V6<-lapply(scan_vcf$`*:*-*`$GENO$AD[,HighBulk],'[[',2)
QTL_Rice$V8<-lapply(scan_vcf$`*:*-*`$GENO$AD[,LowBulk],'[[',1)
QTL_Rice$V9<-lapply(scan_vcf$`*:*-*`$GENO$AD[,LowBulk],'[[',2)
QTL_Rice$V2<- as.integer.Array(QTL_Rice$V2)
QTL_Rice$V3<- as.character.Array(QTL_Rice$V3)
QTL_Rice$V4<- as.character.Array(QTL_Rice$V4)
QTL_Rice$V5<- as.integer.Array(QTL_Rice$V5)
QTL_Rice$V6<- as.integer.Array(QTL_Rice$V6)
QTL_Rice$V7<- as.integer.Array(QTL_Rice$V7)
QTL_Rice$V8<- as.integer.Array(QTL_Rice$V8)
QTL_Rice$V9<- as.integer.Array(QTL_Rice$V9)
QTL_Rice$V10<- as.integer.Array(QTL_Rice$V10)

colnames(QTL_Rice) <- table_colnames


  
write.table(QTL_Rice,'QTL_VCF_to_Table/QTL_Table.csv',row.names = FALSE,col.names = TRUE,sep = '\t')
}

QTL_VCF_to_Table_Parser(config$QTL_Parser$VCF_File_Name_or_Path,config$QTL_Parser$Number_of_Chromosomes,config$QTL_Parser$HighBulk,config$QTL_Parser$LowBulk)


