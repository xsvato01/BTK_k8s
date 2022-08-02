#!/usr/bin/env Rscript
library(data.table)
library(tidyr)

# Run Rscript --vanilla rearange_table.R /PATH_TO/annot "run_ID"
args = commandArgs(trailingOnly=TRUE)
print(args[1])
print(args[2])

# path for testing

final_list<-list()
# rearange table
for(file in list.files(args[1], pattern = "*normVEP.txt")){
  
  # load file
  table_orig <- read.table(paste(args[1],"/", file, sep = ""), header=T)
  
  # sample name
  table_orig<-cbind(Sample_name = gsub(".norm.merged.annot.normVEP.txt", "", file), table_orig)
  
  # number of callers
  for(line in 1:length(table_orig$Sample_name)){
    if(table_orig$set[line] == paste0("vardict-varscan_SNV") || table_orig$set[line] == paste0("vardict-varscan_INDEL")){
      table_orig$VC_DETECTION_COUNT[line] <- "2"
    }else{
      table_orig$VC_DETECTION_COUNT[line] <- "1"
    }
  }
  
  table_orig1 <-table_orig[grep("NM_002661", table_orig$VEP_Feature),]
  table_orig2 <-table_orig[grep("NM_000061.2", table_orig$VEP_Feature),]
  table_orig3 <-table_orig[grep("NM_000546", table_orig$VEP_Feature),]
  # DNMT3A 
  table_orig4 <-table_orig[grep("NM_022552", table_orig$VEP_Feature),]
  # JAK2
  table_orig5 <-table_orig[grep("NM_004972", table_orig$VEP_Feature),]
  # TET2
  table_orig6 <-table_orig[grep("NM_017628", table_orig$VEP_Feature),]
  table_orig <- rbind(table_orig1, table_orig2, table_orig3, table_orig4, table_orig5, table_orig6)
  table_orig <- table_orig[!duplicated(table_orig$VEP_HGVSc),]  
  # final list
  final_list[[file]] <- table_orig
}

final_df <- rbindlist(final_list)

#check variant occurence and count median of AF in all samples
final_df[,Occurence_in_samples := .N, by = VEP_HGVSc]
# for(i in 1:length(final_df$Sample_name)){
#   ## occurence ##
#   tmp1 <-final_df[VEP_HGVSc == final_df$VEP_HGVSc[i]]
#   final_df$Occurence_in_samples[i]<-length(tmp1$VEP_HGVSc)
# }

# remove redundant columns NEW -> not important et all
final2_df<-subset(final_df, select=-c(vardict.FREQ,vardict.GQ,vardict.GT,
                                      vardict.PVAL,vardict.RBQ,vardict.RD,vardict.RDF,vardict.RDR,
                                      vardict.SDP,vardict.VD,ID,QUAL,FILTER,AC,ADJAF,ADP,AF,AMPFLAG,AN,
                                      BIAS,CSQ,GDAMP,HET,HICNT,HOM,LSEQ,MQ,MSI,MSILEN,NC,NCAMP,
                                      NM,ODDRATIO,PMEAN,PSTD,QSTD,QUAL,RSEQ,SAMPLE,SBF,SHIFT3,SN,
                                      TLAMP,TYPE,VARBIAS,VD,VEP_Allele,VEP_BIOTYPE,VEP_DISTANCE,
                                      VEP_Existing_variation,VEP_FLAGS,VEP_Feature_type,VEP_Gene,VEP_HGNC_ID,
                                      VEP_HGVS_OFFSET,VEP_IMPACT,VEP_PolyPhen,VEP_Protein_position,VEP_REFSEQ_MATCH,
                                      VEP_SIFT,VEP_STRAND,VEP_SYMBOL_SOURCE,VEP_cDNA_position,WT,vardict.ABQ))

# put NA to varscan columns in vardict case
for(col_name in grep("varscan", names(final2_df),value = T)){
  final2_df[set == "vardict",(col_name) := NA]
}

# for(j in 1:length(final2_df$Sample_name)){
#   if(final2_df$set[j] == "vardict"){
#     final2_df[j,grep("varscan", names(final2_df))] <- "NA"
#   }  
# }

# join varscan allelic FREQ
final2_df$varscan_INDEL.FREQ <- gsub("^.$", "F", final2_df$varscan_INDEL.FREQ)
final3_df<-unite(final2_df, varscan.FREQ, c(varscan_INDEL.FREQ, varscan_SNV.FREQ), 
                 remove=TRUE, sep = "")

final3_df$varscan.FREQ<-gsub("%.", "%", final3_df$varscan.FREQ)
final3_df$varscan.FREQ<-gsub("F", "", final3_df$varscan.FREQ)  

#### MEDIAN ####
final3_df[varscan.FREQ != "NANA" & varscan.FREQ != ".",final_AF := varscan.FREQ]
final3_df[is.na(final_AF), final_AF := as.numeric(as.character(vardict.AF))*100]
final3_df[,final_AF := as.numeric(gsub("%.*", "", final_AF))]
final3_df[,variant_median := median(final_AF),by = VEP_HGVSc]

#for(p in 1:length(final3_df$Sample_name)){
# 
#   ## median ##
#   tmp <- final3_df[VEP_HGVSc == final3_df$VEP_HGVSc[p]]
#   print("MEDIAN COUNT..")
#   #print(tmp$varscan.FREQ)
# 
#   AFs <- c()
#   # create tmp df with same mutation around sample
#   for(g in 1:length(tmp$Sample_name)){
# 
#     if(tmp$set[g] == paste0("varscan_SNV") || tmp$set[g] == paste0("varscan_INDEL") ||
#        tmp$set[g] == paste0("vardict-varscan_INDEL") || tmp$set[g] == paste0("vardict-varscan_SNV") ||
#        tmp$set[g] == paste0("FilteredInAll") || tmp$set[g] == paste0("filterInvardict-varscan_SNV"))
#     {
#       AFs[g] <- as.numeric(as.character(gsub("%", "", tmp$varscan.FREQ[g])))
#     }
#     else{
#       AFs[g] <- as.numeric(as.character(tmp$vardict.AF[g]))
#     }
# 
#   }
#   final3_df$variant_median[p] <- median(AFs, na.rm=T)
#}

for(p in 1:length(final3_df$Sample_name)){
  
  if(final3_df$varscan.FREQ[p] == "NANA"){
    final3_df$mutation_description[p] <- paste0(gsub("NM_000546.5:","", final3_df$VEP_HGVSc[p]), " ", gsub("NP_000537.3:", "", final3_df$VEP_HGVSp[p]), " ", 
                                                gsub("\\.",",",as.numeric(as.character(final3_df$vardict.AF[p])) * 100), "%") 
  } else {
    final3_df$mutation_description[p] <- paste0(gsub("NM_000546.5:","", final3_df$VEP_HGVSc[p]), " ", gsub("NP_000537.3:", "", final3_df$VEP_HGVSp[p]), " ", 
                                                gsub("\\.", ",", final3_df$varscan.FREQ[p]))
  }
}

### count strand bias for varscan ### NEW
# varFwd/varCoverageFwd
# ---------------------
# varRev/varCoverageRev
# [varscan_SNV.ADF/(varscan_SNV.ADF+varscan_SNV.RDF)]/[varscan_SNV.ADR/(varscan_SNV.ADR+varscan_SNV.RDR)]
print("STRAND BIAS COUNT..")
for(r in 1:length(final3_df$Sample_name)){
  #print(r)
  # all Varscan variants
  if(final3_df$set[r] == paste0("varscan_SNV") || final3_df$set[r] == paste0("varscan_INDEL") || final3_df$set[r] == paste0("vardict-varscan_INDEL") ||
     final3_df$set[r] == paste0("vardict-varscan_SNV") || final3_df$set[r] == paste0("filterInvardict-varscan_SNV")){
    #print(final3_df$set[r])
    
    # all Varscan SNV  
    if(final3_df$set[r] == paste0("varscan_SNV") || final3_df$set[r] == paste0("vardict-varscan_SNV") || final3_df$set[r] == paste0("filterInvardict-varscan_SNV")){
      #print("SNV")
      citatel <- as.numeric(as.character(final3_df$varscan_SNV.ADF[r]))/(as.numeric(as.character(final3_df$varscan_SNV.ADF[r]))+as.numeric(as.character(final3_df$varscan_SNV.RDF[r])))
      jmenovatel <- as.numeric(as.character(final3_df$varscan_SNV.ADR[r]))/(as.numeric(as.character(final3_df$varscan_SNV.ADR[r]))+as.numeric(as.character(final3_df$varscan_SNV.RDR[r])))
      final3_df$varscan.SBIAS[r] <- round(citatel/jmenovatel, digits = 2)
    }else{
      # all Varscan INDEL
      #print("INDEL")
      citatel <- as.numeric(as.character(final3_df$varscan_INDEL.ADF[r]))/(as.numeric(as.character(final3_df$varscan_INDEL.ADF[r]))+as.numeric(as.character(final3_df$varscan_INDEL.RDF[r])))
      jmenovatel <- as.numeric(as.character(final3_df$varscan_INDEL.ADR[r]))/(as.numeric(as.character(final3_df$varscan_INDEL.ADR[r]))+as.numeric(as.character(final3_df$varscan_INDEL.RDR[r])))
      final3_df$varscan.SBIAS[r] <- round(citatel/jmenovatel, digits = 2)
    }
  }else{
    final3_df$varscan.SBIAS[r] <- "NA"
  }
}

### split table to individual samples ### NEW
final_4_df_list <- list() 
for(sample in levels(final3_df$Sample_name)){
  final4_df <-final3_df[grep(paste0("^",sample,"$"), final3_df$Sample_name),]
  # change formats of columns
  final4_df$HIAF <- gsub("\\.",",",final4_df$HIAF)
  final4_df$VEP_INTRON <- gsub("/10", "", final4_df$VEP_INTRON)
  final4_df$VEP_EXON <- gsub("/11", "", final4_df$VEP_EXON)
  final4_df$vardict.AD <- gsub(",",";", final4_df$vardict.AD)
  final4_df$vardict.AF <- paste0(gsub("\\.",",",as.numeric(as.character(final4_df$vardict.AF)) * 100))
  colnames(final4_df)[colnames(final4_df) == 'vardict.AF'] <- 'vardict.AF_percentage'
  final4_df$vardict.ALD <- gsub(",",";",final4_df$vardict.ALD)
  final4_df$varscan.FREQ <- gsub("\\.", ",", final4_df$varscan.FREQ)
  final4_df$varscan_INDEL.PVAL <- gsub("\\.", ",", final4_df$varscan_INDEL.PVAL)
  final4_df$varscan_SNV.PVAL <- gsub("\\.", ",", final4_df$varscan_SNV.PVAL)
  final4_df$variant_median <- gsub("\\.", ",", final4_df$variant_median)
  final4_df$varscan.SBIAS <- gsub("Inf", "10000",gsub("\\.", ",", final4_df$varscan.SBIAS))
  final_4_df_list[[sample]] <- final4_df
  
  # not important in per sample table, but occuring in overal table for all samples
  final5_df<-subset(final4_df, select=-c(DP, HIAF, HICOV, QUAL.1, REFBIAS, vardict.DP, vardict.ADF, vardict.ADR, varscan_INDEL.ABQ, varscan_INDEL.AD,
                                         varscan_INDEL.ADF, varscan_INDEL.ADR, varscan_INDEL.AF, varscan_INDEL.ALD,
                                         varscan_INDEL.GQ, varscan_INDEL.GT, varscan_INDEL.PVAL, varscan_INDEL.RBQ,
                                         varscan_INDEL.RD, varscan_INDEL.RDF, varscan_INDEL.RDR, varscan_INDEL.SDP,
                                         varscan_INDEL.VD, varscan_SNV.ABQ, varscan_SNV.AD, varscan_SNV.ADF, varscan_SNV.ADR,
                                         varscan_SNV.AF, varscan_SNV.ALD, varscan_SNV.GQ, varscan_SNV.GT, varscan_SNV.PVAL,
                                         varscan_SNV.RBQ, varscan_SNV.RD, varscan_SNV.RDF, varscan_SNV.RDR, varscan_SNV.SDP,
                                         varscan_SNV.VD, VC_DETECTION_COUNT, final_AF))
  # save per sample table
  write.table(final5_df, file=paste0(args[1],"/", sample, ".sample.merged.anot.txt"), sep="\t", quote = F, row.names = F)
}

final_big_df<-rbindlist(final_4_df_list)
# save final big table
DIR_file_name <- paste0(args[1],"/",args[2],".allsamples.merged.anot.txt")
print(DIR_file_name)

write.table(final_big_df,file=DIR_file_name, sep="\t", quote = F, row.names = F)
