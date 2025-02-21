rm(list=ls())
library(TCGAbiolinks)
library(GEOquery)
library(dplyr)
library(SummarizedExperiment)
library(pbapply)
library(stringr)
library(sva)
library(parallel)
library(org.Hs.eg.db)
library(lubridate)

# Create a parallel socket cluster (if you have any other option, that could be used) 
num_cluster<-4
clu<-makeCluster(num_cluster)

# Data preparation for the comprehensive NSCLC patient transcriptome analysis
# Terms representing elements not available were defined in advance of the analysis
na.terms<-c(NA, "[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Not Evaluated]", "[Unknown]", "MX", "TX", "NX")

# 1. TCGA data (LUAD)
target_cancer <- "TCGA-LUAD"

###### Clinical data (for recurrence)
query<-GDCquery(
  project = target_cancer,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  file.type = ".txt",
)

GDCdownload(
  query = query
)

clinical<-GDCprepare(query)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})

####### Mutation data
query <- GDCquery(
   project = "TCGA-LUAD", 
   data.category = "Simple Nucleotide Variation", 
   access = "open", 
   legacy = FALSE, 
   data.type = "Masked Somatic Mutation", 
 )
 GDCdownload(query)
 maf <- GDCprepare(query)

####### Expression data
query <- GDCquery(
  project = target_cancer,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

clinical_df<-data.frame(data@colData@listData[c('barcode', 'sample_type', 'ajcc_pathologic_stage', 'age_at_index',
                                                "days_to_last_follow_up", "ajcc_pathologic_t", "ajcc_pathologic_n",
                                                "ajcc_pathologic_m", "ajcc_staging_system_edition", 'race', 'gender',
                                                "vital_status", "days_to_death", "paper_expression_subtype")])

# Sample IDs
clinical_df_TCGA_LUAD<-data.frame(barcode=clinical_df$barcode,
                                  ID=sapply(clinical_df$barcode,FUN=function(s){
                                    paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
                                  }))

# Histology
clinical_df_TCGA_LUAD$Histology<-sapply(clinical_df$sample_type, FUN=function(x){
  if(x=="Solid Tissue Normal"){
    return("Control")
  }else if(x=="Primary Tumor"){
    return("ADC")
  }else{
    return("[Not Available]")
  }
})

# Age: age
clinical_df_TCGA_LUAD$Age<-clinical_df$age
clinical_df_TCGA_LUAD$Age<-sapply(clinical_df_TCGA_LUAD$Age, FUN=function(x){
  if(is.na(x)){
    "[Not Available]"
  }else{
    x
  }
})

# survival_time: survival_time
clinical_df_TCGA_LUAD$survival_time<-sapply(c(1:length(clinical_df$vital_status)), FUN=function(x){
  t=setdiff(c(clinical_df$days_to_last_follow_up[x], clinical_df$days_to_death[x]), na.terms)
  if(length(t)==0){
    return("[Not Available]")
  }
  return(max(as.numeric(t)))
})

# vital_status: vital_status
clinical_df_TCGA_LUAD$vital_status<-clinical_df$vital_status

# Disease-free survival: DFS
clinical_df_TCGA_LUAD$DFS<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  followup.nte<-clinical$clinical_follow_up_v1.0_luad %>%
                   filter(bcr_patient_barcode %in% barcode.cut) %>% 
                   filter(new_tumor_event_dx_indicator %in% "YES")
  
  if(dim(followup.nte)[1]>0){
    DFS<-min(followup.nte$new_tumor_event_dx_days_to)
  }
  
  nofollowup.nte<-clinical$clinical_nte_luad %>% 
                    filter(bcr_patient_barcode %in% barcode.cut)
  
  if(dim(nofollowup.nte)[1]>0){
    DFS<-nofollowup.nte$new_tumor_event_dx_days_to
  }
  
  if(!exists("DFS")){
    DFS<-(clinical_df_TCGA_LUAD %>% filter(barcode == s))$survival_time
  }
  return(DFS)
})

# NTE-status
clinical_df_TCGA_LUAD$NTE_status<-sapply(clinical_df_TCGA_LUAD$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  surv.time<-clinical_df_TCGA_LUAD$survival_time[match(s, clinical_df_TCGA_LUAD$barcode)]
  DFS.time<-clinical_df_TCGA_LUAD$DFS[match(s, clinical_df_TCGA_LUAD$barcode)]
  
  if(surv.time!=DFS.time){
    return("YES")
  }
  
  if(surv.time %in% na.terms){
    return("[Not Available]")
  }
  
  if(!is.na(match(barcode.cut, clinical$clinical_nte_luad$bcr_patient_barcode))){
    return("YES")
  }
  
  indicator<-clinical$clinical_follow_up_v1.0_luad$new_tumor_event_dx_indicator[
    match(barcode.cut, clinical$clinical_follow_up_v1.0_luad$bcr_patient_barcode)]
  
  if("YES" %in% indicator){
    return("YES")
  }
  
  if("NO" %in% indicator){
    return("NO")
  }
  
  return("[Not Available]")
})

# ajcc_pathologic_t
clinical_df_TCGA_LUAD$ajcc_pathologic_t<-clinical_df$ajcc_pathologic_t

# ajcc_pathologic_n
clinical_df_TCGA_LUAD$ajcc_pathologic_n<-clinical_df$ajcc_pathologic_n

# ajcc_pathologic_m
clinical_df_TCGA_LUAD$ajcc_pathologic_m<-clinical_df$ajcc_pathologic_m

# ajcc_pathologic_stage
clinical_df_TCGA_LUAD$ajcc_pathologic_stage<-clinical_df$ajcc_pathologic_stage

# Race: race
clinical_df_TCGA_LUAD$Race<-clinical_df$race

# Gender: gender
clinical_df_TCGA_LUAD$Gender<-clinical_df$gender

# Adjuvant_chemo
clinical_df_TCGA_LUAD$Adjuvant_chemo<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Pharmaceutical Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

# Adjuvant_radiation
clinical_df_TCGA_LUAD$Adjuvant_radiation<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Radiation Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

#  Smoking_status
clinical_df_TCGA_LUAD$Smoking_status<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  ind<-match(barcode.cut, clinical$clinical_patient_luad$bcr_patient_barcode)
  smoking.indicator<-clinical$clinical_patient_luad$tobacco_smoking_history_indicator[ind]
  if(smoking.indicator == "1"){
    res<-"Never smoker"
  }else if (smoking.indicator %in% c("2", "3", "4", "5")){
    res<-"Ever smoker"
  }else{
    res<-"[Not Available]"
  }
  return(res)
})

# Predominant_subtype
clinical_df_TCGA_LUAD$Predominant_subtype<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, clinical$clinical_patient_luad$bcr_patient_barcode)
  clinical$clinical_patient_luad$histologic_diagnosis...68[ind]
})

## Mutation statuses
# EGFR mutation
mut<-maf %>% filter(Hugo_Symbol == "EGFR")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUAD$EGFR_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"EGFR-WT"
  else
    res<-"EGFR-Mut"
  return(res)
})

# KRAS mutation
mut<-maf %>% filter(Hugo_Symbol == "KRAS")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUAD$KRAS_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"KRAS-WT"
  else
    res<-"KRAS-Mut"
  return(res)
})

# TP53 mutation
mut<-maf %>% filter(Hugo_Symbol == "TP53")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUAD$TP53_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"TP53-WT"
  else
    res<-"TP53-Mut"
  return(res)
})

# STK11 mutation
mut<-maf %>% filter(Hugo_Symbol == "STK11")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUAD$STK11_mut <-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"STK11-WT"
  else
    res<-"STK11-Mut"
  return(res)
})

##### Expression data processing
data.luad <- assay(data, "fpkm_unstrand")
data.luad.row.info<-rowRanges(data)
protein_coding_gene.ind<-which(data.luad.row.info@elementMetadata@listData$gene_type=="protein_coding")

data.luad.protein_coding<-data.luad[protein_coding_gene.ind, ]

expression_df_luad<-t(data.luad.protein_coding)
colnames(expression_df_luad)<-row.names(data.luad.protein_coding)
expression_df_luad<-log(expression_df_luad+1, base=2)

gene_info<-data.frame(gene_id=data.luad.row.info@elementMetadata@listData$gene_id, 
                      gene_symbol=data.luad.row.info@elementMetadata@listData$gene_name) %>% 
  filter(gene_id %in% colnames(expression_df_luad))

gene_info$ENTREZID<-pbsapply(gene_info$gene_symbol, FUN=function(x){
  library(org.Hs.eg.db)
  tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                 keys=x, 
                                 columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")$ENTREZID[1], 
           error= function(e){ 
             return(NA)
           })
}, cl=clu)

clusterEvalQ(clu, {
  symbol.check<-read.csv("hgnc-symbol-check-tcga.csv")
  NULL
})

gene_info$gene_symbol<-pbsapply(gene_info$gene_symbol, FUN=function(x){
  library(dplyr)
  chk<-symbol.check%>% filter(Input %in% x)
  
  if(length(which(chk$Match.type=="Approved symbol"))>0){
    return(x)
  }else{
    return(chk$Approved.symbol[1])
  }
}, cl=clu)

gene_info$ENTREZID[which(is.na(gene_info$ENTREZID))]<-
  pbsapply(gene_info$gene_symbol[which(is.na(gene_info$ENTREZID))],FUN=function(x){
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                   keys=x, 
                                   columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")$ENTREZID[1], 
             error= function(e){ 
               return(NA)
             })
  }, cl=clu)

# remove the genes that has no gene symbol or gene id
remove.ind<-c(intersect(which(gene_info$gene_symbol==""), which(is.na(gene_info$ENTREZID))),
# remove the genes which are duplicated              
              which(duplicated(gene_info$gene_symbol)))
gene_info<-gene_info[-remove.ind,]



expression_df_luad<-expression_df_luad[,gene_info$gene_id]

# 1. TCGA data (LUSC)
target_cancer <- "TCGA-LUSC"
###### Clinical data (for recurrence)
query<-GDCquery(
  project = target_cancer,
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  file.type = ".txt",
)

GDCdownload(
  query = query
)

clinical<-GDCprepare(query)
clinical<-lapply(clinical, FUN=function(x){
  return(x[-c(1,2),])
})

###### Mutation data
query <- GDCquery(
  project = target_cancer, 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
)
GDCdownload(query)
maf <- GDCprepare(query)

####### Expression data
query <- GDCquery(
  project = target_cancer,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

clinical_df<-data.frame(data@colData@listData[c('barcode', 'sample_type', 'ajcc_pathologic_stage', 'age_at_index',
                                                "days_to_last_follow_up", "ajcc_pathologic_t", "ajcc_pathologic_n",
                                                "ajcc_pathologic_m", "ajcc_staging_system_edition", 'race', 'gender',
                                                "vital_status", "days_to_death")])

# Sample IDs
clinical_df_TCGA_LUSC<-data.frame(barcode=clinical_df$barcode,
                                  ID=sapply(clinical_df$barcode,FUN=function(s){
                                    paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
                                  }))

# Histology
clinical_df_TCGA_LUSC$Histology<-sapply(clinical_df$sample_type, FUN=function(x){
  if(x=="Solid Tissue Normal"){
    return("Control")
  }else if(x=="Primary Tumor"){
    return("SCC")
  }else{
    return("[Not Available]")
  }
})

# Age: age
clinical_df_TCGA_LUSC$Age<-sapply(clinical_df$age, FUN=function(x){
  if(is.na(x)){
    "[Not Available]"
  }else{
    x
  }
})
# survival_time: survival_time
clinical_df_TCGA_LUSC$survival_time<-sapply(c(1:length(clinical_df$vital_status)), FUN=function(x){
  t=setdiff(c(clinical_df$days_to_last_follow_up[x], clinical_df$days_to_death[x]), na.terms)
  if(length(t)==0){
    return("[Not Available]")
  }
  return(max(as.numeric(t)))
})

# vital_status: vital_status
clinical_df_TCGA_LUSC$vital_status<-clinical_df$vital_status

# Disease-free survival: DFS
clinical_df_TCGA_LUSC$DFS<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  followup.nte<-clinical$clinical_follow_up_v1.0_lusc %>%
    filter(bcr_patient_barcode %in% barcode.cut) %>% 
    filter(new_tumor_event_dx_indicator %in% "YES")
  
  if(dim(followup.nte)[1]>0){
    DFS<-min(followup.nte$new_tumor_event_dx_days_to)
  }
  
  nofollowup.nte<-clinical$clinical_nte_lusc %>% 
    filter(bcr_patient_barcode %in% barcode.cut)
  
  if(dim(nofollowup.nte)[1]>0){
    DFS<-nofollowup.nte$new_tumor_event_dx_days_to
  }
  
  if(!exists("DFS")){
    DFS<-(clinical_df_TCGA_LUSC %>% filter(barcode == s))$survival_time
  }
  return(DFS)
})

# NTE-status
clinical_df_TCGA_LUSC$NTE_status<-sapply(clinical_df_TCGA_LUSC$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  surv.time<-clinical_df_TCGA_LUSC$survival_time[match(s, clinical_df_TCGA_LUSC$barcode)]
  DFS.time<-clinical_df_TCGA_LUSC$DFS[match(s, clinical_df_TCGA_LUSC$barcode)]
  
  if(surv.time!=DFS.time){
    return("YES")
  }
  
  if(surv.time %in% na.terms){
    return("[Not Available]")
  }
  
  if(!is.na(match(barcode.cut, clinical$clinical_nte_lusc$bcr_patient_barcode))){
    return("YES")
  }
  
  indicator<-clinical$clinical_follow_up_v1.0_lusc$new_tumor_event_dx_indicator[
    match(barcode.cut, clinical$clinical_follow_up_v1.0_lusc$bcr_patient_barcode)]
  
  if("YES" %in% indicator){
    return("YES")
  }
  
  if("NO" %in% indicator){
    return("NO")
  }
  
  return("[Not Available]")
})

# ajcc_pathologic_t
clinical_df_TCGA_LUSC$ajcc_pathologic_t<-clinical_df$ajcc_pathologic_t
# ajcc_pathologic_n
clinical_df_TCGA_LUSC$ajcc_pathologic_n<-clinical_df$ajcc_pathologic_n
# ajcc_pathologic_m
clinical_df_TCGA_LUSC$ajcc_pathologic_m<-clinical_df$ajcc_pathologic_m
# ajcc_pathologic_stage
clinical_df_TCGA_LUSC$ajcc_pathologic_stage<-clinical_df$ajcc_pathologic_stage
# Race: race
clinical_df_TCGA_LUSC$Race<-clinical_df$race
# Gender: gender
clinical_df_TCGA_LUSC$Gender<-clinical_df$gender
# Adjuvant_chemo
clinical_df_TCGA_LUSC$Adjuvant_chemo<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Pharmaceutical Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})
# Adjuvant_radiation
clinical_df_TCGA_LUSC$Adjuvant_radiation<-sapply(data@colData@listData$treatments, FUN=function(x){
  res<-x %>% filter(treatment_type=="Radiation Therapy, NOS") %>% dplyr::select(treatment_or_therapy) %>% pull
  if(res=="not reported"){
    res<-"[Not Available]"
  }
  return(res)
})

#  Smoking_status
clinical_df_TCGA_LUSC$Smoking_status<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  
  ind<-match(barcode.cut, clinical$clinical_patient_lusc$bcr_patient_barcode)
  smoking.indicator<-clinical$clinical_patient_lusc$tobacco_smoking_history_indicator[ind]
  if(smoking.indicator == "1"){
    res<-"Never smoker"
  }else if (smoking.indicator %in% c("2", "3", "4", "5")){
    res<-"Ever smoker"
  }else{
    res<-"[Not Available]"
  }
  return(res)
})

# Predominant_subtype
clinical_df_TCGA_LUSC$Predominant_subtype<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, clinical$clinical_patient_lusc$bcr_patient_barcode)
  clinical$clinical_patient_lusc$histologic_diagnosis...19[ind]
})

## Mutation statuses
# EGFR mutation
mut<-maf %>% filter(Hugo_Symbol == "EGFR")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUSC$EGFR_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"EGFR-WT"
  else
    res<-"EGFR-Mut"
  return(res)
})

# KRAS mutation
mut<-maf %>% filter(Hugo_Symbol == "KRAS")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUSC$KRAS_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"KRAS-WT"
  else
    res<-"KRAS-Mut"
  return(res)
})

# TP53 mutation
mut<-maf %>% filter(Hugo_Symbol == "TP53")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUSC$TP53_mut<-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"TP53-WT"
  else
    res<-"TP53-Mut"
  return(res)
})

# STK11 mutation
mut<-maf %>% filter(Hugo_Symbol == "STK11")
mut$barcode.short<-sapply(mut$Tumor_Sample_Barcode, FUN=function(x){
  paste0(unlist(strsplit(x, "-", fixed=T))[1:3], collapse = "-")
})

clinical_df_TCGA_LUSC$STK11_mut <-sapply(clinical_df$barcode, FUN=function(s){
  barcode.cut<-paste0(unlist(strsplit(s, "-"))[1:3], collapse = "-")
  ind<-match(barcode.cut, mut$barcode.short)
  if(is.na(ind))
    res<-"STK11-WT"
  else
    res<-"STK11-Mut"
  return(res)
})

##### Expression data processing
data.lusc <- assay(data, "fpkm_unstrand")
data.lusc.row.info<-rowRanges(data)
protein_coding_gene.ind<-which(data.lusc.row.info@elementMetadata@listData$gene_type=="protein_coding")

data.lusc.protein_coding<-data.lusc[protein_coding_gene.ind, ]

expression_df_lusc<-t(data.lusc.protein_coding)
colnames(expression_df_lusc)<-row.names(data.lusc.protein_coding)
expression_df_lusc<-log(expression_df_lusc+1, base=2)
expression_df_lusc<-expression_df_lusc[,gene_info$gene_id]

clinical_TCGA<-rbind(clinical_df_TCGA_LUAD, clinical_df_TCGA_LUSC)
expression_TCGA<-rbind(expression_df_luad, expression_df_lusc)
colnames(expression_TCGA)<-gene_info$gene_symbol

rm(clinical);rm(clinical_df); rm(clinical_df_TCGA_LUAD); rm(clinical_df_TCGA_LUSC)
rm(expression_df_luad);rm(expression_df_lusc)
rm(data.luad); rm(data.luad.protein_coding); rm(data.luad.row.info)
rm(data.lusc); rm(data.lusc.protein_coding); rm(data.lusc.row.info)
rm(data);rm(maf);rm(query);rm(mut)

# 2. GSE3141
data <- getGEO(filename="GEO/GSE3141_series_matrix.txt.gz")
clinical_df <- data@phenoData@data
labels<-lapply(clinical_df$characteristics_ch1, FUN=function(x){
  strsplit(x, split = ";", fixed=T) %>% unlist
})

# Sample IDs
clinical_GSE3141<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE3141$Histology<-sapply(labels, FUN=function(x){
  unlist(strsplit(x[1], split=": ", fixed=T))[3]
})
# Age: age
clinical_GSE3141$Age<-rep("[Not Available]", length(clinical_df$geo_accession))

# survival_time: survival_time
clinical_GSE3141$survival_time<-sapply(labels, FUN=function(x){
  if(identical(x, labels[[106]])){
    return(72)
  }
  as.numeric(unlist(strsplit(x[2], split=": ", fixed=T))[2])
})
# vital_status: vital_status
clinical_GSE3141$vital_status<-sapply(labels, FUN=function(x){
  unlist(strsplit(x[3], split=": ", fixed=T))[2]
})
# Disease-free survival: DFS
clinical_GSE3141$DFS<-
# NTE-status
clinical_GSE3141$NTE_status<-
# ajcc_pathologic_t
clinical_GSE3141$ajcc_pathologic_t<-
# ajcc_pathologic_n
clinical_GSE3141$ajcc_pathologic_n<-
# ajcc_pathologic_m
clinical_GSE3141$ajcc_pathologic_m<-
# ajcc_pathologic_stage
clinical_GSE3141$ajcc_pathologic_stage<-
# Race: race
clinical_GSE3141$Race<-
# Gender: gender
clinical_GSE3141$Gender<-
# Adjuvant_chemo
clinical_GSE3141$Adjuvant_chemo<-
# Adjuvant_radiation
clinical_GSE3141$Adjuvant_radiation<-
#  Smoking_status
clinical_GSE3141$Smoking_status<-
# Predominant_subtype
clinical_GSE3141$Predominant_subtype<-
## Mutation statuses
# EGFR mutation
clinical_GSE3141$EGFR_mut<-
# KRAS mutation
clinical_GSE3141$KRAS_mut<-
# TP53 mutation
clinical_GSE3141$TP53_mut<-
# STK11 mutation
clinical_GSE3141$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
data.MA<-log2(data.MA)

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE3141<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE3141)<-gene_info$gene_symbol
rownames(expression_GSE3141)<-colnames(data.MA)

# 3. GSE8894
data <- getGEO(filename="GEO/GSE8894_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE8894<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE8894$Histology<-clinical_df$`cell type:ch1`
# Age: age
clinical_GSE8894$Age<-sapply(clinical_df$`age:ch1`, FUN=function(x){
  if(is.na(x)){
    "[Not Available]"
  }else{
    x
  }
})
# Disease-free survival: DFS
clinical_GSE8894$DFS<-clinical_df$`recurrence free survival time (month):ch1`
# NTE-status
clinical_GSE8894$NTE_status<-clinical_df$`status (1=recurrence, 0=non-recurrence):ch1`
# Gender: gender
clinical_GSE8894$Gender<-clinical_df$`gender:ch1`
# survival_time: survival_time
clinical_GSE8894$survival_time<-
  # vital_status: vital_status
  clinical_GSE8894$vital_status<-
  # ajcc_pathologic_t
  clinical_GSE8894$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE8894$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE8894$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE8894$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE8894$Race<-
  # Adjuvant_chemo
  clinical_GSE8894$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE8894$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE8894$Smoking_status<-
  # Predominant_subtype
  clinical_GSE8894$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE8894$EGFR_mut<-
  # KRAS mutation
  clinical_GSE8894$KRAS_mut<-
  # TP53 mutation
  clinical_GSE8894$TP53_mut<-
  # STK11 mutation
  clinical_GSE8894$STK11_mut<-
  rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE8894<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE8894)<-gene_info$gene_symbol
rownames(expression_GSE8894)<-colnames(data.MA)

# 4. GSE10072
data <- getGEO(filename="GEO/GSE10072_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE10072<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE10072$Histology<-clinical_df$source_name_ch1
# Age: age
clinical_GSE10072$Age<-clinical_df$`Age at Diagnosis:ch1`
# ajcc_pathologic_stage
clinical_GSE10072$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
#  Smoking_status
clinical_GSE10072$Smoking_status<-clinical_df$`Cigarette Smoking Status:ch1`
# Gender: gender
clinical_GSE10072$Gender<-clinical_df$`Gender:ch1`
# survival_time: survival_time
clinical_GSE10072$survival_time<-
  # vital_status: vital_status
  clinical_GSE10072$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE10072$DFS<-
  # NTE-status
  clinical_GSE10072$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE10072$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE10072$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE10072$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE10072$Race<-
  # Adjuvant_chemo
  clinical_GSE10072$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE10072$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE10072$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE10072$EGFR_mut<-
  # KRAS mutation
  clinical_GSE10072$KRAS_mut<-
  # TP53 mutation
  clinical_GSE10072$TP53_mut<-
  # STK11 mutation
  clinical_GSE10072$STK11_mut<- rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE10072<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE10072)<-gene_info$gene_symbol
rownames(expression_GSE10072)<-colnames(data.MA)

# 5. GSE10245 
data <- getGEO(filename="GEO/GSE10245_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE10245<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE10245$Histology<-clinical_df$`disease state:ch1`
# Age: age
clinical_GSE10245$Age<-
  # survival_time: survival_time
  clinical_GSE10245$survival_time<-
  # vital_status: vital_status
  clinical_GSE10245$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE10245$DFS<-
  # NTE-status
  clinical_GSE10245$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE10245$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE10245$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE10245$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE10245$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE10245$Race<-
  # Gender: gender
  clinical_GSE10245$Gender<-
  # Adjuvant_chemo
  clinical_GSE10245$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE10245$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE10245$Smoking_status<-
  # Predominant_subtype
  clinical_GSE10245$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE10245$EGFR_mut<-
  # KRAS mutation
  clinical_GSE10245$KRAS_mut<-
  # TP53 mutation
  clinical_GSE10245$TP53_mut<-
  # STK11 mutation
  clinical_GSE10245$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE10245<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE10245)<-gene_info$gene_symbol
rownames(expression_GSE10245)<-colnames(data.MA)

# 6. GSE10445
data <- getGEO(filename="GEO/GSE10445_series_matrix.txt.gz")
clinical_df <- data@phenoData@data
labels<-lapply(clinical_df$characteristics_ch1, FUN=function(x){
  strsplit(x, split = ", ", fixed=T) %>% unlist
})

# Sample IDs
clinical_GSE10445<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE10445$Histology<-rep("ADC or LCC", length(clinical_df$geo_accession))
# Age: age
clinical_GSE10445$Age<-sapply(labels, FUN=function(x) {substr(x[3], nchar(x[3]) - 1, nchar(x[3]))})
# ajcc_pathologic_t
clinical_GSE10445$ajcc_pathologic_t<-rep("T2", length(clinical_df$geo_accession))
# ajcc_pathologic_n
clinical_GSE10445$ajcc_pathologic_n<-rep("N0", length(clinical_df$geo_accession))  
# Gender: gender
clinical_GSE10445$Gender<-sapply(labels, FUN=function(x){x[2]})
# survival_time: survival_time
clinical_GSE10445$survival_time<-
  # vital_status: vital_status
  clinical_GSE10445$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE10445$DFS<-
  # NTE-status
  clinical_GSE10445$NTE_status<-
  # ajcc_pathologic_m
  clinical_GSE10445$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE10445$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE10445$Race<-
  # Adjuvant_chemo
  clinical_GSE10445$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE10445$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE10445$Smoking_status<-
  # Predominant_subtype
  clinical_GSE10445$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE10445$EGFR_mut<-
  # KRAS mutation
  clinical_GSE10445$KRAS_mut<-
  # TP53 mutation
  clinical_GSE10445$TP53_mut<-
  # STK11 mutation
  clinical_GSE10445$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}


clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE10445<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE10445)<-gene_info$gene_symbol
rownames(expression_GSE10445)<-colnames(data.MA)

# 7. GSE10799
data <- getGEO(filename="GEO/GSE10799_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE10799<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE10799$Histology<-clinical_df$source_name_ch1
# Gender: gender
clinical_GSE10799$Gender<-sapply(clinical_df$characteristics_ch1, FUN=function(x){unlist(strsplit(x, " "))[1]})
# Age: age
clinical_GSE10799$Age<-
  # survival_time: survival_time
  clinical_GSE10799$survival_time<-
  # vital_status: vital_status
  clinical_GSE10799$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE10799$DFS<-
  # NTE-status
  clinical_GSE10799$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE10799$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE10799$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE10799$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE10799$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE10799$Race<-
  # Adjuvant_chemo
  clinical_GSE10799$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE10799$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE10799$Smoking_status<-
  # Predominant_subtype
  clinical_GSE10799$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE10799$EGFR_mut<-
  # KRAS mutation
  clinical_GSE10799$KRAS_mut<-
  # TP53 mutation
  clinical_GSE10799$TP53_mut<-
  # STK11 mutation
  clinical_GSE10799$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE10799<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE10799)<-gene_info$gene_symbol
rownames(expression_GSE10799)<-colnames(data.MA)

# 8. GSE11969
data <- getGEO(filename="GEO/GSE11969_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE11969<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE11969$Histology<-clinical_df$`Histology:ch1`
# Age: age
clinical_GSE11969$Age<-sapply(clinical_df$`Age:ch1`, FUN=function(x){
  if(x=="NA"){
    "[Not Available]"
  }else{
    x
  }
})
# survival_time: survival_time
clinical_GSE11969$survival_time<-sapply(clinical_df$`Survival (days):ch1`, FUN=function(x){
  if(x=="NA"){
    "[Not Available]"
  }else{
    x
  }
})
# vital_status: vital_status
clinical_GSE11969$vital_status<-sapply(clinical_df$`Status:ch1`, FUN=function(x){
  if(x=="NA"){
    "[Not Available]"
  }else{
    x
  }
})
# ajcc_pathologic_t
clinical_GSE11969$ajcc_pathologic_t<-sapply(str_sub(clinical_df$`pTNM:ch1`,1,2), FUN=function(x){
  if(x %in% c("T1","T2", "T3", "T4", "TX")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_n
clinical_GSE11969$ajcc_pathologic_n<-sapply(str_sub(clinical_df$`pTNM:ch1`,3,4), FUN=function(x){
  if(x %in% c("NX", "N0", "N1", "N2", "N3")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_m
clinical_GSE11969$ajcc_pathologic_m<-sapply(str_sub(clinical_df$`pTNM:ch1`,5,6), FUN=function(x){
  if(x %in% c("M0", "M1", "MX")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_stage
clinical_GSE11969$ajcc_pathologic_stage<-sapply(clinical_df$`pStage:ch1`, FUN=function(x){
  if(x=="Stage: 1"){
    "Stage I"
  }else if(x=="IA"){
    "Stage IA"
  }else if(x=="IB"){
    "Stage IB"
  }else if(x=="IIA"){
    "Stage IIA"
  }else if(x=="IIB"){
    "Stage IIB"
  }else if(x=="IIIA"){
    "Stage IIIA"
  }else if(x=="IIIB"){
    "Stage IIIB"
  }else if(x=="IV"){
    "Stage IV"
  }else{
    "[Not Available]"
  }
})
# Gender: gender
clinical_GSE11969$Gender<-sapply(clinical_df$characteristics_ch1.3, FUN=function(x){
  if(x=="Sex: F"){
    "female"
  }else if(x=="Sex: M"){
    "male"
  }else{
    "[Not Available]"
  }
})
# NTE-status
clinical_GSE11969$NTE_status<-clinical_df$characteristics_ch1.10
#  Smoking_status
clinical_GSE11969$Smoking_status<-sapply(clinical_df$`Smoking (BI):ch1`, FUN=function(x){
  if(x=="0"){
    "Never smoker"
  }else if(x=="NA"){
    "[Not Available]"
  }else{
    "Ever smoker"
  }
})
## Mutation statuses
# EGFR mutation
clinical_GSE11969$EGFR_mut<-sapply(clinical_df$`EGFR status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "EGFR-WT"
  }else if(x=="Mut"){
    "EGFR-Mut"
  }else{
    "[Not Available]"
  }
})
# KRAS mutation
clinical_GSE11969$KRAS_mut<-sapply(clinical_df$`K-ras Status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "KRAS-WT"
  }else if (x=="Mut"){
    "KRAS-Mut"
  }else{
    "[Not Available]"
  }
})
# TP53 mutation
clinical_GSE11969$TP53_mut<-sapply(clinical_df$`p53 Status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "p53-WT"
  }else if(x=="Mut"){
    "p53-Mut"
  }else{
    "[Not Available]"
  }
})

# Disease-free survival: DFS
clinical_GSE11969$DFS<-
  # Race: race
  clinical_GSE11969$Race<-
  # Adjuvant_chemo
  clinical_GSE11969$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE11969$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE11969$Predominant_subtype<-
  # STK11 mutation
  clinical_GSE11969$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$`Gene symbol`[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[mapping.table$ID,]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)

write.csv(rownames(data.MA), 'gene_GSE11969.csv', quote = F, row.names = F)
clusterEvalQ(clu, {
  symbol.check<-read.csv("gene_GSE11969_update.csv", sep = ",")
  NULL
})

rownames(data.MA)<-pbsapply(rownames(data.MA), FUN=function(x){
  chk<-symbol.check%>% filter(Input %in% x)
  
  if(length(which(chk$Match.type=="Approved symbol"))>0){
    return(x)
  }else{
    return(chk$Approved.symbol[1])
  }
},cl=clu)

expression_GSE11969<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$gene_symbol[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))

expression_GSE11969<-t(apply(expression_GSE11969, MARGIN=1, FUN=function(x){
  x_<-na.omit(x)
  return((x- min(x_)) /(max(x_)-min(x_)))
})*15)
colnames(expression_GSE11969)<-gene_info$gene_symbol
rownames(expression_GSE11969)<-colnames(data.MA)

# 9. GSE12667
data <- getGEO(filename="GEO/GSE12667_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE12667<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE12667$Histology<-clinical_df$`Pathology:ch1`
# Race: race
clinical_GSE12667$Race<-clinical_df$`Race:ch1`
# Gender: gender
clinical_GSE12667$Gender<-clinical_df$`Gender:ch1`
#  Smoking_status
clinical_GSE12667$Smoking_status<-clinical_df$`Smoking_Status:ch1`
# Age: age
clinical_GSE12667$Age<-
  # survival_time: survival_time
  clinical_GSE12667$survival_time<-
  # vital_status: vital_status
  clinical_GSE12667$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE12667$DFS<-
  # NTE-status
  clinical_GSE12667$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE12667$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE12667$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE12667$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE12667$ajcc_pathologic_stage<-
  # Adjuvant_chemo
  clinical_GSE12667$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE12667$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE12667$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE12667$EGFR_mut<-
  # KRAS mutation
  clinical_GSE12667$KRAS_mut<-
  # TP53 mutation
  clinical_GSE12667$TP53_mut<-
  # STK11 mutation
  clinical_GSE12667$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE12667<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE12667)<-gene_info$gene_symbol
rownames(expression_GSE12667)<-colnames(data.MA)

# 10. GSE13213
data <- getGEO(filename="GEO/GSE13213_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE13213<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE13213$Histology<-clinical_df$`Histology:ch1`
# Age: age
clinical_GSE13213$Age<-clinical_df$`Age:ch1`
# survival_time: survival_time
clinical_GSE13213$survival_time<-clinical_df$`Survival (days):ch1`
# vital_status: vital_status
clinical_GSE13213$vital_status<-clinical_df$`Status:ch1`
# NTE-status
clinical_GSE13213$NTE_status<-clinical_df$`Evidence of relapse:ch1`
# ajcc_pathologic_t
clinical_GSE13213$ajcc_pathologic_t<-sapply(str_sub(clinical_df$`TNM (Pathological):ch1`,1,2), FUN=function(x){
  if(x %in% c("T1","T2", "T3", "T4", "TX")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_n
clinical_GSE13213$ajcc_pathologic_n<-sapply(str_sub(clinical_df$`TNM (Pathological):ch1`,3,4), FUN=function(x){
  if(x %in% c("NX", "N0", "N1", "N2", "N3")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_m
clinical_GSE13213$ajcc_pathologic_m<-sapply(str_sub(clinical_df$`TNM (Pathological):ch1`,5,6), FUN=function(x){
  if(x %in% c("M0", "M1", "MX")){
    x
  }else{
    "[Not Available]"
  }
})
# ajcc_pathologic_stage
clinical_GSE13213$ajcc_pathologic_stage<-sapply(clinical_df$`Stage (Pathological ):ch1`, FUN=function(x){
  if(x=="IA"){
    "Stage IA"
  }else if(x=="IB"){
    "Stage IB"
  }else if(x=="IIA"){
    "Stage IIA"
  }else if(x=="IIB"){
    "Stage IIB"
  }else if(x=="IIIA"){
    "Stage IIIA"
  }else if(x=="IIIB"){
    "Stage IIIB"
  }else if(x=="IV"){
    "Stage IV"
  }else{
    "[Not Available]"
  }
})
# Gender: gender
clinical_GSE13213$Gender<-clinical_df$`Sex:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE13213$EGFR_mut<-sapply(clinical_df$`EGFR status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "EGFR-WT"
  }else if(x=="Mut"){
    "EGFR-Mut"
  }else{
    "[Not Available]"
  }
})
# KRAS mutation
clinical_GSE13213$KRAS_mut<-sapply(clinical_df$`K-ras Status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "KRAS-WT"
  }else if (x=="Mut"){
    "KRAS-Mut"
  }else{
    "[Not Available]"
  }
})
# TP53 mutation
clinical_GSE13213$TP53_mut<-sapply(clinical_df$`p53 Status:ch1`, FUN=function(x){
  if(x=="Wt"){
    "p53-WT"
  }else if(x=="Mut"){
    "p53-Mut"
  }else{
    "[Not Available]"
  }
})
#  Smoking_status
clinical_GSE13213$Smoking_status<-sapply(clinical_df$`Smoking (BI):ch1`, FUN=function(x){
  if(x=="0"){
    "Never smoker"
  }else if(x=="NA"){
    "[Not Available]"
  }else{
    "Ever smoker"
  }
})
# Disease-free survival: DFS
clinical_GSE13213$DFS<-
  # Race: race
  clinical_GSE13213$Race<-
  # Adjuvant_chemo
  clinical_GSE13213$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE13213$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE13213$Predominant_subtype<-
  # STK11 mutation
  clinical_GSE13213$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Removing probes which contains NA values
na.ind<-which(apply(data.MA,MARGIN=1, FUN=function(x){
  length(which(is.na(x)))
})>0)
data.MA<-data.MA[-na.ind,]
data.Platform<-data.Platform[-na.ind,]

# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$GENE[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(na.omit(data.MA[entrez.to.id[[k]],])))
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE13213<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
expression_GSE13213<-t(apply(expression_GSE13213, MARGIN=1, FUN=function(x){
  x_<-na.omit(x)
  return((x- min(x_)) /(max(x_)-min(x_)))
})*15)
colnames(expression_GSE13213)<-gene_info$gene_symbol
rownames(expression_GSE13213)<-colnames(data.MA)

# 11. GSE14814
data <- getGEO(filename="GEO/GSE14814_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE14814<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE14814$Histology<-clinical_df$`Histology type:ch1`
# Age: age
clinical_GSE14814$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE14814$survival_time<-clinical_df$`OS time:ch1`
# vital_status: vital_status
clinical_GSE14814$vital_status<-clinical_df$`OS status:ch1`
# ajcc_pathologic_stage
clinical_GSE14814$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Gender: gender
clinical_GSE14814$Gender<-clinical_df$`Sex:ch1`
# Adjuvant_chemo
clinical_GSE14814$Adjuvant_chemo<-sapply(clinical_df$`Post Surgical Treatment:ch1`, FUN=function(x){
  if(x=="OBS"){
    "no"
  }else if(x=="ACT"){
    "yes"
  }else{
    "[Not Available]"
  }
})
# Predominant_subtype
clinical_GSE14814$Predominant_subtype<-clinical_df$`predominant subtype:ch1`
# Disease-free survival: DFS
clinical_GSE14814$DFS<-
  # NTE-status
  clinical_GSE14814$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE14814$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE14814$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE14814$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE14814$Race<-
  # Adjuvant_radiation
  clinical_GSE14814$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE14814$Smoking_status<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE14814$EGFR_mut<-
  # KRAS mutation
  clinical_GSE14814$KRAS_mut<-
  # TP53 mutation
  clinical_GSE14814$TP53_mut<-
  # STK11 mutation
  clinical_GSE14814$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
na.ind<-which(apply(data.MA, MARGIN=1, FUN=function(x){length(which(is.na(x)))})>0)
data.MA<-data.MA[-na.ind,]
mapping.table<-mapping.table[-na.ind, ]

if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE14814<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE14814)<-gene_info$gene_symbol
rownames(expression_GSE14814)<-colnames(data.MA)

# 12 GSE18842
data <- getGEO(filename="GEO/GSE18842_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE18842<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE18842$Histology<-clinical_df$`sample type:ch1`
# Age: age
clinical_GSE18842$Age<-
  # survival_time: survival_time
  clinical_GSE18842$survival_time<-
  # vital_status: vital_status
  clinical_GSE18842$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE18842$DFS<-
  # NTE-status
  clinical_GSE18842$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE18842$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE18842$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE18842$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE18842$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE18842$Race<-
  # Gender: gender
  clinical_GSE18842$Gender<-
  # Adjuvant_chemo
  clinical_GSE18842$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE18842$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE18842$Smoking_status<-
  # Predominant_subtype
  clinical_GSE18842$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE18842$EGFR_mut<-
  # KRAS mutation
  clinical_GSE18842$KRAS_mut<-
  # TP53 mutation
  clinical_GSE18842$TP53_mut<-
  # STK11 mutation
  clinical_GSE18842$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE18842<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE18842)<-gene_info$gene_symbol
rownames(expression_GSE18842)<-colnames(data.MA)

# 13. GSE19188
data <- getGEO(filename="GEO/GSE19188_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE19188<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE19188$Histology<-clinical_df$`cell type:ch1`
# survival_time: survival_time
clinical_GSE19188$survival_time<-sapply(clinical_df$`overall survival:ch1`, FUN=function(x){
  if(x=="Not available"){
    return("[Not Available]")
  }else{
    return(x)
  }
})
# vital_status: vital_status
clinical_GSE19188$vital_status<-clinical_df$`status:ch1`
# Gender: gender
clinical_GSE19188$Gender<-clinical_df$`gender:ch1`
# Age: age
clinical_GSE19188$Age<-
  # Disease-free survival: DFS
  clinical_GSE19188$DFS<-
  # NTE-status
  clinical_GSE19188$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE19188$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE19188$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE19188$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE19188$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE19188$Race<-
  # Adjuvant_chemo
  clinical_GSE19188$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE19188$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE19188$Smoking_status<-
  # Predominant_subtype
  clinical_GSE19188$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE19188$EGFR_mut<-
  # KRAS mutation
  clinical_GSE19188$KRAS_mut<-
  # TP53 mutation
  clinical_GSE19188$TP53_mut<-
  # STK11 mutation
  clinical_GSE19188$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA<-readLines('GSE19188.txt')
samples<- strsplit(data.MA[1], "\t")[[1]]
samples <- grep("^\"GSM", samples, value = TRUE)
samples <- gsub("\"", "", samples)

data.MA<-pblapply(data.MA[-1], FUN=function(line){
  numbers <- str_extract_all(line, "-?\\d+\\.?\\d*")[[1]]  
  # 숫자를 numeric 타입으로 변환
  numbers <- as.numeric(numbers)
  
  list(gene=numbers[1], expr=numbers[-1])
})

genes<-sapply(data.MA, FUN=function(x){x$gene})
data.MA<-do.call(rbind, lapply(data.MA, FUN=function(x){x$expr}))
rownames(data.MA)<-genes
colnames(data.MA)<-samples

clinical_GSE19188<-clinical_GSE19188[-2,]
data.MA<-do.call(cbind, lapply(clinical_GSE19188$ID, FUN=function(sample){
  return(data.MA[,match(sample, colnames(data.MA))])
}))
colnames(data.MA)<-clinical_GSE19188$ID

expression_GSE19188<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE19188)<-gene_info$gene_symbol
rownames(expression_GSE19188)<-colnames(data.MA)

# 14. GSE19804
data <- getGEO(filename="GEO/GSE19804_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE19804<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE19804$Histology<-clinical_df$`tissue:ch1`
# Age: age
clinical_GSE19804$Age<-clinical_df$`age:ch1`
# ajcc_pathologic_stage
clinical_GSE19804$ajcc_pathologic_stage<-clinical_df$`stage:ch1`
# Gender: gender
clinical_GSE19804$Gender<-clinical_df$`gender:ch1`

# survival_time: survival_time
clinical_GSE19804$survival_time<-
  # vital_status: vital_status
  clinical_GSE19804$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE19804$DFS<-
  # NTE-status
  clinical_GSE19804$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE19804$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE19804$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE19804$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE19804$Race<-
  # Adjuvant_chemo
  clinical_GSE19804$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE19804$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE19804$Smoking_status<-
  # Predominant_subtype
  clinical_GSE19804$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE19804$EGFR_mut<-
  # KRAS mutation
  clinical_GSE19804$KRAS_mut<-
  # TP53 mutation
  clinical_GSE19804$TP53_mut<-
  # STK11 mutation
  clinical_GSE19804$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))


data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE19804<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE19804)<-gene_info$gene_symbol
rownames(expression_GSE19804)<-colnames(data.MA)

# 15. GSE26939
data <- getGEO(filename="GEO/GSE26939_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE26939<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE26939$Histology<-rep("ADC", length(clinical_df$geo_accession))
# Age: age
clinical_GSE26939$Age<-clinical_df$`age (90=greater than or equal to 90):ch1`
# survival_time: survival_time
clinical_GSE26939$survival_time<-sapply(clinical_df$`survival_months:ch1`, FUN=function(x){
  as.numeric(x)
})
# vital_status: vital_status
clinical_GSE26939$vital_status<-clinical_df$`survival_status(0='alive',1='dead'):ch1`
# ajcc_pathologic_stage
clinical_GSE26939$ajcc_pathologic_stage<-sapply(c(1:dim(clinical_df)[1]), FUN=function(k){
  x=setdiff(c(clinical_df$`stage:ch1`[k], clinical_df$`Stage:ch1`[k]), NA)
  if(length(x)==0){
    return("[Not Available]")
  }
  if(x=="IA"){
    "Stage IA"
  }else if(x=="IB"){
    "Stage IB"
  }else if(x=="IIA"){
    "Stage IIA"
  }else if(x=="IIB"){
    "Stage IIB"
  }else if(x=="IIIA"){
    "Stage IIIA"
  }else if(x=="IIIB"){
    "Stage IIIB"
  }else if(x=="IV"){
    "Stage IV"
  }else{
    "[Not Available]"
  }
})
# Gender: gender
clinical_GSE26939$Gender<-sapply(c(1:dim(clinical_df)[1]), FUN=function(k){
  x=setdiff(c(clinical_df$`sex:ch1`[k], clinical_df$`Sex:ch1`[k]), NA)
  if(x=="F"){
    "female"
  }else if(x=="M"){
    "male"
  }else{
    "[Not Available]"
  }
})
#  Smoking_status
clinical_GSE26939$Smoking_status<-sapply(clinical_df$`smoking_status(0=nonsmoker,1=smoker):ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }
  if(x=="Smoker"){
    "Ever smoker"
  }else if(x=="NeverSmoker"){
    "Never smoker"
  }else{
    "[Not Available]"
  }
})
# Predominant_subtype
clinical_GSE26939$Predominant_subtype<-clinical_df$`subtype:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE26939$EGFR_mut<-clinical_df$`egfr (0='wt',1='mutated'):ch1`
# KRAS mutation
clinical_GSE26939$KRAS_mut<-clinical_df$`kras (0='wt',1='mutated'):ch1`
# TP53 mutation
clinical_GSE26939$TP53_mut<-clinical_df$`tp53 (0='wt',1='mutated'):ch1`
# STK11 mutation
clinical_GSE26939$STK11_mut<-clinical_df$`stk11 (0='wt',1='mutated'):ch1`
# Disease-free survival: DFS
clinical_GSE26939$DFS<-
  # NTE-status
  clinical_GSE26939$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE26939$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE26939$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE26939$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE26939$Race<-
  # Adjuvant_chemo
  clinical_GSE26939$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE26939$Adjuvant_radiation<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Symbol
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Symbol=unlist(strsplit(data.Platform$ORF[k], " ",fixed = T))
  if(length(Symbol)!=0){
    df<-data.frame(ID, Symbol)  
  }else{
    df<-data.frame(ID=ID, Symbol=NA)
  }
  return(df)
})) %>% filter(!is.na(Symbol))

data.MA<-data.MA[match(mapping.table$ID, rownames(data.MA)),]

symbol.to.id<-pblapply(unique(mapping.table$Symbol), FUN=function(geneid){
  indices=which(mapping.table$Symbol %in% geneid)
})
names(symbol.to.id)<-unique(mapping.table$Symbol)

data.MA<-do.call(rbind, pblapply(c(1:length(symbol.to.id)), FUN=function(k){
  if(length(symbol.to.id[[k]])>1){
    return(colMeans(data.MA[symbol.to.id[[k]],]))  
  }else{
    return(data.MA[symbol.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(symbol.to.id)

write.csv(rownames(data.MA), 'gene_GSE26939.csv', quote = F, row.names = F)

clusterEvalQ(clu, {
  symbol.check<-read.csv("gene_GSE26939_update.csv", sep = ",")
  NULL
})

rownames(data.MA)<-pbsapply(rownames(data.MA), FUN=function(x){
  chk<-symbol.check%>% filter(Input %in% x)
  
  if(length(which(chk$Match.type=="Approved symbol"))>0){
    return(x)
  }else{
    return(chk$Approved.symbol[1])
  }
},cl=clu)

expression_GSE26939<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$gene_symbol[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))

expression_GSE26939<-t(apply(expression_GSE26939, MARGIN=1, FUN=function(x){
  x_<-na.omit(x)
  return((x- min(x_)) /(max(x_)-min(x_)))
})*15)
colnames(expression_GSE26939)<-gene_info$gene_symbol
rownames(expression_GSE26939)<-colnames(data.MA)

# 16. GSE28571
data <- getGEO(filename="GEO/GSE28571_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE28571<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE28571$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE28571$Age<-
  # survival_time: survival_time
  clinical_GSE28571$survival_time<-
  # vital_status: vital_status
  clinical_GSE28571$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE28571$DFS<-
  # NTE-status
  clinical_GSE28571$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE28571$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE28571$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE28571$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE28571$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE28571$Race<-
  # Gender: gender
  clinical_GSE28571$Gender<-
  # Adjuvant_chemo
  clinical_GSE28571$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE28571$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE28571$Smoking_status<-
  # Predominant_subtype
  clinical_GSE28571$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE28571$EGFR_mut<-
  # KRAS mutation
  clinical_GSE28571$KRAS_mut<-
  # TP53 mutation
  clinical_GSE28571$TP53_mut<-
  # STK11 mutation
  clinical_GSE28571$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE28571<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE28571)<-gene_info$gene_symbol
rownames(expression_GSE28571)<-colnames(data.MA)

# 17. GSE29016
data <- getGEO(filename="GEO/GSE29016_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE29016<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE29016$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE29016$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE29016$survival_time<-sapply(clinical_df$`os:ch1`, FUN=function(x){
  if(x=="NA"){
    return("[Not Available]")
  }else{
    return(as.numeric(x)*365.25)
  }
})
# vital_status: vital_status
clinical_GSE29016$vital_status<-sapply(clinical_df$`osbin:ch1`, FUN=function(x){
  if(x=="NA"){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# ajcc_pathologic_t
clinical_GSE29016$ajcc_pathologic_t<-sapply(clinical_df$`tnm:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 1, 3), substr(tnm, 1, 2)))
})
# ajcc_pathologic_n
clinical_GSE29016$ajcc_pathologic_n<-sapply(clinical_df$`tnm:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 4, 5), substr(tnm, 3, 4)))
})
# ajcc_pathologic_m
clinical_GSE29016$ajcc_pathologic_m<-sapply(clinical_df$`tnm:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 6, 7), substr(tnm, 5, 6)))
})
# ajcc_pathologic_stage
clinical_GSE29016$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Gender: gender
clinical_GSE29016$Gender<-clinical_df$`Sex:ch1`
#  Smoking_status
clinical_GSE29016$Smoking_status<-clinical_df$`smoking:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE29016$EGFR_mut<-clinical_df$`egfr_mut:ch1`
# KRAS mutation
clinical_GSE29016$KRAS_mut<-clinical_df$`kras_mut:ch1`

# Disease-free survival: DFS
clinical_GSE29016$DFS<-
  # NTE-status
  clinical_GSE29016$NTE_status<-
  # Race: race
  clinical_GSE29016$Race<-
  # Adjuvant_chemo
  clinical_GSE29016$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE29016$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE29016$Predominant_subtype<-
  # TP53 mutation
  clinical_GSE29016$TP53_mut<-
  # STK11 mutation
  clinical_GSE29016$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE29016<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE29016)<-gene_info$gene_symbol
rownames(expression_GSE29016)<-colnames(data.MA)

# 18. GSE30219
data <- getGEO(filename="GEO/GSE30219_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE30219<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE30219$Histology<-clinical_df$`histology:ch1`
clinical_GSE30219$Histology[which(clinical_GSE30219$Histology=="SCC")]<-"SCLC"
# Age: age
clinical_GSE30219$Age<-sapply(clinical_df$`age at surgery:ch1`, FUN=function(x){
  if(x%in%c("ND", "NTL")){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# survival_time: survival_time
clinical_GSE30219$survival_time<-sapply(clinical_df$`follow-up time (months):ch1`, FUN=function(x){
  if(x%in%c("ND", "NTL")){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# vital_status: vital_status
clinical_GSE30219$vital_status<-clinical_df$`status:ch1`
# Disease-free survival: DFS
clinical_GSE30219$DFS<-sapply(clinical_df$`disease free survival in months:ch1`, FUN=function(x){
  if(x%in%c("na")||is.na(x)){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# NTE-status
clinical_GSE30219$NTE_status<-clinical_df$`relapse (event=1; no event=0):ch1`
# ajcc_pathologic_t
clinical_GSE30219$ajcc_pathologic_t<-clinical_df$`pt stage:ch1`
# ajcc_pathologic_n
clinical_GSE30219$ajcc_pathologic_n<-clinical_df$`pn stage:ch1`
# ajcc_pathologic_m
clinical_GSE30219$ajcc_pathologic_m<-clinical_df$`pm stage:ch1`
# Gender: gender
clinical_GSE30219$Gender<-clinical_df$`gender:ch1`
# ajcc_pathologic_stage
clinical_GSE30219$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE30219$Race<-
  # Adjuvant_chemo
  clinical_GSE30219$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE30219$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE30219$Smoking_status<-
  # Predominant_subtype
  clinical_GSE30219$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE30219$EGFR_mut<-
  # KRAS mutation
  clinical_GSE30219$KRAS_mut<-
  # TP53 mutation
  clinical_GSE30219$TP53_mut<-
  # STK11 mutation
  clinical_GSE30219$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE30219<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE30219)<-gene_info$gene_symbol
rownames(expression_GSE30219)<-colnames(data.MA)

# 19. GSE31210
data <- getGEO(filename="GEO/GSE31210_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE31210<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE31210$Histology<-sapply(clinical_df$`tissue:ch1`, FUN=function(x){
  if(x=="primary lung tumor"){
    "ADC"
  }else{
    "Control"
  }
})
# Age: age
clinical_GSE31210$Age<-sapply(1:length(clinical_df$`age (years):ch1`), FUN=function(k){
  if(is.na(clinical_df$`age (years):ch1`[k])){
    return(clinical_df$`age:ch1`[k])
  }else{
    return(clinical_df$`age (years):ch1`[k])
  }
})
# survival_time: survival_time
clinical_GSE31210$survival_time<-sapply(clinical_df$`days before death/censor:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# vital_status: vital_status
clinical_GSE31210$vital_status<-clinical_df$`death:ch1` 
# Disease-free survival: DFS
clinical_GSE31210$DFS<-sapply(clinical_df$`days before relapse/censor:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# NTE-status
clinical_GSE31210$NTE_status<-clinical_df$`relapse:ch1`
# ajcc_pathologic_stage
clinical_GSE31210$ajcc_pathologic_stage<-clinical_df$`pathological stage:ch1`
# Gender: gender
clinical_GSE31210$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE31210$Smoking_status<-clinical_df$`smoking status:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE31210$EGFR_mut<-sapply(clinical_df$`gene alteration status:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }
  if(x=="EGFR mutation +"){
    "EGFR-Mut"
  }else{
    "EGFR-WT"
  }
})
# KRAS mutation
clinical_GSE31210$KRAS_mut<-sapply(clinical_df$`gene alteration status:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }
  if(x=="KRAS mutation +"){
    "KRAS-Mut"
  }else{
    "KRAS-WT"
  }
})
# ajcc_pathologic_t
clinical_GSE31210$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE31210$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE31210$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE31210$Race<-
  # Adjuvant_chemo
  clinical_GSE31210$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE31210$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE31210$Predominant_subtype<-
  # TP53 mutation
  clinical_GSE31210$TP53_mut<-
  # STK11 mutation
  clinical_GSE31210$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE31210<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE31210)<-gene_info$gene_symbol
rownames(expression_GSE31210)<-colnames(data.MA)

# 20. GSE32863
data <- getGEO(filename="GEO/GSE32863_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE32863<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE32863$Histology<-clinical_df$`type:ch1`
# Age: age
clinical_GSE32863$Age<-clinical_df$`age:ch1`
# NTE-status
clinical_GSE32863$NTE_status<-clinical_df$`recurrence:ch1`
# ajcc_pathologic_stage
clinical_GSE32863$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Race: race
clinical_GSE32863$Race<-clinical_df$`ethnicity:ch1`
# Gender: gender
clinical_GSE32863$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE32863$Smoking_status<-clinical_df$`smoking status:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE32863$EGFR_mut<-clinical_df$`egfr mutation status:ch1`
# KRAS mutation
clinical_GSE32863$KRAS_mut<-clinical_df$`kras mutation status:ch1`

# survival_time: survival_time
clinical_GSE32863$survival_time<-
  # vital_status: vital_status
  clinical_GSE32863$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE32863$DFS<-
  # ajcc_pathologic_t
  clinical_GSE32863$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE32863$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE32863$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE32863$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE32863$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE32863$Predominant_subtype<-
  # TP53 mutation
  clinical_GSE32863$TP53_mut<-
  # STK11 mutation
  clinical_GSE32863$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE32863<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE32863)<-gene_info$gene_symbol
rownames(expression_GSE32863)<-colnames(data.MA)

# 21. GSE37745
data <- getGEO(filename="GEO/GSE37745_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE37745<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE37745$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE37745$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE37745$survival_time<-clinical_df$`days to determined death status:ch1`
# vital_status: vital_status
clinical_GSE37745$vital_status<-clinical_df$`dead:ch1`
# Disease-free survival: DFS
clinical_GSE37745$DFS<-sapply(clinical_df$`days to recurrence / to last visit:ch1`, FUN=function(x){
  if(x=="not known"){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# NTE-status
clinical_GSE37745$NTE_status<-clinical_df$`recurrence:ch1`
# ajcc_pathologic_stage
clinical_GSE37745$ajcc_pathologic_stage<-clinical_df$`tumor stage:ch1`
# Adjuvant_chemo
clinical_GSE37745$Adjuvant_chemo<-clinical_df$`adjuvant treatment:ch1`
# Gender: gender
clinical_GSE37745$Gender<-clinical_df$`gender:ch1`

# ajcc_pathologic_t
clinical_GSE37745$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE37745$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE37745$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE37745$Race<-
  # Adjuvant_radiation
  clinical_GSE37745$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE37745$Smoking_status<-
  # Predominant_subtype
  clinical_GSE37745$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE37745$EGFR_mut<-
  # KRAS mutation
  clinical_GSE37745$KRAS_mut<-
  # TP53 mutation
  clinical_GSE37745$TP53_mut<-
  # STK11 mutation
  clinical_GSE37745$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE37745<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE37745)<-gene_info$gene_symbol
rownames(expression_GSE37745)<-colnames(data.MA)

# 22.GSE40419
data <- getGEO(filename="GEO/GSE40419_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE40419<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE40419$Histology<-sapply(clinical_df$source_name_ch1, FUN=function(x){
  if(x=="Lung cancer cells"){
    return("ADC")
  }else{
    return("Control")
  }
})
# Age: age
clinical_GSE40419$Age<-clinical_df$`age_at_diagnosis:ch1`
# ajcc_pathologic_stage
clinical_GSE40419$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Gender: gender
clinical_GSE40419$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE40419$Smoking_status<-clinical_df$`smoking_status:ch1`

# survival_time: survival_time
clinical_GSE40419$survival_time<-
  # vital_status: vital_status
  clinical_GSE40419$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE40419$DFS<-
  # NTE-status
  clinical_GSE40419$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE40419$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE40419$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE40419$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE40419$Race<-
  # Adjuvant_chemo
  clinical_GSE40419$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE40419$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE40419$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE40419$EGFR_mut<-
  # KRAS mutation
  clinical_GSE40419$KRAS_mut<-
  # TP53 mutation
  clinical_GSE40419$TP53_mut<-
  # STK11 mutation
  clinical_GSE40419$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA<-read.table("GEO/GSE40419_LC-87_RPKM_expression.txt", sep="\t", header=T)
genes<-data.MA$gene
samples<-colnames(data.MA)[-c(1:6)]
data.MA<-data.MA[,-c(1:6)]
data.MA<-do.call(cbind, lapply(1:dim(data.MA)[2], FUN=function(k){
  log2(as.numeric(data.MA[,k])+1)
}))

clusterExport(clu, "genes")
entrez.to.id<-pblapply(unique(genes), FUN=function(geneid){
  indices=which(genes %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(genes)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

write.csv(unique(genes), "gene_GSE40419.csv",quote = FALSE, row.names = F)

clusterEvalQ(clu, {
  symbol.check<-read.csv("gene_GSE40419_update.csv", sep = ",")
  NULL
})

rownames(data.MA)<-pbsapply(rownames(data.MA), FUN=function(x){
  chk<-symbol.check%>% filter(Input %in% x)
  
  if(length(which(chk$Match.type=="Approved symbol"))>0){
    return(x)
  }else{
    return(chk$Approved.symbol[1])
  }
},cl=clu)

expression_GSE40419<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$gene_symbol[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))

expression_GSE40419<-do.call(rbind, lapply(clinical_df$title, FUN=function(x){
  ind<-match(x, samples)
  expression_GSE40419[ind,]
}))
colnames(expression_GSE40419)<-gene_info$gene_symbol
rownames(expression_GSE40419)<-clinical_df$geo_accession


# 23. GSE41271
data <- getGEO(filename="GEO/GSE41271_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE41271<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE41271$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE41271$Age<-sapply(c(1:dim(clinical_df)[1]), FUN=function(k){
  surgery.date <- as.Date(clinical_df$`date of surgery:ch1`[k], format = "%Y-%m-%d")
  birth.date.raw <- clinical_df$`date of birth:ch1`[k]
  if (grepl("/", birth.date.raw)) {
    birth.date <- as.Date(birth.date.raw, format = "%m/%d/%Y")
  } else {
    birth.date <- as.Date(birth.date.raw, format = "%Y-%m-%d")
  }
  # 생존 일수 계산
  survival.days <- as.numeric(surgery.date - birth.date)
  # NA인 경우 처리
  if (is.na(survival.days)) {
    return("[Not Available]")
  } else {
    return(round(survival.days / 365.25))  # 생존 기간을 연도로 변환
  }
})
# survival_time: survival_time
clinical_GSE41271$survival_time<-sapply(c(1:dim(clinical_df)[1]), FUN=function(k){
  surgery.date<-as_date(clinical_df$`date of surgery:ch1`[k])
  last.followup.date<-as_date(clinical_df$`last follow-up survival:ch1`)[k]
  survival.days<-as.numeric(last.followup.date-surgery.date)
  if(is.na(survival.days)){
    return("[Not Available]")
  }else{
    return(survival.days)  
  }
})
# vital_status: vital_status
clinical_GSE41271$vital_status<-clinical_df$`vital statistics:ch1`
# Disease-free survival: DFS
clinical_GSE41271$DFS<-sapply(c(1:dim(clinical_df)[1]), FUN=function(k){
  surgery.date<-as_date(clinical_df$`date of surgery:ch1`[k])
  last.followup.date<-as_date(clinical_df$`last follow-up recurrence:ch1`)[k]
  survival.days<-as.numeric(last.followup.date-surgery.date)
  if(is.na(survival.days)){
    return("[Not Available]")
  }else{
    return(survival.days)  
  }
})
# NTE-status
clinical_GSE41271$NTE_status<-clinical_df$`recurrence:ch1`
# ajcc_pathologic_stage
clinical_GSE41271$ajcc_pathologic_stage<-clinical_df$`final patient stage:ch1`
# Race: race
clinical_GSE41271$Race<-clinical_df$`race:ch1`
# Gender: gender
clinical_GSE41271$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE41271$Smoking_status<-clinical_df$`tobacco history:ch1`

# ajcc_pathologic_t
clinical_GSE41271$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE41271$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE41271$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE41271$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE41271$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE41271$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE41271$EGFR_mut<-
  # KRAS mutation
  clinical_GSE41271$KRAS_mut<-
  # TP53 mutation
  clinical_GSE41271$TP53_mut<-
  # STK11 mutation
  clinical_GSE41271$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE41271<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE41271)<-gene_info$gene_symbol
rownames(expression_GSE41271)<-colnames(data.MA)

# 24. GSE42127
data <- getGEO(filename="GEO/GSE42127_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE42127<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE42127$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE42127$Age<-clinical_df$`age at surgery:ch1`
# survival_time: survival_time
clinical_GSE42127$survival_time<-clinical_df$`overall survival months:ch1`
# vital_status: vital_status
clinical_GSE42127$vital_status<-clinical_df$`survival status:ch1`
# ajcc_pathologic_stage
clinical_GSE42127$ajcc_pathologic_stage<-clinical_df$`final.pat.stage:ch1`
# Gender: gender
clinical_GSE42127$Gender<-clinical_df$`gender:ch1`
# Adjuvant_chemo
clinical_GSE42127$Adjuvant_chemo<-clinical_df$`had_adjuvant_chemo:ch1`

# Disease-free survival: DFS
clinical_GSE42127$DFS<-
  # NTE-status
  clinical_GSE42127$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE42127$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE42127$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE42127$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE42127$Race<-
  # Adjuvant_radiation
  clinical_GSE42127$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE42127$Smoking_status<-
  # Predominant_subtype
  clinical_GSE42127$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE42127$EGFR_mut<-
  # KRAS mutation
  clinical_GSE42127$KRAS_mut<-
  # TP53 mutation
  clinical_GSE42127$TP53_mut<-
  # STK11 mutation
  clinical_GSE42127$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE42127<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE42127)<-gene_info$gene_symbol
rownames(expression_GSE42127)<-colnames(data.MA)

# 25. GSE43458
data <- getGEO(filename="GEO/GSE43458_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE43458<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE43458$Histology<-clinical_df$`histology:ch1`
#  Smoking_status
clinical_GSE43458$Smoking_status<-clinical_df$`smoking status:ch1`
# Age: age
clinical_GSE43458$Age<-
  # survival_time: survival_time
  clinical_GSE43458$survival_time<-
  # vital_status: vital_status
  clinical_GSE43458$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE43458$DFS<-
  # NTE-status
  clinical_GSE43458$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE43458$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE43458$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE43458$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE43458$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE43458$Race<-
  # Gender: gender
  clinical_GSE43458$Gender<-
  # Adjuvant_chemo
  clinical_GSE43458$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE43458$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE43458$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE43458$EGFR_mut<-
  # KRAS mutation
  clinical_GSE43458$KRAS_mut<-
  # TP53 mutation
  clinical_GSE43458$TP53_mut<-
  # STK11 mutation
  clinical_GSE43458$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

genes<-pbsapply(data.Platform$gene_assignment, FUN=function(x){
  unlist(strsplit(x, split=" // "))[2]
}, cl=clu)

clusterExport(clu, "genes")
entrez.to.id<-pblapply(unique(genes), FUN=function(geneid){
  indices=which(genes %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(genes)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE43458<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$gene_symbol[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE43458)<-gene_info$gene_symbol
rownames(expression_GSE43458)<-colnames(data.MA)

# 26. GSE43580
data <- getGEO(filename="GEO/GSE43580_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE43580<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE43580$Histology<-sapply(clinical_df$source_name_ch1, FUN=function(x){
  unlist(strsplit(x, ","))[1]
})
# ajcc_pathologic_t
clinical_GSE43580$ajcc_pathologic_t<-sapply(clinical_df$`tnm stage:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 1, 3), substr(tnm, 1, 2)))
})
# ajcc_pathologic_n
clinical_GSE43580$ajcc_pathologic_n<-sapply(clinical_df$`tnm stage:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 4, 5), substr(tnm, 3, 4)))
})
# ajcc_pathologic_m
clinical_GSE43580$ajcc_pathologic_m<-sapply(clinical_df$`tnm stage:ch1`, FUN=function(tnm){
  ifelse(tnm == "neo", "[Not Available]", 
         ifelse(nchar(tnm) == 7, substr(tnm, 6, 7), substr(tnm, 5, 6)))
})  
# ajcc_pathologic_stage
clinical_GSE43580$ajcc_pathologic_stage<-clinical_df$`ajcc uicc stage:ch1`
# Race: race
clinical_GSE43580$Race<-clinical_df$`ethnicity:ch1`
# Gender: gender
clinical_GSE43580$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE43580$Smoking_status<-clinical_df$`smoking status:ch1` 
# Age: age
clinical_GSE43580$Age<-
  # survival_time: survival_time
  clinical_GSE43580$survival_time<-
  # vital_status: vital_status
  clinical_GSE43580$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE43580$DFS<-
  # NTE-status
  clinical_GSE43580$NTE_status<-
  # Adjuvant_chemo
  clinical_GSE43580$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE43580$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE43580$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE43580$EGFR_mut<-
  # KRAS mutation
  clinical_GSE43580$KRAS_mut<-
  # TP53 mutation
  clinical_GSE43580$TP53_mut<-
  # STK11 mutation
  clinical_GSE43580$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE43580<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE43580)<-gene_info$gene_symbol
rownames(expression_GSE43580)<-colnames(data.MA)

# GSE46539
data <- getGEO(filename="GEO/GSE46539-GPL6883_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE46539<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE46539$Histology<-sapply(clinical_df$source_name_ch1, FUN=function(x){
  if(x=="Normal lung tissue adjacent to lung adenocarcinoma"){
    return("Control")
  }else{
    "ADC"
  }
})
# Age: age
clinical_GSE46539$Age<-clinical_df$`age:ch1`
#  Smoking_status
clinical_GSE46539$Smoking_status<-clinical_df$`smoking status:ch1`
# Gender: gender
clinical_GSE46539$Gender<-clinical_df$`gender:ch1`

# survival_time: survival_time
clinical_GSE46539$survival_time<-
  # vital_status: vital_status
  clinical_GSE46539$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE46539$DFS<-
  # NTE-status
  clinical_GSE46539$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE46539$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE46539$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE46539$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE46539$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE46539$Race<-
  # Adjuvant_chemo
  clinical_GSE46539$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE46539$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE46539$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE46539$EGFR_mut<-
  # KRAS mutation
  clinical_GSE46539$KRAS_mut<-
  # TP53 mutation
  clinical_GSE46539$TP53_mut<-
  # STK11 mutation
  clinical_GSE46539$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE46539<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
expression_GSE46539<-t(apply(expression_GSE46539, MARGIN=1, FUN=function(x){
  x_<-na.omit(x)
  return((x- min(x_)) /(max(x_)-min(x_)))
})*15)
colnames(expression_GSE46539)<-gene_info$gene_symbol
rownames(expression_GSE46539)<-colnames(data.MA)

# 28. GSE50081
data <- getGEO(filename="GEO/GSE50081_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE50081<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE50081$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE50081$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE50081$survival_time<-as.numeric(clinical_df$`survival time:ch1`)*365.25
# vital_status: vital_status
clinical_GSE50081$vital_status<-clinical_df$`status:ch1`
# Disease-free survival: DFS
clinical_GSE50081$DFS<-as.numeric(clinical_df$`disease-free survival time:ch1`)*365.25
# NTE-status
clinical_GSE50081$NTE_status<-clinical_df$`recurrence:ch1`
# ajcc_pathologic_t
clinical_GSE50081$ajcc_pathologic_t<-clinical_df$`t-stage:ch1`
# ajcc_pathologic_n
clinical_GSE50081$ajcc_pathologic_n<-clinical_df$`n-stage:ch1`
# ajcc_pathologic_m
clinical_GSE50081$ajcc_pathologic_m<-clinical_df$`m-stage:ch1`
# ajcc_pathologic_stage
clinical_GSE50081$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Gender: gender
clinical_GSE50081$Gender<-clinical_df$`Sex:ch1`
#  Smoking_status
clinical_GSE50081$Smoking_status<-clinical_df$`smoking:ch1`
# Race: race
clinical_GSE50081$Race<-
  # Adjuvant_chemo
  clinical_GSE50081$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE50081$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE50081$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE50081$EGFR_mut<-
  # KRAS mutation
  clinical_GSE50081$KRAS_mut<-
  # TP53 mutation
  clinical_GSE50081$TP53_mut<-
  # STK11 mutation
  clinical_GSE50081$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE50081<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE50081)<-gene_info$gene_symbol
rownames(expression_GSE50081)<-colnames(data.MA)

# 29. GSE60644
data <- getGEO(filename="GEO/GSE60644_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE60644<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE60644$Histology<-clinical_df$`histological subtype:ch1`
  # Age: age
  clinical_GSE60644$Age<-
  # survival_time: survival_time
  clinical_GSE60644$survival_time<-
  # vital_status: vital_status
  clinical_GSE60644$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE60644$DFS<-
  # NTE-status
  clinical_GSE60644$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE60644$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE60644$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE60644$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE60644$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE60644$Race<-
  # Gender: gender
  clinical_GSE60644$Gender<-
  # Adjuvant_chemo
  clinical_GSE60644$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE60644$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE60644$Smoking_status<-
  # Predominant_subtype
  clinical_GSE60644$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE60644$EGFR_mut<-
  # KRAS mutation
  clinical_GSE60644$KRAS_mut<-
  # TP53 mutation
  clinical_GSE60644$TP53_mut<-
  # STK11 mutation
  clinical_GSE60644$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE60644<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE60644)<-gene_info$gene_symbol
rownames(expression_GSE60644)<-colnames(data.MA)

# 30. GSE63459
data <- getGEO(filename="GEO/GSE63459_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE63459<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE63459$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE63459$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE63459$survival_time<-clinical_df$`survival after surgery (days):ch1`
# vital_status: vital_status
clinical_GSE63459$vital_status<-clinical_df$`death due to cancer:ch1`
# ajcc_pathologic_stage
clinical_GSE63459$ajcc_pathologic_stage<-clinical_df$`7th tnm stage:ch1`
# Race: race
clinical_GSE63459$Race<-clinical_df$`race:ch1`
# Gender: gender
clinical_GSE63459$Gender<-clinical_df$`Sex:ch1`
#  Smoking_status
clinical_GSE63459$Smoking_status<-clinical_df$`smoking status:ch1`
# KRAS mutation
clinical_GSE63459$KRAS_mut<-clinical_df$`kras status:ch1` 
# TP53 mutation
clinical_GSE63459$TP53_mut<-clinical_df$`tp53 mutation:ch1`

# Disease-free survival: DFS
clinical_GSE63459$DFS<-
  # NTE-status
  clinical_GSE63459$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE63459$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE63459$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE63459$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE63459$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE63459$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE63459$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE63459$EGFR_mut<-
  # STK11 mutation
  clinical_GSE63459$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE63459<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE63459)<-gene_info$gene_symbol
rownames(expression_GSE63459)<-colnames(data.MA)

# 31. GSE68465
data <- getGEO(filename="GEO/GSE68465_series_matrix.txt.gz")
clinical_df <- read.csv("./LungCancer_validation_data/GSE68465_clinical.CSV")

# Sample IDs
clinical_GSE68465<-data.frame(ID=clinical_df$accession)
# Histology
clinical_GSE68465$Histology<-clinical_df$histology
# Age: age
clinical_GSE68465$Age<-clinical_df$age
# survival_time: survival_time
clinical_GSE68465$survival_time<-clinical_df$last_contact_days_to
# vital_status: vital_status
clinical_GSE68465$vital_status<-clinical_df$vital_status
# NTE-status
clinical_GSE68465$NTE_status<-clinical_df$new_tumor_event_dx_indicator
# Disease-free survival: DFS
clinical_GSE68465$DFS<-clinical_df$new_tumor_event_dx_days_to
# ajcc_pathologic_t
clinical_GSE68465$ajcc_pathologic_t<-clinical_df$ajcc_pathologic_t
# ajcc_pathologic_n
clinical_GSE68465$ajcc_pathologic_n<-clinical_df$ajcc_pathologic_n
# ajcc_pathologic_stage
clinical_GSE68465$ajcc_pathologic_stage<-clinical_df$ajcc_pathologic_stage
# Race: race
clinical_GSE68465$Race<-clinical_df$race
# Gender: gender
clinical_GSE68465$Gender<-clinical_df$gender
# Adjuvant_chemo
clinical_GSE68465$Adjuvant_chemo<-clinical_df$pharmaceutical_tx_adjuvant
# Adjuvant_radiation
clinical_GSE68465$Adjuvant_radiation<-clinical_df$radiation_treatment_adjuvant
#  Smoking_status
clinical_GSE68465$Smoking_status<-clinical_df$tobacco_smoking_history_indicator
# ajcc_pathologic_m
clinical_GSE68465$ajcc_pathologic_m<-
  # Predominant_subtype
  clinical_GSE68465$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE68465$EGFR_mut<-
  # KRAS mutation
  clinical_GSE68465$KRAS_mut<-
  # TP53 mutation
  clinical_GSE68465$TP53_mut<-
  # STK11 mutation
  clinical_GSE68465$STK11_mut<-rep("[Not Available]", length(clinical_df$accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data
# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
data.MA<-log2(data.MA)

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE68465<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE68465)<-gene_info$gene_symbol
rownames(expression_GSE68465)<-colnames(data.MA)

# 32. GSE72094
data <- getGEO(filename="GEO/GSE72094_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE72094<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE72094$Histology<-clinical_df$source_name_ch1
# Age: age
clinical_GSE72094$Age<-sapply(clinical_df$`age_at_diagnosis:ch1`, FUN=function(x){
  if(x=="NA"){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# survival_time: survival_time
clinical_GSE72094$survival_time<-sapply(clinical_df$`survival_time_in_days:ch1`, FUN=function(x){
  if(x=="NA"){
    "[Not Available]"
  }else{
    as.numeric(x)
  }
})
# vital_status: vital_status
clinical_GSE72094$vital_status<-clinical_df$`vital_status:ch1`
# ajcc_pathologic_stage
clinical_GSE72094$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Race: race
clinical_GSE72094$Race<-clinical_df$`race:ch1`
# Gender: gender
clinical_GSE72094$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE72094$Smoking_status<-clinical_df$`smoking_status:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE72094$EGFR_mut<-clinical_df$`egfr_status:ch1`
# KRAS mutation
clinical_GSE72094$KRAS_mut<-clinical_df$`kras_status:ch1`
# TP53 mutation
clinical_GSE72094$TP53_mut<-clinical_df$`tp53_status:ch1`
# STK11 mutation
clinical_GSE72094$STK11_mut<-clinical_df$`stk11_status:ch1`

# Disease-free survival: DFS
clinical_GSE72094$DFS<-
  # NTE-status
  clinical_GSE72094$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE72094$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE72094$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE72094$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE72094$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE72094$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE72094$Predominant_subtype<-
  rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$EntrezGeneID[k], " ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE72094<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE72094)<-gene_info$gene_symbol
rownames(expression_GSE72094)<-colnames(data.MA)

# 33. GSE75037
data <- getGEO(filename="GEO/GSE75037_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE75037<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE75037$Histology<-clinical_df$`histology:ch1`
# Age: age
clinical_GSE75037$Age<-clinical_df$`age (yrs):ch1`
# ajcc_pathologic_stage
clinical_GSE75037$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Race: race
clinical_GSE75037$Race<-clinical_df$`race:ch1`
# Gender: gender
clinical_GSE75037$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE75037$Smoking_status<-clinical_df$`smoker:ch1`
## Mutation statuses
# EGFR mutation
clinical_GSE75037$EGFR_mut<-sapply(clinical_df$`egfr:ch1`, FUN=function(x){
  if(is.na(x)){
    "[Not Available]"
  }else{
    x
  }
})
# KRAS mutation
clinical_GSE75037$KRAS_mut<-sapply(clinical_df$`kras:ch1`, FUN=function(x){
  if(is.na(x)){
    "[Not Available]"
  }else{
    x
  }
})

# survival_time: survival_time
clinical_GSE75037$survival_time<-
  # vital_status: vital_status
  clinical_GSE75037$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE75037$DFS<-
  # NTE-status
  clinical_GSE75037$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE75037$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE75037$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE75037$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE75037$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE75037$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE75037$Predominant_subtype<-
  # TP53 mutation
  clinical_GSE75037$TP53_mut<-
  # STK11 mutation
  clinical_GSE75037$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=data.Platform$Entrez_Gene_ID[k]
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]

entrez.to.id<-lapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
})
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE75037<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE75037)<-gene_info$gene_symbol
rownames(expression_GSE75037)<-colnames(data.MA)

# 34. GSE81089
data <- getGEO(filename="GEO/GSE81089_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE81089<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE81089$Histology<-sapply(clinical_df$`histology:ch1`, FUN=function(x){
  if(is.na(x)){
    return("Control") 
  }
  if(x=='1'){
    return("SCC")
  }else if(x=='2'){
    return("ADC")
  }else if(x=='3'){
    return("LCC")
  }
})
# Age: age
clinical_GSE81089$Age<-sapply(clinical_df$`age:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }else{
    return(as.numeric(x))
  }
})
# survival_time: survival_time
clinical_GSE81089$survival_time<-sapply(1:nrow(clinical_df), FUN=function(k){
  if(is.na(clinical_df$`surgery date:ch1`[k])){
    return("[Not Available]")
  }
  
  surgery.date<-as_date(clinical_df$`surgery date:ch1`[k])
  
  if(clinical_df$`vital date:ch1`[k]=="n/a"){
    return("[Not Available]")
  }
  vital_date<-as_date(clinical_df$`vital date:ch1`[k])
  
  survival.days<-as.numeric(vital_date-surgery.date)
  if(is.na(survival.days)){
    return("[Not Available]")
  }else{
    return(survival.days)  
  }
})
# vital_status: vital_status
clinical_GSE81089$vital_status<-clinical_df$`dead:ch1`
# ajcc_pathologic_stage
clinical_GSE81089$ajcc_pathologic_stage<-sapply(clinical_df$`stage tnm:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }
  if(x=='1'){
    return("1a")
  }else if(x=="2"){
    return("1b")
  }else if(x=="3"){
    return("2a")
  }else if(x=="4"){
    return("2b")
  }else if(x=="5"){
    return("3a")
  }else if(x=="6"){
    return("3b")
  }else if(x=="7"){
    return("IV")
  }
})
# Gender: gender
clinical_GSE81089$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE81089$Smoking_status<-sapply(clinical_df$`smoking:ch1`, FUN=function(x){
  if(is.na(x)){
    return("[Not Available]")
  }
  if(x %in% c("1", "2")){
    "Ever smoker"
  }else{
    "Never smoker"
  }
})
# Disease-free survival: DFS
clinical_GSE81089$DFS<-
  # NTE-status
  clinical_GSE81089$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE81089$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE81089$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE81089$ajcc_pathologic_m<-
  # Race: race
  clinical_GSE81089$Race<-
  # Adjuvant_chemo
  clinical_GSE81089$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE81089$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE81089$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE81089$EGFR_mut<-
  # KRAS mutation
  clinical_GSE81089$KRAS_mut<-
  # TP53 mutation
  clinical_GSE81089$TP53_mut<-
  # STK11 mutation
  clinical_GSE81089$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA<-read.table("GEO/GSE81089_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header = T)
genes<-data.MA$GeneID
samples<-colnames(data.MA)[-1]
data.MA<-data.MA[,-1]
data.MA<-do.call(cbind, lapply(1:dim(data.MA)[2], FUN=function(k){
  log2(as.numeric(data.MA[,k])+1)
}))
colnames(data.MA)<-samples
rownames(data.MA)<-genes

expression_GSE81089<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE81089)<-gene_info$gene_symbol
rownames(expression_GSE81089)<-colnames(data.MA)

# 35. GSE101929
data <- getGEO(filename="GEO/GSE101929_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE101929<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE101929$Histology<-clinical_df$`tumor_normal status:ch1`
# Age: age
clinical_GSE101929$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE101929$survival_time<-clinical_df$`survival (days):ch1`
# vital_status: vital_status
clinical_GSE101929$vital_status<-clinical_df$`death due to lung cancer (all years):ch1`
# ajcc_pathologic_stage
clinical_GSE101929$ajcc_pathologic_stage<-clinical_df$`Stage:ch1`
# Race: race
clinical_GSE101929$Race<-clinical_df$`race:ch1`
# Gender: gender
clinical_GSE101929$Gender<-clinical_df$`gender:ch1`
#  Smoking_status
clinical_GSE101929$Smoking_status<-sapply(clinical_df$`smoking pack years:ch1`, FUN=function(x){
  if(x=="0"){
    "Never smoker"
  }else{
    "Ever smoker"
  }
})

# Disease-free survival: DFS
clinical_GSE101929$DFS<-
  # NTE-status
  clinical_GSE101929$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE101929$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE101929$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE101929$ajcc_pathologic_m<-
  # Adjuvant_chemo
  clinical_GSE101929$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE101929$Adjuvant_radiation<-
  # Predominant_subtype
  clinical_GSE101929$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE101929$EGFR_mut<-
  # KRAS mutation
  clinical_GSE101929$KRAS_mut<-
  # TP53 mutation
  clinical_GSE101929$TP53_mut<-
  # STK11 mutation
  clinical_GSE101929$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA=data@assayData$exprs
data.Platform<-data@featureData@data

# Mapping Probe ID to Entrez Gene ID
mapping.table<-do.call(rbind, pblapply(c(1:dim(data.Platform)[1]), FUN=function(k){
  ID=data.Platform$ID[k]
  Entrez=unlist(strsplit(data.Platform$ENTREZ_GENE_ID[k], " /// ",fixed = T))
  if(length(Entrez)!=0){
    df<-data.frame(ID, Entrez)  
  }else{
    df<-data.frame(ID=ID, Entrez=NA)
  }
  return(df)
})) %>% filter(!is.na(Entrez))

data.MA<-data.MA[mapping.table$ID,]
if(max(data.MA)>1000){
  data.MA<-log2(data.MA)
}


clusterExport(clu, "mapping.table")
entrez.to.id<-pblapply(unique(mapping.table$Entrez), FUN=function(geneid){
  indices=which(mapping.table$Entrez %in% geneid)
}, cl=clu)
names(entrez.to.id)<-unique(mapping.table$Entrez)

data.MA<-do.call(rbind, pblapply(c(1:length(entrez.to.id)), FUN=function(k){
  if(length(entrez.to.id[[k]])>1){
    return(colMeans(data.MA[entrez.to.id[[k]],]))  
  }else{
    return(data.MA[entrez.to.id[[k]],])
  }
}))
rownames(data.MA)<-names(entrez.to.id)

expression_GSE101929<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE101929)<-gene_info$gene_symbol
rownames(expression_GSE101929)<-colnames(data.MA)

# 36. GSE101929
data <- getGEO(filename="GEO/GSE140343_series_matrix.txt.gz")
clinical_df <- data@phenoData@data

# Sample IDs
clinical_GSE140343<-data.frame(ID=clinical_df$geo_accession)
# Histology
clinical_GSE140343$Histology<-sapply(clinical_df$`disease state:ch1`, FUN=function(x){
  if(is.na(x)){
    "ADC"
  }else{
    "Control"
  }
})
# Age: age
clinical_GSE140343$Age<-clinical_df$`age:ch1`
# survival_time: survival_time
clinical_GSE140343$survival_time<-
  # vital_status: vital_status
  clinical_GSE140343$vital_status<-
  # Disease-free survival: DFS
  clinical_GSE140343$DFS<-
  # NTE-status
  clinical_GSE140343$NTE_status<-
  # ajcc_pathologic_t
  clinical_GSE140343$ajcc_pathologic_t<-
  # ajcc_pathologic_n
  clinical_GSE140343$ajcc_pathologic_n<-
  # ajcc_pathologic_m
  clinical_GSE140343$ajcc_pathologic_m<-
  # ajcc_pathologic_stage
  clinical_GSE140343$ajcc_pathologic_stage<-
  # Race: race
  clinical_GSE140343$Race<-
  # Gender: gender
  clinical_GSE140343$Gender<-
  # Adjuvant_chemo
  clinical_GSE140343$Adjuvant_chemo<-
  # Adjuvant_radiation
  clinical_GSE140343$Adjuvant_radiation<-
  #  Smoking_status
  clinical_GSE140343$Smoking_status<-
  # Predominant_subtype
  clinical_GSE140343$Predominant_subtype<-
  ## Mutation statuses
  # EGFR mutation
  clinical_GSE140343$EGFR_mut<-
  # KRAS mutation
  clinical_GSE140343$KRAS_mut<-
  # TP53 mutation
  clinical_GSE140343$TP53_mut<-
  # STK11 mutation
  clinical_GSE140343$STK11_mut<-rep("[Not Available]", length(clinical_df$geo_accession))

data.MA<-read.table("GEO/GSE140343_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header = T)
genes<-data.MA$GeneID
samples<-colnames(data.MA)[-1]
data.MA<-data.MA[,-1]
data.MA<-do.call(cbind, lapply(1:dim(data.MA)[2], FUN=function(k){
  log2(as.numeric(data.MA[,k])+1)
}))
colnames(data.MA)<-samples
rownames(data.MA)<-genes

expression_GSE140343<-do.call(cbind, lapply(1:dim(gene_info)[1], FUN=function(i){
  ind<-match(gene_info$ENTREZID[i], rownames(data.MA))
  if(is.na(ind)){
    return(rep(NA, dim(data.MA)[2]))
  }else{
    return(data.MA[ind,])
  }
}))
colnames(expression_GSE140343)<-gene_info$gene_symbol
rownames(expression_GSE140343)<-colnames(data.MA)

# Unifying clinical labels
clinical.data.list.total<-list(clinical_TCGA, clinical_GSE3141, clinical_GSE8894, clinical_GSE10072,
                         clinical_GSE10245, clinical_GSE10445, clinical_GSE10799, clinical_GSE11969,
                         clinical_GSE12667, clinical_GSE13213, clinical_GSE14814, clinical_GSE18842,
                         clinical_GSE19188, clinical_GSE19804, clinical_GSE26939, clinical_GSE28571,
                         clinical_GSE29016, clinical_GSE30219, clinical_GSE31210, clinical_GSE32863,
                         clinical_GSE37745, clinical_GSE40419, clinical_GSE41271, clinical_GSE42127,
                         clinical_GSE43458, clinical_GSE43580, clinical_GSE46539, clinical_GSE50081,
                         clinical_GSE60644, clinical_GSE63459, clinical_GSE68465, clinical_GSE72094,
                         clinical_GSE75037, clinical_GSE81089, clinical_GSE101929, clinical_GSE140343)
clinical.data.list.total[[1]]<-clinical.data.list.total[[1]][,-2]
colnames(clinical.data.list.total[[1]])[1]<-"ID"

for(k in 1:length(clinical.data.list.total)){
  time<-clinical.data.list.total[[k]]$survival_time
  if(length(setdiff(time, na.terms))<=5){
    next()
  }
  time.max<-max(as.numeric(setdiff(time, na.terms)))
  if(time.max<10){
    time<-sapply(time, FUN=function(k){
      if(k %in% na.terms){
        return(k)
      }else{
        return(as.numeric(k)*365.25)
      }
    })
    clinical.data.list.total[[k]]$survival_time<-time
    next()
  }
  
  if(time.max<300){
    time<-sapply(time, FUN=function(k){
      if(k %in% na.terms){
        return(k)
      }else{
        return(as.numeric(k)*(365.25/12))
      }
    })
  }
  
  clinical.data.list.total[[k]]$survival_time<-time
  
  time<-clinical.data.list.total[[k]]$DFS
  if(length(setdiff(time, na.terms))<=5){
    next()
  }
  time.max<-max(as.numeric(setdiff(time, na.terms)))
  if(time.max<10){
    time<-sapply(time, FUN=function(k){
      if(k %in% na.terms){
        return(k)
      }else{
        return(as.numeric(k)*365.25)
      }
    })
    clinical.data.list.total[[k]]$DFS<-time
    next()
  }  
  if(time.max<300){
    time<-sapply(time, FUN=function(k){
      if(k %in% na.terms){
        return(k)
      }else{
        return(as.numeric(k)*(365.25/12))
      }
    })
  }  
  clinical.data.list.total[[k]]$DFS<-time
}

category<-c("Histology", "Gender", "Race", "Smoking_status", "vital_status", "NTE_status",
            "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
            "Adjuvant_chemo", "Adjuvant_radiation","Predominant_subtype",
            "EGFR_mut", "KRAS_mut", "TP53_mut", "STK11_mut")

label.df<-do.call(rbind, lapply(category, FUN=function(cat){
  labels<-unlist(lapply(clinical.data.list.total, FUN=function(clinical_df){
    clinical_df[,match(cat, colnames(clinical_df))]
  })) %>% unique %>% sort
  data.frame(Category=cat, label=labels)
}))

write.csv(label.df, "clinical_labels_origin.csv", row.names = F)
label.df<-read.csv("clinical_labels_refined.csv")
label.df$label[which(is.na(label.df$label))]<-"NA"

# NA to [Not Available]
for(i in 1:length(clinical.data.list.total)){
  for(j in 1:ncol(clinical.data.list.total[[i]])){
    clinical.data.list.total[[i]][,j][which(is.na(clinical.data.list.total[[i]][,j]))]<-"[Not Available]"  
  }
}

# label update
for(cat in unique(label.df$Category)){
  labs<-label.df %>% filter(Category==cat)
  for(i in 1:nrow(labs)){
    for(k in 1:length(clinical.data.list.total)){
      label.origin<-clinical.data.list.total[[k]][,cat]
      label.ind<-which(label.origin==labs$label[i])
      clinical.data.list.total[[k]][label.ind,cat]<-labs$Label_refined[i]
    }
  }
}

# indexing samples per dataset
for(i in 1:length(clinical.data.list.total)){
  clinical.data.list.total[[i]]$Dataset<-rep(i, nrow(clinical.data.list.total[[i]]))
  clinical.data.list.total[[i]]$sample_ind<-1:nrow(clinical.data.list.total[[i]])
}

expression.data.list.total<-list(expression_TCGA, expression_GSE3141, expression_GSE8894, expression_GSE10072,
                           expression_GSE10245, expression_GSE10445, expression_GSE10799, expression_GSE11969,
                           expression_GSE12667, expression_GSE13213, expression_GSE14814, expression_GSE18842,
                           expression_GSE19188, expression_GSE19804, expression_GSE26939, expression_GSE28571,
                           expression_GSE29016, expression_GSE30219, expression_GSE31210, expression_GSE32863,
                           expression_GSE37745, expression_GSE40419, expression_GSE41271, expression_GSE42127,
                           expression_GSE43458, expression_GSE43580, expression_GSE46539, expression_GSE50081,
                           expression_GSE60644, expression_GSE63459, expression_GSE68465, expression_GSE72094,
                           expression_GSE75037, expression_GSE81089, expression_GSE101929, expression_GSE140343)

dataset.name<-c("TCGA", "GSE3141", "GSE8894", "GSE10072", "GSE10245", 
                "GSE10445", "GSE10799", "GSE11969", "GSE12667", 
                "GSE13213", "GSE14814", "GSE18842", "GSE19188", 
                "GSE19804", "GSE26939", "GSE28571", "GSE29016", 
                "GSE30219", "GSE31210", "GSE32863", "GSE37745", 
                "GSE40419", "GSE41271", "GSE42127", "GSE43458", 
                "GSE43580", "GSE46539", "GSE50081", "GSE60644",
                "GSE63459", "GSE68465", "GSE72094", "GSE75037",
                "GSE81089", "GSE101929", "GSE140343")

clinical.annotation.total<-do.call(rbind, clinical.data.list.total)
rm(list=c('clinical_TCGA', 'clinical_GSE3141', 'clinical_GSE8894', 'clinical_GSE10072',
          'clinical_GSE10245', 'clinical_GSE10445', 'clinical_GSE10799', 'clinical_GSE11969',
          'clinical_GSE12667', 'clinical_GSE13213', 'clinical_GSE14814', 'clinical_GSE18842',
          'clinical_GSE19188', 'clinical_GSE19804', 'clinical_GSE26939', 'clinical_GSE28571',
          'clinical_GSE29016', 'clinical_GSE30219', 'clinical_GSE31210', 'clinical_GSE32863',
          'clinical_GSE37745', 'clinical_GSE40419', 'clinical_GSE41271', 'clinical_GSE42127',
          'clinical_GSE43458', 'clinical_GSE43580', 'clinical_GSE46539', 'clinical_GSE50081',
          'clinical_GSE60644', 'clinical_GSE63459', 'clinical_GSE68465', 'clinical_GSE72094',
          'clinical_GSE75037', 'clinical_GSE81089', 'clinical_GSE101929', 'clinical_GSE140343',
          'expression_TCGA', 'expression_GSE3141', 'expression_GSE8894', 'expression_GSE10072',
          'expression_GSE10245', 'expression_GSE10445', 'expression_GSE10799', 'expression_GSE11969',
          'expression_GSE12667', 'expression_GSE13213', 'expression_GSE14814', 'expression_GSE18842',
          'expression_GSE19188', 'expression_GSE19804', 'expression_GSE26939', 'expression_GSE28571',
          'expression_GSE29016', 'expression_GSE30219', 'expression_GSE31210', 'expression_GSE32863',
          'expression_GSE37745', 'expression_GSE40419', 'expression_GSE41271', 'expression_GSE42127',
          'expression_GSE43458', 'expression_GSE43580', 'expression_GSE46539', 'expression_GSE50081',
          'expression_GSE60644', 'expression_GSE63459', 'expression_GSE68465', 'expression_GSE72094',
          'expression_GSE75037', 'expression_GSE81089', 'expression_GSE101929', 'expression_GSE140343'))
rm(list=c('ptm', 'symbol.to.id', 'labs', 'labels', 'label.df', 'entrez.to.id', 'symbol.to.id', 'data', 'data_long', 
     'data.MA','data.Platform', 'na.ind', 'remove.ind', 'samples', 'tb', 'time', 'time.max',
     'label.ind', 'label.origin'))
