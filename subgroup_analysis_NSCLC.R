### Selection of samples to be analyzed
# Target sample: Non-small cell lung cancer (Adenocarcinoma + squamous cell carcinoma + large cell carcinoma)
target.samples<-clinical.annotation.total %>% filter(Histology %in% c("ADC", "LCC", "SCC")) %>% 
  filter(!vital_status %in% na.terms) %>% filter(!survival_time %in% na.terms) %>% filter(survival_time>0)
target.samples <- target.samples %>% filter(Dataset %in% names(which(target.samples$Dataset %>% table>10)))

expression.data.list<-lapply(target.samples$Dataset %>% unique, FUN=function(ind){
  sample.ind=(target.samples %>% filter(Dataset==ind))$sample_ind
  expression.data.list.total[[ind]][sample.ind,]
})

dat<-do.call(rbind, expression.data.list)
na.ind<-which(is.na(apply(dat, MARGIN=2, sum)))
dat<-dat[,-na.ind]

# Isolation of RNA-seq datasets
rnaseq.ind<-which(target.samples$Dataset %in% c(1, 22, 34, 36))
dat_rnaseq<-dat[rnaseq.ind,]
target.samples_rnaseq<-target.samples[rnaseq.ind,]

datalist.rnaseq<-lapply(unique(target.samples_rnaseq$Dataset), FUN=function(k){
  dat_rnaseq[which(target.samples_rnaseq$Dataset==k),]
})

dat<-dat[-rnaseq.ind,]
target.samples<-target.samples[-rnaseq.ind,]

merged.data<-ComBat(dat=t(dat),batch = target.samples$Dataset) %>% t

# FSQN to RNA-seq datasets
datalist.rnaseq<-lapply(datalist.rnaseq, FUN=function(x){
  quantileNormalizeByFeature(matrix_to_normalize=x, target_distribution_matrix=merged.data)
})

datalist.rnaseq<-do.call(rbind, datalist.rnaseq)

merged.data<-rbind(merged.data, datalist.rnaseq)
target.samples<-rbind(target.samples, target.samples_rnaseq)

rm(datalist.rnaseq);rm(target.samples_rnaseq)

expression.data.list<-lapply(1:length(unique(target.samples$Dataset)), FUN=function(i){
  dataset.ind<-sort(unique(target.samples$Dataset))[i]
  merged.data[which(target.samples$Dataset==dataset.ind),]
})

clinical.data.list<-lapply(1:length(unique(target.samples$Dataset)), FUN=function(i){
  dataset.ind<-sort(unique(target.samples$Dataset))[i]
  target.samples[which(target.samples$Dataset==dataset.ind),]
})
target.dataset.name<-dataset.name[sort(unique(target.samples$Dataset))]

# Final merged expression data and clinical annotations
merged.data<-do.call(rbind, expression.data.list)
target.samples<-do.call(rbind, clinical.data.list)

# Processing time-to-event data
surv_df<-target.samples %>% dplyr::select(survival_time, vital_status) 
y<-data.frame(time=as.numeric(surv_df$survival_time),
              event=as.numeric(surv_df$vital_status))

clusterExport(clu, "y")
clusterExport(clu, "merged.data")
clusterEvalQ(clu, {
  f<-readRDS("f_processed.RDS")
  f2g<-readRDS("f2g_processed.RDS")
  NULL
})

### Mechanism list selection
# Significant gene candidates from randomly sampled patient subgroup for n times
n = 1000
sig_genes_subsample<-pblapply(1:n, FUN=function(n){
  sample_ind<-unique(sort(sample(1:nrow(y), nrow(y), replace=T)))
  
  sig.gene.res<-do.call(rbind, pblapply(c(1:ncol(merged.data)), FUN=function(g.ind){
    cox_ph_model<-coxph(Surv(as.numeric(y$time)[sample_ind], as.numeric(y$event)[sample_ind])~
                          merged.data[sample_ind,g.ind])
    model.summary<-summary(cox_ph_model)
    HR<-model.summary$conf.int[1, 1]
    lower<-model.summary$conf.int[1, 3]
    upper<-model.summary$conf.int[1, 4]
    p.value<-model.summary$coefficients[1, 5]
    c.ind=model.summary$concordance[1]
    res<-data.frame(gene=colnames(merged.data)[g.ind], HR, lower, upper, p.value, c.ind)
    return(res)
  }))
  sig.gene.res$adjusted_p<-p.adjust(sig.gene.res$p.value)
  
  return((sig.gene.res %>% filter(adjusted_p<=0.05))$gene %>% sort)
  
}, cl=clu)


# Selecting total significant gene candidates which contained in significant gene candidates above proportion p1
p1 = 0.1
# Selecting significant mechanism candidates which contain significant gene candidates above proportion p2
p2 = 0.3
sig_gene_candidates<-names(which(table(unlist(sig_genes_subsample)) >= p*length(sig_genes_subsample)))
clusterExport(clu, "sig_gene_candidates")


mech_gene_length<-pbsapply(unique(f2g$MOD_ID), FUN=function(k){
  length((f2g %>% filter(MOD_ID==k))$Symbol)
}, cl=clu)

mech_sig_gene_list<-pblapply(unique(f2g$MOD_ID), FUN=function(k){
  sort((f2g %>% filter(MOD_ID==k) %>% filter(Symbol %in% sig_gene_candidates))$Symbol)
}, cl=clu)

mechanism_candidates<-tibble(Name=f$Name, Source=f$Source,
                             Gene_list=mech_sig_gene_list,
                             len_mech_gene=mech_gene_length,
                             len_sig_gene=sapply(mech_sig_gene_list, length)) %>% 
  filter(len_sig_gene>=5) %>% mutate(Percentage=len_sig_gene/len_mech_gene)

clusterExport(clu, 'mechanism_candidates')

mechanism_candidates<-do.call(rbind, pblapply(unique(mechanism_candidates$Gene_list), FUN=function(m){
  matching_indices<-which(sapply(1:nrow(mechanism_candidates), FUN=function(i){
    identical(m, mechanism_candidates$Gene_list[[i]])  
  }))
  
  Name<-list(mechanism_candidates$Name[matching_indices])
  Source<-unique(mechanism_candidates$Source[matching_indices])
  Gene_list=list(m)
  len_mech_gene<-unique(mechanism_candidates$len_mech_gene[matching_indices])
  len_sig_gene<-unique(mechanism_candidates$len_sig_gene[matching_indices])
  Percentage<-unique(mechanism_candidates$Percentage[matching_indices])
  
  return(tibble(Name, Gene_list, len_mech_gene, len_sig_gene, Percentage))
}, cl=clu)) %>% filter(Percentage >= p2)


##
set.seed(42)
K=5 # Number of the patient subgroup

# 1. Random sampling of latent patient subgroups
Z<-sample(1:K, size=dim(target.samples)[1], replace=TRUE)

# 2. Random sampling of gene sets
M_Z<-sample(dim(mechanism_candidates)[1], size=K, replace=FALSE)

# 3. Subgroup-specific model optimization
weibull_aft_models<-pblapply(1:K, FUN=function(k){
  # Optimizing Weibull AFT model
  sample_ind<-which(Z==k)
  genes=mechanism_candidates$Gene_list[[M_Z[k]]]
  
  g.ind<-match(genes, colnames(merged.data))
  
  weibull_aft_model <-survreg(
    Surv(as.numeric(y$time)[sample_ind], as.numeric(y$event)[sample_ind])~
      merged.data[sample_ind,g.ind],
    dist="weibull"
  )

  return(weibull_aft_model)  
})
