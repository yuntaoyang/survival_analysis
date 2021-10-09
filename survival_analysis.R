# survival analysis of gene signature low and gene signature high (IL6_JAK_STAT3 pathway)
# TCGA-BLCA data

#---- set up libraries----------------------------------------------------------
library(dplyr)
library(io)
library(SummarizedExperiment)
library(survival)
library(survminer)

#---- set up parameters---------------------------------------------------------
project <- 'tcga_blca_survival'

#---- set up path --------------------------------------------------------------
path_data <- '/Users/yyang18/Google Drive/Data/TCGA_GEO/TCGA_BLCA/'
path_gsea <- '/Users/yyang18/Google Drive/Data/GSEA/'
path_gene_id <- '/Users/yyang18/Google Drive/Data/Gene_ID_Convert/'
path_out <- './out/'
 
#---- input --------------------------------------------------------------------
clinical <- read.csv(paste0(path_data,"tcga-blca_clinical_information.tsv"),sep = '\t')
biospecimen <- read.csv(paste0(path_data,"tcga-blca_biospecimen.tsv"),sep = '\t')
expdat <- read.csv(paste0(path_data,"tcga-blca_expdat_fpkm.csv"),header = TRUE, row.names = 1,check.names = FALSE)
hallmark <- read.csv(paste0(path_gsea,"h.all.v7.2.symbols.gmt"),sep = '\t',header = FALSE, row.names = 1)
gene_id_convert <- read.csv(paste0(path_gene_id,"Ensembl_ID.csv"))[,c('gene_id','gene_name')]

#---- data proprecessing -------------------------------------------------------
survival_data <- clinical[,c('submitter_id','days_to_last_follow_up','vital_status')]
survival_data <- survival_data[!is.na(survival_data$days_to_last_follow_up), ] %>% 
                  filter(vital_status == 'Alive' | vital_status== 'Dead') # keep alive & dead in vital status
primary_cancer <- filter(biospecimen, sample_type == 'Primary Tumor') # get primary cancer samples from biospecimen
IL6_JAK_STAT3 <- as.character(hallmark['HALLMARK_IL6_JAK_STAT3_SIGNALING',])[-c(1)] # select IL6_JAK_STAT3 gene signatrues
IL6_JAK_STAT3 <- IL6_JAK_STAT3[IL6_JAK_STAT3 != ""] # remove NA 
gene_signature <- data.frame('gene_name' = IL6_JAK_STAT3) %>% # convert to gene_id
                  merge(gene_id_convert) %>%
                  mutate(gene_id = substr(gene_id,1,15)) %>%
                  unique()

#---- activity score -----------------------------------------------------------
FPKM <- expdat[gene_signature$gene_id,]
names(FPKM) <- substring(names(FPKM),1,16)
FPKM_primary_cancer <- FPKM[,which(names(FPKM) %in% primary_cancer$submitter_id)] # primary cancer sample
FPKM_primary_cancer_Normalized <- t(scale(t(FPKM_primary_cancer))) %>%
                                  colSums() %>%
                                  sort()
rank <- data.frame(submitter_id = names(FPKM_primary_cancer_Normalized), activity_score = as.numeric(FPKM_primary_cancer_Normalized)) %>%
        mutate(submitter_id = substr(submitter_id,1,12))

#---- input of survival model --------------------------------------------------
survival <- merge(survival_data,rank)[-c(54, 20, 21), ] %>% # remove unusual survival days and duplicates
            arrange(activity_score) %>% # sort based on expression value of gene signature
            mutate(group = c(rep("low",149),rep("high",149))) %>% # set group information
            mutate(vital_status = replace(vital_status, vital_status == 'Alive',1)) %>%
            mutate(vital_status = replace(vital_status, vital_status == 'Dead',2)) %>%
            mutate(group = replace(group, group == 'low',1)) %>%
            mutate(group = replace(group, group == 'high',2)) %>%
            mutate(vital_status = as.numeric(as.character(vital_status))) %>%
            mutate(group = as.numeric(as.character(group))) %>%
            mutate(years_to_last_follow_up = days_to_last_follow_up/365)

#---- survival analysis --------------------------------------------------------
fit <- survfit(Surv(years_to_last_follow_up, vital_status) ~ group, data = survival)
print(fit)
pdf(paste0(path_out,project,'.pdf'),width = 4, height = 4, onefile = FALSE)
ggsurvplot(fit, 
           pval = TRUE, 
           conf.int = FALSE, 
           conf.int.style = "step",
           xlab = "Time in years", 
           break.time.by = 2, 
           ggtheme = theme_classic(),
           risk.table = "absolute", 
           risk.table.y.text.col = TRUE, 
           risk.table.y.text = FALSE,
           ncensor.plot =FALSE, 
           surv.median.line = "hv", 
           fontsize = 2.5,
           xlim = c(0, 10), 
           legend.labs = c("Low", "High"),
           palette = c("#E7B800", "#2E9FDF"))
dev.off()




