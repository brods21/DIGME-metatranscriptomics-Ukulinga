 
# 1/27/2025
# eggnog transdecoder but with salmon 
  # Using ADJUSTED KO AND CAZY counts: divided count by # assignments the eggnog algo gave it
  # adjusting for "uncertainty according to Brandon, Corinn
# WHAT I AM DOING:
  # use relative abundance - multiple wayss
  # use deseq normalization f
  # compare results

rm(list=ls())

library(tidyverse) # not working right now. 
library(ggplot2)
library(dplyr)
library(readr)
library(DESeq2) 
library(broom)
library(lme4)
library(vsn)
library(hexbin)
library(apeglm)


source("/Users/caitlinbroderick/Documents/DIGME/metatranscriptomics/DIGME_RNA_functionsModified.R")


#############################################################################
# BRING in and process FILES !!!
#############################################################################

setwd("/Users/caitlinbroderick/Documents/DIGME/metatranscriptomics")


# FIRST OFF : bring in total number of functionally annotated genes by eggNOG

# bring in counts files
  # these are KO and CAZy adjusted_inclTPM/
ko_counts_init= read.csv("L0/ko_counts_scaled_011725_RNA.csv")
#cazy_counts_init= read.csv("L0/cazy_counts_scaled_011725_RNA.csv")

# total number of annotated reads
total_eggnog_annotations  = read.csv("L0/annotated_eggnog_reads_persample_RNA.csv")

#add coreuniqueID col.
ko_counts_init$CoreUniqueID <- gsub(pattern = "_S[0-9].*", replacement = "", ko_counts_init$SampleID)
#cazy_counts_init$CoreUniqueID <- gsub(pattern = "_S[0-9].*", replacement = "", cazy_counts_init$SampleID)
total_eggnog_annotations$CoreUniqueID <- gsub(pattern = "_S[0-9].*", replacement = "", total_eggnog_annotations$SampleID)

# READ IN METADATA
meta  = read.csv ("L0/metadata.csv" )

# merge with meta
ko_counts_init2 =merge(ko_counts_init,meta_new,by= "CoreUniqueID")
#cazy_counts_init2 =merge(cazy_counts_init,meta_new,by= "CoreUniqueID")
total_eggnog_annotations_wmeta =merge( total_eggnog_annotations, meta,by= "CoreUniqueID")

# remove ribosomal rrna!!
pathway_to_ko_init <- read.csv("L1/Extra files for analysis/KEGG KO to pathway.csv")
names(pathway_to_ko_init)
ribosomal = pathway_to_ko_init %>% filter(Pathway =="map03010")

ko_counts_init3  = ko_counts_init2 %>% filter(!KO %in% ribosomal$KO)

# calculate totalKO and cazy counts per sample
ko_counts_all = ko_counts_init3 %>% group_by(SampleID) %>% summarize(totalKOsample = sum(count))
#cazy_counts_all = cazy_counts_init2 %>% group_by(SampleID) %>% summarize(totalCAZYreadssample = sum(count))

# merge count data with both "sample counts" - database specific. and total annotated. 
ko_counts_init4 = merge (ko_counts_init3 , total_eggnog_annotations, by = c("CoreUniqueID", "SampleID") )
#cazy_counts_init3 = merge (cazy_counts_init2 , total_eggnog_annotations, by = c("CoreUniqueID", "SampleID") )



ko_counts = merge(ko_counts_init4, ko_counts_all, by= "SampleID")
#cazy_counts = merge(cazy_counts_init3, cazy_counts_all, by= "SampleID")

# Count up BOTH relative abundance measurements
  # relabundance out of all annotated eggnog
ko_counts$relabun_eggnog = ko_counts$count / ko_counts$totalsamplereads * 100
#cazy_counts$relabun_eggnog = cazy_counts$count / cazy_counts$totalsamplereads * 100
  
# relative abundance out of all annotated database-specific
ko_counts$relabun_KO = ko_counts$count / ko_counts$totalKOsample * 100
#cazy_counts$relabun_cazy = cazy_counts$count / cazy_counts$totalCAZYreadssample* 100

# a couple plots to look aat how closely correlated these are
plot(ko_counts$relabun_KO,ko_counts$relabun_eggnog) # very CLOse!
plot(ko_counts$relabun_KO,ko_counts$TPM) # less close,,,,,, 
plot(ko_counts$relabun_eggnog,ko_counts$TPM)
#plot(cazy_counts$relabun_cazy, cazy_counts$relabun_eggnog) # interesting, far less close. 
#plot(cazy_counts$relabun_cazy, cazy_counts$TPM) # MUCH less close!!!!! how come. 

#cazy_counts$logWP <- -log(cazy_counts$ActualWP)
############################################################
# Site x species to do deseq normalization. 
############################################################

#### KO
sitexsp_raw_wmeta <-ko_counts %>% dplyr::select(c(SampleID,RainTrt,logWP,CO2_eq,KO,count)) %>% 
  pivot_wider(names_from=KO, values_from = count, values_fill = 0)
names(sitexsp_raw_wmeta)

sitexsp_raw <- as.data.frame (sitexsp_raw_wmeta[,-(1:4)] )
rownames(sitexsp_raw) = sitexsp_raw_wmeta$SampleID

sitexsp_raw_metadata <- sitexsp_raw_wmeta[,(1:4)]
sitexsp_raw_metadata$RainTrt <- as.factor(sitexsp_raw_metadata$RainTrt)
table(sitexsp_raw_metadata$RainTrt )
levels(sitexsp_raw_metadata$RainTrt )

sitexsp_raw_t =t(sitexsp_raw )
sitexsp_raw_t <- round(sitexsp_raw_t)


#write.csv(sitexsp_wmeta, "RNA_SitexSp_KOs.csv",row.names=F)
############################################################
# KO : run deseq2  
  # in order to get filtered, normalized (NOT transformed!) counts for downstresm
# To assess which KOs have a relationship with GWC in A and D treatments.
############################################################

######ALL !
# CHANGING this to pply, what happens????
dds <- DESeqDataSetFromMatrix(countData = sitexsp_raw_t,
                                  colData = sitexsp_raw_metadata,
                                  design = ~RainTrt + RainTrt:poly(logWP,2)) # TESTED THIS !!!! Design does not matter for normalization output.

dds <- estimateSizeFactors(dds) # estimate size factors
length(dds ) # y

plot( sizeFactors(dds), colSums(counts(dds)), # assess them
      ylab = "library sizes", xlab = "size factors", cex = .6 ) #read counts normalized for sequencing depth 

deseq_counts_init = counts(dds,normalized =T) # 

# GOING FORWARD: using filtered-by-abundance, NOT !!! rlog- transformed data.
deseq_counts_t = as.data.frame (t(deseq_counts_init))
deseq_counts= rownames_to_column(deseq_counts_t, var="SampleID")
deseq_counts_long = deseq_counts %>% pivot_longer(!SampleID, names_to = "KO", values_to ="deseq_counts")

write.csv(deseq_counts_long, "L1/deseq_normalized_counts_072125.csv", row.names=F)


# QUICK compaarison of normalizaation methods
allcountmethods  = deseq_counts_long  %>% 
  merge (ko_counts %>% select (KO, relabun_KO, relabun_eggnog,SampleID), by = c ("SampleID", "KO")) # %>% # they are just zeroes
  #pivot_longer(cols = c (deseq_counts, relabun_KO,relabun_eggnog ), names_to = "method", values_to = "count")


relabun_KO = allcountmethods %>% 
  ggplot (aes ( x  =deseq_counts , y =  relabun_KO )) +
  geom_point(alpha = 0.4) +
  theme_bw()  +
    stat_cor(aes(label = after_stat(rr.label)), color = "black", r.digits = 4)


relabun_egg = allcountmethods %>% 
  ggplot (aes ( x  =deseq_counts , y =  relabun_eggnog )) +
  geom_point(alpha = 0.4) +
  theme_bw()  +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", r.digits = 4)

combo = ggarrange (relabun_KO ,relabun_egg ,nrow = 1 )

ggsave ("Figures/Figure_Sup_NormCompare.png", plot = combo, width =6, height = 3, dpi = 300)



########################
# the hypothesis testing. 
########################

# NOW we filter !!!!
counts(dds)
#keep_genes <- rowSums (counts(dds) > 5  ) > 3 # keep genes with at more than 5 genes in at least 4 samples.

32*.2
table(rowSums (counts(dds)!=0) > 5) # must be in at least 6 samples
table(rowSums (counts(dds)!=0) > 3) # must be in at least 4 samples 
table(rowSums (counts(dds)!=0) > 4) 
table(rowSums (counts(dds) > 5  ) > 3 ) # okay so this knocks off like a few? 
# KEEP 7432, toss 4994
# you can filter genes after normalization!!!!!

# this would be keeping genes that are present in at least 4 sampless 
keep_genes = rowSums (counts(dds)!=0) > 3
table(keep_genes)
dds <- dds[ keep_genes, ] 
dim(dds) # 7593   32 # GOING WITH THIS!@@@@
# 7469   32 wgen taking out the ribosomal


dds= DESeq(dds)

names(dds)

# results : RainTrt
res.trt <- results(dds, name= "RainTrt_Drought_vs_Ambient") # get results. You can customize this
hist(res.trt$padj)

# shrink log fold changes association with condition:
res_lfcShrink.trt <- lfcShrink(dds, coef=2, type = "apeglm")
hist(res.trt$log2FoldChange)  # fucked
hist(res_lfcShrink.trt$log2FoldChange) # thsssee aare incredibly shrunken now....

plotMA(res_lfcShrink.trt ) # sigs have very low LFC

res_lfcShrink.trt %>% as.data.frame %>% filter(padj<0.05) %>% nrow() # 98 significant p value, but very low LFC
res_lfcShrink.trt %>% as.data.frame %>% filter(padj<0.05 & abs(log2FoldChange) >=1) %>% nrow()

#res_lfcShrink_sig.trt <- res_lfcShrink.trt %>%  as.data.frame() %>% filter(padj<0.05  & abs(log2FoldChange) >= 1) \
#nrow(res_lfcShrink_sig.trt) # only one gene with LFC > 1....
res_lfcShrink_sig.trt <- res_lfcShrink.trt %>%  as.data.frame() %>% filter(padj<0.05 &  abs(log2FoldChange) >= 0.5)
hist(res_lfcShrink_sig.trt$log2FoldChange, breaks=200) # these are useless. 
#write.csv(res_lfcShrink_sig.trt, "L2/A_vs_D_deseq_sig_028925.csv")

#write.csv(res_lfcShrink.trt, "L2/A_vs_D_deseq_all_028925.csv")

###########
# results : WP under ambient - LINEAR
resultsNames(dds)
res.amb.lin <- results(dds, name= "RainTrtAmbient.poly.logWP..2.1") # get results. You can customize this
hist(res.amb.lin$padj) 

res_lfcShrink_amb_lin <- lfcShrink(dds, coef=3, type = "apeglm")
hist(res.amb.lin$log2FoldChange)  # 
hist(res_lfcShrink_amb_lin$log2FoldChange) # better!

plotMA(res_lfcShrink_amb_lin ) # good.


res_lfcShrink_sig_amb_lin <- res_lfcShrink_amb_lin %>%  as.data.frame() %>% filter(padj<0.05 &  abs(log2FoldChange) >= 0.5)



###########
# results : WP under drought - LINEAR
resultsNames(dds)
res.drt.lin <- results(dds, name= "RainTrtDrought.poly.logWP..2.1") # get results. You can customize this
res.drt.lin
hist(res.drt.lin$pvalue)
hist(res.drt.lin$padj) # lots


res_lfcShrink_drt_lin <- lfcShrink(dds, coef=4, type = "apeglm")
#apeglm was what the tutorial recommended but apeglm not working
hist(res.drt.lin$log2FoldChange)  # fucked
hist(res_lfcShrink_drt_lin$log2FoldChange) # better!

plotMA(res_lfcShrink_drt_lin ) # some sig have very low LFCs. 

res.drt.lin %>% as.data.frame %>% filter(padj<0.05) %>% nrow() # 1141 significant
res_lfcShrink_drt_lin %>%  as.data.frame() %>% filter(padj<0.05  & abs(log2FoldChange) >= 2)  %>% nrow() # only 100 > 2
res_lfcShrink_drt_lin %>%  as.data.frame() %>% filter(padj<0.05  & abs(log2FoldChange) >= 1)  %>% nrow() # 438 > 1
res_lfcShrink_drt_lin %>%  as.data.frame() %>% filter(padj<0.01) %>% nrow() #


#res_Amb_lfcShrink_sig <- res_Amb_lfcShrink %>%  as.data.frame() %>% filter(padj<0.05) # & ab s(log2FoldChange) > 2) 
res_lfcShrink_sig_drt_lin <- res_lfcShrink_drt_lin %>%  as.data.frame() %>% filter(padj<0.05 &  abs(log2FoldChange) >= 0.5)
dim(res_lfcShrink_sig_drt_lin) # damn, 1186. 
# 438 if 1, 100 if 2

###########
# results : WP under ambient - quadratic
resultsNames(dds)
res.amb.quad <- results(dds, name= "RainTrtAmbient.poly.logWP..2.2") # get results. You can customize this

res_lfcShrink_amb_quad <- lfcShrink(dds, coef=5, type = "apeglm") # sone riws diudnit cinverge



plotMA(res_lfcShrink_amb_quad ) # good.

res_lfcShrink_sig_amb_quad <- res_lfcShrink_amb_quad %>%  as.data.frame() %>% filter(padj<0.05 &  abs(log2FoldChange) >= 0.5)



###########
# results : WP under drought - QUAD
resultsNames(dds)
res.drt.quad <- results(dds, name= "RainTrtDrought.poly.logWP..2.2") # get results. You can customize this
res.drt.quad
hist(res.drt.quad$pvalue)
hist(res.drt.quad$padj) # lots


res_lfcShrink_drt_quad <- lfcShrink(dds, coef=6, type = "apeglm")

hist(res.drt.quad$log2FoldChange)  # fucked
hist(res_lfcShrink_drt_quad$log2FoldChange) # better!

plotMA(res_lfcShrink_drt_quad ) # some sig have very low LFCs. 

res.drt.quad %>% as.data.frame %>% filter(padj<0.05) %>% nrow() # 1141 significant
res_lfcShrink_drt_quad %>%  as.data.frame() %>% filter(padj<0.05  & abs(log2FoldChange) >= 2)  %>% nrow() # only 100 > 2
res_lfcShrink_drt_quad %>%  as.data.frame() %>% filter(padj<0.05  & abs(log2FoldChange) >= 1)  %>% nrow() # 438 > 1
res_lfcShrink_drt_quad %>%  as.data.frame() %>% filter(padj<0.01) %>% nrow() #

res_lfcShrink_sig_drt_quad <- res_lfcShrink_drt_quad %>%  
  as.data.frame() %>% filter(padj<0.05 &  abs(log2FoldChange) >= 0.5)
dim(res_lfcShrink_sig_drt_quad) # damn, 1186. 





# combine the four
# SIG ONLY:
res_lfcShrink_sig_amb_lin$TrtSig = "Ambient"
res_lfcShrink_sig_drt_lin$TrtSig = "Drought"
res_lfcShrink_sig_amb_quad$TrtSig = "Ambient"
res_lfcShrink_sig_drt_quad$TrtSig = "Drought"

res_lfcShrink_sig_amb_lin$Term = "Linear"
res_lfcShrink_sig_drt_lin$Term = "Linear"
res_lfcShrink_sig_amb_quad$Term = "Quadratic"
res_lfcShrink_sig_drt_quad$Term = "Quadratic"

res_lfcShrink_sig_amb_lin <- rownames_to_column(res_lfcShrink_sig_amb_lin, var = "KO")
res_lfcShrink_sig_drt_lin <- rownames_to_column(res_lfcShrink_sig_drt_lin, var = "KO")
res_lfcShrink_sig_amb_quad<- rownames_to_column(res_lfcShrink_sig_amb_quad, var = "KO")
res_lfcShrink_sig_drt_quad <- rownames_to_column(res_lfcShrink_sig_drt_quad, var = "KO")


res_lfcShrink_sig_all = rbind(res_lfcShrink_sig_amb_lin, res_lfcShrink_sig_drt_lin,
                              res_lfcShrink_sig_amb_quad, res_lfcShrink_sig_drt_quad)
head(res_lfcShrink_sig_all)
dim(res_lfcShrink_sig_all)
70+24+146+43
write.csv(res_lfcShrink_sig_all, "L2/DESeq2_sig_genes_031325.csv", row.names=F)

min(abs(res_lfcShrink_sig_all$log2FoldChange))
max(res_lfcShrink_sig_all$pvalue )

# NOW all genes including nonsig
res_lfcShrink_amb_lin$TrtSig = "Ambient"
res_lfcShrink_drt_lin$TrtSig = "Drought"
res_lfcShrink_amb_quad$TrtSig = "Ambient"
res_lfcShrink_drt_quad$TrtSig = "Drought"

res_lfcShrink_amb_lin$Term = "Linear"
res_lfcShrink_drt_lin$Term = "Linear"
res_lfcShrink_amb_quad$Term = "Quadratic"
res_lfcShrink_drt_quad$Term = "Quadratic"

res_lfcShrink_amb_lin <- rownames_to_column(as.data.frame(res_lfcShrink_amb_lin), var = "KO")
res_lfcShrink_drt_lin <- rownames_to_column(as.data.frame(res_lfcShrink_drt_lin), var = "KO")
res_lfcShrink_amb_quad<- rownames_to_column(as.data.frame(res_lfcShrink_amb_quad), var = "KO")
res_lfcShrink_drt_quad <- rownames_to_column(as.data.frame(res_lfcShrink_drt_quad), var = "KO")
#intersect(res_lfcShrink_amb$KO, res_lfcShrink_drt$KO)
# 224 genes are shared . (238)( now 236 with no ribosomes)

res_lfcShrink_all = rbind(res_lfcShrink_amb_lin, res_lfcShrink_drt_lin,
                              res_lfcShrink_amb_quad, res_lfcShrink_drt_quad)
head(res_lfcShrink_all)
dim(res_lfcShrink_all)

write.csv(res_lfcShrink_all, "L2/DESeq2_genes_ALL_00721.csv", row.names=F)



########### PLOT shape of curve  for significant genes. 


ko_desc <- read.csv("L1/Extra files for analysis/ko_descriptions.csv")


##### AMBIENT  - lienar
significant_lin_amb  = deseq_counts_long  %>% 
  mutate (CoreUniqueID = gsub( "_.*" , "", SampleID)  ) %>% 
  merge (meta_new, by = "CoreUniqueID") %>% 
  filter (RainTrt == "Ambient") %>%
  merge (res_lfcShrink_sig_amb_lin %>% select(KO, log2FoldChange) , by = "KO") %>%   # merge with ambient. 
  merge (ko_desc, by = "KO", all.x = T) %>% 
  mutate(short1 = gsub("\\;.*", "", ko_description)) %>%  # add short label
  mutate(short =  gsub ("E[0-9]\\..*" , "" , gsub ("[A-Za-z0-9]{6,}", "", short1)) ) %>% 
  filter (!KO %in% res_lfcShrink_sig_amb_quad$KO ) %>%  # remove if in quad 
  mutate(lfc_posorneg = log2FoldChange>0)
  
significant_lin_amb_fig = significant_lin_amb %>% 
  ggplot(aes (x = ActualWP, y = deseq_counts )) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_x_continuous(trans = ggforce::trans_reverser('log10'), labels = label_trueminus , expand = c (0,0)) +
  geom_text(
    aes(x = Inf, y = Inf, label = KO , color = lfc_posorneg) ,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
    size = 4) +
    geom_text(
      data = significant_lin_amb %>% select (KO,deseq_counts,short,lfc_posorneg) %>% group_by(short,KO,lfc_posorneg) %>% summarize (max = max(deseq_counts) ), 
      aes(x = Inf, y = max*.85, label = short , color = lfc_posorneg) ,
      hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
      vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
      size = 4) +
  geom_text(
    aes(x = Inf, y = Inf, label = short , color = lfc_posorneg) ,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = -10 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
    size = 4) +
  facet_wrap(~KO, scales = "free") + theme_bw() + 
  scale_color_manual(values = c ("red", "blue")) +
  theme( #panel.spacing.x =   unit(-1, "lines"), # 5
       # panel.spacing.y =   unit(-0.5, "lines"), # 3
         axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
         axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_blank() ,
        legend.position = "none") 

significant_lin_amb_fig

ggsave ("Figures/Figure_Sup_SigGenes_A_L.png", plot = significant_lin_amb_fig, width =12, height = 12, dpi = 300)

##### AMBIENT  - quad

significant_quad_amb = deseq_counts_long  %>% 
  mutate (CoreUniqueID = gsub( "_.*" , "", SampleID)  ) %>% 
        merge (meta_new, by = "CoreUniqueID") %>% 
  filter (RainTrt == "Ambient") %>%
  merge (res_lfcShrink_sig_amb_quad %>% select(KO, log2FoldChange) , by = "KO") %>% 
  merge (ko_desc, by = "KO", all.x = T) %>% 
  mutate(short1 = gsub("\\;.*", "", ko_description)) %>%  # add short label
  mutate(short =  gsub ("E[0-9]\\..*" , "" , gsub ("[A-Za-z0-9]{6,}", "", short1)) ) %>% 
  mutate(lfc_posorneg = log2FoldChange>0) 

significant_quad_amb_fig = significant_quad_amb %>% 
  ggplot(aes (x = ActualWP, y = deseq_counts )) +
  geom_point(size = 0.5, alpha = 0.5) +
 # geom_point(data = predicted_drt_w_lfc_wquad_long, aes (x=ActualWP, y = predicted_counts ), size = 0.5, alpha = 0.5 , color= "blue",
       #      inherit.aes = F) +
  #scale_x_log10() +
  scale_x_continuous(trans = ggforce::trans_reverser('log10'), labels = label_trueminus) +
    geom_text(
      data = significant_quad_amb %>% select (KO,deseq_counts,short,lfc_posorneg) %>% group_by(short,KO,lfc_posorneg) %>% summarize (max = max(deseq_counts) ), 
      aes(x = Inf, y = max*.85, label = short , color = lfc_posorneg) ,
      hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
      vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
      size = 4) +
  geom_text(
    aes(x = Inf, y = Inf, label = KO , color = lfc_posorneg) ,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = 1.1   # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
  ) +
  facet_wrap(~KO, scales = "free") + theme_bw() + 
  scale_color_manual(values = c ("red", "blue")) +
  theme(panel.spacing.x =   unit(0, "lines"),
        panel.spacing.y =   unit(0, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1 , size = 7),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_blank() ,
        legend.position = "none")




ggsave ("Figures/Figure_Sup_SigGenes_A_Q.png", plot = significant_quad_amb_fig, width =8, height = 8, dpi = 300)



# DROUght - ambient 
significant_lin_drt = deseq_counts_long  %>% 
  mutate (CoreUniqueID = gsub( "_.*" , "", SampleID)  ) %>% 
  merge (meta_new, by = "CoreUniqueID") %>% 
  filter (RainTrt == "Drought") %>%
  merge (res_lfcShrink_sig_drt_lin %>% select(KO, log2FoldChange) , by = "KO") %>%   # merge with ambient. 
  filter (!KO %in% res_lfcShrink_sig_drt_quad$KO ) %>%  # remove if in quad 
  merge (ko_desc, by = "KO", all.x = T) %>% 
  mutate(short1 = gsub("\\;.*", "", ko_description)) %>%  # add short label
  mutate(short =  gsub ("E[0-9]\\..*" , "" , gsub ("[A-Za-z0-9]{6,}", "", short1)) ) %>% 
  mutate(lfc_posorneg = log2FoldChange>0)

significant_lin_drt_fig = significant_lin_drt %>% 
  ggplot(aes (x = ActualWP, y = deseq_counts )) +
  geom_point(size = 0.5, alpha = 0.5) +
  # geom_point(data = predicted_drt_w_lfc_wquad_long, aes (x=ActualWP, y = predicted_counts ), size = 0.5, alpha = 0.5 , color= "blue",
  #      inherit.aes = F) +
  #scale_x_log10() +
  scale_x_continuous(trans = ggforce::trans_reverser('log10'), labels = label_trueminus , expand = c (0,0)) +
    geom_text(
      data = significant_lin_drt %>% select (KO,deseq_counts,short,lfc_posorneg) %>% group_by(short,KO,lfc_posorneg) %>% summarize (max = max(deseq_counts) ), 
      aes(x = Inf, y = max*.85, label = short , color = lfc_posorneg) ,
      hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
      vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
      size = 4) +
  geom_text(
    aes(x = Inf, y = Inf, label = KO , color = lfc_posorneg) ,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
    size = 4) +
  facet_wrap(~KO, scales = "free") + theme_bw() + 
  scale_color_manual(values = c ("red", "blue")) +
  theme( #panel.spacing.x =   unit(-1, "lines"), # 5
    # panel.spacing.y =   unit(-0.5, "lines"), # 3
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_blank() ,
    legend.position = "none") 

#significant_lin_drt_fig

ggsave ("Figures/Figure_Sup_SigGenes_D_L.png", plot = significant_lin_amb_fig, width =12, height = 12, dpi = 300)



# DROUGHT!!  QUAD !!!!!!!!!!!
significant_quad_drt =deseq_counts_long  %>% 
  mutate (CoreUniqueID = gsub( "_.*" , "", SampleID)  ) %>% 
  merge (meta_new, by = "CoreUniqueID") %>% 
  filter (RainTrt == "Drought") %>%
  merge (res_lfcShrink_sig_drt_quad %>% select(KO, log2FoldChange) , by = "KO") %>% 
  merge (ko_desc, by = "KO", all.x = T) %>% 
  mutate(short1 = gsub("\\;.*", "", ko_description)) %>%  # add short label
  mutate(short =  gsub ("E[0-9]\\..*" , "" , gsub ("[A-Za-z0-9]{6,}", "", short1)) )  %>%
  mutate(lfc_posorneg = log2FoldChange>0)

significant_quad_drt_fig = significant_quad_drt %>% 
  ggplot(aes (x = ActualWP, y = deseq_counts )) +
  geom_point(size = 0.5, alpha = 0.5) +
  # geom_point(data = predicted_drt_w_lfc_wquad_long, aes (x=ActualWP, y = predicted_counts ), size = 0.5, alpha = 0.5 , color= "blue",
  #      inherit.aes = F) +
  #scale_x_log10() +
  scale_x_continuous(trans = ggforce::trans_reverser('log10'), labels = label_trueminus) +
  geom_text(
    data = significant_quad_drt %>% select (KO,deseq_counts,short,lfc_posorneg) %>% group_by(short,KO,lfc_posorneg) %>% summarize (max = max(deseq_counts) ), 
    aes(x = Inf, y = max*.85, label = short , color = lfc_posorneg) ,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = 1.1 ,  # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
    size = 4) +
  geom_text(
    aes(x = Inf, y = Inf, label = KO , color = lfc_posorneg) ,  size = 7,
    hjust = -0.1  , # Adjusts horizontal alignment (1 is flush right, 1.1 adds padding)
    vjust = 1.1   # Adjusts vertical alignment (1 is flush top, 1.5 adds padding)
  ) +
  facet_wrap(~KO, scales = "free") + theme_bw() + 
  scale_color_manual(values = c ("red", "blue")) +
  theme(panel.spacing.x =   unit(0, "lines"),
        panel.spacing.y =   unit(0, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1 , size = 7),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_blank() ,
        legend.position = "none")

ggsave ("Figures/Figure_Sup_SigGenes_D_Q.png", plot = significant_quad_drt_fig, width =11, height = 11, dpi = 300)



