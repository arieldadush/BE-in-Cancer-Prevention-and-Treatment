
########################################################################
# This file contains all the code for analyzing the germline component of the paper

# paper:The context-specific role of germline pathogenicity in tumorigenesis
# https://www.nature.com/articles/s41588-021-00949-1
# website for MSK genes and mutations:  https://www.signaldb.org/
# upload the msk table to "msk_entire_info"
# upload (main_clinvar_table.tsv) to "clinvar_pathogenic_mutations" 

two_color_pallate = c("#ff968a","#393b45")
library(readxl)
library(tidyr)
 
# upload the msk_entire_info dataframe. 
# FS del = 139, FS ins=47, IF del=3, Missense=88, Nonsense=119, Nonstart=3, Other=1, Splice=50 --> SNVs
msk_entire_info = msk_entire_info %>% filter(mismatch!="del") # 1 mutation was mistakenly included ---> 239 SNVs
msk_entire_info  = msk_entire_info %>% 
  mutate(ADAR=case_when(
    ((mismatch=="G>A" & strand == "+") |(mismatch=="C>T"& strand == "-"))~TRUE,TRUE~FALSE)) 

#142 SNVs can be edited, here they are with all added calcs
msk_entire_info_BE  = msk_entire_info %>% filter(ADAR)

#Add the AlphaMissense analysis - bystander edits analysis  
#!!! Warning: very heavy run: 
dat.big <- data.table::fread("/Users/merdlerrabinor2/Downloads/AlphaMissense_hg19.tsv.gz", verbose = T)
saveRDS(dat.big, file = "AlphaMissense_hg19", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
AlphaMissense = readRDS("AlphaMissense_hg19", refhook = NULL)
remove (dat.big)
AlphaMissense = AlphaMissense %>% filter(am_class=="likely_pathogenic"| am_class=="pathogenic")
AlphaMissense = AlphaMissense %>% mutate(indexer = paste(.$`#CHROM`,POS,REF,ALT,sep="_")) 
AlphaMissense = AlphaMissense %>% select (indexer,am_class) 
AlphaMissense <- AlphaMissense[!duplicated(AlphaMissense), ] # few of the lines are duplicated -> take them off
print (head (AlphaMissense))

mutations_positive_strand_for_DNA = msk_entire_info_BE %>% filter(strand=="+" )
mutations_negative_strand_for_DNA = msk_entire_info_BE %>% filter(strand=="-" ) 

# for the positive strand - 87 mutations
mutations_pos_for_intersection <- data.frame(
  indexer = character(0),   
  clinVarIndexer = character(0)
)

for (mutation_row in 1:nrow(mutations_positive_strand_for_DNA)) {
  
  chor <- mutations_positive_strand_for_DNA$chrom[mutation_row]
  Ref <- substr(mutations_positive_strand_for_DNA$mismatch[mutation_row], 1, 1)
  Alt <- substr(mutations_positive_strand_for_DNA$mismatch[mutation_row], 3, 3)
  ClinVarIndexer <- mutations_positive_strand_for_DNA$HGVSc[mutation_row]
  sequence <- mutations_positive_strand_for_DNA$X21_dna_seq[mutation_row]
  Pos <- mutations_positive_strand_for_DNA$pos[mutation_row]
  
  running_pos = Pos - 10
  for (i in 1:21)
  {
    if (substr(sequence, i, i)==Alt)
    {
      new_row <- data.frame(indexer = paste(chor, running_pos, Alt, Ref, sep = "_"), clinVarIndexer = ClinVarIndexer)
      mutations_pos_for_intersection = rbind (mutations_pos_for_intersection,new_row)
    }
    running_pos = running_pos + 1
  }
}

# for the negative strand - 55 mutations
mutations_neg_for_intersection <- data.frame(
  indexer = character(0),   
  clinVarIndexer = character(0)
)

for (mutation_row in 1:nrow(mutations_negative_strand_for_DNA)) {
  chor <- mutations_negative_strand_for_DNA$chrom[mutation_row]
  Ref <- substr(mutations_negative_strand_for_DNA$mismatch[mutation_row], 1, 1)
  if (Ref=="A" ) {Real_Ref="T"}
  if (Ref=="T" ) {Real_Ref="A"}
  if (Ref=="C" ) {Real_Ref="G"}
  if (Ref=="G" ) {Real_Ref="C"}
  Alt <- substr(mutations_negative_strand_for_DNA$mismatch[mutation_row], 3, 3)
  if (Alt=="A" ) {Real_Alt="T"}
  if (Alt=="T" ) {Real_Alt="A"}
  if (Alt=="C" ) {Real_Alt="G"}
  if (Alt=="G" ) {Real_Alt="C"}
  ClinVarIndexer <- mutations_negative_strand_for_DNA$HGVSc[mutation_row]
  sequence <- mutations_negative_strand_for_DNA$X21_dna_seq[mutation_row]
  Pos <- mutations_negative_strand_for_DNA$pos[mutation_row]
  
  running_pos = Pos + 10
  for (i in 1:21)
  {
    if (substr(sequence, i, i)==Alt)
    {
      new_row <- data.frame(indexer = paste(chor, running_pos, Real_Alt, Real_Ref, sep = "_"), clinVarIndexer = ClinVarIndexer)
      mutations_neg_for_intersection = rbind (mutations_neg_for_intersection,new_row)
    }
    running_pos = running_pos - 1
  }
}

#intersect with alpha-missense annotations
merged_positive <- merge(mutations_pos_for_intersection, AlphaMissense, by = "indexer")
merged_negative <- merge(mutations_neg_for_intersection, AlphaMissense, by = "indexer")
merged_all = rbind (merged_positive,merged_negative)
merged_all <- merged_all %>%
  group_by(clinVarIndexer) %>%
  summarise(count_pathogenic_bystander_edits = n())

#  how many local-of-targets are possible
library(stringr)
msk_entire_info_BE = msk_entire_info_BE %>%
  mutate(count_RNA_bystander_edits = str_count(X21_dna_seq, str_sub(mismatch, start = 3, end = 3)))  

#add the results to a new column called "count_pathogenic_bystander_edits"
msk_entire_info_BE = merge(msk_entire_info_BE,merged_all, by.x = "HGVSc", by.y = "clinVarIndexer", all.x = TRUE)
#replace with 0 if none was found
msk_entire_info_BE <- msk_entire_info_BE %>%
  mutate(count_pathogenic_bystander_edits = ifelse(is.na(count_pathogenic_bystander_edits), 0, count_pathogenic_bystander_edits)) %>%
  rename(count_RNA_pathogenic_bystander_edits = count_pathogenic_bystander_edits)

#local off-target :
msk_entire_info_BE = msk_entire_info_BE%>% 
  mutate (count_RNA_bystander_edits = as.numeric(as.factor(count_RNA_bystander_edits)),
          count_RNA_pathogenic_bystander_edits = as.numeric(as.factor(count_RNA_pathogenic_bystander_edits)))

msk_entire_info_endogenous_ADAR = msk_entire_info_BE %>% filter(non_cas_editing)
msk_entire_info_gnomAD = msk_entire_info_BE %>% filter (!is.na(gnomAD))

mean (msk_entire_info_endogenous_ADAR$count_RNA_bystander_edits)
mean (msk_entire_info_endogenous_ADAR$count_RNA_pathogenic_bystander_edits)
table (msk_entire_info_endogenous_ADAR$count_RNA_pathogenic_bystander_edits)

# calculate to motif of each in msk_entire_info_endogenous_ADAR
msk_entire_info_endogenous_ADAR   = msk_entire_info_endogenous_ADAR %>% mutate(motifBase = ifelse (strand == '+',substr(X21_dna_seq, 10, 10),substr(X21_dna_seq, 12, 12)))
motif = nrow(subset(msk_entire_info_endogenous_ADAR, (motifBase != "G")))
print (motif)  # The ADAR motif is seen in 28 (58%) of the Endogenous-ADAR-compatible variants

# graphs:

#split mutations in which more than one cancer type
msk_single_cancers_per_mutation= msk_entire_info %>% filter (!grepl(",",Cancer_type))
msk_multiple_cancers_per_mutation = msk_entire_info %>% filter (grepl(",",Cancer_type))

for (row_index in 1:nrow(msk_multiple_cancers_per_mutation))
{
  single_mutation = msk_multiple_cancers_per_mutation[row_index,]
  multiple_cancer_per_single_mutation = trimws(unlist(strsplit(single_mutation$Cancer_type,",")))
  for (cancer_index in 1:length(multiple_cancer_per_single_mutation))
  {
    new_raw = single_mutation
    new_raw$Cancer_type = multiple_cancer_per_single_mutation[cancer_index]
    msk_single_cancers_per_mutation = rbind(msk_single_cancers_per_mutation,new_raw)
  }
} 

msk_analyse = msk_single_cancers_per_mutation %>% count(Gene,Cancer_type)  
 
ADAR = rep(c(0),times=nrow(msk_analyse))
non_editing = rep(c(0),times=nrow(msk_analyse))
msk_analyse <- cbind(msk_analyse, ADAR,non_editing) 

for (ind in 1:nrow(msk_analyse))
{
  my_row = msk_analyse[ind,]
  sub_for_gene_by_tumor = msk_single_cancers_per_mutation %>% filter(Gene==my_row$Gene & Cancer_type == my_row$Cancer_type)
  msk_analyse[ind,'ADAR'] = nrow( sub_for_gene_by_tumor %>% filter(ADAR==TRUE))
  msk_analyse[ind,'non_editing'] = nrow( sub_for_gene_by_tumor %>% filter(ADAR==FALSE))
}    

#install.packages("scatterpie")
#library(scatterpie)
#install.packages("wesanderson")
#install.packages("formattable")
library(formattable)
library(scatterpie)

msk_analyse = msk_analyse %>% 
  mutate (gene_n = as.numeric(as.factor(Gene)),cancer_n = as.numeric(as.factor(Cancer_type)))

cancer_types = sort(unique(msk_analyse$Cancer_type))
genes = unique(msk_analyse$Gene)
# calculate the number of mutations to edit by ADAR
edit_percentage = (nrow(msk_entire_info %>% filter (ADAR==TRUE)))*100/nrow(msk_entire_info)
edit_percentage= formattable(edit_percentage,format="f",digits=1)

msk_analyse = msk_analyse %>% mutate(radius = case_when(n<3 ~ 0.3,(n>2&n<5)~0.5,(n>4&n<9)~0.6,TRUE~0.7))

library(ggpattern)
ggplot() + geom_scatterpie(data=msk_analyse, aes(y=cancer_n, x=gene_n, r=radius),
                           cols = c("ADAR", "non_editing")) +
  scale_x_continuous(breaks=c(1:35), labels=genes) + 
  scale_y_continuous(breaks=c(1:38), labels=cancer_types) + 
  labs(x="Gene", y="Cancer type",title = paste0(edit_percentage,"% of the germline pathogenic high-penetrance SNVs \n amenable by Endo-ADAR")) + 
  coord_fixed()+ 
  scale_fill_manual(values = two_color_pallate)+
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+
  theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black"),
        axis.title = element_blank(),  
        panel.grid.minor = element_line(color = "transparent"))

########################################################################################

# CPG - NCBI gene (screening for kids)
# https://www.ncbi.nlm.nih.gov/gtr/tests/562503/
cpg_ncbi= c("CDKN1C","CEBPA","DIS3L2","ETV6","EXT1","EXT2","SAMD9L","FH","SORT",
            "KIF1B", "DICER1","ALK", "POT1","GATA2","GPC3", "MSH6","APC","HRAS","EPCAM",
            "SMAD4","MAX","MEN1","MLH1","MSH2","NBN","ATM","NF1","NF2","DDX41","SUFU",
            "PMS2","SDHAF2","TMEM127","PRKAR1A","PTCH1","PTEN","RB1","REST","RET","SDHA",
            "SDHB","SDHC","SDHD","BLM","BMPR1A","SMARCA4","SMARCB1","SMARCE1", "SRP72","STK11",
            "TERC","TERT","TP53","TSC1","TSC2","VHL","WT1","CDC73","AXIN2","BAP1","RUNX1","PTCH2","PHOX2B","TRIP13")

clinvar_pathogenic_ncbi = clinvar_pathogenic_mutations %>% 
  filter(gene %in% cpg_ncbi) %>% 
  filter(grepl("C0027672",phenotype))  # 2820 mutations 

# changing into 2 categories : ADAR or non_editing
clinvar_pathogenic_ncbi = clinvar_pathogenic_ncbi %>% mutate(EditingOptionByCas = case_when(
  (mutation=="G>A")~"ADAR",
  TRUE ~ "non_editing"))

## building the tables for the analysis
clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_ncbi %>% filter(EditingOptionByCas == "ADAR" )
clinvar_pathogenic_ncbi_gnomAD_ADAR = clinvar_pathogenic_endo_ADAR %>% filter(!is.na(gnomAD_freq))
mean(clinvar_pathogenic_ncbi_gnomAD_ADAR$gnomAD_freq)

# bystander edits - statistics:
#Endogenous-ADAR: 
clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_endo_ADAR %>% mutate(DNA_21_bases = str_sub(gRNA_41bases,start = 11, end = 31))
mutations_positive_strand_for_DNA = clinvar_pathogenic_endo_ADAR %>% filter(strand=="+" )  
mutations_negative_strand_for_DNA = clinvar_pathogenic_endo_ADAR %>% filter(strand=="-" )  
 
# for the positive strand
mutations_pos_for_intersection <- data.frame(
  indexer = character(0),   
  clinVarIndexer = character(0)
)

for (mutation_row in 1:nrow(mutations_positive_strand_for_DNA)) {
 
  split_location <- strsplit(mutations_positive_strand_for_DNA$location[mutation_row], ":")
  chor <-  trimws(split_location[[1]][1])
  Pos <-  as.numeric(trimws(split_location[[1]][2]))
  sequence <- mutations_positive_strand_for_DNA$DNA_21_bases[mutation_row]
  
  running_pos = Pos - 10
  for (i in 1:21)
  { 
    if (substr(sequence, i, i)=="A")
    {
      new_row <- data.frame(indexer = paste(chor, running_pos, "A", "G", sep = "_"), clinVarIndexer = mutations_positive_strand_for_DNA$orig_name[mutation_row])
      mutations_pos_for_intersection = rbind (mutations_pos_for_intersection,new_row)
    }
    running_pos = running_pos + 1
  }
}

# for the negative strand
mutations_neg_for_intersection <- data.frame(
  indexer = character(0),   
  clinVarIndexer = character(0)
)
 
for (mutation_row in 1:nrow(mutations_negative_strand_for_DNA)) {
  
  split_location <- strsplit(mutations_positive_strand_for_DNA$location[mutation_row], ":")
  chor <-  trimws(split_location[[1]][1])
  Pos <-  as.numeric(trimws(split_location[[1]][2]))
  
  sequence <- mutations_positive_strand_for_DNA$DNA_21_bases[mutation_row]
  
  running_pos = Pos + 10
  for (i in 1:21)
  {
    if (substr(sequence, i, i)=="A")
    {
      new_row <- data.frame(indexer = paste(chor, running_pos, "T", "C", sep = "_"), clinVarIndexer = mutations_positive_strand_for_DNA$orig_name[mutation_row])
      mutations_neg_for_intersection = rbind (mutations_neg_for_intersection,new_row)
    }
    running_pos = running_pos - 1
  }
}


#warning: very heavy run:
dat.big <- data.table::fread(paste0("Desktop/BE genetic project","/tables/AlphaMissense_hg38.tsv.gz"), verbose = T)
saveRDS(dat.big, file = "AlphaMissense_hg38", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
AlphaMissense = readRDS("AlphaMissense_hg38", refhook = NULL)
remove (dat.big)
AlphaMissense = AlphaMissense %>% filter(am_class=="likely_pathogenic"| am_class=="pathogenic")
AlphaMissense = AlphaMissense %>% mutate(indexer = paste(.$`#CHROM`,POS,REF,ALT,sep="_")) 
AlphaMissense = AlphaMissense %>% select (indexer,am_class) 
AlphaMissense <- AlphaMissense[!duplicated(AlphaMissense), ] # few of the lines are duplicated -> take them off

#intersect positive strand with alpha-missense annotations
merged_positive <- merge(mutations_pos_for_intersection, AlphaMissense, by = "indexer")
merged_negative <- merge(mutations_neg_for_intersection, AlphaMissense, by = "indexer")
merged_all = rbind (merged_positive,merged_negative)
merged_all <- merged_all %>%
  group_by(clinVarIndexer) %>%
  summarise(count_pathogenic_bystander_edits = n())

# add the bystander edits analysis to the main table
library(stringr)
clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_endo_ADAR %>% 
  mutate(count_final_DNA_bystander_edits = str_count(DNA_21_bases, str_sub(mutation, start = 3, end = 3))) 

clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_endo_ADAR %>% mutate(count_final_DNA_patho_bystander_edits=0)
for (i in 1:nrow(merged_all)) {
  ClinVarIndexer = merged_all$clinVarIndexer[i]
  number_edits = merged_all$count_pathogenic_bystander_edits[i]
  clinvar_pathogenic_endo_ADAR$count_final_DNA_patho_bystander_edits[clinvar_pathogenic_endo_ADAR$orig_name==ClinVarIndexer] = number_edits
}

# save this table for publish:
clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_endo_ADAR %>% select(-phenotype,-RNA_41n_seq_for_guide,-RNA_off_target_hits_85_min_identity,
                                                                       -count_RNA_bystander_edits,-count_RNA_pathogenic_bystander_edits,-NGG_editing_window,
                                                                       -MIT,-count_DNA_bystander_edits,-count_DNA_pathogenic_bystander_edits,-DNA_21_bases,-EditingOptionByCas)
 
clinvar_pathogenic_endo_ADAR = clinvar_pathogenic_endo_ADAR %>% rename(count_bystander_edits = count_final_DNA_bystander_edits, 
                                                                       count_pathogenic_bystander_edits = count_final_DNA_patho_bystander_edits,
                                                                       off_targets_85_min_identity = DNA_off_target_hits_85_min_identity)

# calculate the bystander edits : 
clinvar_pathogenic_endo_ADAR$count_bystander_edits <- as.numeric(as.character(clinvar_pathogenic_endo_ADAR$count_bystander_edits))
clinvar_pathogenic_endo_ADAR$count_pathogenic_bystander_edits <- as.numeric(as.character(clinvar_pathogenic_endo_ADAR$count_pathogenic_bystander_edits))
mean (clinvar_pathogenic_endo_ADAR$count_bystander_edits)
mean (clinvar_pathogenic_endo_ADAR$count_pathogenic_bystander_edits)
table (clinvar_pathogenic_endo_ADAR$count_pathogenic_bystander_edits)
#calculate distant off-target hits
table (clinvar_pathogenic_endo_ADAR$off_targets_85_min_identity)

# generating the figure:
clinvar_pathogenic_ncbi_for_figure = clinvar_pathogenic_ncbi %>% count (gene,EditingOptionByCas)
 
library(ggforce)
library (ggpattern)
# Create a list to store the pie charts
pie_charts <- vector("list", length = 40)
i=1
# Generate the pie charts and store them in the list
for (myGene in unique(clinvar_pathogenic_ncbi_for_figure$gene)) {
  # Generate each pie chart (replace this with your own code)
  myData = clinvar_pathogenic_ncbi_for_figure %>%filter(gene == myGene)
 
  
  if (nrow(myData %>% filter (EditingOptionByCas=="ADAR")) == 0)
  {
    
    new_row <- data.frame(gene = myGene,EditingOptionByCas = "ADAR",n=0)
    myData = rbind (myData,new_row)
  }
  
  if (nrow(myData %>% filter (EditingOptionByCas=="non_editing")) == 0)
  {
    new_row <- data.frame(gene = myGene,EditingOptionByCas = "non_editing",n=0)
    myData = rbind (myData,new_row)
  }
  
  pie_charts[[i]] <- ggplot(myData, aes(x = "", y = n, fill = EditingOptionByCas)) +
    geom_col(color="black", aes(fill = EditingOptionByCas)) + 
    # scale_pattern_discrete(choices = three_patter_choices) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = two_color_pallate)+
    labs(x = NULL, y = NULL, fill = NULL, title = paste0(myGene," (n=", sum(myData$n),")")) +
    theme_void()+
    theme(plot.title = element_text(size=10,color = "black", face = "bold"))+
    theme(legend.position = "none")
  
  i=i+1
}

# Sort the pie charts by size
pie_charts <- pie_charts[order(sapply(pie_charts, function(x) sum(x$data$n)), decreasing = FALSE)]

# Arrange the pie charts using gridExtra
library(cowplot)
plot_grid(plotlist = pie_charts, ncol = 5,nrow = 8)
