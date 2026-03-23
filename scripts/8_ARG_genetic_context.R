##Project: Wild house mice landscape project
##Aim: ARG genetic context analysis (Assembly-based)
##Author: Víctor Hugo Jarquín-Díaz
##Date: 20.03.2026

##Libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)

##Load required information
##Get sample data
sample_data_all_house<-readRDS(file="data/raw/Mus_sample_data_complete_VHJD.rds")
##Get contigs annotation
rgi.arg.data<- readRDS(file="data/raw/ARG_annotation_contigs_VHJD.rds")
##Get ARG abundance
arg.contigs.norm<- readRDS(file="data/raw/ARG_abundance_contigs_VHJD.rds")

##Color pallets
pal.resm<- c("Efflux"= "#6B4E09", "Target replacement"= "#E4C466",
             "Inactivation"= "#39F2D4", "Target alteration"= "#9EA8DC",
             "Target protection"= "#1DB326","Multiple mechanisms"= "#7A9C99")

pal.abx<-c("Tetracycline"="#994F00", "Glycopeptide"="#006CD1",
           "Aminoglycoside"="#B672B8","Fluoroquinolone"="#FCB71B",
           "Phosphonic acid"="#24D582", "Peptide"="#F3CA95",
           "Diaminopyrimidine"="#748B31", "Phenicol"="#C5B51B",
           "MLS"="#89BDB6", "Beta lactam"="#E651BA",
           "Multidrug"="#A2A2A2","Other"="#e2dfe1")

##Retrain only those very well identified ARG in long contigs
rgi.arg.data%>%
  dplyr::filter(Best_Identities>=80&Length_of_refseq>=80)%>%
  dplyr::select(c(Sample_ID, Best_Hit_ARO, Best_Identities, Length_of_refseq,
                  Drug_Class, Resistance_Mechanism, ARO, AMR_Gene_Family))%>%
  distinct(Sample_ID, ARO)%>%
  dplyr::mutate(present = 1)%>%
  pivot_wider(names_from = Sample_ID,
              values_from = present,
              values_fill = 0)-> arg.contigs

##Get plasmer annotations
plasmer.prediction.1<- readRDS("data/raw/Plasmid_annotation_contigs_1_VHJD.rds")
plasmer.prediction.2<- readRDS("data/raw/Plasmid_annotation_contigs_2_VHJD.rds")

##Check for plasmid contigs with ARGs
rgi.arg.data%>%
  dplyr::filter(Best_Identities>=80&Length_of_refseq>=80)%>%
  dplyr::filter(Batch=="Batch_1")%>%
  dplyr::mutate(Contig = str_replace(Contig, "_[^_]*$", ""))%>%
  left_join(by= c("Contig", "Sample_ID"), plasmer.prediction.1)-> tmp

rgi.arg.data%>%
  dplyr::filter(Best_Identities>=80&Length_of_refseq>=80)%>%
  dplyr::filter(Batch=="Batch_2")%>%
  dplyr::mutate(Contig = str_replace(Contig, "_[^_]*$", ""))%>%
  left_join(by= c("Contig", "Sample_ID"), plasmer.prediction.2)%>%
  rbind(tmp)-> plasmid.arg.data

##Get only those genes that were detected in plasmids
plasmid.arg.data%>%
  dplyr::mutate(Length_of_refseq=case_when(Length_of_refseq>=100~100,
                                           T~Length_of_refseq))%>%
  dplyr::mutate(Prediction= case_when(is.na(Prediction)==T~ "chromosome",
                                      T~Prediction))%>%
  dplyr::mutate(Prediction= case_when(Prediction=="chromosome"~ "Chromosome",
                                      T~"Plasmid"))%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted=="beta lactam antibiotics"~ "Beta lactam",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="mls antibiotics"~ "MLS",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Others"))%>%
  ggplot(aes(x= Best_Identities, y= Length_of_refseq))+
  geom_jitter(alpha = 0.5, width = 0.03, size=3,
              aes(fill= Drug_Class_adjusted),
              shape = 21,
              color= "black")+
  labs(x = "Identity (%)", y = "Coverage (%)",
       tag= "a", fill  = "Drug class")+
  scale_fill_manual(values= pal.abx)+
  theme_minimal()+
  facet_grid(~Prediction)+
  guides(fill=guide_legend(nrow=4, byrow=TRUE, alpha=1))+
  theme(text = element_text(size=16), legend.position="bottom")-> Supp00A


plasmid.arg.data%>%
  dplyr::filter(Prediction=="plasmid")%>%
  dplyr::select(!c(Contig, ORF_lenght, Start, Stop))%>%
  distinct()%>%
  dplyr::group_by(ARO)%>%
  dplyr::mutate(mean_identity = mean(Best_Identities),
                mean_coverage= mean(Length_of_refseq))%>%
  dplyr::select(c(Best_Hit_ARO, ARO, Resistance_Mechanism,
                  Drug_Class_adjusted, mean_identity, mean_coverage))%>%
  ungroup()%>%
  dplyr::mutate(Best_Hit_ARO= case_when(Best_Hit_ARO=="tetO"~ "tet(O)",
                                        T~Best_Hit_ARO))%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted=="beta lactam antibiotics"~ "Beta lactam",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="mls antibiotics"~ "MLS",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Others"))%>%
  distinct()-> tmp

arg.contigs.norm%>%
  dplyr::mutate(total_positive= n_distinct(Sample_ID))%>%
  dplyr::group_by(ARO)%>%
  dplyr::mutate(n_samples = n_distinct(Sample_ID),
                proportion = n_samples/total_positive*100,
                mean_reads= mean(arg_mapped),
                mean_fpkm= mean(FPKM))%>%
  ungroup()%>%
  dplyr::select(c(ARO, n_samples, proportion,
                  mean_reads, mean_fpkm))%>%
  distinct()%>%
  dplyr::filter(ARO%in%tmp$ARO)%>%
  left_join(tmp, by = "ARO")%>%
  dplyr::mutate(Best_Hit_ARO= fct_reorder(Best_Hit_ARO, proportion))%>%
  ggplot(aes(x= proportion, y= Best_Hit_ARO))+
  geom_point(shape= 21, aes(fill= Drug_Class_adjusted, size=mean_identity), color= "black")+
  labs(x = "Occurrence in plasmids (%)",
       tag= "b", fill  = "Drug class", size= "Identity (%)")+
  scale_fill_manual(values= pal.abx)+
  scale_size_continuous(limits = c(70, 100))+
  theme_minimal()+
  guides(fill=guide_legend(override.aes = list(size = 3), ncol = 1, byrow=TRUE))+
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position="left")-> Supp00B.1


arg.contigs.norm%>%
  dplyr::mutate(total_positive= n_distinct(Sample_ID))%>%
  dplyr::group_by(ARO)%>%
  dplyr::mutate(n_samples = n_distinct(Sample_ID),
                proportion = n_samples/total_positive*100,
                mean_reads= mean(arg_mapped),
                mean_fpkm= mean(FPKM))%>%
  ungroup()%>%
  dplyr::select(c(ARO, n_samples, proportion,
                  mean_reads, mean_fpkm))%>%
  distinct()%>%
  dplyr::filter(ARO%in%tmp$ARO)%>%
  left_join(tmp, by = "ARO")%>%
  dplyr::mutate(Best_Hit_ARO= fct_reorder(Best_Hit_ARO, proportion))%>%
  ggplot(aes(x= log10(mean_reads), y= Best_Hit_ARO))+
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  labs(x = "Abundance", tag = "  ")+
  theme_minimal()+
  theme(text = element_text(size=16),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())->Supp00B.2

##Merge two figures
Supp00B<-ggarrange(Supp00B.1, Supp00B.2,  ncol = 2, nrow = 1, widths = c(5, 1))


