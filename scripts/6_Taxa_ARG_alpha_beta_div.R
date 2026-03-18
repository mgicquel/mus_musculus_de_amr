##Project: Wild house mice landscape project
##Aim: Microbiome and ARG alpha and beta diversity analysis
##Author: Víctor Hugo Jarquín-Díaz
##Date: 28.08.2024

##Libraries
library(tidyverse)
library(phyloseq)
library(vroom)
library(rtk)
library(microbiome)

##Load dataset only with those mice that were positive for ARGs
PS.mOTUS.filt <- readRDS("data/raw/PS.mOTUS.filtered.rds")
PS.ARG.raw.filt <- readRDS("data/raw/PS.ARG.filtered.rds")
PS.ARG.FPKM.filt <- readRDS("data/raw/PS.ARG.FPKM.filtered.rds")

##Transform taxa data
PS.mOTUS.filt.clr <- microbiome::transform(PS.mOTUS.filt, "clr") #Centered log ratio transformation
PS.mOTUS.filt.comp<- microbiome::transform(PS.mOTUS.filt, "compositional") ##Total sum transformation
PS.mOTUS.rare<- rarefy_even_depth(PS.mOTUS.filt, rngseed=2020, sample.size=1000, replace=F)

##Estimate alpha diversity
#Rarefy the dataset to an even sequencing depth using rtk
data.rarefied <- rtk(input = PS.mOTUS.filt@otu_table, depth = 1000, repeats = 100,
                     threads = 10, ReturnMatrix = 2)

#Summarize diversity per sample
# Extract summary stats per sample and metric
alpha.rtk <- map_dfr(data.rarefied[["divvs"]], function(sample_entry) {
  #Extract the sample name
  sample_name <- sample_entry[[1]]
  #Extract diversity (skip the sample name)
  metric_names <- names(sample_entry)[-1]
  #Compute mean, sd and 95%ci
  map_dfr(metric_names, function(metric) {
    vals <- sample_entry[[metric]]
    m <- mean(vals)
    s <- sd(vals)
    ci <- qt(0.975, df = length(vals) - 1) * s / sqrt(length(vals))

    tibble(sequencing_name = sample_name,
           Metric = metric,
           Mean = m,
           SD = s,
           CI95_lower = m - ci,
           CI95_upper = m + ci
    )
  })
})

##Get the mean diversity per sample
alpha.rtk%>%
  dplyr::mutate(across(
    c(Mean:CI95_upper), ~ replace_na(.x, 0)))%>%
  dplyr::select(c(sequencing_name, Metric, Mean))%>%
  pivot_wider(names_from = Metric, values_from = Mean)%>%
  dplyr::rename(
    Sample_ID=sequencing_name,
    observed_mOTUs_rtk= richness,
    chao1_mOTUs_rtk= chao1,
    diversity_shannon_mOTUs_rtk= shannon,
    evenness_simpson_mOTUs_rtk=eveness,
    dominance_simpson_mOTUs_rtk=simpson,
    diversity_invsimpson_mOTUs_rtk= invsimpson)->alpha.diversity

##Get alpha diversity from phyloseq
##Raw counts
alpha.raw<- microbiome::alpha(PS.mOTUS.filt, index =c("Observed","Chao1", "Shannon", "Simpson"))

alpha.raw%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::rename("observed_mOTUs_raw" = "observed",
                "chao1_mOTUs_raw" = "chao1",
                "diversity_shannon_mOTUs_raw"= "diversity_shannon",
                "evenness_simpson_mOTUs_raw"= "evenness_simpson",
                "dominance_simpson_mOTUs_raw"= "dominance_simpson")%>%
  left_join(alpha.diversity, by= "Sample_ID")-> alpha.diversity

##Rarefied phyloseq
alpha.rar<- microbiome::alpha(PS.mOTUS.rare, index =c("Observed","Chao1", "Shannon", "Simpson"))

alpha.rar%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::rename("observed_mOTUs_rare" = "observed",
                "chao1_mOTUs_rare" = "chao1",
                "diversity_shannon_mOTUs_rare"= "diversity_shannon",
                "evenness_simpson_mOTUs_rare"= "evenness_simpson",
                "dominance_simpson_mOTUs_rare"= "dominance_simpson")-> alpha.rar

alpha.diversity%>%
  left_join(alpha.rar, by= "Sample_ID")%>%
  dplyr::mutate(across(
    c(observed_mOTUs_rare:dominance_simpson_mOTUs_rare), ~ replace_na(.x, 0)))-> alpha.diversity

rm(alpha.rar, alpha.raw, alpha.rtk, data.rarefied)

##Get alpha diversity from ARGs
alpha.arg <- microbiome::alpha(PS.ARG.FPKM.filt, index =c("Observed","Chao1", "Shannon", "Simpson"))

alpha.arg%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::rename("observed_arg" = "observed",
                "chao1_arg" = "chao1",
                "diversity_shannon_arg"= "diversity_shannon",
                "evenness_simpson_arg"= "evenness_simpson",
                "dominance_simpson_arg"= "dominance_simpson")-> alpha.arg

alpha.diversity%>%
  left_join(alpha.arg, by= "Sample_ID")%>%
  dplyr::mutate(across(
    c(observed_arg:dominance_simpson_arg), ~ replace_na(.x, 0)))-> alpha.diversity

##Rarefied ARG diversity
data.rarefied <- rtk(input = PS.ARG.FPKM.filt@otu_table, depth = 400, repeats = 100,
                     threads = 10, ReturnMatrix = 2)

#Summarize diversity per sample
# Extract summary stats per sample and metric
alpha.rtk <- map_dfr(data.rarefied[["divvs"]], function(sample_entry) {
  #Extract the sample name
  sample_name <- sample_entry[[1]]
  #Extract diversity (skip the sample name)
  metric_names <- names(sample_entry)[-1]
  #Compute mean, sd and 95%ci
  map_dfr(metric_names, function(metric) {
    vals <- sample_entry[[metric]]
    m <- mean(vals)
    s <- sd(vals)
    ci <- qt(0.975, df = length(vals) - 1) * s / sqrt(length(vals))

    tibble(sequencing_name = sample_name,
           Metric = metric,
           Mean = m,
           SD = s,
           CI95_lower = m - ci,
           CI95_upper = m + ci
    )
  })
})

##Get the mean diversity per sample
alpha.rtk%>%
  dplyr::mutate(across(
    c(Mean:CI95_upper), ~ replace_na(.x, 0)))%>%
  dplyr::select(c(sequencing_name, Metric, Mean))%>%
  pivot_wider(names_from = Metric, values_from = Mean)%>%
  dplyr::rename(
    Sample_ID=sequencing_name,
    observed_arg_rtk= richness,
    chao1_arg_rtk= chao1,
    diversity_shannon_arg_rtk= shannon,
    evenness_simpson_arg_rtk=eveness,
    dominance_simpson_arg_rtk=simpson,
    diversity_invsimpson_arg_rtk= invsimpson)%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  dplyr::mutate(across(
    c(observed_arg_rtk:dominance_simpson_arg_rtk), ~ replace_na(.x, 0)))-> alpha.diversity

##Extract the reads mapped to ARGs
##ARG information
arg_data_all<-readRDS(file="data/raw/ARG_annotation_complete_VHJD.rds")

arg_data_all%>%
  dplyr::select(c(Sample_ID,arg_str_mapped))%>%
  distinct()-> arg.mapped

##Get sample data
sample_data_all_house<-readRDS(file="data/raw/Mus_sample_data_complete_VHJD.rds")

##Make some plots
require("ggExtra")
pal.host<- c("House mouse"= "#3E5496", "Other rodent"= "#0A9086")
pal.batch<- c("Batch_1"= "#68855c", "Batch_2"= "#a06177")
pal.transect<- c("MV"="#cc392c", "BB"="#08324e", "BY"="#3f382a")
pal.level<-c("low"="#526a83", "high"="#af6458")
pal.year<- c("2015"= "#E64B35FF", "2016"= "#4DBBD5FF",
             "2017"= "#00A087FF", "2018"= "#F39B7FFF",
             "2019"= "#8491B4FF","2021"="#91D1C2FF","2022"= "#B09C85FF")

##ARG diversity by year (Main Figure 1)
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg= replace_na(observed_arg, 0),
                observed_arg_rtk= replace_na(observed_arg_rtk, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  ggplot(aes(x= Year, y= observed_arg))+
  geom_violin(alpha= 0.5, trim=F)+
  geom_point(shape=21, position=position_jitter(0.2), size=2, aes(fill= Host), color= "black")+
  scale_fill_manual(values = pal.host)+
  geom_hline(yintercept = 100, linetype='dashed', color= "red")+
  xlab("Year")+
  ylab("Number of ARGs")+
  labs(tag= "b")+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> B

##Statistical analysis
require("rstatix")
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  wilcox_test(observed_arg ~ Year, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.arg.year

sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  wilcox_effsize(observed_arg ~ Year, alternative = "two.sided")%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.arg.year)-> stats.test.arg.year ##Table 1

##Year comparisons within batches
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::group_by(Seq_batch)%>%
  wilcox_test(observed_arg ~ Year, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.arg.year.batch

sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::group_by(Seq_batch)%>%
  wilcox_effsize(observed_arg ~ Year, alternative = "two.sided")%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.arg.year.batch)-> stats.test.arg.year.batch ##Supplement table S6

##Over all batch and sequencing effect
require("ggh4x")
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                Filtered_reads= Filtered_reads+1,
                Tax_mapped= Tax_mapped+1,
                observed_arg_rtk= replace_na(observed_arg_rtk, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::select(c(Total_reads, Filtered_reads, Seq_batch, observed_mOTUs_rare, observed_arg, observed_arg_rtk))%>%
  dplyr::rename("Total Reads" = "Total_reads",
                "Filtered Reads" = "Filtered_reads",
                "Taxa (mOTUs)" = "observed_mOTUs_rare",
                "ARGs"= "observed_arg",
                "ARGs rare"= "observed_arg_rtk")%>%
  pivot_longer(!Seq_batch, names_to = "measurment", values_to = "value")%>%
  dplyr::mutate(measurment= fct_relevel(measurment, "Total Reads", "Filtered Reads",
                                        "Taxa (mOTUs)", "ARGs", "ARGs rare"))%>%
  ggplot(aes(x= Seq_batch, y= value))+
  geom_violin(alpha= 0.5, trim=F)+
  geom_point(shape=21, position=position_jitter(0.2), size=1.5, alpha=0.1, aes(fill= Seq_batch), color= "black")+
  scale_fill_manual(values = pal.batch)+
  labs(tag= "a")+
  theme_minimal()+
  facet_wrap(~measurment,  ncol=5, scales= "free")+
  facetted_pos_scales(
    y = list(measurment %in% c("Total Reads", "Filtered Reads") ~
               scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x)),
                             limits = c(10000, 500000000),
                             guide  = "axis_logticks"),
             measurment %in% c("Taxa (mOTUs)", "ARGs", "ARGs rare") ~
               scale_y_continuous(labels = scales::label_number(accuracy = 1),
                                  limits = c(0, 250))))+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")-> Supp0A

##Add the median
Supp0A + stat_summary(fun = median, geom = "point",
    size = 3, shape = 21,  fill = "black", color = "black",
    position = position_dodge(width = 0.6))-> Supp0A

##Statistical comparison
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                Filtered_reads= Filtered_reads+1,
                Tax_mapped= Tax_mapped+1,
                observed_arg_rtk= replace_na(observed_arg_rtk, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::select(c(Total_reads, Filtered_reads, Seq_batch, observed_mOTUs_rare, observed_arg, observed_arg_rtk))%>%
  dplyr::rename("Total Reads" = "Total_reads",
                "Filtered Reads" = "Filtered_reads",
                "Taxa (mOTUs)" = "observed_mOTUs_rare",
                "ARGs"= "observed_arg",
                "ARGs rare"= "observed_arg_rtk")%>%
  pivot_longer(!Seq_batch, names_to = "measurment", values_to = "value")%>%
  dplyr::mutate(measurment= fct_relevel(measurment, "Total Reads", "Filtered Reads",
                                        "Taxa (mOTUs)", "ARGs", "ARGs rare"))%>%
  dplyr::group_by(measurment)%>%
  wilcox_test(value ~ Seq_batch, alternative = "two.sided")%>%
  add_significance()-> stats.test.batch

##Effect size
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                Filtered_reads= Filtered_reads+1,
                Tax_mapped= Tax_mapped+1,
                observed_arg_rtk= replace_na(observed_arg_rtk, 0),
                observed_arg= replace_na(observed_arg, 0))%>%
  dplyr::select(c(Total_reads, Filtered_reads, Seq_batch, observed_mOTUs_rare, observed_arg, observed_arg_rtk))%>%
  dplyr::rename("Total Reads" = "Total_reads",
                "Filtered Reads" = "Filtered_reads",
                "Taxa (mOTUs)" = "observed_mOTUs_rare",
                "ARGs"= "observed_arg",
                "ARGs rare"= "observed_arg_rtk")%>%
  pivot_longer(!Seq_batch, names_to = "measurment", values_to = "value")%>%
  dplyr::mutate(measurment= fct_relevel(measurment, "Total Reads", "Filtered Reads",
                                        "Taxa (mOTUs)", "ARGs", "ARGs rare"))%>%
  dplyr::group_by(measurment)%>%
  wilcox_effsize(value ~ Seq_batch, alternative = "two.sided")%>%
  dplyr::select(c(measurment, effsize, magnitude))%>%
  left_join(stats.test.batch, by= "measurment")%>%
  relocate(c(effsize, magnitude), .after = p.signif)%>%
  dplyr::select(!`.y.`)->stats.test.batch

##Temporal change using rarefied ARG
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg= replace_na(observed_arg, 0),
                observed_arg_rtk= replace_na(observed_arg_rtk, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  ggplot(aes(x= Year, y= observed_arg_rtk))+
  geom_violin(alpha= 0.5, trim=F)+
  geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.25, aes(fill= Seq_batch), color= "black")+
  scale_fill_manual(values = pal.batch)+
  geom_hline(yintercept = 100, linetype='dashed', color= "red")+
  xlab("Year")+
  ylab("Number of ARGs")+
  labs(tag= "b")+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> Supp0B

##Add the median
Supp0B + stat_summary(fun = median, geom = "point",
                      size = 3, shape = 21,  fill = "black", color = "black",
                      position = position_dodge(width = 0.6))-> Supp0B

##Stats
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                Filtered_reads= Filtered_reads+1,
                Tax_mapped= Tax_mapped+1,
                observed_arg_rtk= replace_na(observed_arg_rtk, 0),
                observed_arg= replace_na(observed_arg, 0),
                Year= as.factor(Year),
                Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::select(c(Total_reads, Filtered_reads, Seq_batch, observed_mOTUs_rare, observed_arg_rtk, Year))%>%
  dplyr::rename("Total Reads" = "Total_reads",
                "Filtered Reads" = "Filtered_reads",
                "Taxa (mOTUs)" = "observed_mOTUs_rare",
                "ARG rare" = "observed_arg_rtk")%>%
  pivot_longer(!c(Seq_batch, Year), names_to = "measurment", values_to = "value")%>%
  dplyr::mutate(measurment= fct_relevel(measurment, "Total Reads",
                                        "Filtered Reads", "Taxa (mOTUs)", "ARG rare"))%>%
  dplyr::group_by(measurment)%>%
  wilcox_test(value ~ Year, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.year ##Supplement table S7

##Effect size
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                Filtered_reads= Filtered_reads+1,
                Tax_mapped= Tax_mapped+1,
                observed_arg_rtk= replace_na(observed_arg_rtk, 0),
                observed_arg= replace_na(observed_arg, 0),
                Year= as.factor(Year),
                Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::select(c(Total_reads, Filtered_reads, Seq_batch, observed_mOTUs_rare, observed_arg_rtk, Year))%>%
  dplyr::rename("Total Reads" = "Total_reads",
                "Filtered Reads" = "Filtered_reads",
                "Taxa (mOTUs)" = "observed_mOTUs_rare",
                "ARG rare"= "observed_arg_rtk")%>%
  pivot_longer(!c(Seq_batch, Year), names_to = "measurment", values_to = "value")%>%
  dplyr::mutate(measurment= fct_relevel(measurment, "Total Reads",
                                        "Filtered Reads", "Taxa (mOTUs)","ARG rare"))%>%
  dplyr::group_by(measurment)%>%
  wilcox_effsize(value ~ Year, alternative = "two.sided")%>%
  dplyr::select(c(measurment, group1, group2, effsize, magnitude))%>%
  left_join(stats.test.year)%>%
  relocate(c(effsize, magnitude), .after = p.adj.signif)%>%
  dplyr::select(!`.y.`)->stats.test.year ##Supplement table S7

##Now for non rarefied ARGs
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  wilcox_test(observed_arg ~ Year, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.arg.year

sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  wilcox_effsize(observed_arg ~ Year, alternative = "two.sided")%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.arg.year)-> stats.test.arg.year

##Year comparisons within batches
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::group_by(Seq_batch)%>%
  wilcox_test(observed_arg ~ Year, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.arg.year.batch

sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::group_by(Seq_batch)%>%
  wilcox_effsize(observed_arg ~ Year, alternative = "two.sided")%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.arg.year)-> stats.test.arg.year.batch

##correlation between rarefied and non rarefied ARG richness
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0),
                observed_arg_rtk= replace_na(observed_arg_rtk, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::mutate(arg_level=
                  case_when(observed_arg>=100 ~ "high",
                            TRUE ~ "low"))%>%
  dplyr::mutate(arg_level= fct_relevel(arg_level,
                                       "high", "low"))%>%
  ggplot(aes(x= observed_arg, y= observed_arg_rtk, color=Seq_batch))+
  geom_point(shape=21, alpha= 0.2, position=position_jitter(0.2), size=3, aes(fill= Seq_batch), color= "black")+
  labs(tag= "c", x= "ARGs (unrarefied)", y= "ARGs (rarefied)")+
  scale_fill_manual(values = pal.batch)+
  scale_color_manual(values = pal.batch)+
  coord_cartesian(ylim = c(0, NA))+
  geom_smooth(method = lm, se = T)+
  geom_hline(yintercept = 100, linetype='dotted')+
  geom_vline(xintercept = 100, linetype='dotted')+
  scale_y_continuous(limits = c(0, 150))+
  scale_x_continuous(limits = c(0, 150))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  ggpubr::stat_cor(label.y = c(120, 120),  label.x = 50, method = "spearman",
                   aes(label= paste("rho","'='", after_stat(r), after_stat(p.label), sep= "~` `~")))+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="bl")+
  theme_minimal()+
  facet_wrap(~Year,  ncol=2, scales= "free")+
  theme(text = element_text(size=16),
        legend.position = "none")-> Supp0C

###Check Geographical correlation (only for presentations)
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::mutate(arg_level=
                  case_when(observed_arg>=100 ~ "high",
                            TRUE ~ "low"))%>%
  dplyr::mutate(arg_level= fct_relevel(arg_level,
                                       "high", "low"))%>%
  dplyr::mutate(transect=case_when(Latitude>=53.5 ~ "MV",
                                   Latitude<=50.0 ~ "BY",
                                   TRUE ~ "BB"))%>%
  ggplot(aes(y= observed_arg, x= Latitude, color=transect))+
  geom_point(shape=21, alpha= 0.5, position=position_jitter(0.2), size=3, aes(fill= transect), color= "black")+
  scale_fill_manual(values = pal.transect)+
  scale_color_manual(values = pal.transect)+
  geom_hline(yintercept = 100, linetype='dotted', alpha=0.5)+
  labs(tag= "d", y= "ARGs richness", x= "Latitude")+
  scale_y_continuous(limits = c(0, 210))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  geom_rug(col=rgb(.5,0,0, alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "none")-> tmp.fig

##Number of farms
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::mutate(arg_level=
                  case_when(observed_arg>=100 ~ "high",
                            TRUE ~ "low"))%>%
  dplyr::mutate(arg_level= fct_relevel(arg_level,
                                       "high", "low"))%>%
  dplyr::mutate(transect=case_when(Latitude>=53.5 ~ "MV",
                                   Latitude<=50.0 ~ "BY",
                                   TRUE ~ "BB"))%>%
  dplyr::select(c(transect, Longitude,Latitude))%>%
  dplyr::mutate(farmid= paste(Longitude, "_", Latitude))%>%
  distinct()

##Relationship between Taxa composition and ARG richness level
##Beta diversity
#Compositional (Aitchison distance)
aitch_dist<- vegan::vegdist(t(PS.mOTUS.filt@otu_table),
                            method="aitchison", pseudocount=1)
##For PCA with Aitchitson
ordination.clr <- phyloseq::ordinate(PS.mOTUS.filt.clr, "RDA") #principal components analysis

##Define strata
sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.diversity, by= "Sample_ID")%>%
  left_join(arg.mapped, by= "Sample_ID")%>%
  dplyr::mutate(arg_str_mapped= replace_na(arg_str_mapped, 0))%>%
  dplyr::mutate(Year= as.factor(Year))%>%
  dplyr::mutate(Year= fct_relevel(Year,
                                  "2015", "2016",
                                  "2017", "2018",
                                  "2019", "2021","2022"))%>%
  dplyr::mutate(arg_level=
                  case_when(observed_arg>=100 ~ "high",
                            TRUE ~ "low"))%>%
  dplyr::mutate(arg_level= fct_relevel(arg_level,
                                       "high", "low"))%>%
  dplyr::mutate(transect=case_when(Latitude>=53.5 ~ "MV",
                                   Latitude<=50.0 ~ "BY",
                                   TRUE ~ "BB"))%>%
  dplyr::filter(Sample_ID%in%sample_names(PS.mOTUS.filt))%>%
  column_to_rownames("Sample_ID")-> ord.data

##Save for the farm comparisons only the data from ARG positive mice
saveRDS(ord.data, file="data/raw/Mus_sample_data_farm_VHJD.rds")

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(ord.data, Seq_batch)

mOTU.adonis.ait<- vegan::adonis2(aitch_dist~ Year + Sex + Longitude + Latitude + transect +
                                   arg_level + Tax_mapped,
                                 data = ord.data, na.action = na.fail,
                                 permutations = perm, by="margin")

#PCA Aitchitson Select Axis 1 and 2
seg.data<-data.frame(ordination.clr$CA$u[,1:2])
seg.data%>%
  rownames_to_column("Sample_ID")-> seg.data

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::left_join(seg.data, by="Sample_ID")->seg.data

##Just to have an overview
ggplot() +
  geom_point(data=seg.data, aes(x=PC1,y=PC2, fill= Year, shape=arg_level), size=2) +
  stat_ellipse(data=seg.data, aes(x=PC1,y=PC2, color= Year))+
  guides(fill = guide_legend(override.aes=list(shape=c(21), size= 3)), color= "none")+
  labs(tag= "a", fill  = "Year", shape= "ARG level")+
  theme_minimal()+
  scale_fill_manual(values = pal.year)+
  scale_shape_manual(values = c(24, 25), labels= c("High", "Low"))+
  scale_color_manual(values = pal.year)+
  theme(text = element_text(size=16))+
  annotate("text", x = 0.05, y = 0.090, label= "Permanova (Year)", size= 3)+
  annotate("text", x = 0.05, y = 0.080, label= "Aitchison distance", size= 3)+
  annotate("text", x = 0.05, y = 0.070, label= paste0(label = "R²= ", round(mOTU.adonis.ait$R2[1], digits = 3),
                                                      ", p = ", mOTU.adonis.ait$`Pr(>F)`[1]), color = "black", size= 3)+
  xlab(paste0("PC 1 [", round(ordination.clr$CA$eig[1] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(ordination.clr$CA$eig[2] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))-> Supp1A

##Taxonomic composition
##Barplot by sample
PS.mOTUS.filt%>%
  microbiome::aggregate_taxa(level = "Family") %>%
  microbiome::transform(transform = "compositional")-> PS.mOTUS.filter.Fam

###Get the abundances for all families before agglomerate
data.frame(t(abundances(PS.mOTUS.filter.Fam)))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "Family", values_to = "abundance")%>%
  dplyr::mutate(Family = gsub("\\.", " ", Family))%>%
  dplyr::mutate(Family = gsub("fam ", "fam.", Family))->abund.family

# merge all taxa that are detected rare
PS.mOTUS.filter.Fam <- microbiome::aggregate_rare(PS.mOTUS.filter.Fam, level="Family",
                                                     detection = 0.05, prevalence = 0.1)

otu <- as(otu_table(PS.mOTUS.filter.Fam), "matrix")
otu <- t(otu)
otu <- otu_table(otu, taxa_are_rows = F)
tax <- tax_table(PS.mOTUS.filter.Fam)
# Coerce to data.frame
n <- as.data.frame(tax)
n%>%
  rownames_to_column()%>%
  dplyr::rename(tax = rowname)%>%
  dplyr::mutate(tax= str_replace_all(tax, "NA ", ""))%>%
  dplyr::mutate(Family= str_replace_all(Family, "NA ", ""))%>%
  dplyr::mutate(unique= str_replace_all(unique, "NA ", ""))-> n

j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
j2 <- lapply(j1,'[[',"x") # select for Names

m <- data.frame(unlist(j2))

m%>%
  rownames_to_column()%>%
  dplyr::mutate(rowname = gsub("\\.", "_", rowname))%>%
  dplyr::mutate(rowname= str_replace_all(rowname, "NA ", ""))%>%
  dplyr::mutate(rowname= str_replace_all(rowname, "fam_", "fam\\."))%>%
  separate(rowname, sep = "_", c("Sample_ID","tax"))%>%
  dplyr::mutate(tax= str_replace_all(tax, "NA ", ""))%>%
  dplyr::group_by(Sample_ID)%>%
  dplyr::rename(Abundance = unlist.j2.)%>%
  left_join(n, by="tax")%>%
  arrange(Sample_ID, desc(tax))->m

sample_data_all_house%>%
  rownames_to_column("Sample_ID")%>%
  left_join(m, by="Sample_ID")-> fam.hmhz

rm(j1, j2, m, n)

#plot
fam.hmhz%>%
  dplyr::filter(!is.na(Family))%>%
  dplyr::mutate(tax= fct_relevel(tax,
                                 "Bacteroidaceae",
                                 "Bacteroidales fam. incertae sedis",
                                 "Clostridiaceae",
                                 "Clostridiales fam. incertae sedis",
                                 "Helicobacteraceae",
                                 "Lachnospiraceae",
                                 "Lactobacillaceae",
                                 "Muribaculaceae",
                                 "Odoribacteraceae",
                                 "Oscillospiraceae",
                                 "Prevotellaceae",
                                 "Rikenellaceae",
                                 "Other",
                                 "unknown"))%>%
  ggplot(aes(x=Sample_ID, y=Abundance*100, fill=tax))+
  geom_bar(aes(), stat="identity", position="stack", width=.75) +
  scale_fill_manual(values=c("Bacteroidaceae"="#008ece",
                             "Bacteroidales fam. incertae sedis"="#59c7eb",
                             "Clostridiaceae"="#e0607e",
                             "Clostridiales fam. incertae sedis"="#eca0b2",
                             "Helicobacteraceae"="#077187",
                             "Lachnospiraceae"="#6aaab7",
                             "Lactobacillaceae"="#8e2043",
                             "Muribaculaceae"="#bc7a8f",
                             "Odoribacteraceae"= "#0a9086",
                             "Oscillospiraceae"="#54bfb7",
                             "Prevotellaceae"="#fea090",
                             "Rikenellaceae"="#3e5496",
                             "Other"="#E1E2E5",
                             "unknown"="#B8BCC1")) +
  theme_minimal()+
  labs(tag= "b",
       fill= "Family")+
  ylab("Relative abundance")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=5))+
  theme(text = element_text(size=16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())-> Supp1B

##Get summary for report purposes
fam.hmhz%>%
  dplyr::filter(!is.na(Family))%>%
  dplyr::mutate(tax= fct_relevel(tax,
                                 "Bacteroidaceae",
                                 "Bacteroidales fam. incertae sedis",
                                 "Clostridiaceae",
                                 "Clostridiales fam. incertae sedis",
                                 "Helicobacteraceae",
                                 "Lachnospiraceae",
                                 "Lactobacillaceae",
                                 "Muribaculaceae",
                                 "Odoribacteraceae",
                                 "Oscillospiraceae",
                                 "Prevotellaceae",
                                 "Rikenellaceae",
                                 "Other",
                                 "unknown"))%>%
  dplyr::group_by(tax)%>%
  summarise(mean_abundance= mean(Abundance, na.rm = TRUE),
            sd_abundance= sd(Abundance, na.rm = TRUE),
            n_abundance= length(Abundance),
            se_abundance= sd_abundance/sqrt(n_abundance),
            lower_ci= mean_abundance - qt(1 - (0.05 / 2),
                                          n_abundance - 1) * se_abundance,
            upper_ci= mean_abundance + qt(1 - (0.05 / 2),
                                          n_abundance - 1) * se_abundance)->summarized_taxa_family

write.csv(summarized_taxa_family, "Tables/Taxa_summary_family_Mus.csv")

##Get mOTU prevalence
data.frame(t(prevalence(PS.mOTUS.filt)))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "mOTU", values_to = "prevalence")%>%
  dplyr::select(!row)%>%
  dplyr::left_join(data.frame(PS.mOTUS.filt@tax_table)%>%
                     rownames_to_column("mOTU"), by="mOTU")->prev.motus

##Get mOTU mean abundance
data.frame(t(abundances(PS.mOTUS.filt.comp)))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "mOTU", values_to = "abundance")%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(mOTU) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")%>%
  left_join(prev.motus, by="mOTU")-> prev.abund.motus

##Get mOTU abundance per sample
data.frame(t(abundances(PS.mOTUS.filt.comp)))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "mOTU", values_to = "abundance")%>%
  left_join(prev.motus, by="mOTU")-> prev.abund.all.motus

##Check relationship between Taxa and ARG Richness
ord.data%>%
  rownames_to_column("Sample_ID")%>%
  left_join(abund.family, by="Sample_ID")%>%
  dplyr::filter(Family=="Enterobacteriaceae")%>%
  ggplot(aes(x= log2(abundance), y= observed_arg, color=arg_level))+
  geom_point(shape=21, alpha= 0.2, size=3, aes(fill= arg_level), color= "black")+
  labs(tag= "c", x= "Enterobacteriaceae abundance", y= "ARG richness")+
  scale_fill_manual(values = pal.level)+
  scale_color_manual(values = pal.level)+
  coord_cartesian(ylim = c(0, NA))+
  geom_smooth(method = lm, se = T)+
  geom_hline(yintercept = 100, linetype='dotted')+
  scale_y_continuous(limits = c(0, 210))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  ggpubr::stat_cor(label.y = c(200, 10), label.x = log2(0.0001), method = "spearman",
                   aes(label= paste("rho","'='", after_stat(r), after_stat(p.label), sep= "~` `~")))+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "none")-> Supp1Ca

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  left_join(abund.family, by="Sample_ID")%>%
  dplyr::filter(Family=="Enterococcaceae")%>%
  ggplot(aes(x= log2(abundance), y= observed_arg, color=arg_level))+
  geom_point(shape=21, alpha= 0.2, size=3, aes(fill= arg_level), color= "black")+
  labs(x= "Enterococcaceae abundance", y= "ARG richness")+
  scale_fill_manual(values = pal.level)+
  scale_color_manual(values = pal.level)+
  coord_cartesian(ylim = c(0, NA))+
  geom_smooth(method = lm, se = T)+
  geom_hline(yintercept = 100, linetype='dotted')+
  scale_y_continuous(limits = c(0, 210))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  ggpubr::stat_cor(label.y = c(200, 10), label.x = log2(0.0001), method = "spearman",
                   aes(label= paste("rho","'='", after_stat(r), after_stat(p.label), sep= "~` `~")))+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "none")-> Supp1Cb

##Now let's go to the general description of ARGs composition
##Determine ARG abundance in each single metagenome
##Overall prevalence
data.frame(t(prevalence(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)->prev.aro

data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(prev.aro, by="ARO")-> prev.abund.aro

##Get mean abundance by ARO
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(ARO) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")%>%
  left_join(prev.aro, by="ARO")-> prev.mabund.aro


##Try to determine which genes have a greater change from lower to high ARG level
###Differential abundance
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(ord.data%>%rownames_to_column("Sample_ID"), by="Sample_ID")%>%
  dplyr::group_by(ARO)%>%
  nest()%>%
  dplyr::mutate(model = map(data, ~ glm(log1p(abundance) ~ arg_level, data = .x, family = gaussian())),
                results = map(model, tidy))%>%
  select(ARO, results)%>%
  unnest(results)%>%
  dplyr::filter(term=="arg_levellow")%>%
  dplyr::mutate(padj = p.adjust(p.value, method = "BH"))%>%
  dplyr::mutate(adj.sig = case_when(padj<0.001  ~ "***",
                                    padj<0.01  ~ "**",
                                    padj<0.05  ~ "*",
                                    padj>0.05  ~ "NS"))-> diff.abund.arg

##Using wilcoxon
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(ord.data%>%rownames_to_column("Sample_ID"), by="Sample_ID")%>%
  dplyr::group_by(ARO)%>%
  nest()%>%
  dplyr::mutate(test = map(data, ~ wilcox.test(abundance ~ arg_level, data = .x)),
                esize = map(data, ~ rstatix::wilcox_effsize(abundance ~ arg_level, data = .x)),
                p.value = map_dbl(test, ~ .x$p.value),
                effsize= map_dbl(esize, ~ .x$effsize),
                magnitude= map_dbl(esize, ~ .x$magnitude))%>%
  dplyr::mutate(padj = p.adjust(p.value, method = "BH")) %>%
  select(ARO, p.value, padj, effsize, magnitude)%>%
  dplyr::mutate(adj.sig = case_when(padj<0.001  ~ "***",
                                    padj<0.01  ~ "**",
                                    padj<0.05  ~ "*",
                                    padj>0.05  ~ "NS"))%>%
  dplyr::rename(p.value.wc = p.value,
                padj.wc = padj,
                adj.sig.wc = adj.sig,
                effsize.wc =effsize,
                magnitude.wc = magnitude)-> diff.abund.wilcox.arg

##Generate a volcano plot
##Get annotation for genes
ARG.annot<-readRDS(file="data/raw/ARG_annotation_VHJD.rds")

##Merge diff abundance results
diff.abund.arg%>%
  left_join(diff.abund.wilcox.arg, by= "ARO")%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by= "ARO")-> diff.abund.arg.all

pal.abx<-c("tetracycline antibiotic"="#994F00", "glycopeptide antibiotic"="#006CD1",
           "aminoglycoside antibiotic"="#B672B8","fluoroquinolone antibiotic"="#FCB71B",
           "phosphonic acid antibiotic"="#24D582", "peptide antibiotic"="#F3CA95",
           "diaminopyrimidine antibiotic"="#748B31", "phenicol antibiotic"="#C5B51B",
           "mls antibiotics"="#89BDB6", "beta lactam antibiotics"="#E651BA",
           "multidrug"="#A2A2A2","other"="#e2dfe1")

##Differential abundance
require("ggrepel")
diff.abund.arg.all%>%
  dplyr::filter(padj>0)%>%
  dplyr::mutate(level = case_when(adj.sig=="NS" ~ 0,
                                  T ~ 1))%>%
  dplyr::mutate(label= case_when(level==1&!Drug_Class_adjusted %in% c("other", "multidrug")~ ARG_CARD_Short,
                                 T~ ""))%>%
  dplyr::mutate(label= case_when(Drug_Class_adjusted!="glycopeptide antibiotic"~ label,
                                 T~ ""))%>%
  ggplot(aes(x = estimate, y = -log10(padj), label= label)) +
  geom_point(aes(fill= Drug_Class_adjusted, size=prevalence*100, alpha=level),
             color="black", shape= 21,) +
  scale_fill_manual(values = pal.abx)+
  labs(x= expression(High~ARG~level~~phantom(0)~Effect~size~(beta~estimate)~~phantom(0)~Lower~ARG~level),
       y= expression(-log[10]~(Q)),
       fill= "Drug class",
       size= "Overall prevalence (%)",
       shape= "Status",
       tag = "a")+
  guides(fill=  guide_legend(override.aes=list(shape=21, size= 3)),
         alpha= "none")+
  geom_text_repel(max.overlaps = 25)+
  geom_hline(yintercept=-log10(0.1), col="black", linetype= "dashed") +
  geom_vline(xintercept=0, col="black", linetype= "dashed") +
  theme_minimal()+
  theme(text = element_text(size=16))-> Supp0a

##Composition plot by Mechanism, Drug class and genes by year
ord.data%>%
  rownames_to_column("Sample_ID")-> tmp3

pal.abx.2<-c("Tetracycline"="#994F00", "Glycopeptide"="#006CD1",
             "Aminoglycoside"="#B672B8","Fluoroquinolone"="#FCB71B",
             "Phosphonic acid"="#24D582", "Peptide"="#F3CA95",
             "Diaminopyrimidine"="#748B31", "Phenicol"="#C5B51B",
             "MLS"="#89BDB6", "Beta lactam"="#E651BA",
             "Multidrug"="#A2A2A2","Other"="#e2dfe1")

data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted%in%c("beta lactam antibiotics", "Beta lactam")~ "Beta lactam",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="mls antibiotics" ~ "MLS"  ,
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="phenicol antibiotic"~ "Phenicol",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Other"))%>%
  dplyr::mutate(Drug_Class_adjusted = fct_relevel(Drug_Class_adjusted,
                                                  "Aminoglycoside","Beta lactam",
                                                  "Diaminopyrimidine", "Fluoroquinolone",
                                                  "Glycopeptide", "MLS", "Peptide",
                                                  "Phenicol", "Phosphonic acid",
                                                  "Tetracycline", "Multidrug", "Other"))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(Year, Drug_Class_adjusted) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  ungroup()%>%
  group_by(Year) %>%
  mutate(rel_abundance = (abundance / sum(abundance))*100) %>%
  ungroup()->summarized_dc_year

##Get a relationship between specific ARGs and taxa
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted%in%c("beta lactam antibiotics", "Beta lactam")~ "Beta lactam",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="mls antibiotics" ~ "MLS"  ,
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="phenicol antibiotic"~ "Phenicol",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Other"))%>%
  dplyr::mutate(Drug_Class_adjusted = fct_relevel(Drug_Class_adjusted,
                                                  "Aminoglycoside","Beta lactam",
                                                  "Diaminopyrimidine", "Fluoroquinolone",
                                                  "Glycopeptide", "MLS", "Peptide",
                                                  "Phenicol", "Phosphonic acid",
                                                  "Tetracycline", "Multidrug", "Other"))%>%
  dplyr::filter(ARG_CARD_Short%in%c("CblA-1","tet(W)","tet(40)","tet(Q)",
                                    "tet(32)","tet(O)","tet(W/N/W)"))%>%
  dplyr::select(c(Sample_ID, abundance, ARG_CARD_Short, Drug_Class_adjusted,
                  prevalence, observed_arg, arg_level))%>%
  dplyr::rename(abundance_aro= abundance,
                prevalence_aro= prevalence)-> high.prev.abund.aro

prev.abund.all.motus%>%
  dplyr::filter(!Kingdom%in%c("unknown"))-> high.prev.abund.motus

##Get a relationship between Bacteroides and Cbl1 beta lactamase
high.prev.abund.aro%>%
  dplyr::select(c(Sample_ID, abundance_aro, ARG_CARD_Short, observed_arg, arg_level))%>%
  pivot_wider(id_cols= c(Sample_ID, observed_arg, arg_level),
              names_from = "ARG_CARD_Short",
              values_from = "abundance_aro")-> tmp

high.prev.abund.motus%>%
  dplyr::filter(Species%in%c("Escherichia coli",
                             "Enterococcus villorum",
                             "Enterococcus ratti",
                             "Enterococcus faecalis",
                             "Enterococcus hirae",
                             "Enterococcus faecium",
                             "Bacteroides caecimuris",
                             "Bacteroides acidifaciens",
                             "Bacteroides dorei/vulgatus"))%>%
  dplyr::select(c(Sample_ID, abundance, Species))%>%
  pivot_wider(id_cols= Sample_ID,
              names_from = "Species",
              values_from = "abundance")%>%
  left_join(tmp, by= "Sample_ID")-> aro.motus

##Plot now
aro.motus%>%
  ggplot(aes(x= log2(`Escherichia coli`), y= observed_arg))+
  geom_point(shape=21, alpha= 0.2, size=3, aes(fill= arg_level), color= "black")+
  labs(x= "E. coli abundance", y= "ARG richness")+
  scale_fill_manual(values = pal.level)+
  coord_cartesian(ylim = c(0, NA))+
  geom_smooth(method = lm, se = T, color= "black")+
  geom_hline(yintercept = 100, linetype='dotted')+
  scale_y_continuous(limits = c(0, 210))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  ggpubr::stat_cor(label.y = 200, label.x = log2(0.0001), method = "spearman",
                   aes(label= paste("rho","'='", after_stat(r), after_stat(p.label), sep= "~` `~")))+
  geom_rug(col=rgb(.5,0,0,alpha=.2), sides="l")+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "none")-> Supp1Cc

###Bacterides
##Plot now
aro.motus%>%
  dplyr::select(c(Sample_ID, `Bacteroides acidifaciens`,`Bacteroides caecimuris`,
                  `Bacteroides dorei/vulgatus`, `CblA-1`))%>%
  pivot_longer(!c(Sample_ID,`CblA-1`),
               names_to = "Species",
               values_to = "abundance")%>%
  ggplot(aes(x= log2(abundance), y= log2(`CblA-1`), color=Species))+
  geom_point(shape=21, alpha= 0.2, size=3, aes(fill= Species), color= "black")+
  labs(tag= "d", x= "mOTU abundance", y= "CblA-1 abundance")+
  scale_fill_manual(values = c("Bacteroides acidifaciens"="#88CCEE",
                               "Bacteroides caecimuris"="#882255",
                               "Bacteroides dorei/vulgatus"= "#117733"))+
  scale_color_manual(values = c("Bacteroides acidifaciens"="#88CCEE",
                                "Bacteroides caecimuris"="#882255",
                                "Bacteroides dorei/vulgatus"= "#117733"))+
  coord_cartesian(ylim = c(0, NA))+
  geom_smooth(method = lm, se = T)+
  scale_y_continuous(limits = c(0, 30))+
  guides(color= "none", size= "none",
         fill= guide_legend(override.aes=list(shape=c(21), size= 3)))+
  ggpubr::stat_cor(label.y = log2(20000), label.x = log2(0.0001), method = "spearman",
                   aes(label= paste("rho","'='", after_stat(r), after_stat(p.label), sep= "~` `~")))+
  facet_grid(~Species)+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "none")-> Supp1D

### Get a general prevalence-abundance plot
prev.mabund.aro%>%
  left_join(ARG.annot, by= "ARO")%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted%in%c("beta lactam antibiotics", "Beta lactam")~ "Beta lactam",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="mls antibiotics" ~ "MLS"  ,
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="phenicol antibiotic"~ "Phenicol",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Other"))%>%
  dplyr::mutate(Drug_Class_adjusted = fct_relevel(Drug_Class_adjusted,
                                                  "Aminoglycoside","Beta lactam",
                                                  "Diaminopyrimidine", "Fluoroquinolone",
                                                  "Glycopeptide", "MLS", "Peptide",
                                                  "Phenicol", "Phosphonic acid",
                                                  "Tetracycline", "Multidrug", "Other"))%>%
  dplyr::mutate(level = case_when((prevalence>=0.50)&!Drug_Class_adjusted%in%c("Multidrug", "Other") ~ 1,
                                  T ~ 0))%>%
  dplyr::mutate(label= case_when(level==1~ ARG_CARD_Short,
                                 T~ ""))%>%
  ggplot(aes(x= prevalence*100, y= log2(mean_abundance+1), label= label,
             fill= Drug_Class_adjusted))+
  annotate(geom = "segment", x = 50, xend = 50, y = -Inf, yend = Inf, linetype = "dashed")+
  geom_point(shape=21, color="black", size= 3)+
  scale_fill_manual(values = pal.abx.2)+
  scale_x_continuous("Overall prevalence (% mice with ARG)", seq(0, 100, by = 30),
                     limits = c(0, 100)) +
  ylab("Mean ARG abundance \n(log2(FPKM))")+
  labs(tag= "a", fill= "Drug class")+
  geom_text_repel(max.overlaps = 100)+
  guides(fill = guide_legend(nrow = 4, ncol = 3), color= "none")+
  theme_minimal()+
  theme(text = element_text(size=16),
        legend.position = "bottom",
        panel.border = element_blank())-> Supp2A

##Resistance mechanism
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Resistance_Mechanism= if_else(Resistance_Mechanism=="antibiotic target alteration",
                                              "Target alteration", "Other mechanisms"))%>%
  dplyr::mutate(Resistance_Mechanism = fct_relevel(Resistance_Mechanism,
                                                   "Other mechanisms", "Target alteration"))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(Resistance_Mechanism) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(rel_abundance = (abundance / sum(abundance))*100) %>%
  ungroup()%>%
  ggplot(aes(x= 1, y= rel_abundance, fill= Resistance_Mechanism)) +
  scale_fill_manual(name = "Resistance_Mechanism", values = c("Target alteration"= "#9EA8DC",
                                                              "Other mechanisms"="#e2dfe1")) +
  geom_bar(aes(), stat="identity", position="stack") +
  theme_minimal()+
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
        text = element_text(size=16), axis.line.y = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        strip.background = element_blank(), legend.position = "none",
        panel.border = element_blank(), panel.spacing = unit(0, "lines")) +
  labs(y= "Relative abundance \n(FPKM)", tag = "a")+
  guides(fill= guide_legend(nrow = 13, ncol = 1))-> F2A

##Drug class general
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Drug_Class_adjusted= if_else(Drug_Class_adjusted=="multidrug", "Multidrug", "Other"))%>%
  dplyr::mutate(Drug_Class_adjusted = fct_relevel(Drug_Class_adjusted,
                                                  "Other", "Multidrug"))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(Drug_Class_adjusted) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(rel_abundance = (abundance / sum(abundance))*100) %>%
  ungroup()%>%
  ggplot(aes(x= 1, y= rel_abundance, fill= Drug_Class_adjusted)) +
  scale_fill_manual(name = "Drug class", values = c("Multidrug"="#A2A2A2","Other"="#e2dfe1")) +
  geom_bar(aes(), stat="identity", position="stack") +
  theme_minimal()+
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
        text = element_text(size=16), axis.line.y = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        strip.background = element_blank(), legend.position = "none",
        axis.text.y =  element_blank(),  axis.title.y =  element_blank(),
        panel.border = element_blank(), panel.spacing = unit(0, "lines")) +
  labs(y= "Relative abundance \n(Gene Coverage Per Million)")-> F2B

##Now based on resistance mechanism
pal.resm.2<- c("Efflux"= "#df8ac7", "Target replacement"= "#1b8b93",
               "Inactivation"= "#f8eb40", "Target alteration"= "#9EA8DC",
               "Target protection"= "#003a00","Reduced permeability"="#f29d51",
               "Multiple mechanisms"= "#e2dfe1")

data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Resistance_Mechanism= case_when(
    Resistance_Mechanism=="antibiotic inactivation"~ "Inactivation",
    Resistance_Mechanism=="antibiotic target protection"~ "Target protection",
    Resistance_Mechanism=="antibiotic target alteration"~ "Target alteration",
    Resistance_Mechanism=="antibiotic target replacement"~ "Target replacement",
    Resistance_Mechanism=="reduced permeability to antibiotic"~ "Reduced permeability",
    Resistance_Mechanism=="antibiotic efflux" ~ "Efflux"  ,
    Resistance_Mechanism=="multiple mechanisms"~ "Multiple mechanisms"))%>%
  dplyr::mutate(Resistance_Mechanism = fct_relevel(Resistance_Mechanism,
                                                   "Efflux","Inactivation",
                                                   "Target alteration",
                                                   "Target protection",
                                                   "Target replacement",
                                                   "Reduced permeability",
                                                   "Multiple mechanisms"))%>%
  dplyr::filter(!is.na(abundance))%>%
  dplyr::filter(!Drug_Class_adjusted%in%c("multidrug", "other"))%>%
  group_by(Year, Resistance_Mechanism) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  ungroup()%>%
  group_by(Year) %>%
  dplyr::mutate(rel_abundance = abundance / sum(abundance)) %>%
  ungroup()->summarized_rm_year

##Get per mouse plot
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Resistance_Mechanism= case_when(
    Resistance_Mechanism=="antibiotic inactivation"~ "Inactivation",
    Resistance_Mechanism=="antibiotic target protection"~ "Target protection",
    Resistance_Mechanism=="antibiotic target alteration"~ "Target alteration",
    Resistance_Mechanism=="antibiotic target replacement"~ "Target replacement",
    Resistance_Mechanism=="reduced permeability to antibiotic"~ "Reduced permeability",
    Resistance_Mechanism=="antibiotic efflux" ~ "Efflux"  ,
    Resistance_Mechanism=="multiple mechanisms"~ "Multiple mechanisms"))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(Sample_ID) %>%
  dplyr::mutate(rel_abundance = (abundance / sum(abundance))*100) %>%
  dplyr::mutate(Resistance_Mechanism = fct_relevel(Resistance_Mechanism,
                                                   "Target protection",
                                                   "Efflux",
                                                   "Target alteration",
                                                   "Inactivation",
                                                   "Reduced permeability",
                                                   "Target replacement",
                                                   "Multiple mechanisms"))%>%
  ggplot(aes(x= Sample_ID, y= rel_abundance, fill= Resistance_Mechanism)) +
  scale_fill_manual(name = "Resistance_Mechanism", values = pal.resm.2) +
  geom_bar(aes(), stat="identity", position="stack") +
  theme_minimal()+
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
        text = element_text(size=16), axis.line.y = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        strip.background = element_blank(), legend.position = "bottom",
        panel.border = element_blank(), panel.spacing = unit(0, "lines")) +
  labs(y= "Relative abundance \n(FPKM)", tag = "b")+
  guides(fill= guide_legend(nrow = 4, ncol = 3))-> Supp2B

##Alluvial plot
require(ggalluvial)
data.frame(t(abundances(transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x)))))%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "abundance")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  left_join(ARG.annot, by= "ARO")%>%
  left_join(prev.aro, by="ARO")%>%
  left_join(tmp3, by="Sample_ID")%>%
  dplyr::mutate(Drug_Class_adjusted= case_when(
    Drug_Class_adjusted=="aminoglycoside antibiotic"~ "Aminoglycoside",
    Drug_Class_adjusted%in%c("beta lactam antibiotics", "Beta lactam")~ "Beta lactam",
    Drug_Class_adjusted=="diaminopyrimidine antibiotic"~ "Diaminopyrimidine",
    Drug_Class_adjusted=="fluoroquinolone antibiotic"~ "Fluoroquinolone",
    Drug_Class_adjusted=="glycopeptide antibiotic"~ "Glycopeptide",
    Drug_Class_adjusted=="mls antibiotics" ~ "MLS"  ,
    Drug_Class_adjusted=="peptide antibiotic"~ "Peptide",
    Drug_Class_adjusted=="phenicol antibiotic"~ "Phenicol",
    Drug_Class_adjusted=="phosphonic acid antibiotic"~ "Phosphonic acid",
    Drug_Class_adjusted=="tetracycline antibiotic"~ "Tetracycline",
    Drug_Class_adjusted=="multidrug"~ "Multidrug",
    T~ "Other"))%>%
  dplyr::mutate(Drug_Class_adjusted = fct_relevel(Drug_Class_adjusted,
                                                  "Aminoglycoside","Beta lactam",
                                                  "Diaminopyrimidine", "Fluoroquinolone",
                                                  "Glycopeptide", "MLS", "Peptide",
                                                  "Phenicol", "Phosphonic acid",
                                                  "Tetracycline", "Multidrug", "Other"))%>%
  dplyr::mutate(Resistance_Mechanism= case_when(
    Resistance_Mechanism=="antibiotic inactivation"~ "Inactivation",
    Resistance_Mechanism=="antibiotic target protection"~ "Target protection",
    Resistance_Mechanism=="antibiotic target alteration"~ "Target alteration",
    Resistance_Mechanism=="antibiotic target replacement"~ "Target replacement",
    Resistance_Mechanism=="reduced permeability to antibiotic"~ "Reduced permeability",
    Resistance_Mechanism=="antibiotic efflux" ~ "Efflux"  ,
    Resistance_Mechanism=="multiple mechanisms"~ "Multiple mechanisms"))%>%
  dplyr::mutate(Resistance_Mechanism = fct_relevel(Resistance_Mechanism,
                                                   "Multiple mechanisms",
                                                   "Efflux",
                                                   "Reduced permeability","Inactivation",
                                                   "Target alteration",
                                                   "Target replacement",
                                                   "Target protection"))%>%
  dplyr::filter(!is.na(abundance))%>%
  group_by(Year, Drug_Class_adjusted, Resistance_Mechanism) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")%>%
  ungroup()%>%
  group_by(Year) %>%
  mutate(rel_abundance = (abundance / sum(abundance))*100) %>%
  ungroup()%>%
  dplyr::filter(!Drug_Class_adjusted%in%c("Multidrug", "Other"))%>%
  ggplot(aes(axis1 = Year, axis2 = Drug_Class_adjusted, axis3 = Resistance_Mechanism, y = rel_abundance)) +
  geom_alluvium(aes(fill = Drug_Class_adjusted), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/50, fill = "grey95", color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Year", "Drug class", "Resistance mechanism"), expand = c(.05, .05)) +
  scale_fill_manual(values = pal.abx.2) +
  labs(y = "Relative ARG abundance", x = NULL, fill = "Drug class") +
  guides(fill= guide_legend(override.aes=list(shape=c(22), size= 3)))+
  theme_minimal() +
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")-> F1C

##Get the PCA for ARGs and color by year
##Get dataset at same level
PS.ARG.scaled.HM<- transform_sample_counts(PS.ARG.FPKM.filt, function(x) 1E9 * x/sum(x))

#Compositional (Aitchison distance)
aitch_dist<- vegan::vegdist(t(PS.ARG.scaled.HM@otu_table),
                            method="aitchison", pseudocount=1)
##For PCA with Aitchitson
ordination.clr <- phyloseq::ordinate(PS.ARG.scaled.HM, "RDA") #principal components analysis

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(ord.data, Seq_batch)

arg.adonis.ait<- vegan::adonis2(aitch_dist~ Year + Sex + Longitude + Latitude + transect +
                                  arg_level + Tax_mapped,
                                data = ord.data, na.action = na.fail,
                                permutations = perm, by="margin")

#PCA Aitchitson Select Axis 1 and 2
seg.data<-data.frame(ordination.clr$CA$u[,1:2])
seg.data%>%
  rownames_to_column("Sample_ID")-> seg.data

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::left_join(seg.data, by="Sample_ID")->seg.data

##Just to have an overview
ggplot() +
  geom_point(data=seg.data, aes(x=PC1,y=PC2, fill= Year), shape= 21, size=2) +
  stat_ellipse(data=seg.data, aes(x=PC1,y=PC2, color= Year))+
  guides(fill = guide_legend(override.aes=list(shape=c(21), size= 3)), color= "none")+
  labs(tag= "c", fill  = "Year")+
  theme_minimal()+
  scale_fill_manual(values = pal.year)+
  scale_color_manual(values = pal.year)+
  theme(text = element_text(size=16))+
  annotate("text", x = 0.15, y = 0.100, label= "Permanova (Year level)", size= 3)+
  annotate("text", x = 0.15, y = 0.090, label= "Aitchison distance", size= 3)+
  annotate("text", x = 0.15, y = 0.080, label= paste0(label = "R²= ", round(arg.adonis.ait$R2[1], digits = 3),
                                                      ", p = ", arg.adonis.ait$`Pr(>F)`[1]), color = "black", size= 3)+
  xlab(paste0("PC 1 [", round(ordination.clr$CA$eig[1] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(ordination.clr$CA$eig[2] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))-> C

##Craft the Final Figure 2
BC<-ggarrange(B, C,  ncol = 2, nrow = 1)
A<- ggarrange(F2A, F2B, F1C, ncol = 3, nrow = 1)
Fig2<-ggarrange(A, BC,  ncol = 1, nrow = 2)


##Craft the Final Supplement 1
Supp1AB<-ggarrange(Supp1A, Supp1B,  ncol = 2, nrow = 1)
Supp1C<-ggarrange(Supp1Ca, Supp1Cb, Supp1Cc, ncol = 3, nrow = 1)
Supp1<-ggarrange(Supp1AB, Supp1C, Supp1D, ncol = 1, nrow = 3)
