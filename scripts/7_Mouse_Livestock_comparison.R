##Project: Wild HMHZ project
##Aim: House mice ARG comparison to farm resistomes
##Author: Víctor Hugo Jarquín-Díaz
##Date: 02.09.2025 (Update filtering)

##Libraries
library(tidyverse)
library(phyloseq)
library(vroom)
library(microbiome)

##Load house mice data that were positive for ARGs
PS.ARG.FPKM.filt <- readRDS("data/raw/PS.ARG.FPKM.filtered.rds")
sample_data_house<-readRDS(file="data/raw/Mus_sample_data_farm_VHJD.rds")

##Load farm data
PS.ARG.FPKM.farm<- readRDS(file = "data/raw/PS.ARG.FPKM.farm.rds")

##Pallet
pal.host<- c("mouse"= "#3E5496", "pig"="#a06177",
             "chicken"="#d9af6b", "cattle"="#8c785d", "turkey"="#68855c")

pal.country<- c("Germany"="#ffcc00",
                "Denmark"="#356c4e",
                "Poland"="#dc143c",
                "Netherlands"="#e8640c",
                "Belgium"="#a47d6f",
                "France"="#24185d",
                "Italy"="#0f5224",
                "Spain"="#ffd691",
                "Bulgaria"="#009663",
                "Ghana"="#000000")

##Merge farm and house mice
##Get a uniform sample data for both
##We need: Host, Country, Year, arg alpha diversity
sample_data_house%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::select(c(Sample_ID,
                  Mouse_ID,
                  Host,
                  Year,
                  observed_arg:dominance_simpson_arg,
                  transect))%>%
  dplyr::rename(collection_date = Year,
                sample_title= Mouse_ID,
                geo_loc_name = transect,
                env_biome = Host)%>%
  dplyr::mutate(geo_loc_name = as.factor("Germany"),
                env_biome = as.factor("mouse"))%>%
  dplyr::mutate(bioproject = "PRJEBXXXXX")%>%
  relocate(sample_title)%>%
  relocate(collection_date, .after = sample_title)%>%
  relocate(geo_loc_name, .after = collection_date)%>%
  relocate(env_biome, .after = geo_loc_name)%>%
  column_to_rownames("Sample_ID")-> data_merge_HM

## Alpha diversity farm
##Get alpha diversity from ARGs
alpha.arg <- microbiome::alpha(PS.ARG.FPKM.farm, index =c("Observed","Chao1", "Shannon", "Simpson"))

alpha.arg%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::rename("observed_arg" = "observed",
                "chao1_arg" = "chao1",
                "diversity_shannon_arg"= "diversity_shannon",
                "evenness_simpson_arg"= "evenness_simpson",
                "dominance_simpson_arg"= "dominance_simpson")-> alpha.arg

data.frame(PS.ARG.FPKM.farm@sam_data)%>%
  rownames_to_column("Sample_ID")%>%
  left_join(alpha.arg, by= "Sample_ID")%>%
  dplyr::mutate(across(
    c(observed_arg:dominance_simpson_arg), ~ replace_na(.x, 0)))%>%
  dplyr::select(c(Sample_ID,
                  sample_title,
                  env_biome,
                  collection_date,
                  observed_arg:dominance_simpson_arg,
                  geo_loc_name))%>%
  relocate(sample_title)%>%
  relocate(collection_date, .after = sample_title)%>%
  relocate(geo_loc_name, .after = collection_date)%>%
  relocate(env_biome, .after = geo_loc_name)%>%
  dplyr::mutate(bioproject = case_when(
    geo_loc_name!="Ghana"&env_biome%in%c("pig","poultry")~ "PRJEB22062",
    geo_loc_name!="Ghana"&env_biome%in%c("cattle","turkey")~ "PRJEB39685",
    T~ "PRJEB62878"),
    env_biome= case_when(env_biome=="poultry"~ "chicken",
                         T~ env_biome))%>%
  column_to_rownames("Sample_ID")-> data_merge_farm

##Now change the sample data for each dataset
PS.ARG.FPKM.filt@sam_data<- sample_data(data_merge_HM)
PS.ARG.FPKM.farm@sam_data<- sample_data(data_merge_farm)

##Filter samples with zero counts from farm
PS.ARG.FPKM.farm<- phyloseq::prune_samples(sample_sums(PS.ARG.FPKM.farm)>1, PS.ARG.FPKM.farm)

##Merge into a single dataset for comparison
PS.ARG.FPKM.all<- merge_phyloseq(PS.ARG.FPKM.filt, PS.ARG.FPKM.farm)

##Get data for ordination
PS.ARG.FPKM.all@sam_data-> ord.data
data.frame(ord.data)-> ord.data
##Beta diversity
#Compositional (Aitchison distance)
aitch_dist<- vegan::vegdist(t(PS.ARG.FPKM.all@otu_table),
                            method="aitchison", pseudocount=1)
##For PCA with Aitchitson
PS.ARG.FPKM.all.clr <- microbiome::transform(PS.ARG.FPKM.all, "clr") #Centered log ratio transformation
ordination.clr <- phyloseq::ordinate(PS.ARG.FPKM.all.clr, "RDA") #principal components analysis

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(ord.data, bioproject)

arg.adonis.ait<- vegan::adonis2(aitch_dist~ env_biome + geo_loc_name,
                                data = ord.data, na.action = na.fail,
                                permutations = perm, by="margin")

#PCA Aitchitson Select Axis 1, 2 and 3
seg.data<-data.frame(ordination.clr$CA$u[,1:3])
seg.data%>%
  rownames_to_column("Sample_ID")-> seg.data

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::left_join(seg.data, by="Sample_ID")->seg.data

##PCA plot
ggplot() +
  geom_point(data=seg.data, aes(x=PC1,y=PC2, fill= env_biome), shape= 21, size=2) +
  guides(fill = guide_legend(override.aes=list(shape=c(21), size= 3)), color= "none")+
  labs(tag= "c", fill  = "Host")+
  theme_bw()+
  scale_fill_manual(values = pal.host)+
  scale_color_manual(values = pal.host)+
  theme(text = element_text(size=16))+
  annotate(geom = "segment", x = 0, xend = 0, y = -Inf, yend = Inf, linetype = "dashed")+
  annotate(geom = "segment", x = -Inf, xend = Inf, y = 0, yend = 0, linetype = "dashed")+
  annotate("text", x = -0.01, y = 0.050, label= "Permanova (Host)", size= 4)+
  annotate("text", x = -0.01, y = 0.045, label= "Aitchison distance", size= 4)+
  annotate("text", x = -0.01, y = 0.040, label= paste0(label = "R²= ", round(arg.adonis.ait$R2[1], digits = 3),
                                                       ", p = ", arg.adonis.ait$`Pr(>F)`[1]), color = "black", size= 4)+
  xlab(paste0("PC 1 [", round(ordination.clr$CA$eig[1] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(ordination.clr$CA$eig[2] / sum(ordination.clr$CA$eig)*100, digits = 2), "%]"))-> F3C

##Scale and check with other distances
PS.ARG.FPKM.comp<- transform_sample_counts(PS.ARG.FPKM.all, function(x) 1E9 * x/sum(x))
PS.ARG.FPKM.comp<-microbiome::transform(PS.ARG.FPKM.comp, "compositional") ##Total sum transformation

#1) Presence/absence (Jaccard index)
jaccard_dist <- phyloseq::distance(PS.ARG.FPKM.comp, method="jaccard",
                                   type="samples", binary=T)

arg.adonis.jac<- vegan::adonis2(jaccard_dist~ env_biome + geo_loc_name,
                                data = ord.data, na.action = na.fail,
                                permutations = perm, by="margin")

#2) Abundance (Bray-Curtis index)
bray_dist<- phyloseq::distance(PS.ARG.FPKM.comp,
                               method="bray", weighted=T)

ordination<- ordinate(PS.ARG.FPKM.comp,
                      method="PCoA", distance="bray")

arg.adonis.bc<- vegan::adonis2(bray_dist~ env_biome + geo_loc_name,
                               data = ord.data, na.action = na.fail,
                               permutations = perm, by="margin")

##PCoA
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, ord.data$env_biome, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2
seg.data<-cbind(vectors[,1:4],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:4])
names(seg.data)<-c("env_biome","v.PCoA1","v.PCoA2","v.PCoA3","PCoA1","PCoA2","PCoA3")

##Add sample data
ord.data%>%
  dplyr::select(!c(env_biome))%>%
  cbind(seg.data)-> seg.data

##Just to have an overview
ggplot() +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= env_biome), shape=21, size=2) +
  scale_fill_manual(values =pal.host, limits = c("mouse", "pig", "chicken",
                                                 "cattle", "turkey"))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(tag= "a", fill  = "Host")+
  theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom")+
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps , group=grps ),size=4, shape= 4) +
  scale_color_manual(values = pal.host)+
  annotate(geom = "segment", x = 0, xend = 0, y = -Inf, yend = Inf, linetype = "dashed")+
  annotate(geom = "segment", x = -Inf, xend = Inf, y = 0, yend = 0, linetype = "dashed")+
  annotate("text", x = -0.2, y = -0.49, label= "Permanova (Host)", size= 4)+
  annotate("text", x = -0.2, y = -0.53, label= "Bray-Curtis distance", size= 4)+
  annotate("text", x = -0.2, y = -0.57, label= paste0(label = "R²= ", round(arg.adonis.bc$R2[1], digits = 3),
                                                      ", p = ", arg.adonis.bc$`Pr(>F)`[1]), color = "black", size= 4)+

  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> SuppA

ggplot() +
  geom_point(data=seg.data, aes(x=v.PCoA2,y=v.PCoA3, fill= env_biome), shape=21, size=2) +
  scale_fill_manual(values =pal.host, limits = c("mouse", "pig", "chicken",
                                                 "cattle", "turkey"))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(tag= "b", fill  = "Host")+
  theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom")+
  geom_point(data=centroids[,1:5], aes(x=PCoA2,y=PCoA3, color= grps , group=grps ),size=4, shape= 4) +
  scale_color_manual(values = pal.host)+
  annotate(geom = "segment", x = 0, xend = 0, y = -Inf, yend = Inf, linetype = "dashed")+
  annotate(geom = "segment", x = -Inf, xend = Inf, y = 0, yend = 0, linetype = "dashed")+
  annotate("text", x = -0.2, y = -0.49, label= "Permanova (Host)", size= 4)+
  annotate("text", x = -0.2, y = -0.53, label= "Bray-Curtis distance", size= 4)+
  annotate("text", x = -0.2, y = -0.57, label= paste0(label = "R²= ", round(arg.adonis.bc$R2[1], digits = 3),
                                                      ", p = ", arg.adonis.bc$`Pr(>F)`[1]), color = "black", size= 4)+

  xlab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 3 [", round(ordination$values[3,2]*100, digits = 2), "%]"))-> SuppB

##Let's check only presence absence distance
mvd<- vegan::betadisper(jaccard_dist, ord.data$env_biome, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2
seg.data<-cbind(vectors[,1:4],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:4])
names(seg.data)<-c("env_biome","v.PCoA1","v.PCoA2","v.PCoA3","PCoA1","PCoA2","PCoA3")

##Add sample data
ord.data%>%
  dplyr::select(!c(env_biome))%>%
  cbind(seg.data)-> seg.data

ordination<- ordinate(PS.ARG.FPKM.comp,
                      method="PCoA", distance="jaccard", binary = TRUE)

##Just to have an overview
ggplot() +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= env_biome), shape=21, size=2) +
  scale_fill_manual(values =pal.host, limits = c("mouse", "pig", "chicken",
                                                 "cattle", "turkey"))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(tag= "c", fill  = "Host")+
  theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom")+
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps , group=grps ),size=4, shape= 4) +
  scale_color_manual(values = pal.host)+
  annotate(geom = "segment", x = 0, xend = 0, y = -Inf, yend = Inf, linetype = "dashed")+
  annotate(geom = "segment", x = -Inf, xend = Inf, y = 0, yend = 0, linetype = "dashed")+
  annotate("text", x = -0.2, y = -0.49, label= "Permanova (Host)", size= 4)+
  annotate("text", x = -0.2, y = -0.53, label= "Jaccard distance", size= 4)+
  annotate("text", x = -0.2, y = -0.57, label= paste0(label = "R²= ", round(arg.adonis.jac$R2[1], digits = 3),
                                                      ", p = ", arg.adonis.jac$`Pr(>F)`[1]), color = "black", size= 4)+

  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> SuppC

ggplot() +
  geom_point(data=seg.data, aes(x=v.PCoA2,y=v.PCoA3, fill= env_biome), shape=21, size=2) +
  scale_fill_manual(values =pal.host, limits = c("mouse", "pig", "chicken",
                                                 "cattle", "turkey"))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(tag= "c", fill  = "Host")+
  theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom")+
  geom_point(data=centroids[,1:5], aes(x=PCoA2,y=PCoA3, color= grps , group=grps ),size=4, shape= 4) +
  scale_color_manual(values = pal.host)+
  annotate(geom = "segment", x = 0, xend = 0, y = -Inf, yend = Inf, linetype = "dashed")+
  annotate(geom = "segment", x = -Inf, xend = Inf, y = 0, yend = 0, linetype = "dashed")+
  annotate("text", x = -0.2, y = -0.49, label= "Permanova (Host)", size= 4)+
  annotate("text", x = -0.2, y = -0.53, label= "Jaccard distance", size= 4)+
  annotate("text", x = -0.2, y = -0.57, label= paste0(label = "R²= ", round(arg.adonis.jac$R2[1], digits = 3),
                                                      ", p = ", arg.adonis.jac$`Pr(>F)`[1]), color = "black", size= 4)+

  xlab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 3 [", round(ordination$values[3,2]*100, digits = 2), "%]"))-> SuppD

##
##ARG matrices
bc.mat<- as.matrix(bray_dist)
jac.mat<- as.matrix(jaccard_dist)
ait.mat<- as.matrix(aitch_dist)

##Combine distance matrices
comb_dist<- function(feature_matix_1, feature_matix_2) {

  pairs <- expand.grid(colnames(feature_matix_1), colnames(feature_matix_1))
  pairs <- data.frame(id_1 = pairs[,1], id_2 = pairs[,2])
  pairs <- pairs[order(match(pairs$id_1, colnames(feature_matix_1))),] # Keep the original order

  combi.dist<- data.frame(pairs,
                          type_matrix_A= as.vector(feature_matix_1),
                          type_matrix_B= as.vector(feature_matix_2))

  ##Just take the lower part of the combinations
  rowCol <- expand.grid(rownames(feature_matix_1), colnames(feature_matix_1))
  labs <- rowCol[as.vector(lower.tri(feature_matix_1, diag = T)),]
  labs%>%
    dplyr::mutate(sample_pair= paste0(Var1, "_", Var2))%>%
    dplyr::select(c(sample_pair))-> labs

  labs$sample_pair->labs

  # Filter combi.dist based on sample pairs
  combi.dist <- combi.dist %>%
    mutate(sample_pair = paste0(id_1, "_", id_2)) %>%
    filter(sample_pair %in% labs)%>%
    select(sample_pair, type_matrix_A, type_matrix_B)

  return(combi.dist)
}

# Generate all pairwise combinations of column names while keeping the order
##Bray-Jacc
distances.arg.df<- comb_dist(bc.mat, jac.mat)
##Jacc-Gi
tmp<- comb_dist(bc.mat, ait.mat)
tmp%>%
  dplyr::select(!type_matrix_A)%>%
  dplyr::rename("Ait_arg_all"= type_matrix_B)-> tmp1

distances.arg.df%>%
  dplyr::rename("Bray_arg_all"= type_matrix_A,
                "Jac_arg_all"= type_matrix_B)%>%
  dplyr::left_join(tmp1, by= "sample_pair")-> distances.arg.df

##Normalize to 0 to 1
require(scales)
distances.arg.df%>%
  dplyr::mutate(across(all_of("Ait_arg_all"), ~ rescale(.x, to = c(0, 1))))-> distances.arg.df

##Add metadata required for the comparisons
##Get vectors of samples (Needs to improve this would be hell for more localities or categories)
##Host
ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(env_biome=="mouse")%>%
  dplyr::pull(Sample_ID)-> mouse

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(env_biome=="pig")%>%
  dplyr::pull(Sample_ID)-> pig

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(env_biome=="chicken")%>%
  dplyr::pull(Sample_ID)-> chicken

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(env_biome=="cattle")%>%
  dplyr::pull(Sample_ID)-> cattle

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(env_biome=="turkey")%>%
  dplyr::pull(Sample_ID)-> turkey

##Country
ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Germany")%>%
  dplyr::pull(Sample_ID)-> Germany

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Belgium")%>%
  dplyr::pull(Sample_ID)-> Belgium

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Bulgaria")%>%
  dplyr::pull(Sample_ID)-> Bulgaria

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Denmark")%>%
  dplyr::pull(Sample_ID)-> Denmark

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="France")%>%
  dplyr::pull(Sample_ID)-> France

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Ghana")%>%
  dplyr::pull(Sample_ID)-> Ghana

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Italy")%>%
  dplyr::pull(Sample_ID)-> Italy

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Netherlands")%>%
  dplyr::pull(Sample_ID)-> Netherlands

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Poland")%>%
  dplyr::pull(Sample_ID)-> Poland

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(geo_loc_name=="Spain")%>%
  dplyr::pull(Sample_ID)-> Spain

distances.arg.df%>%
  tidyr::separate(sample_pair, c("Sample_A", "Sample_B"), sep= "_",remove = F)%>%
  dplyr::mutate(env_biome_A = case_when(Sample_A%in%mouse  ~ "mouse",
                                        Sample_A%in%pig  ~ "pig",
                                        Sample_A%in%chicken  ~ "chicken",
                                        Sample_A%in%cattle  ~ "cattle",
                                        Sample_A%in%turkey  ~ "turkey"))%>%
  dplyr::mutate(env_biome_B = case_when(Sample_B%in%mouse  ~ "mouse",
                                        Sample_B%in%pig  ~ "pig",
                                        Sample_B%in%chicken  ~ "chicken",
                                        Sample_B%in%cattle  ~ "cattle",
                                        Sample_B%in%turkey  ~ "turkey"))%>%
  dplyr::mutate(same_env_biome = case_when(env_biome_A == env_biome_B  ~ 1,
                                           env_biome_A != env_biome_B ~ 0))%>%
  dplyr::mutate(geo_loc_name_A = case_when(Sample_A%in%Germany  ~ "Germany",
                                           Sample_A%in%Ghana  ~ "Ghana",
                                           Sample_A%in%Belgium  ~ "Belgium",
                                           Sample_A%in%Bulgaria  ~ "Bulgaria",
                                           Sample_A%in%Denmark  ~ "Denmark",
                                           Sample_A%in%France  ~ "France",
                                           Sample_A%in%Italy  ~ "Italy",
                                           Sample_A%in%Poland  ~ "Poland",
                                           Sample_A%in%Spain  ~ "Spain",
                                           Sample_A%in%Netherlands  ~ "Netherlands"))%>%
  dplyr::mutate(geo_loc_name_B = case_when(Sample_B%in%Germany  ~ "Germany",
                                           Sample_B%in%Ghana  ~ "Ghana",
                                           Sample_B%in%Belgium  ~ "Belgium",
                                           Sample_B%in%Bulgaria  ~ "Bulgaria",
                                           Sample_B%in%Denmark  ~ "Denmark",
                                           Sample_B%in%France  ~ "France",
                                           Sample_B%in%Italy  ~ "Italy",
                                           Sample_B%in%Poland  ~ "Poland",
                                           Sample_B%in%Spain  ~ "Spain",
                                           Sample_B%in%Netherlands  ~ "Netherlands"))%>%
  dplyr::mutate(Same_geo_loc_name = case_when(geo_loc_name_A == geo_loc_name_B  ~ 1,
                                              geo_loc_name_A != geo_loc_name_B ~ 0))%>%
  dplyr::mutate(Same_sample = case_when(Sample_A == Sample_B  ~ 1,
                                        Sample_A != Sample_B ~ 0))%>%
  dplyr::mutate(env_biome_pair = paste0(env_biome_A, "_", env_biome_B))%>%
  dplyr::mutate(geo_loc_name_pair = paste0(geo_loc_name_A, "_", geo_loc_name_B))-> distances.arg.df

##Are mice resistomes more similar to any farm animal
##Check alpha diversity
ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::mutate(env_biome = fct_relevel(env_biome,
                                        "mouse", "pig", "cattle",
                                        "chicken", "turkey"))%>%
  dplyr::group_by(env_biome)%>%
  dplyr::mutate(median_arg= median(observed_arg))%>%
  dplyr::ungroup()%>%
  ggplot()+
  geom_violin(trim=F, aes(x= env_biome, y= observed_arg, color= env_biome, fill = env_biome),
              alpha= 0.5)+
  geom_boxplot(aes(x= env_biome, y= observed_arg), width = 0.01, outlier.shape = NA, color= "black")+
  geom_point(shape=21, size=2, alpha= 0.5,
             aes(fill = env_biome, x= env_biome, y=median_arg),
             color= "black")+
  scale_fill_manual(values = pal.host)+
  scale_color_manual(values = pal.host)+
  geom_hline(yintercept = max(ord.data%>%dplyr::filter(env_biome=="mouse")%>%pull(observed_arg)),
             linetype='dotted')+
  ylab("Number of ARGs per sample")+
  scale_x_discrete(labels=c("mouse" = "Mouse",
                            "pig" = "Pig",
                            "cattle" = "Cattle",
                            "chicken"="Chicken",
                            "turkey"="Turkey"))+
  labs(tag= "a")+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> F3A

##Check alpha diversity by country
ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::mutate(env_biome = fct_relevel(env_biome,
                                        "mouse", "pig", "cattle",
                                        "chicken", "turkey"))%>%
  dplyr::mutate(geo_loc_name = fct_relevel(geo_loc_name,
                                           "Belgium", "Bulgaria", "Denmark",
                                           "France",  "Germany",  "Ghana", "Italy",
                                           "Netherlands", "Poland","Spain"))%>%
  dplyr::group_by(env_biome, geo_loc_name)%>%
  dplyr::mutate(median_arg= median(observed_arg))%>%
  dplyr::ungroup()%>%
  ggplot()+
  geom_violin(trim=F, aes(x= env_biome, y= observed_arg, color= env_biome, fill = env_biome),
              alpha= 0.5)+
  geom_boxplot(aes(x= env_biome, y= observed_arg), width = 0.01, outlier.shape = NA, color= "black")+
  geom_point(shape=21, size=2, alpha= 0.5,
             aes(fill = env_biome, x= env_biome, y=median_arg),
             color= "black")+
  scale_fill_manual(values = pal.host)+
  scale_color_manual(values = pal.host)+
  geom_hline(yintercept = max(ord.data%>%dplyr::filter(env_biome=="mouse")%>%pull(observed_arg)),
             linetype='dotted')+
  ylab("Number of ARGs per sample")+
  scale_x_discrete(labels=c("mouse" = "Mouse",
                            "pig" = "Pig",
                            "cattle" = "Cattle",
                            "chicken"="Chicken",
                            "turkey"="Turkey"))+
  labs(tag= "a")+
  facet_wrap(~geo_loc_name, ncol=5)+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")-> SuppI

##Stats
require(rstatix)
ord.data%>%
  wilcox_test(observed_arg ~ env_biome, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.arg

ord.data%>%
  wilcox_effsize(observed_arg ~ env_biome)%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.arg)-> stats.test.arg

##Determine that the difference between mouse and pig is less than against the other hosts
require(emmeans)
model.richness <- lm(observed_arg ~ env_biome, data = ord.data)
emms.richness <- emmeans(model.richness, ~ env_biome)
pairs.richeness<- pairs(emmeans(model.richness, ~ env_biome), adjust = "bonferroni")

# Compare contrasts (e.g., difference in differences)
contrast_diff <-  contrast(emms.richness, method = "trt.vs.ctrl", ref = "mouse", adjust = "bonferroni") # example
summary(contrast_diff)

##Distances against mouse
##Determine "higher" similarity
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))%>%
  dplyr::group_by(env_biome_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(env_biome_pair, median_ait))%>%
  distinct()%>%
  dplyr::mutate(median_ait= 1-median_ait)%>%
  dplyr::arrange(desc(median_ait))

##Store the median
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))%>%
  dplyr::group_by(env_biome_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(env_biome_A, env_biome_pair, median_ait))%>%
  distinct()-> host.median

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))%>%
  dplyr::group_by(env_biome_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  ggplot()+
  geom_violin(aes(x= env_biome_pair, y= 1-Ait_arg_all, color=  env_biome_A, fill =  env_biome_A),
              alpha= 0.5, trim=F, draw_quantiles = T)+
  geom_boxplot(aes(x= env_biome_pair, y= 1-Ait_arg_all),
               width = 0.01, outlier.shape = NA, color= "black")+
  geom_point(data= host.median, shape=21, size=2, alpha= 1,
             aes(fill =  env_biome_A, x=  env_biome_pair, y=1-median_ait),
             color= "black")+
  scale_fill_manual(values = pal.host)+
  scale_color_manual(values = pal.host)+
  geom_hline(yintercept = 0.5,
             linetype='dotted')+
  scale_x_discrete(labels=c("mouse_mouse" = "Mouse",
                            "pig_mouse" = "Pig",
                            "cattle_mouse" = "Cattle",
                            "chicken_mouse"="Chicken",
                            "turkey_mouse"="Turkey"))+
  ylab("ARG Similarity")+
  labs(tag= "d")+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> F3D

##Determine that the difference between mouse and pig is less than against the other hosts
require(emmeans)
require(lme4)

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))-> lmm.data

model.distances <- lmer (data = lmm.data,
                         1-Ait_arg_all ~ env_biome_pair + geo_loc_name_pair + (1 | Sample_A) + (1 | Sample_B),
                         REML = F)
emms.distances <- emmeans(model.distances, ~ env_biome_pair)
pairs.distances<- pairs(emmeans(model.distances, ~ env_biome_pair), adjust = "bonferroni")

# Compare contrasts (e.g., difference in differences)
contrast.diff.dist <-  contrast(emms.distances, method = "trt.vs.ctrl",
                                ref = "mouse_mouse", adjust = "bonferroni")
summary(contrast.diff.dist)

##Stats
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))%>%
  wilcox_test(Bray_arg_all ~ env_biome_pair, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.dist

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("mouse_mouse", "chicken_mouse", "pig_mouse",
                                    "cattle_mouse", "turkey_mouse"))%>%
  dplyr::mutate(env_biome_pair = fct_relevel(env_biome_pair,
                                             "mouse_mouse", "pig_mouse", "cattle_mouse",
                                             "chicken_mouse", "turkey_mouse"))%>%
  wilcox_effsize(Bray_arg_all ~ env_biome_pair,)%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.dist)-> stats.test.dist

##Check number of genes shared between hosts
##ARG counts per sample
ARG.counts<- data.frame(t(PS.ARG.raw.all@otu_table[rowSums(PS.ARG.raw.all@otu_table) >= 1,]))

ARG.counts%>%
  rownames_to_column("Sample_ID")%>%
  pivot_longer(!Sample_ID, names_to = "ARO", values_to = "count")%>%
  dplyr::filter(count>=1)-> ARG.counts

ord.data%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::select(c(Sample_ID, env_biome))%>%
  dplyr::mutate(env_biome = fct_relevel(env_biome,
                                        "mouse", "pig", "cattle",
                                        "chicken", "turkey"))%>%
  left_join(ARG.counts)%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))-> ARG.counts

##For all samples
###Sample groups
sam.grp<- c( "mouse", "pig", "cattle",
             "chicken", "turkey")
list_core <- c() # an empty object to store information
for (n in sam.grp){ # for each variable n in Group
  tmp<- ARG.counts%>%
    dplyr::filter(env_biome==n)%>%
    dplyr::select(ARO)%>%
    unique()
  core_m <- tmp$ARO
  print(paste0("No. of ARGs in ", n, " : ", length(core_m))) # print otus identified in each group.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

require(eulerr)
plot(eulerr::euler(list_core),
     fills = c("#3E5496","#a06177","#8c785d","#d9af6b","#68855c"),
     alpha= 0.5,
     quantities = TRUE,
     labels = c( "Mouse", "Pig", "Cattle",
                 "Chicken", "Turkey"))-> F3B

cowplot::plot_grid(
  F3B, labels = "b",
  label_fontfamily = "sans",
  label_fontface = "plain",
  label_size = 16)-> F3B

ggsave(file = "Rplots.pdf", plot = F3B,
       width = 5, height =5, dpi = 600)

##Determine the number of intersected genes with Mouse
length(intersect(list_core$mouse, list_core$pig))
length(intersect(list_core$mouse, list_core$cattle))
length(intersect(list_core$mouse, list_core$chicken))
length(intersect(list_core$mouse, list_core$turkey))

##Run a Fisher test to determine the significance of the overlap
envs <- names(list_core)
pairs <- combn(envs, 2, simplify = FALSE)

fisher_results <- lapply(pairs, function(pair) {
  g1 <- list_core[[pair[1]]]
  g2 <- list_core[[pair[2]]]
  all_genes <- unique(unlist(list_core))

  k <- length(intersect(g1, g2))
  a <- length(setdiff(g1, g2))
  b <- length(setdiff(g2, g1))
  N <- length(all_genes)

  pval <- fisher.test(matrix(c(k, a, b, N - (a + b + k)), nrow = 2),
                      alternative = "greater")$p.value

  data.frame(env1 = pair[1], env2 = pair[2], overlap = k, p.value = pval)
})

fisher_results <- do.call(rbind, fisher_results)
fisher_results$p.adj <- p.adjust(fisher_results$p.value, method = "bonferroni")

##Prevalence plot to identify genes prevalent in both
data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="mouse"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_mouse")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)->prev.aro

data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="pig"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_pig")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="cattle"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_cattle")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

##Only "Germans"
data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="pig"&geo_loc_name=="Germany"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_pig_Ger")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="cattle"&geo_loc_name=="Germany"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_cattle_Ger")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

##
data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="chicken"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_chicken")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome=="turkey"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_turkey")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

data.frame(t(prevalence(subset_samples(PS.ARG.FPKM.comp, env_biome!="mouse"))))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_livestock")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

data.frame(t(prevalence(PS.ARG.FPKM.comp,)))%>%
  rownames_to_column("row")%>%
  pivot_longer(!row, names_to = "ARO", values_to = "prevalence_overall")%>%
  dplyr::mutate(ARO= gsub("\\.", "\\:", ARO))%>%
  dplyr::select(!row)%>%
  dplyr::left_join(prev.aro, by="ARO")->prev.aro

##add ARO information
data.frame(PS.ARG.FPKM.comp@tax_table)%>%
  rownames_to_column("ARO")%>%
  dplyr::select(c("ARO","Variant", "Resistance_Mechanism", "ARG_CARD_Short", "Drug_Class_adjusted"))%>%
  dplyr::left_join(prev.aro, by="ARO")-> prev.aro

##Mobility condition
arg.mobility<- read_csv("data/raw/arg_mobility.csv")

arg.mobility%>%
  dplyr::select(c(ARO, localization, mobility))-> arg.mobility

pal.abx.2<-c("Tetracycline"="#994F00", "Glycopeptide"="#006CD1",
             "Aminoglycoside"="#B672B8","Fluoroquinolone"="#FCB71B",
             "Phosphonic acid"="#24D582", "Peptide"="#F3CA95",
             "Diaminopyrimidine"="#748B31", "Phenicol"="#C5B51B",
             "MLS"="#89BDB6", "Beta lactam"="#E651BA",
             "Multidrug"="#A2A2A2","Other"="#e2dfe1")

##Now plot it!
require("ggrepel")
prev.aro%>%
  dplyr::filter(prevalence_mouse>0)%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::filter(!Drug_Class_adjusted %in% c("Other"))%>%
  dplyr::mutate(level = case_when((prevalence_mouse>=0.15|localization=="integron")&prevalence_pig>=0.05 ~ 1,
                                  T ~ 0))%>%
  dplyr::mutate(label= case_when(level==1~ ARG_CARD_Short,
                                 T~ ""))%>%
  ggplot(aes(x= prevalence_mouse*100, y= prevalence_pig*100, alpha=level,
             fill=  Drug_Class_adjusted, label= label, shape= mobility))+
  geom_hline(yintercept = 10,
             linetype='dotted')+
  geom_vline(xintercept = 10,
             linetype='dotted')+
  geom_point(size=3, color= "black")+
  scale_shape_manual(values = c(21, 24))+
  geom_text_repel(max.overlaps = 100)+
  scale_fill_manual(values = pal.abx.2)+
  labs(tag= "e", x= "Mouse gut (prevalence %)", y= "Pig manure (prevalence %)",
       fill= "Drug class", shape= "Mobility potential")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         alpha="none")+
  theme_bw()+
  theme(text = element_text(size=16))-> F3E

##Number of promoted, no promoted and co promoted genes
prev.aro%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::mutate(prom_class = case_when(prevalence_mouse>=0.10&prevalence_pig>=0.10 ~ "co-promoted",
                                       prevalence_mouse>=0.10&prevalence_pig<0.10 ~ "mouse-promoted",
                                       prevalence_mouse<0.10&prevalence_pig>=0.10 ~ "pig-promoted",
                                       T ~ "not-promoted"))%>%
  dplyr::select(!c(prevalence_overall,prevalence_livestock,prevalence_turkey,prevalence_chicken,
                   prevalence_cattle_Ger,prevalence_pig_Ger,prevalence_cattle))-> prev.mus.pig

##Cattle
prev.aro%>%
  dplyr::filter(prevalence_mouse>0)%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::filter(!Drug_Class_adjusted %in% c("Other"))%>%
  dplyr::mutate(level = case_when((prevalence_mouse>=0.15|localization=="integron")&prevalence_cattle>=0.05 ~ 1,
                                  T ~ 0))%>%
  dplyr::mutate(label= case_when(level==1~ ARG_CARD_Short,
                                 T~ ""))%>%
  ggplot(aes(x= prevalence_mouse*100, y= prevalence_cattle*100, alpha=level,
             fill=  Drug_Class_adjusted, label= label, shape= mobility))+
  geom_hline(yintercept = 10,
             linetype='dotted')+
  geom_vline(xintercept = 10,
             linetype='dotted')+
  geom_point(size=3, color= "black")+
  scale_shape_manual(values = c(21, 24))+
  geom_text_repel(max.overlaps = 100)+
  scale_fill_manual(values = pal.abx.2)+
  labs(tag= "f", x= "Mouse gut (prevalence %)", y= "Cattle manure (prevalence %)",
       fill= "Drug class", shape= "Mobility potential")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         alpha="none")+
  theme_bw()+
  theme(text = element_text(size=16))-> F3F

prev.aro%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::mutate(prom_class = case_when(prevalence_mouse>=0.10&prevalence_cattle>=0.10 ~ "co-promoted",
                                       prevalence_mouse>=0.10&prevalence_cattle<0.10 ~ "mouse-promoted",
                                       prevalence_mouse<0.10&prevalence_cattle>=0.10 ~ "cattle-promoted",
                                       T ~ "not-promoted"))%>%
  dplyr::select(!c(prevalence_overall,prevalence_livestock,prevalence_turkey,prevalence_chicken,
                   prevalence_cattle_Ger,prevalence_pig_Ger,prevalence_pig))-> prev.mus.cattle

##Poultry and turkey
prev.aro%>%
  dplyr::filter(prevalence_mouse>0)%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::filter(!Drug_Class_adjusted %in% c("Other"))%>%
  dplyr::mutate(level = case_when((prevalence_mouse>=0.15|localization=="integron")&prevalence_chicken>=0.05 ~ 1,
                                  T ~ 0))%>%
  dplyr::mutate(label= case_when(level==1~ ARG_CARD_Short,
                                 T~ ""))%>%
  ggplot(aes(x= prevalence_mouse*100, y= prevalence_chicken*100, alpha=level,
             fill=  Drug_Class_adjusted, label= label, shape= mobility))+
  geom_hline(yintercept = 10,
             linetype='dotted')+
  geom_vline(xintercept = 10,
             linetype='dotted')+
  geom_point(size=3, color= "black")+
  scale_shape_manual(values = c(21, 24))+
  geom_text_repel(max.overlaps = 100)+
  scale_fill_manual(values = pal.abx.2)+
  labs(tag= "a", x= "Mouse gut (prevalence %)", y= "Chicken manure (prevalence %)",
       fill= "Drug class", shape= "Mobility potential")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         alpha="none")+
  theme_bw()+
  theme(text = element_text(size=16))-> SuppE

##Turkey
prev.aro%>%
  dplyr::filter(prevalence_mouse>0)%>%
  dplyr::left_join(arg.mobility, by="ARO")%>%
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
  dplyr::filter(!Drug_Class_adjusted %in% c("Other"))%>%
  dplyr::mutate(level = case_when((prevalence_mouse>=0.15|localization=="integron")&prevalence_turkey>=0.05 ~ 1,
                                  T ~ 0))%>%
  dplyr::mutate(label= case_when(level==1~ ARG_CARD_Short,
                                 T~ ""))%>%
  ggplot(aes(x= prevalence_mouse*100, y= prevalence_turkey*100, alpha=level,
             fill=  Drug_Class_adjusted, label= label, shape= mobility))+
  geom_hline(yintercept = 10,
             linetype='dotted')+
  geom_vline(xintercept = 10,
             linetype='dotted')+
  geom_point(size=3, color= "black")+
  scale_shape_manual(values = c(21, 24))+
  geom_text_repel(max.overlaps = 100)+
  scale_fill_manual(values = pal.abx.2)+
  labs(tag= "b", x= "Mouse gut (prevalence %)", y= "Turkey manure (prevalence %)",
       fill= "Drug class", shape= "Mobility potential")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         alpha="none")+
  theme_bw()+
  theme(text = element_text(size=16))-> SuppF

###Arrange Figures
require("ggpubr")
Figure.5<-ggarrange(F3A, F3B, F3C, F3D, F3E, F3F, ncol = 2, nrow = 3)

##Check for Pigs
##Which genes are shared with pigs from different countries
##Determine "higher" similarity order
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("pig_mouse"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(geo_loc_name_pair, median_ait))%>%
  distinct()%>%
  dplyr::mutate(median_ait= 1-median_ait)%>%
  dplyr::arrange(desc(median_ait))

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("pig_mouse"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(geo_loc_name_A, geo_loc_name_pair, median_ait))%>%
  distinct()-> country.median

##Now plot it!
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("pig_mouse"))%>%
  dplyr::mutate(geo_loc_name_pair = fct_relevel(geo_loc_name_pair,
                                                "Denmark_Germany", "Netherlands_Germany", "France_Germany",
                                                "Belgium_Germany", "Germany_Germany", "Poland_Germany",
                                                "Italy_Germany", "Spain_Germany", "Bulgaria_Germany",
                                                "Ghana_Germany"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_bc= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  ggplot()+
  geom_violin(aes(x= geo_loc_name_pair, y= 1-Ait_arg_all, color=  geo_loc_name_A, fill =  geo_loc_name_A),
              alpha= 0.5, trim=F, draw_quantiles = T)+
  geom_boxplot(aes(x= geo_loc_name_pair, y= 1-Ait_arg_all),
               width = 0.01, outlier.shape = NA, color= "black")+
  geom_point(data=country.median, shape=21, size=2, alpha= 1,
             aes(fill =  geo_loc_name_A, x=  geo_loc_name_pair, y=1-median_ait),
             color= "black")+
  scale_fill_manual(values = pal.country)+
  scale_color_manual(values = pal.country)+
  geom_hline(yintercept = 0.5,
             linetype='dotted')+
  scale_x_discrete(labels=c("Germany_Germany"= "Germany",
                            "Denmark_Germany"= "Denmark",
                            "Poland_Germany"="Poland",
                            "Netherlands_Germany"="Netherlands",
                            "Belgium_Germany"="Belgium",
                            "France_Germany"="France",
                            "Italy_Germany"="Italy",
                            "Spain_Germany"="Spain",
                            "Bulgaria_Germany"="Bulgaria",
                            "Ghana_Germany"="Ghana"))+
  ylab("Mouse-Pig ARG Similarity")+
  labs(tag= "g")+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> SuppG

ggsave(file = "Rplots.pdf", plot = SuppG,
       width = 10, height =6, dpi = 600)

##Statistics
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("pig_mouse"))%>%
  dplyr::mutate(geo_loc_name_pair = fct_relevel(geo_loc_name_pair,
                                                "Germany_Germany","Denmark_Germany","Poland_Germany",
                                                "Netherlands_Germany", "Belgium_Germany", "France_Germany",
                                                "Italy_Germany", "Spain_Germany", "Bulgaria_Germany",
                                                "Ghana_Germany"))%>%
  wilcox_test(Ait_arg_all ~ geo_loc_name_pair, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.country

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("pig_mouse"))%>%
  dplyr::mutate(geo_loc_name_pair = fct_relevel(geo_loc_name_pair,
                                                "Germany_Germany","Denmark_Germany","Poland_Germany",
                                                "Netherlands_Germany", "Belgium_Germany", "France_Germany",
                                                "Italy_Germany", "Spain_Germany", "Bulgaria_Germany",
                                                "Ghana_Germany"))%>%
  wilcox_effsize(Ait_arg_all ~ geo_loc_name_pair)%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.country)-> stats.test.country

##Check for Cattle
##Which genes are shared with pigs from different countries
##Determine "higher" similarity order
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("cattle_mouse"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(geo_loc_name_pair, median_ait))%>%
  distinct()%>%
  dplyr::mutate(median_ait= 1-median_ait)%>%
  dplyr::arrange(desc(median_ait))

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("cattle_mouse"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(geo_loc_name_A, geo_loc_name_pair, median_ait))%>%
  distinct()-> country.median

##Now plot it!
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("cattle_mouse"))%>%
  dplyr::mutate(geo_loc_name_pair = fct_relevel(geo_loc_name_pair,
                                                "Germany_Germany",
                                                "Netherlands_Germany",
                                                "France_Germany"))%>%
  dplyr::group_by(geo_loc_name_pair)%>%
  dplyr::mutate(median_ait= median(Ait_arg_all))%>%
  dplyr::ungroup()%>%
  ggplot()+
  geom_violin(aes(x= geo_loc_name_pair, y= 1-Ait_arg_all, color=  geo_loc_name_A, fill =  geo_loc_name_A),
              alpha= 0.5, trim=F, draw_quantiles = T)+
  geom_boxplot(aes(x= geo_loc_name_pair, y= 1-Ait_arg_all),
               width = 0.01, outlier.shape = NA, fill="black",color= "black")+
  geom_point(data=country.median, shape=21, size=2, alpha= 1,
             aes(fill =  geo_loc_name_A, x=  geo_loc_name_pair, y=1-median_ait),
             color= "black")+
  scale_fill_manual(values = pal.country)+
  scale_color_manual(values = pal.country)+
  geom_hline(yintercept = 0.5,
             linetype='dotted')+
  scale_x_discrete(labels=c("Germany_Germany"= "Germany",
                            "Netherlands_Germany"="Netherlands",
                            "France_Germany"="France"))+
  ylab("Mouse-Cattle ARG Similarity")+
  labs(tag= "h")+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.title.x=element_blank(),
        legend.position = "none")-> SuppH

##Statistics
distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("cattle_mouse"))%>%
  wilcox_test(Ait_arg_all ~ geo_loc_name_pair, alternative = "two.sided")%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()-> stats.test.cattle

distances.arg.df%>%
  rownames_to_column("Sample_ID")%>%
  dplyr::filter(Same_sample==0)%>%
  dplyr::filter(env_biome_pair%in%c("cattle_mouse"))%>%
  wilcox_effsize(Ait_arg_all ~ geo_loc_name_pair)%>%
  dplyr::select(!c(n1, n2))%>%
  left_join(stats.test.cattle)-> stats.test.cattle


##Put all supplements
Supp1AB<-ggarrange(SuppA, SuppB, ncol = 2, nrow = 1, legend = "bottom", common.legend = T)
Supp1CD<-ggarrange(SuppC, SuppD, ncol = 2, nrow = 1, legend = "bottom", common.legend = T)
Supp1<-ggarrange(Supp1AB, Supp1CD, ncol = 1, nrow = 2, legend = "bottom", common.legend = T)

##Put all supplements
Supp2<-ggarrange(SuppE, SuppF, ncol = 2, nrow = 1, legend = "right", common.legend = T)

##Put all supplements
Supp3<-ggarrange(SuppG, SuppH, ncol = 2, nrow = 1)

##Which genes are shared with cattle and pigs
data.frame(PS.ARG.FPKM.all@tax_table)%>%
  rownames_to_column("ARO")%>%
  dplyr::filter(ARO%in%intersect(list_core$mouse, list_core$pig))-> pig.mouse.arg

data.frame(PS.ARG.FPKM.all@tax_table)%>%
  rownames_to_column("ARO")%>%
  dplyr::filter(ARO%in%setdiff(list_core$mouse, list_core$cattle))-> cattle.mouse.arg
