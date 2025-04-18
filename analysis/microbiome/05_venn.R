############### relative abundance at each site 1/9/23

## tutorial: https://github.com/surh/scip_barplot/blob/a775e0d97b47e8dd40498b6713016bbb7e045003/extended_example.Rmd#L285-L301

##what we want: relative abundance bar charts with:
## a facet for each bee genus at each site
## y axis is rel abund
## x axis is individual bee
## so 7 sites and four bee genera (rows will be site and cols will be genus)
## colored by ASV (family probs)
## maybe will have to subet to each genus then facet only by site but we will see
#rm(list=ls())
wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands'
setwd(wdpath)

library(stringr)
library(tidyverse)
library(ggVennDiagram)
library(gt)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(viridis)
library(UpSetR)
library(ComplexUpset)


source("analysis/microbiome/src/vennFunctions.R")

load("data/spec_RBCL_16s.Rdata")

spec16s <- spec.net %>%
  filter(Apidae == 1) %>%
  select(UniqueID, Site, Genus, GenusSpecies, starts_with('16s')) %>%
  na.omit() %>%
  pivot_longer(-c(UniqueID, Site, Genus, GenusSpecies), names_to = 'Bacteria', values_to = 'Abundance')

gensp_summary_table <- spec.net %>%
  filter(Apidae == 1) %>%
  select(GenusSpecies) %>%
  group_by(GenusSpecies) %>%
  summarize(n = n()) %>%
  arrange(desc(n))



apis_subset <- genus_filter(spec16s, 'Apis')
bombus_subset <- genus_filter(spec16s, 'Bombus')
anthophora_subset <- genus_filter(spec16s, 'Anthophora')
megachile_subset <- genus_filter(spec16s, 'Megachile')
melissodes_subset <- genus_filter(spec16s, 'Melissodes')


venn_data <- list(Anthophora = anthophora_subset,
                  Apis = apis_subset,
                  Bombus = bombus_subset,
                  Megachile = megachile_subset,
                  Melissodes = melissodes_subset)


ggVennDiagram(venn_data,
              label_alpha = 0.5) + 
  scale_color_manual(values = plasma(5)) +
  scale_fill_gradient(low="white",high = "black") +
  scale_x_continuous(expand = expansion(mult = .2))

upset_genus <- upset(fromList(venn_data),
                     colnames(fromList(venn_data)),
                     queries=query_by_degree(
                       fromList(venn_data),
                       colnames(fromList(venn_data)),
                       params_by_degree=list(
                         '0'=list(color='black', fill='black'),
                         '1'=list(color='black', fill='black'),
                         '2'=list(color='black', fill='black'),
                         '3'=list(color='black', fill='black'),
                         '4'=list(color='black', fill='black'),
                         '5'=list(color='red', fill='red')
                       )))
upset_genus

#genus intersections
# intersection_vals <- process_region_data(Venn(venn_data))
# 
# shared_otus_all <- unlist(intersection_vals$item[intersection_vals$id == '12345'])
# 
# bombus_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Bombus'])
# 
# apis_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Apis'])
# 
# megachile_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Megachile'])
# 
# melissodes_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Melissodes'])
# 
# anthophora_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Anthophora'])
# 



#now making a nice table to show the shared OTUs between all genera
# shared_table_data <- spec16s %>%
#   filter(Bacteria %in% shared_otus_all) %>%
#   select(UniqueID, Bacteria, Abundance, GenusSpecies) %>%
#   filter(Abundance > 0) %>%
#   mutate(Feature.ID = row_number()) %>%
#   mutate(Taxon = Bacteria) %>%
#   select(-Bacteria) %>%
#   mutate(Taxon = str_replace_all(Taxon, '16s:', '')) %>%
#   na.omit() %>%
#   filter(grepl("Bifidobacterium_indicum|Lactobacillus_apis|Fructobacillus_tropaeoli|Bartonella_apis|Rickettsia_bellii|Frischella_perrara|Acinetobacter_nectaris", Taxon)) %>%
#   group_by(GenusSpecies, Taxon) %>%
#   mutate(num_indiv = n()) %>%
#   filter(GenusSpecies %in% c("Apis mellifera","Bombus huntii","Melissodes confusus","Bombus centralis","Bombus bifarius")) %>%
#   ungroup() %>%
#   group_by(Taxon)
# 
# shared_table_data$Taxon = str_remove(shared_table_data$Taxon, '.*(?=s__)')
# 
# shared_boxplots <- ggplot(shared_table_data, 
#                          aes(x=GenusSpecies, y=Abundance)) +
#                   geom_jitter(alpha=0.2) +
#                   geom_boxplot() + facet_wrap(~Taxon)
# shared_boxplots
# gt_shared <- gt(shared_table_data) %>%
#   sub_missing(columns=1:7, missing_text = '') %>%
#   tab_row_group(label = 'Present in all genera',
#                 rows=1:31)
# 
# 
# gt_shared

## what I want:
## 1. choose several strains
## 2. make violin plot for each strain 
##      One panel for each strain 
##      5 boxes for relative abundance in each species




### subsetting different intersections based on how many individuals we have of each
### will do one for ones >50, >20, >10, and > 4


## >= 50 indiv
Amellifera_subset <- species_filter(spec16s, 'Apis mellifera') 
Bcentralis_subset <- species_filter(spec16s, 'Bombus centralis')
Bhuntii_subset <- species_filter(spec16s, 'Bombus huntii')
Bbifarius_subset <- species_filter(spec16s, 'Bombus bifarius')
Mconfusus_subset <- species_filter(spec16s, 'Melissodes confusus')
## >= 20 indiv
Amontana_subset <- species_filter(spec16s, 'Anthophora montana')
Brufocinctus_subset <- species_filter(spec16s, 'Bombus rufocinctus')
Bflavifrons_subset <- species_filter(spec16s, 'Bombus flavifrons')
Mfrigida_subset <- species_filter(spec16s, 'Megachile frigida')
## >= 10 indiv
Bnevadensis_subset <- species_filter(spec16s, "Bombus nevadensis")
Mcomata_subset <- species_filter(spec16s, "Megachile comata")
## >= 4 indiv
Bfervidus_subset <- species_filter(spec16s, "Bombus fervidus")
Aurbana_subset <- species_filter(spec16s, "Anthophora urbana")
Bmixtus_subset <- species_filter(spec16s, "Bombus mixtus")
Bmorrisoni_subset <- species_filter(spec16s, "Bombus morrisoni")
Bsonorus_subset <- species_filter(spec16s, "Bombus sonorus")
Dmaura_subset <- species_filter(spec16s, "Dufourea maura")
Mrelativa_subset <- species_filter(spec16s, "Megachile relativa")



venn_data_50 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset)

venn_50 <- ggVennDiagram(venn_data_50,
                         label_alpha = 0.5) + 
  scale_color_manual(values = plasma(5)) +
  scale_fill_gradient(low="white",high = "black") +
  scale_x_continuous(expand = expansion(mult = .2))
venn_50

upset_50 <- upset(fromList(venn_data_50),
                  colnames(fromList(venn_data_50)),
                  queries=query_by_degree(
                    fromList(venn_data_50),
                    colnames(fromList(venn_data_50)),
                    params_by_degree=list(
                      '0'=list(color='black', fill='black'),
                      '1'=list(color='black', fill='black'),
                      '2'=list(color='black', fill='black'),
                      '3'=list(color='black', fill='black'),
                      '4'=list(color='black', fill='black'),
                      '5'=list(color='red', fill='red')
                    )))
upset_50


## determining shared 
intersection_vals <- process_region_data(Venn(venn_data_50))

shared_otus_all <- unlist(intersection_vals$item[intersection_vals$id == '12345'])

#now making a nice table to show the shared OTUs between all genera
shared_table_data <- spec16s %>%
  filter(Bacteria %in% shared_otus_all) %>%
  select(UniqueID, Bacteria, Abundance) %>%
  filter(Abundance > 0) %>%
  filter(!duplicated(Bacteria)) %>%
  select(-Abundance) %>%
  mutate(Feature.ID = row_number()) %>%
  mutate(Taxon = Bacteria) %>%
  select(-c(Bacteria, UniqueID)) %>%
  mutate(Taxon = str_replace_all(Taxon, '16s:', '')) %>%
  parse_taxonomy() %>%
  mutate(Kingdom = str_replace_all(Kingdom, 'd__', '')) %>%
  arrange(Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(Genus = str_replace_all(Genus, '_', ' ')) %>%
  mutate(Species = str_replace_all(Species, '_', ' ')) %>%
  #mutate(which_genera = 'Shared Among All') %>%
  as_tibble()

boxplot_data <- spec16s %>%
  filter(Bacteria %in% shared_otus_all) %>%
  select(UniqueID, Bacteria, Abundance, GenusSpecies) %>%
  filter(Abundance > 0) %>%
  mutate(Feature.ID = row_number()) %>%
  mutate(Taxon = Bacteria) %>%
  select(-Bacteria) %>%
  mutate(Taxon = str_replace_all(Taxon, '16s:', '')) %>%
  na.omit() %>%
  filter(grepl("Yersiniaceae|Erwiniaceae|Enterobacteriaceae|Acetobacteraceae|Streptococcaceae|Leuconostocaceae|Lactobacillaceae|Leuconostocaceae|Hafniaceae|Orbaceae|Neisseriaceae|Bifidobacteriaceae|Moraxellaceae", Taxon)) %>%
  group_by(GenusSpecies, Taxon) %>%
  mutate(num_indiv = n()) %>%
  filter(GenusSpecies %in% c("Apis mellifera","Bombus huntii","Melissodes confusus","Bombus centralis","Bombus bifarius")) %>%
  ungroup() 

boxplot_data$Taxon <- str_replace_all(boxplot_data$Taxon, '; s__.*', "")
boxplot_data$Taxon <- str_replace_all(boxplot_data$Taxon, '; g__.*', "")
boxplot_data$Taxon <- str_remove(boxplot_data$Taxon, '.*(?=f__)') 


boxplot_data <- boxplot_data %>%
  group_by(Taxon, GenusSpecies) %>%
  filter(!Taxon == 'g__uncultured')



shared_boxplots <- ggplot(boxplot_data,
                          aes(x=GenusSpecies, y=Abundance)) +
  geom_jitter(aes(fill=GenusSpecies), shape=21) +
  geom_boxplot(aes(fill=GenusSpecies),alpha=0.2) +
  facet_wrap(~Taxon) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)
shared_boxplots

### now do overlap between just the species that are in > 20

venn_data_20 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset,
                     Anthophora_montana = Amontana_subset,
                     Bombus_rufocinctus = Brufocinctus_subset,
                     Bombus_flavifrons = Bflavifrons_subset,
                     Megachile_frigida = Mfrigida_subset)

upset_20 <- upset(fromList(venn_data_20),
                  colnames(fromList(venn_data_20)),
                  n_intersections=15)
# queries=query_by_degree(
#   fromList(venn_data_20),
#   colnames(fromList(venn_data_20)),
#   params_by_degree=list(
#     '0'=list(color='black', fill='black'),
#     '1'=list(color='black', fill='black'),
#     '2'=list(color='black', fill='black'),
#     '3'=list(color='black', fill='black'),
#     '4'=list(color='black', fill='black'),
#     '5'=list(color='black', fill='black'),
#     '6'=list(color='black', fill='black'),
#     '7'=list(color='black', fill='black'),
#     '8'=list(color='black', fill='black'),
#     '9'=list(color='red', fill='red')
#   )))
upset_20

## now do 10
venn_data_10 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset,
                     Anthophora_montana = Amontana_subset,
                     Bombus_rufocinctus = Brufocinctus_subset,
                     Bombus_flavifrons = Bflavifrons_subset,
                     Megachile_frigida = Mfrigida_subset,
                     Bombus_nevadensis = Bnevadensis_subset,
                     Megachile_comata = Mcomata_subset)

upset_10 <- upset(fromList(venn_data_10), 
                  colnames(fromList(venn_data_10)),
                  n_intersections=15)
upset_10
## now do 4
venn_data_4 <- list(Apis_mellifera = Amellifera_subset,
                    Bombus_centralis = Bcentralis_subset,
                    Bombus_huntii = Bhuntii_subset,
                    Bombus_bifarius = Bbifarius_subset,
                    Melissodes_confusus = Mconfusus_subset,
                    Anthophora_montana = Amontana_subset,
                    Bombus_rufocinctus = Brufocinctus_subset,
                    Bombus_flavifrons = Bflavifrons_subset,
                    Megachile_frigida = Mfrigida_subset,
                    Bombus_nevadensis = Bnevadensis_subset,
                    Megachile_comata = Mcomata_subset,
                    Bombus_fervidus = Bfervidus_subset,
                    Anthophora_urbana = Aurbana_subset,
                    Bombus_mixtus = Bmixtus_subset,
                    Bombus_morrisoni = Bmorrisoni_subset,
                    Bombus_sonorus = Bsonorus_subset,
                    Dufourea_maura = Dmaura_subset)

upset_4 <- upset(fromList(venn_data_4),
                 colnames(fromList(venn_data_4)))
upset_4