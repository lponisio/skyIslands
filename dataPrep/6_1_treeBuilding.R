
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html

## packages
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(tibble)
library(pals)

## working dir

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

## Data imports

spec16s <- read.csv('spec_RBCL_16s.csv')

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

meta <- spec16s %>%
  select(all_of(meta_cols), Apidae) %>%
  filter(Apidae == 1)

##########################

##match unique IDs
match_shared_ID <- function(first_df, second_df){
  shared <- first_df$UniqueID[first_df$UniqueID %in% second_df$UniqueID]
  matched_df <- first_df %>%
  filter(UniqueID %in% shared)
  matched_df
}

match_shared_tiplabels <- function(tree, pres_abs_table){
  tree_tips <- tree$tip.label
  #browser()
  all_cols <- pres_abs_table %>%
    select(-UniqueID) 
  #browser()
  match_cols <- all_cols[colnames(all_cols) %in% tree_tips]
  #browser()
  match_cols$UniqueID <- pres_abs_table$UniqueID
  match_cols
}

matched_presabs <- match_shared_tiplabels(tree.16s, comm_presabs)
matched_pres_meta <- match_shared_ID(matched_presabs, meta)


comm_presabs <- as.data.frame(indiv.comm.16sR0) 
comm_presabs[comm_presabs > 0] <- 1 #converts from abundance to P/A -- check if needs to actualy do this

comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID")
matched_pres_meta <- match_shared_ID(comm_presabs, meta)
matched_id <- matched_pres_meta$UniqueID
row.names(matched_pres_meta) <- matched_id

#making presence/abs columns in metadata
meta_match_sites <- match_shared_ID(meta, matched_pres_meta) %>%
  select(UniqueID, Site) %>%
  mutate(Site = as.factor(Site)) %>%
  group_by(UniqueID, Site) %>%
  count() %>%
  pivot_wider(UniqueID, 
              names_from=Site, 
              values_from = n, 
              values_fill=0,
              names_expand = TRUE,
              id_expand=TRUE) %>%
  pivot_longer(cols=CH:SM,
               names_to='Site',
               values_to='Site_present')


meta_match_genus <- match_shared_ID(meta, matched_pres_meta) %>%
  select(UniqueID, Genus) %>%
  mutate(Site = as.factor(Genus)) %>%
  group_by(UniqueID, Genus) %>%
  count() %>%
  pivot_wider(UniqueID, 
              names_from=Genus, 
              values_from = n, 
              values_fill=0,
              names_expand = TRUE,
              id_expand=TRUE) %>%
  pivot_longer(cols=Agapostemon:Megachile,
               names_to='Genus',
               values_to='Genus_present')

## now need to figure out how to incorporate
## the tip labels so that the features are 
## tied to the metadata -- features as rows
## columns for metadata

##matched_pres_meta is a table with unique ID 
## as rows and features as columns while both meta_match
## have unique ID as rows and metadata as columns
## need to create an item that has features as rows
## and metadata as columns



features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
  right_join(meta_match_sites, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>% 
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))



features_genus_metadata <- match_shared_ID(matched_pres_meta, meta_match_genus) %>%
  right_join(meta_match_genus, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>% 
  mutate(Genus = factor(Genus, levels=c("Apis", ## ordered by relative (eyeballed) PCA microbiome similarity
                                      "Bombus",
                                      "Anthophora",
                                      "Megachile",
                                      "Agapostemon")))
  

## probs want to use trimmed tree for final
## visualization but right now will retain


#p <- ggtree(tree.16s, layout='rectangular')

p <- ggtree(tree.16s, layout="fan", open.angle=10, size=0.5)


p2 <- p + 
  geom_fruit(
    data=features_site_metadata,
    geom=geom_tile,
    mapping=aes(y=bacteria, 
                x=Site, 
                alpha=as.factor(Site_present),
                fill=Site_present)) 

p3 <- p + 
  geom_fruit(
    data=features_genus_metadata,
    geom=geom_tile,
    mapping=aes(y=bacteria, 
                x=Genus, 
                alpha=as.factor(Genus_present),
                fill=Genus)) 
p3

############### relative abundance at each site 1/9/23

## tutorial: https://github.com/surh/scip_barplot/blob/a775e0d97b47e8dd40498b6713016bbb7e045003/extended_example.Rmd#L285-L301

##what we want: relative abundance bar charts with:
## a facet for each bee genus at each site
## y axis is rel abund
## x axis is individual bee
## so 7 sites and four bee genera (rows will be site and cols will be genus)
## colored by ASV (family probs)
## maybe will have to subet to each genus then facet only by site but we will see




relabund.dat <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(UniqueID, Site, Genus, starts_with('X16s')) %>%
  na.omit() %>%
  pivot_longer(-c(UniqueID, Site, Genus), names_to = 'Bacteria', values_to = 'Abundance') %>% 
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))

# #change the column names to just include the bacteria family
relabund.dat$Bacteria <- gsub(".*D_4__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
relabund.dat$Bacteria <- gsub('\\.D_5__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

# #change the column names to just include the bacteria order
# relabund.dat$Bacteria <- gsub(".*D_3__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
# relabund.dat$Bacteria <- gsub('\\.D_4__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

## some issues with resolution: for now to get around this i am filtering out samples that included anything that was not resolved to family

library(randomcoloR)

relabund.dat.clean <- relabund.dat %>%
  filter(!grepl('X16s', Bacteria)) %>%
  filter(Genus != 'Agapostemon') %>%
  filter(Abundance > 0.01) ## too many groups -- decide what is the cutoff to show on relabund bars



n <- length(unique(relabund.dat.clean$Bacteria))
palette <- distinctColorPalette(n)
strains <- unique(relabund.dat.clean$Bacteria)

color_dict <- data.frame(palette, strains)


## good enough for now probs

## now to make prelim plot

# relabund.dat.clean %>%
#   ggplot(aes(x=UniqueID, y=Abundance)) +
#   geom_bar(aes(fill=Bacteria), stat='identity', position='fill') +
#   facet_grid(rows=vars(relabund.dat.clean$Site), cols=vars(relabund.dat.clean$Genus))

## yuck :( okay needs to do each genus separately

plot_genus_by_site_relabund <- function(data, genus){

  genus_subset <- data %>%
    filter(Genus == genus) %>%
    arrange(Bacteria)

  these_colors <- color_dict %>%
    filter(strains %in% genus_subset$Bacteria) %>%
    arrange(strains)

  ggplot(data=genus_subset, aes(x=UniqueID, y=Abundance)) +
    geom_bar(aes(fill=Bacteria), stat='identity', position='fill', width=1) +
    facet_grid(~Site, scales='free_x', space='free_x') +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.text.y = element_text(color = "black")) +
    scale_fill_manual(values=as.vector(these_colors$palette), name = "Bacteria Family") + 
    labs(x='Individuals', y='Relative Abundance')
}

apis_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Apis')
apis_plot

bombus_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Bombus')
bombus_plot

anthophora_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Anthophora')
anthophora_plot

megachile_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Megachile')
megachile_plot

## lets see how they look combined 
library(ggpubr)

ggarrange(
  apis_plot, 
  bombus_plot,
  anthophora_plot,
  megachile_plot,
  labels = c("Ap", "B", "An", "M"),
  common.legend = FALSE, legend = "bottom"
)