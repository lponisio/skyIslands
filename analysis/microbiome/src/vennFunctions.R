genus_filter <- function(dataframe, genus){
  new_df <- dataframe %>%
    filter(Genus == genus) %>%
    filter(Abundance > 0) %>%
    select(Bacteria)
  unique(new_df$Bacteria)
}

species_filter <- function(dataframe, species){
  new_df <- dataframe %>%
    filter(GenusSpecies == species) %>%
    filter(Abundance > 0) %>%
    select(Bacteria)
  unique(new_df$Bacteria)
}