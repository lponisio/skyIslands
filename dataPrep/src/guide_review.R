
guild.plot.data <- spec.net %>%
  dplyr::mutate(
    Year=as.character(Year),
    GenusLabel=trimws(as.character(Genus)),
    GenusLabel=ifelse(is.na(GenusLabel) | !nzchar(GenusLabel),
                      'Unidentified genus', GenusLabel),
    SpeciesLabel=trimws(as.character(Species)),
    IDResolution=dplyr::case_when(
      GenusLabel == 'Unidentified genus' ~ 'Unidentified genus',
      is.na(SpeciesLabel) | !nzchar(SpeciesLabel) |
        tolower(SpeciesLabel) %in% c('sp', 'sp.', 'spp', 'spp.',
                                     'unknown', 'unidentified') ~ 'Genus only',
      TRUE ~ 'Species / morphospecies'
    )
  )

## Guild membership is assigned by Family, the All bees panel is where
## to look for unidentified bees.
plot.guilds <- dplyr::bind_rows(
  guild.plot.data %>% dplyr::filter(Family %in% bee.families) %>%
    dplyr::mutate(Guild='All bees'),
  guild.plot.data %>% dplyr::filter(GenusSpecies == 'Apis mellifera') %>%
    dplyr::mutate(Guild='Honey bees'),
  guild.plot.data %>% dplyr::filter(Genus %in% 'Bombus') %>%
    dplyr::mutate(Guild='Bumble bees'),
  guild.plot.data %>%
    dplyr::filter(Family %in% butterfly.hawk.families) %>%
    dplyr::mutate(Guild='Butterflies + hawk moths'),
  guild.plot.data %>% dplyr::filter(Family %in% 'Syrphidae') %>%
    dplyr::mutate(Guild='Syrphids'),
  guild.plot.data %>% dplyr::filter(Family %in% 'Bombyliidae') %>%
    dplyr::mutate(Guild='Bee flies'),
  guild.plot.data %>%
    dplyr::filter(Family %in% focal.pollinator.families) %>%
    dplyr::mutate(Guild='Bees + syrphids + butterflies/hawks'),
  guild.plot.data %>%
    dplyr::filter(Family %in% all.pollinator.families) %>%
    dplyr::mutate(Guild='Bees + syrphids + bee flies + butterflies/hawks')
)
guild.order <- c('All bees', 'Honey bees', 'Bumble bees',
                 'Butterflies + hawk moths', 'Syrphids', 'Bee flies',
                 'Bees + syrphids + butterflies/hawks',
                 'Bees + syrphids + bee flies + butterflies/hawks')
plot.guilds$Guild <- factor(plot.guilds$Guild, levels=guild.order)

## Full genus counts and per-year identification status.
## These tables make years with blank genus IDs apparent rather than silently
## dropping them from barplots. Individual genera can be inspected in CSV.
guild.genus.year <- plot.guilds %>%
  dplyr::count(Guild, Year, GenusLabel, name='SpecimenCount')
guild.id.year <- plot.guilds %>%
  dplyr::count(Guild, Year, IDResolution, name='SpecimenCount') %>%
  dplyr::group_by(Guild, Year) %>%
  dplyr::mutate(TotalIndividuals=sum(SpecimenCount),
                Proportion=SpecimenCount/TotalIndividuals) %>%
  dplyr::ungroup()

write.csv(guild.genus.year, file.path(dir.guild.audit,
          'net_genus_counts_by_guild_year.csv'),
          row.names=FALSE)
write.csv(guild.id.year, file.path(dir.guild.audit,
          'net_id_resolution_by_guild_year.csv'),
          row.names=FALSE)

##  Count distinct identified genera per guild-year.
## The separate UnidentifiedCount makes low taxonomic resolution visible.
guild.genus.richness <- guild.genus.year %>%
  dplyr::group_by(Guild, Year) %>%
  dplyr::summarise(NGenera=sum(GenusLabel != 'Unidentified genus'),
                   UnidentifiedCount=sum(SpecimenCount[
                     GenusLabel == 'Unidentified genus']),
                   TotalIndividuals=sum(SpecimenCount), .groups='drop')
write.csv(guild.genus.richness, file.path(dir.guild.audit,
          'net_number_of_genera_by_guild_year.csv'), row.names=FALSE)

## Some raw Lepidoptera rows lack Family and cannot be assigned to butterflies
## versus other moths. Review those rather than assuming they are butterflies.
guild.unassigned <- guild.plot.data %>%
  dplyr::filter(is.na(Family) | !nzchar(trimws(as.character(Family))) |
                (Order %in% 'Lepidoptera' &
                 !(Family %in% butterfly.hawk.families)) |
                Family %in% c('Apidaae')) %>%
  dplyr::count(Year, Order, Family, GenusLabel, name='SpecimenCount')
write.csv(guild.unassigned,
          file.path(dir.guild.audit, 'net_guild_unassigned_taxa_by_year.csv'),
          row.names=FALSE)

##  Faceted stacked genus barplot by year. The diagnostic
## CSV above is the source for exact counts; many genera can make the legend
## large, so a second per-guild plot is written below with genus on the axis.
genus.overview <- ggplot2::ggplot(guild.genus.year,
    ggplot2::aes(x=Year, y=SpecimenCount, fill=GenusLabel)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~Guild, ncol=2, scales='free_y') +
  ggplot2::labs(x='Year', y='Netted individuals', fill='Genus',
    title='Netted pollinators by genus and year',
    subtitle='Guilds overlap; unidentified genera are included where Family is known') +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                 legend.position='none')
ggplot2::ggsave(file.path(dir.guild.audit,
  'net_genus_stacked_by_guild_year.png'),
  genus.overview, width=13, height=11, dpi=180)

## One barplot per guild, grouped by year within genus.
## This is easier for finding individual genera with missing identification
## coverage in particular years than the all-guild stacked overview.
plot.slugs <- c('all_bees', 'honey_bees', 'bumble_bees',
                'butterflies_hawkmoths', 'syrphids', 'bee_flies',
                'combined', 'combined_with_bee_flies')
for (i in seq_along(guild.order)) {
  this.guild <- guild.genus.year[
    as.character(guild.genus.year$Guild) == guild.order[i], ,drop=FALSE]
  if (!nrow(this.guild)) next
  this.guild$GenusLabel <- factor(this.guild$GenusLabel,
    levels=rev(sort(unique(as.character(this.guild$GenusLabel)))))
  p <- ggplot2::ggplot(this.guild,
      ggplot2::aes(x=GenusLabel, y=SpecimenCount, fill=Year)) +
    ggplot2::geom_col(position='dodge') + ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x='Genus', y='Netted individuals', fill='Year',
                  title=paste0('Genus counts by year: ', guild.order[i]))
  plot.height <- max(5, min(28, 2.5 + 0.32*nlevels(this.guild$GenusLabel)))
  ggplot2::ggsave(file.path(dir.guild.audit, paste0(
                'net_genera_', plot.slugs[i], '_by_year.png')),
                  p, width=11, height=plot.height,
                  dpi=180, limitsize=FALSE)


## Number of distinct identified genera per year.
## This complements the abundance bars, which count individuals per genus.
genus.number.plot <- ggplot2::ggplot(guild.genus.richness,
    ggplot2::aes(x=Year, y=NGenera)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~Guild, ncol=2, scales='free_y') +
  ggplot2::theme_bw() +
  ggplot2::labs(x='Year', y='Identified genera',
                title='Number of identified genera by guild and year') +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
ggplot2::ggsave(file.path(dir.guild.audit,
                'net_number_of_genera_by_guild_year.png'),
                genus.number.plot, width=13, height=10, dpi=180)

## Separate ID-resolution barplot by year and guild.
## Species/morphospecies, genus-only, and unidentified genus are separate.
## The percentages answer how identification completeness changes by year.
id.plot <- ggplot2::ggplot(guild.id.year,
    ggplot2::aes(x=Year, y=Proportion, fill=IDResolution)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~Guild, ncol=2) +
  ggplot2::scale_y_continuous(labels=function(x) paste0(round(x*100), '%')) +
  ggplot2::theme_bw() +
  ggplot2::labs(x='Year', y='Proportion of netted individuals',
                fill='Identification',
                title='Taxonomic identification resolution by guild and year') +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
ggplot2::ggsave(file.path(dir.guild.audit,
                'net_guild_ID_completeness_by_year.png'),
                id.plot, width=13, height=10, dpi=180)
}
