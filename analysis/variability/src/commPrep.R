
calcSiteBeta <- function(x, species.type, spec, species.type.int,
                         date.cut.off = 2,
                         observation.cut.off = 2,
                         type.net){
    ## count interactions
    prep.comm <- aggregate(list(Abund=spec[, species.type]),
                           list(GenusSpecies=spec[, species.type],
                                InterGenusSpecies=
                                    spec[, species.type.int],
                                Site=spec$Site,
                                Date=spec$SampleRound,
                                Year=spec$Year),
                           length)

    ## subset to a single year/site
    prep.comm <-  prep.comm[prep.comm[, type.net] == x,]

    if(type.net == "Year"){
        ## take the average across sample rounds
        prep.comm <- aggregate(list(Abund=prep.comm$Abund),
                               list(GenusSpecies=prep.comm$GenusSpecies,
                                    InterGenusSpecies=
                                        prep.comm$InterGenusSpecies,
                                    Site=prep.comm$Site,
                                    Year=prep.comm$Year),
                               mean)

    }
    prep.comm$SiteYear <- paste(prep.comm$Site,
                                prep.comm$Year,
                                ## prep.comm$Date,
                                sep=":")
    prep.comm <- prep.comm[!prep.comm$InterGenusSpecies == "",]
    by.species <- split(prep.comm, prep.comm$GenusSpecies)

    num.observations <- sapply(by.species, function(x) sum(x$Abund))
    ## subset to species seen in > date.cut.off years and at least
    ## observation.cut.off times
    by.species <- by.species[num.observations >= observation.cut.off]
    num.dates <- sapply(by.species, function(x) length(unique(x$SiteYear)))
    by.species <- by.species[num.dates >= date.cut.off]
    by.species <- by.species[!sapply(by.species, is.null)]

    ## year plant combinations
    empty.matrix <- matrix(0, nrow=length(unique(prep.comm$SiteYear)),
                           ncol=length(unique(prep.comm$InterGenusSpecies)))

    rownames(empty.matrix) <- sort(unique(prep.comm$SiteYear))
    colnames(empty.matrix) <-
        sort(unique(prep.comm$InterGenusSpecies))

    if(length(by.species) != 0){
        comm <- vector("list", length=length(by.species))
        for(i in 1:length(by.species)){
            comm[[i]] <- empty.matrix
            this.by.species <- by.species[[i]]
            for(j in 1:nrow(this.by.species)){
                this.row <- this.by.species[j,]
                comm[[i]][match(this.row["SiteYear"], rownames(comm[[i]])),
                          match(this.row["InterGenusSpecies"],
                                colnames(comm[[i]]))] <-
                    as.numeric(this.row[["Abund"]])
            }
            comm[[i]] <- comm[[i]][rowSums(comm[[i]]) > 0,]
            if(!is.matrix(comm[[i]])) comm[[i]] <- NA
        }
        names(comm) <- names(by.species)
        comm <- comm[!sapply(comm, function(x) all(is.na(x)))]
        return(list(comm=comm, year=x))
    }
}

makePretty <- function(comms, spec, net.type){
    ## year is a placeholder for whatever the networks are split by
    ## (year, site)
    year <- sapply(comms, function(x) x$year)
    comms <- lapply(comms, function(x) x$comm)
    comms <- comms[!sapply(year, is.null)]
    year <- unlist(year[!sapply(year, is.null)])
    names(comms) <- year
    site.date <- lapply(comms, function(x) sapply(x, rownames))

    if(net.type == "Year"){
        comm.pp <- list(comm=comms,
                        site.date=site.date,
                        years= rep(names(comms),
                                   sapply(comms, length)))
    }

    if(net.type == "Site"){
        comm.pp <- list(comm=comms,
                        site.date=site.date,
                        sites= rep(names(comms),
                                   sapply(comms, length)))
    }

    return(comm.pp)
}
