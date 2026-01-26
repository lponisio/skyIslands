


calcSiteBeta <- function(x, species.type, spec, species.type.int,
                         turnover.cut.off = 3,
                         type.net){
    if(type.net == "Year"){
        col.agg <- "Site"
    } else if(type.net == "Site"){
        col.agg <- "Year"
    } else if(type.net == "SiteYear"){
        col.agg <- "Date"
    }

    ## count interactions
    prep.comm <- aggregate(list(Abund=spec[, species.type]),
                           list(GenusSpecies=spec[, species.type],
                                InterGenusSpecies=
                                    spec[, species.type.int],
                                Site=spec$Site,
                                Date=spec$SampleRound,
                                Year=spec$Year),
                           length)

    prep.comm$SiteYear <- paste(prep.comm$Site,
                                prep.comm$Year,
                                sep=":")

    ## subset to a single year/site
    prep.comm <-  prep.comm[prep.comm[, type.net] == x,]

    if(type.net == "Year" | type.net == "Site"){
        ## take the average across sample rounds
        prep.comm <- aggregate(list(Abund=prep.comm$Abund),
                               list(GenusSpecies=prep.comm$GenusSpecies,
                                    InterGenusSpecies=
                                        prep.comm$InterGenusSpecies,
                                    Site=prep.comm$Site,
                                    Year=prep.comm$Year),
                               mean)
    }
    prep.comm <- prep.comm[!prep.comm$InterGenusSpecies == "",]
    by.species <- split(prep.comm, prep.comm$GenusSpecies)
    ## subset to species seen in > turnover cut off (needs to be
    ## enough unique observations to use for beta-diversity
    ## calculation)
    num.4.turnover <- sapply(by.species, function(x)
        length(unique(x[,col.agg])))

    by.species <- by.species[num.4.turnover >= turnover.cut.off]
    by.species <- by.species[!sapply(by.species, is.null)]

    ## year/site plant combinations
    empty.matrix <- matrix(0, nrow=length(unique(prep.comm[, col.agg])),
                           ncol=length(unique(prep.comm$InterGenusSpecies)))
    rownames(empty.matrix) <- sort(unique(prep.comm[,col.agg]))
    colnames(empty.matrix) <-
        sort(unique(prep.comm$InterGenusSpecies))

    if(length(by.species) != 0){
        comm <- vector("list", length=length(by.species))
        for(i in 1:length(by.species)){
            comm[[i]] <- empty.matrix
            this.by.species <- by.species[[i]]
            for(j in 1:nrow(this.by.species)){
                this.row <- this.by.species[j,]
                comm[[i]][match(as.character(this.row[,col.agg]),
                                rownames(comm[[i]])),
                          match(this.row["InterGenusSpecies"],
                                colnames(comm[[i]]))] <-
                    as.numeric(this.row[["Abund"]])
            }
            comm[[i]] <- comm[[i]][rowSums(comm[[i]]) > 0,]
            comm[[i]] <- comm[[i]][,colSums(comm[[i]]) > 0]
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

    if(net.type == "SiteYear"){
        comm.pp <- list(comm=comms,
                        site.date=site.date,
                        sites= sapply(strsplit(names(site.date), ":"),
                                      function(x) x[1]),
                        years= sapply(strsplit(names(site.date), ":"),
                                      function(x) x[2])
                        )
    }
    return(comm.pp)
}
