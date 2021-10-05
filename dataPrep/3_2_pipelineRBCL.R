

setwd('/Volumes/bombus/Dropbox (University of Oregon)')
setwd('skyIslands_saved/SI_pipeline/merged/RBCL/split')

rm(list=ls())
files.taxa <- file.path("f1_taxa_barplots",
                        list.files("f1_taxa_barplots"))

taxa <- lapply(files.taxa, read.csv)

dna.controls <- c("DNActrl1")

ill.controls <- c("RBCLIlluminaCtrl1", "RBCLIlluminaCtrl2",
                  "RBCLIlluminaCtrl3",
                  "RBCLIlluminaCtrl4")

primer.names <- c("forwardbarcode", "revbarcode", "barcodesequence",
                  "forwardgenomicprimer", "revgenomicprimer")

getControls <- function(x, control.names){
    y <- x[x$index %in% control.names, !colnames(x)
      %in% primer.names]
    return(y)
}

dna.control.bact <- lapply(taxa, getControls, dna.controls)

ill.control.bact <- lapply(taxa, getControls, ill.controls)
ill.control.bact[[1]] <- NULL

countControlBact <- function(control.bact){
    rownames(control.bact) <- control.bact$index
    control.bact$index <- NULL
    control.bact.mat <- as.matrix(control.bact)
    bact.count <- colSums(control.bact.mat)
    bacteria.in.controls <- names(bact.count[bact.count > 0])
    return(bacteria.in.controls)
}


dna.bact <- lapply(dna.control.bact, countControlBact)

ill.bact <- lapply(ill.control.bact, countControlBact)


dna.bact.cleaned  <- lapply(dna.bact, function(x){
    gsub("\\.", ";", x)
})
names(dna.bact.cleaned) <- paste0("R", 0)

## names(dna.bact.cleaned) <- paste0("R", 0:4)

ill.bact.cleaned  <- lapply(ill.bact, function(x){
    gsub("\\.", ";", x)
})

names(ill.bact.cleaned) <- paste0("R", 0)
## names(ill.bact.cleaned) <- paste0("R", 1:4)


## uncombined
mapply(function(x, y){
    write.table(x,
                row.names=FALSE,
 file=sprintf("f1_taxa_barplots/DNAcontrolBact_%s.txt",
                             y))
},
x=dna.bact.cleaned,
y=names(dna.bact.cleaned)
)

mapply(function(x, y){
    write.table(x,
                row.names=FALSE,
              file=sprintf("f1_taxa_barplots/ILLcontrolBact_%s.txt",
                           y))
    },
    x=ill.bact.cleaned,
    y=names(ill.bact.cleaned)
)
