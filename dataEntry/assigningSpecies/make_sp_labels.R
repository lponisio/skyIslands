## load list of species labels

sp.labels <- unlist(lapply(sp.ids, function(x){
	rep(paste(x$genus, x$species, x$subspecies), length(x$temp.id))
}))


