bombus.labs <- rep(names(bombus.morpho), sapply(bombus.morpho, length))

write.csv(bombus.labs, file= "bombus_species.csv", row.names=FALSE)