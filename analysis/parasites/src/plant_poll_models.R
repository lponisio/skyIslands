remove_subset_formula <- function(form){
    char.form <- as.character(form)
    no.sub <-  gsub(" \\| subset\\(Weights\\)",  "", char.form[2])
    form.out <- formula(paste(no.sub, "~", char.form[3]))
    return(form.out)
}

if(site.or.lat ==  "lat"){

    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity | subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                          Year + (1|Site) + Lat 
                                  )
    formula.flower.div.nosub <- remove_subset_formula(formula.flower.div)

    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance | subset(Weights) ~
                                        Year +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            (1|Site) + Lat 
                                    )
    formula.flower.abund.nosub <- remove_subset_formula(formula.flower.abund)

    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity | subset(Weights)~
                                   MeanFloralDiversity +
                                       SRDoyPoly1 + SRDoyPoly2 +
                                       Year + (1|Site) + Lat  
                               )
    formula.bee.div.nosub <-
        remove_subset_formula(formula.bee.div)

    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance | subset(Weights)~
                                        MeanFloralAbundance +
                                            MeanFloralDiversity +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Year + (1|Site) + Lat 
                                    )
    formula.bombus.abund.nosub <-
        remove_subset_formula(formula.bombus.abund)

    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance | subset(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Year + (1|Site) + Lat 
                                )
    formula.HB.abund.nosub <-
        remove_subset_formula(formula.HB.abund)

    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance | subset(Weights)~
                                     MeanFloralAbundance +
                                         MeanFloralDiversity +
                                         SRDoyPoly1 + SRDoyPoly2 +
                                         Year + (1|Site) + Lat 
                                 )
    formula.bee.abund.nosub <-
        remove_subset_formula(formula.bee.abund)

} else{
    
    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity | subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                          Year + Site 
                                  )
    formula.flower.div.nosub <- remove_subset_formula(formula.flower.div)

    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance | subset(Weights) ~
                                        Year +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Site
                                    )
    formula.flower.abund.nosub <- remove_subset_formula(formula.flower.abund)

    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity | subset(Weights)~
                                   MeanFloralDiversity +
                                       SRDoyPoly1 + SRDoyPoly2 +
                                       Year + Site 
                               )
    formula.bee.div.nosub <-
        remove_subset_formula(formula.bee.div)

    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance | subset(Weights)~
                                        MeanFloralAbundance +
                                            MeanFloralDiversity +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Year + Site
                                    )
    formula.bombus.abund.nosub <-
        remove_subset_formula(formula.bombus.abund)

    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance | subset(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Year + Site
                                )
    formula.HB.abund.nosub <-
        remove_subset_formula(formula.HB.abund)

    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance | subset(Weights)~
                                     MeanFloralAbundance +
                                         MeanFloralDiversity +
                                         SRDoyPoly1 + SRDoyPoly2 +
                                         Year + Site
                                 )
    formula.bee.abund.nosub <-
        remove_subset_formula(formula.bee.abund)

}

## **********************************************************
## convert formulas to brms forma
## **********************************************************
bf.fabund <- bf(formula.flower.abund, family="student")
bf.fdiv <- bf(formula.flower.div, family="student")
bf.babund <- bf(formula.bee.abund, family="student")
bf.bombusabund <- bf(formula.bombus.abund, family="student")
bf.HBabund <- bf(formula.HB.abund, family="student")
bf.bdiv <- bf(formula.bee.div)


bf.fabund.nosub <- bf(formula.flower.abund.nosub, "student")
bf.fdiv.nosub <- bf(formula.flower.div.nosub, "student")
bf.babund.nosub <- bf(formula.bee.abund.nosub, family="student")
bf.bombusabund.nosub <- bf(formula.bombus.abund.nosub, family="student")
bf.HBabund.nosub <- bf(formula.HB.abund.nosub, family="student")
bf.bdiv.nosub <- bf(formula.bee.div.nosub)
