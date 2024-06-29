
remove_subset_formula <- function(form){
    char.form <- as.character(form)
    no.sub <-
        gsub("\\| subset\\(Weights[:alpha:⁠]*\\)",
             "", char.form[2])
    form.out <- formula(paste(no.sub, "~", char.form[3]))
    return(form.out)
}


if(site.or.lat ==  "lat"){

    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity | subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                          Year + (1|Site) + Lat 
                                  )
    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance | subset(Weights) ~
                                        Year +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            (1|Site) + Lat 
                                    )
    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity | subset(Weights)~
                                   MeanFloralDiversity +
                                       SRDoyPoly1 + SRDoyPoly2 +
                                       Year + (1|Site) + Lat  
                               )
    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance | subset(Weights)~
                                        MeanFloralAbundance +
                                            MeanFloralDiversity +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Year + (1|Site) + Lat 
                                    )
    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance | subset(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Year + (1|Site) + Lat 
                                )
    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance | subset(Weights)~
                                     MeanFloralAbundance +
                                         MeanFloralDiversity +
                                         SRDoyPoly1 + SRDoyPoly2 +
                                         Year + (1|Site) + Lat 
                                 )
} else{
    
    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity | subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                          Year + Site 
                                  )
    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance | subset(Weights) ~
                                        Year +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Site
                                    )
    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity | subset(Weights)~
                                   MeanFloralDiversity +
                                       SRDoyPoly1 + SRDoyPoly2 +
                                       Year + Site 
                               )
    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance | subset(Weights)~
                                        MeanFloralAbundance +
                                            MeanFloralDiversity +
                                            SRDoyPoly1 + SRDoyPoly2 +
                                            Year + Site
                                    )
    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance | subset(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Year + Site
                                )
    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance | subset(Weights)~
                                     MeanFloralAbundance +
                                         MeanFloralDiversity +
                                         SRDoyPoly1 + SRDoyPoly2 +
                                         Year + Site
                                 )

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
