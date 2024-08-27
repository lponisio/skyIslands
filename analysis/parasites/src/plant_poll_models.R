
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
    formula.flower.div <- formula(MeanFloralDiversity |
                                  subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                      Year + Lat +
                                      Area +
                                      (1|Site) 
                                  )
    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance |
                                    subset(Weights) ~
                                        Year +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        (1|Site)
                                    )
    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity |
                               subset(Weights)~
                                   MeanFloralDiversity +
                                   SRDoyPoly1 + SRDoyPoly2 +
                                   Year +
                                   Lat + 
                                   Area +
                                   (1|Site)  
                               )
    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance |
                                    subset(Weights)~
                                        MeanFloralAbundance +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Year + Lat +
                                        Area + 
                                        (1|Site)
                                    )
    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance |
                                subset(Weights)~
                                    MeanFloralAbundance +
                                    SRDoyPoly1 + SRDoyPoly2 +
                                    Year + Lat +
                                    Area + 
                                    (1|Site)
                                )
    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance |
                                 subset(Weights) ~
                                     MeanFloralAbundance +
                                     MeanFloralDiversity +
                                     SRDoyPoly1 + SRDoyPoly2 +
                                     Year +
                                     Area + 
                                     (1|Site)
                                 )
} else{
    
    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity |
                                  subset(Weights) ~
                                      SRDoyPoly1 + SRDoyPoly2 +
                                      Site 
                                  )
    ## flower abund
    formula.flower.abund <- formula(MeanFloralAbundance |
                                    subset(Weights) ~
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Site
                                    )
    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity | 
                               subset(Weights) ~
                                   MeanFloralDiversity +
                                   SRDoyPoly1 + SRDoyPoly2 +
                                   Site 
                               )
    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance |
                                    subset(Weights)~
                                        MeanFloralAbundance +
                                        MeanFloralDiversity +
                                        SRDoyPoly1 + SRDoyPoly2 +
                                        Site
                                    )
    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance | 
                                subset(Weights) ~
                                    Site
                                )
    ## bee abund
    formula.bee.abund <- formula(Net_BeeAbundance | 
                                 subset(Weights) ~
                                     MeanFloralAbundance +
                                     MeanFloralDiversity +
                                     SRDoyPoly1 + SRDoyPoly2 +
                                     Site
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
