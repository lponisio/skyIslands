## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                SRDoy + I(SRDoy)^2 +
                                Year + Site
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                  Year + SRDoy + I(SRDoy)^2 +
                                  Site
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(Net_BeeDiversity | weights(Weights)~
                             MeanFloralDiversity + SRDoy + (SRDoy)^2 +
                             Year + Site 
                           )
## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance | weights(Weights)~
                                  MeanFloralAbundance + SRDoy + (SRDoy)^2 +
                                    Year + Site
                                )
## HB abund
formula.HB.abund <- formula(Net_HBAbundance | weights(Weights)~
                              MeanFloralAbundance + SRDoy + (SRDoy)^2 +
                                Year + Site
                            )
## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance | weights(Weights)~
                               MeanFloralAbundance + SRDoy + (SRDoy)^2 +
                                 Year + Site
                             )

## **********************************************************
## convert formulas to brms forma
## **********************************************************
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund, family="student")
bf.bombusabund <- bf(formula.bombus.abund, family="student")
bf.HBabund <- bf(formula.HB.abund, family="student")
bf.bdiv <- bf(formula.bee.div)
