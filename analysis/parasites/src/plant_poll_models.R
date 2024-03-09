## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat +
                                  Year +
                                  SRDoy + I(SRDoy^2) 
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                  Year+ Lat + 
                                    SRDoy + I(SRDoy^2) 
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(Net_BeeDiversity | weights(Weights)~
                             MeanFloralDiversity +
                             Lat  + Year+ SRDoy + I(SRDoy^2) 
                           )
## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance | weights(Weights)~
                                  MeanFloralAbundance + 
                                    SRDoy + I(SRDoy^2) +
                                    Lat
                                )
## HB abund
formula.HB.abund <- formula(Net_HBAbundance | weights(Weights)~
                              MeanFloralAbundance +  
                                SRDoy + I(SRDoy^2)+
                                Lat 
                            )
## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance | weights(Weights)~
                               MeanFloralAbundance +  
                                 SRDoy + I(SRDoy^2)+
                                 Lat 
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
