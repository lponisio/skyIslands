## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat +
                                  Year +
                                  SRDoy + I(SRDoy^2) +
                                  (1|Site)
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                  Year+ Lat + 
                                    SRDoy + I(SRDoy^2) +
                                    (1|Site)
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(Net_BeeDiversity | weights(Weights)~
                             MeanFloralDiversity +
                               Lat + Year +
                               SRDoy + I(SRDoy^2) +
                               (1|Site)
                           )
## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance | weights(Weights)~
                                  MeanFloralAbundance + Year +
                                    SRDoy + I(SRDoy^2) +
                                    Lat + 
                                    (1|Site)
                                )
## HB abund
formula.HB.abund <- formula(Net_HBAbundance | weights(Weights)~
                              MeanFloralAbundance +  Year +
                                SRDoy + I(SRDoy^2) +
                                Lat +
                                (1|Site)
                            )
## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance | weights(Weights)~
                               MeanFloralAbundance +  Year +
                                 SRDoy + I(SRDoy^2) +
                                 Lat +
                                 (1|Site)
                             )

## **********************************************************
## convert formulas to brms forma
## **********************************************************
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund, family="gaussian")
bf.bombusabund <- bf(formula.bombus.abund, family="gaussian")
bf.HBabund <- bf(formula.HB.abund, family="gaussian")
bf.bdiv <- bf(formula.bee.div)
