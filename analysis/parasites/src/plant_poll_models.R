## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat +
                                  SRDoy*Year + I(SRDoy^2)*Year +
                                  SRDoy*Lat + I(SRDoy^2)*Lat +
                                  (1+SRDoy|Site)
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                  SRDoy*Year + I(SRDoy^2)*Year +
                                    SRDoy*Lat + I(SRDoy^2)*Lat +
                                    (1+SRDoy|Site)
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(Net_BeeDiversity | weights(Weights)~
                             MeanFloralDiversity +
                               Lat + Year +
                               SRDoy*Year + I(SRDoy^2)*Year +
                               SRDoy*Lat + I(SRDoy^2)*Lat +
                               (1+SRDoy|Site)
                           )
## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance | weights(Weights)~
                                  MeanFloralAbundance +
                                    MeanFloralDiversity + Year +
                                    SRDoy*Year + I(SRDoy^2)*Year +
                                    SRDoy*Lat + I(SRDoy^2)*Lat +
                                    Lat + 
                                    (1|Site)
                                )
## HB abund
formula.HB.abund <- formula(Net_HBAbundance | weights(Weights)~
                              MeanFloralAbundance +
                                MeanFloralDiversity+  Year +
                                SRDoy*Year + I(SRDoy^2)*Year +
                                SRDoy*Lat + I(SRDoy^2)*Lat +
                                (Site)
                            )
## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance | weights(Weights)~
                               MeanFloralAbundance +
                                 MeanFloralDiversity +  Year +
                                 SRDoy*Year + I(SRDoy^2)*Year +
                                 SRDoy*Lat + I(SRDoy^2)*Lat +
                                 Lat +
                                 (1+SRDoy|Site)
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
