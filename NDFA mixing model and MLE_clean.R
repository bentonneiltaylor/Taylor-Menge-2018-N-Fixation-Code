setwd("C:/Users/Benton/Desktop/Work Files/Research/Costa Rica/Greenhouse Competition Study/Taylor & Menge_Light Regulates Fixation/Data and Analysis Files")

library(bbmle) # for mle2
library(Hmisc) # for errbar
library(MASS)
source("C:\\RCode\\Homegrown Functions.R")

iso<-read.csv("Plant data for Isotope-enriched pots.csv")
act<-iso[iso$Fixer.Type=='a',]#A dataframe of just the Pentaclethra N-fixers
act<-act[c("Plant.Number","X15Nat","Light.Treatment","Nitrogen.Treatment","nodule.biomass..g.","seed.coat.biomass",
           "belowground.biomass","Env.Trt","dry.seed.mass","Seed.used","Total.biomass","seed.N","Plant.N","prop.N.seed")]
sdl<-read.csv("Seedling Nodule Sampling Data_2017.csv")#read in field-sampled seedling data
#########################################################################################################
#Calculating the 15N Soil Value based on what we added
#########################################################################################################
lit.N.u<-109.7857#mean contribution of N from litterfall for tropical forests from Brookshire 2011 supplement
lit.N.sd<-39.3039#standard deviation of litterfall N from Brookshire 2011 supplement
asym.N<-7.8#contribution of N from asymbiotic fixation taken from Reed et al. 2011
#Getting the soil N mineralization estimate (we'll create a gamma distribution of N mineralization values based on literature estimates)
min.u<-lit.N.u+asym.N
min.sd<-lit.N.sd
iters <- 10000 # iterations
minN<-rnorm(iters,mean=min.u,sd=min.sd)
min.shape <- min.u^2/min.sd^2
min.rate <- min.u/min.sd^2
minN <- rgamma(iters,shape=min.shape,rate=min.rate)

minpot<-((((minN/2)/10)/2)*0.03141593)#gives mineralized N in a pot over the experiment duration (6 months).
#divide by 2 to go from 1yr to 6 month, divide by 10 for kg/ha to g/m2, multiply by .03 to scale to size of pot, divide by 2 to factor in sand in each pot

isoN.pot<-(.266667/45)*.35#Amount of isotopic-enriched fertilizer that went into each pot (.35 converts g NH4NO3 to gN)
#### First for the low N treatment. The only factors are the expected mineralized N and the ######################
#isotopically enriched fertilizer
frac.iso<-isoN.pot/(isoN.pot+minpot)#fraction of nitrogen in the pot that came from the fertilizer
frac.min<-minpot/(isoN.pot+minpot)#fraction of nitrogen in the pot that came from soil N mineralization
lowfert.15Nexp<-data.frame("X15Nsoil"=100*((.003663*frac.min)+(.98*frac.iso)))
#### Now for the middle fertilizer treatment. This is the same as for the low N treatment, ######################
#but has the added N of the non-isotopically enriched fertilizer
mid.fertN.pot<-0.3141593#gN of non-iso fertilizer per pot over entire experiment duration
frac.iso.mid<-isoN.pot/(isoN.pot+minpot+mid.fertN.pot)
frac.min.mid<-minpot/(isoN.pot+minpot+mid.fertN.pot)
frac.fert.mid<-mid.fertN.pot/(isoN.pot+minpot+mid.fertN.pot)
midfert.15Nexp<-data.frame("X15Nsoil"=100*((.003663*frac.min.mid)+(.98*frac.iso.mid)+(.003663*frac.fert.mid)))#0.008792522 (0.879% 15N)
##### Now for the High fertilizer treatment. This is the same as for the low N treatment, ######################
#but has the added N of the non-isotopically enriched fertilizer
hi.fertN.pot<-0.6283185#gN of non-iso fertilizer per pot over entire experiment duration
frac.iso.hi<-isoN.pot/(isoN.pot+minpot+hi.fertN.pot)
frac.min.hi<-minpot/(isoN.pot+minpot+hi.fertN.pot)
frac.fert.hi<-hi.fertN.pot/(isoN.pot+minpot+hi.fertN.pot)
hifert.15Nexp<-data.frame("X15Nsoil"=100*((.003663*frac.min.hi)+(.98*frac.iso.hi)+(.003663*frac.fert.hi)))#0.006519405 (0.65194% 15N)

#######################################################################################################
###### Calculating "Used Seed Mass" and "Proportion of Plant N Derived from Seed"#####
#Factoring in decay for seed coat numbers
act$Seed.used<-act$dry.seed.mass-act$seed.coat.biomass
act$seed.N.2<-act$Seed.used*.03359
act$prop.N.seed<-with(act, (seed.N.2/Plant.N))
act$prop.N.seed<-ifelse(act$prop.N.seed>1, .9, no=act$prop.N.seed)

###### Storing mean and sd of prop.N.seed for active fixers in each light treatment ####### (for use in Ndfa calcs below)
lo.s.mn<-mean(with(act[act$Light.Treatment=="L",], prop.N.seed), na.rm=T)
lo.s.sd<-max(.05,sd(with(act[act$Light.Treatment=="L",], prop.N.seed), na.rm=T))
md.s.mn<-mean(with(act[act$Light.Treatment=="M",], prop.N.seed), na.rm=T)
md.s.sd<-max(.05,sd(with(act[act$Light.Treatment=="M",], prop.N.seed), na.rm=T))
hi.s.mn<-mean(with(act[act$Light.Treatment=="H",], prop.N.seed), na.rm=T)
hi.s.sd<-sd(with(act[act$Light.Treatment=="H",], prop.N.seed), na.rm=T)

lo.a<-(lo.s.mn/lo.s.sd)^2
lo.b<-(lo.s.sd^2)/lo.s.mn
md.a<-(md.s.mn/md.s.sd)^2
md.b<-(md.s.sd^2)/md.s.mn
hi.a<-(hi.s.mn/hi.s.sd)^2
hi.b<-(hi.s.sd^2)/hi.s.mn
#########################################################################################################
### Adding in the SD of N15% for each light treatment
lo.15n.sd<-sd(with(act[act$Light.Treatment=="L",], X15Nat), na.rm=T)
md.15n.sd<-sd(with(act[act$Light.Treatment=="M",], X15Nat), na.rm=T)
hi.15n.sd<-sd(with(act[act$Light.Treatment=="H",], X15Nat), na.rm=T)

act$X15Nat_sd<-ifelse(act$Light.Treatment=="L", yes=lo.15n.sd, no=md.15n.sd)
act$X15Nat_sd<-ifelse(act$Light.Treatment=="H", yes=hi.15n.sd, no=act$X15Nat_sd)
#########################################################################################################
#########################################################################################################
#Calculating Ndfa
#########################################################################################################
act2<-act[!is.na(act$prop.N.seed),]
Plant.Number<-unique(act2$Plant.Number)
Ndfa<-NULL
pass<-NULL
for (i in 1:length(Plant.Number)){
  print(i)
  temp<-act2[act2$Plant.Number==Plant.Number[i],]
  pass<-as.data.frame(ifelse(temp$Nitrogen.Treatment=="L", lowfert.15Nexp, ifelse(temp$Nitrogen.Treatment=="M", midfert.15Nexp, hifert.15Nexp)))
  colnames(pass)<-"X15Nsoil"
  pass$X15Nat_fix<-rnorm(n=iters, mean=0.3663, sd=.05)
  pass$prop.N.seed<-sample(seq(.01,.99,length.out=iters),size=iters,replace=T,
                           prob=dnorm(seq(.1,.99,length.out=iters),
                                      mean=temp$prop.N.seed, 
                                      sd=ifelse(temp$Light.Treatment=="L",lo.s.sd,ifelse(temp$Light.Treatment=="M",md.s.sd,hi.s.sd))))
  temp$X15Nsoil<-mean(pass$X15Nsoil)
  temp$X15Nsoil_sd<-sd(pass$X15Nsoil)
  
  Ndfadist <- 100*((temp$X15Nat-(temp$prop.N.seed*pass$X15Nat_fix)-((1-temp$prop.N.seed)*pass$X15Nsoil))/
                     ((pass$X15Nat_fix-pass$X15Nsoil)*(1-temp$prop.N.seed)))
  Ndfadist <- 100*((temp$X15Nat-(pass$prop.N.seed*pass$X15Nat_fix)-((1-pass$prop.N.seed)*pass$X15Nsoil))/
                     ((pass$X15Nat_fix-pass$X15Nsoil)))
  # The previous line defines Ndfa as the percent of total N, not newly acquired N
  temp$Ndfa <- mean(Ndfadist)
  temp$Ndfa_sd <- sd (Ndfadist)
  temp$Ndfa_LCI <- quantile(Ndfadist,0.025)
  temp$Ndfa_UCI <- quantile(Ndfadist,0.975)
  Ndfa<-rbind(Ndfa, temp)
}

# Now calculating Nfixd as the amount of total plant N that comes from fixation
Ndfa$Nfixd<-with(Ndfa, Plant.N*Ndfa/100)
Ndfa$Nfixd_LCI<-with(Ndfa, Plant.N*Ndfa_LCI/100)
Ndfa$Nfixd_UCI<-with(Ndfa, Plant.N*Ndfa_UCI/100)

Ndfa$nod.prop.tot<-with(Ndfa, (nodule.biomass..g./Total.biomass)*100)
Ndfa$nod.prop.bg<-with(Ndfa, (nodule.biomass..g./belowground.biomass)*100)

#write.csv(Ndfa, file="Ndfa mixing model output.csv", row.names=F)
with(Ndfa, aggregate(Ndfa, by=list(Light.Treatment, Nitrogen.Treatment),FUN=mean))

summary(Ndfa$Ndfa)
with(Ndfa, aggregate(Nfixd, by=list(Light.Treatment, Nitrogen.Treatment),FUN=mean))
ndfa<-Ndfa

##############################################################################################################################
##############################################################################################################################
###############################Calculating N in each pot in mg N per kg soil ##################################################
##############################################################################################################################
#Estimating N in each pot's sand/soil mix. We'll take the average N content of the forest soil samples
#as the amount of N from soil, but divide it by half because half the pot was sand
forestN<-mean(sdl$tot.n)#mean concentration from forest soil in mg N per kg soil
potforN<-forestN*(9.8/2)#gives the amount of N in mg for the pot mixture before N additions (9.8 is kg's of potting mix which is halved for the sand fraction)
lowpotN<-(potforN+(isoN.pot*1000))/9.8

midadd<-mid.fertN.pot+isoN.pot
midpotN<-(potforN+(midadd*1000))/9.8

hiadd<-hi.fertN.pot+isoN.pot
hipotN<-(potforN+(hiadd*1000))/9.8
##############################################################################################################################
#########################################################################################################################
##############################################################################################################################
### Now the Maximum Likelihood Functions to put Statistics to these Ndfa numbers #########################3

#We'll run MLE for %Ndfa first

# Adding columns of 0s and 1s for each treatment; multiplying light by nitrogen
# treatments gives unique treatments, and makes it easier to code things in by vector
ndfa$lL <- ndfa$mL <- ndfa$hL <- ndfa$lN <- ndfa$mN <- ndfa$hN <- 0
ndfa[ndfa$Light.Treatment=="L",]$lL <- 1
ndfa[ndfa$Light.Treatment=="M",]$mL <- 1
ndfa[ndfa$Light.Treatment=="H",]$hL <- 1
ndfa[ndfa$Nitrogen.Treatment=="L",]$lN <- 1
ndfa[ndfa$Nitrogen.Treatment=="M",]$mN <- 1
ndfa[ndfa$Nitrogen.Treatment=="H",]$hN <- 1

# All the same (formerly null)

normNLL_Ndfa_aaaaaaaaa <- function(X15N_samp_sd,Ndfa_u,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd){
  
  Ndfa_u_vec <- Ndfa_u
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.aaaaaaaaa <- mle2(normNLL_Ndfa_aaaaaaaaa,start=list(X15N_samp_sd=0.1,Ndfa_u=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005))
summary(ndfafit.aaaaaaaaa)

# All different (formerly nitlit)

normNLL_Ndfa_abcdefghi <- function(X15N_samp_sd,Ndfa_u_LL,Ndfa_u_LM,Ndfa_u_LH,Ndfa_u_ML,
                                   Ndfa_u_MM,Ndfa_u_MH,Ndfa_u_HL,Ndfa_u_HM,Ndfa_u_HH,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_LL*lL*lN + Ndfa_u_LM*lL*mN + Ndfa_u_LH*lL*hN + 
    Ndfa_u_ML*mL*lN + Ndfa_u_MM*mL*mN + Ndfa_u_MH*mL*hN + 
    Ndfa_u_HL*hL*lN + Ndfa_u_HM*hL*mN + Ndfa_u_HH*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.abcdefghi <- mle2(normNLL_Ndfa_abcdefghi,start=list(X15N_samp_sd=0.1,Ndfa_u_LL=50,Ndfa_u_LM=50,Ndfa_u_LH=50,
                                                            Ndfa_u_ML=50,Ndfa_u_MM=50,Ndfa_u_MH=50,Ndfa_u_HL=50,Ndfa_u_HM=50,Ndfa_u_HH=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.abcdefghi)

# light treatments are different: aaa bbb ccc

normNLL_Ndfa_aaabbbccc <- function(X15N_samp_sd,Ndfa_u_a,Ndfa_u_b,Ndfa_u_c,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_a*lL*lN + Ndfa_u_a*lL*mN + Ndfa_u_a*lL*hN + 
    Ndfa_u_b*mL*lN + Ndfa_u_b*mL*mN + Ndfa_u_b*mL*hN + 
    Ndfa_u_c*hL*lN + Ndfa_u_c*hL*mN + Ndfa_u_c*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.aaabbbccc <- mle2(normNLL_Ndfa_aaabbbccc,start=list(X15N_samp_sd=0.1,Ndfa_u_a=50,Ndfa_u_b=50, Ndfa_u_c=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.aaabbbccc)

# Nitrogen treatments are different: abc abc abc

normNLL_Ndfa_abcabcabc <- function(X15N_samp_sd,Ndfa_u_a,Ndfa_u_b,Ndfa_u_c,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_a*lL*lN + Ndfa_u_b*lL*mN + Ndfa_u_c*lL*hN + 
    Ndfa_u_a*mL*lN + Ndfa_u_b*mL*mN + Ndfa_u_c*mL*hN + 
    Ndfa_u_a*hL*lN + Ndfa_u_b*hL*mN + Ndfa_u_c*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.abcabcabc <- mle2(normNLL_Ndfa_abcabcabc,start=list(X15N_samp_sd=0.1,Ndfa_u_a=50,Ndfa_u_b=50, Ndfa_u_c=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.abcabcabc)


# High light different: aaa aaa bbb

normNLL_Ndfa_aaaaaabbb <- function(X15N_samp_sd,Ndfa_u_a,Ndfa_u_b,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_a*lL*lN + Ndfa_u_a*lL*mN + Ndfa_u_a*lL*hN + 
    Ndfa_u_a*mL*lN + Ndfa_u_a*mL*mN + Ndfa_u_a*mL*hN + 
    Ndfa_u_b*hL*lN + Ndfa_u_b*hL*mN + Ndfa_u_b*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.aaaaaabbb <- mle2(normNLL_Ndfa_aaaaaabbb,start=list(X15N_samp_sd=0.1,Ndfa_u_a=50,Ndfa_u_b=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.aaaaaabbb)

# aaa aaa bcc

normNLL_Ndfa_aaaaaabcc <- function(X15N_samp_sd,Ndfa_u_a,Ndfa_u_b,Ndfa_u_c,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_a*lL*lN + Ndfa_u_a*lL*mN + Ndfa_u_a*lL*hN + 
    Ndfa_u_a*mL*lN + Ndfa_u_a*mL*mN + Ndfa_u_a*mL*hN + 
    Ndfa_u_b*hL*lN + Ndfa_u_c*hL*mN + Ndfa_u_c*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.aaaaaabcc <- mle2(normNLL_Ndfa_aaaaaabcc,start=list(X15N_samp_sd=0.1,Ndfa_u_a=50,Ndfa_u_b=50,Ndfa_u_c=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.aaaaaabcc)

# aaa aaa bcd

normNLL_Ndfa_aaaaaabcd <- function(X15N_samp_sd,Ndfa_u_a,Ndfa_u_b,Ndfa_u_c,Ndfa_u_d,
                                   X15N_dat,X15N_ref_u,X15N_ref_sd,X15N_fix_u,N.seed,X15N_fix_sd,
                                   lL,mL,hL,lN,mN,hN){
  
  Ndfa_u_vec <- Ndfa_u_a*lL*lN + Ndfa_u_a*lL*mN + Ndfa_u_a*lL*hN + 
    Ndfa_u_a*mL*lN + Ndfa_u_a*mL*mN + Ndfa_u_a*mL*hN + 
    Ndfa_u_b*hL*lN + Ndfa_u_c*hL*mN + Ndfa_u_d*hL*hN
  X15N_u_vec <- (Ndfa_u_vec/100 + N.seed)*X15N_fix_u + (1-Ndfa_u_vec/100-N.seed)*X15N_ref_u
  sdtot <- sqrt((max((1-Ndfa_u_vec/100-N.seed),0)*X15N_ref_sd)^2 + (Ndfa_u_vec/100*X15N_fix_sd)^2 + 
                  N.seed*X15N_fix_sd^2 + X15N_samp_sd^2)
  -sum(dnorm(X15N_dat,mean=X15N_u_vec,sd=sdtot,log=TRUE))
}

ndfafit.aaaaaabcd <- mle2(normNLL_Ndfa_aaaaaabcd,start=list(X15N_samp_sd=0.1,Ndfa_u_a=50,Ndfa_u_b=50,Ndfa_u_c=50,Ndfa_u_d=50),
                          data=list(X15N_dat=ndfa$X15Nat,X15N_ref_u=ndfa$X15Nsoil,X15N_ref_sd=ndfa$X15Nsoil_sd,
                                    X15N_fix_u=0.3663,N.seed=ndfa$prop.N.seed, X15N_fix_sd=.005,
                                    lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,mN=ndfa$mN,hN=ndfa$hN))
summary(ndfafit.aaaaaabcd)

# Compare models

deltaAICs <- function(AICvec){
  best <- 0
  for(i in 1:length(AICvec)){
    print(AICvec[i]-min(AICvec))
    if(AICvec[i]==min(AICvec)) best <- i
  }
  best
}

deltaAICs(c(AIC(ndfafit.aaaaaaaaa),AIC(ndfafit.aaabbbccc),AIC(ndfafit.abcabcabc),AIC(ndfafit.aaaaaabbb),AIC(ndfafit.aaaaaabcc),AIC(ndfafit.aaaaaabcd),AIC(ndfafit.abcdefghi)))
# Looks like aaa aaa bcc is best for the high light. The abcdefghi says that none of the other ones are significantly > 0 (see the Pr(z)
# in the last column for those). We're treating the negative means as 0s for display and biological interpretation

##################################################################################################################
##### Now MLE Models for Total Plant Biomass (g)#############################################################################################################
##################################################################################################################

# All the same 

normNLL_bio_aaaaaaaaa <- function(bio_u, bio_sd){
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}

biofit.aaaaaaaaa <- mle2(normNLL_bio_aaaaaaaaa,start=list(bio_u=5, bio_sd=.05),
                          data=list(biomass=ndfa$Total.biomass))
summary(biofit.aaaaaaaaa)

# All different
normNLL_bio_abcdefghi <- function(bio_u_LL, bio_u_LM, bio_u_LH,
                                  bio_u_ML, bio_u_MM, bio_u_MH,
                                  bio_u_HL, bio_u_HM, bio_u_HH, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_LL*lL*lN + bio_u_LM*lL*mN + bio_u_LH*lL*hN + 
  bio_u_ML*mL*lN + bio_u_MM*mL*mN + bio_u_MH*mL*hN + 
  bio_u_HL*hL*lN + bio_u_HM*hL*mN + bio_u_HH*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.abcdefghi <- mle2(normNLL_bio_abcdefghi,start=list(bio_u_LL=5, bio_u_LM=5, bio_u_LH=5,
                                                          bio_u_ML=5, bio_u_MM=5, bio_u_MH=5,
                                                          bio_u_HL=5, bio_u_HM=5, bio_u_HH=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.abcdefghi)

#All light treatments different aaa bbb ccc

normNLL_bio_aaabbbccc <- function(bio_u_a, bio_u_b, bio_u_c, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_b*mL*lN + bio_u_b*mL*mN + bio_u_b*mL*hN + 
    bio_u_c*hL*lN + bio_u_c*hL*mN + bio_u_c*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaabbbccc <- mle2(normNLL_bio_aaabbbccc,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaabbbccc)

#All Nitrogen treatments different abc abc abc

normNLL_bio_abcabcabc <- function(bio_u_a, bio_u_b, bio_u_c, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_b*lL*mN + bio_u_c*lL*hN + 
    bio_u_a*mL*lN + bio_u_b*mL*mN + bio_u_c*mL*hN + 
    bio_u_a*hL*lN + bio_u_b*hL*mN + bio_u_c*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.abcabcabc <- mle2(normNLL_bio_abcabcabc,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.abcabcabc)

#High light different aaa aaa bbb

normNLL_bio_aaaaaabbb <- function(bio_u_a, bio_u_b,  
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_a*mL*lN + bio_u_a*mL*mN + bio_u_a*mL*hN + 
    bio_u_b*hL*lN + bio_u_b*hL*mN + bio_u_b*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaaaaabbb <- mle2(normNLL_bio_aaaaaabbb,start=list(bio_u_a=5, bio_u_b=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaaaaabbb)

# BCC in High light different aaa aaa bcc

normNLL_bio_aaaaaabcc <- function(bio_u_a, bio_u_b, bio_u_c, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_a*mL*lN + bio_u_a*mL*mN + bio_u_a*mL*hN + 
    bio_u_b*hL*lN + bio_u_c*hL*mN + bio_u_c*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaaaaabcc <- mle2(normNLL_bio_aaaaaabcc,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaaaaabcc)

# Medium and High light different aaa bbb ccc

normNLL_bio_aaabbbccc <- function(bio_u_a, bio_u_b, bio_u_c, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_b*mL*lN + bio_u_b*mL*mN + bio_u_b*mL*hN + 
    bio_u_c*hL*lN + bio_u_c*hL*mN + bio_u_c*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaabbbccc <- mle2(normNLL_bio_aaabbbccc,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaabbbccc)

# Medium and High light different with N different in high-light aaa bbb cdd

normNLL_bio_aaabbbcdd <- function(bio_u_a, bio_u_b, bio_u_c, bio_u_d,
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_b*mL*lN + bio_u_b*mL*mN + bio_u_b*mL*hN + 
    bio_u_c*hL*lN + bio_u_d*hL*mN + bio_u_d*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaabbbcdd <- mle2(normNLL_bio_aaabbbcdd,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5,bio_u_d=5, bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaabbbcdd)

# BCD in High light different aaa aaa bcd

normNLL_bio_aaaaaabcd <- function(bio_u_a, bio_u_b, bio_u_c, bio_u_d, 
                                  bio_sd,lL,mL,hL,lN,mN,hN){
  bio_u <- bio_u_a*lL*lN + bio_u_a*lL*mN + bio_u_a*lL*hN + 
    bio_u_a*mL*lN + bio_u_a*mL*mN + bio_u_a*mL*hN + 
    bio_u_b*hL*lN + bio_u_c*hL*mN + bio_u_d*hL*hN
  -sum(dnorm(biomass,mean=bio_u,sd=bio_sd,log=TRUE))
}
biofit.aaaaaabcd <- mle2(normNLL_bio_aaaaaabcd,start=list(bio_u_a=5, bio_u_b=5, bio_u_c=5, bio_u_d=5,bio_sd=.05),
                         data=list(biomass=ndfa$Total.biomass,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(biofit.aaaaaabcd)

deltaAICs(c(AIC(biofit.aaaaaaaaa),AIC(biofit.aaabbbccc),AIC(biofit.abcabcabc),AIC(biofit.aaaaaabbb),AIC(biofit.aaaaaabcc),AIC(biofit.aaaaaabcd),
            AIC(biofit.aaabbbccc),AIC(biofit.aaabbbcdd),AIC(biofit.abcdefghi)))
#Looks like, similar to the Ndfa models, the aaa aaa bcc model is the best fit by quite a bit

##################################################################################################################
##### Now MLE Models for Nodule Biomass (percent of belowground biomass allocated to nodules)#############################################################################################################
##################################################################################################################

# All the same 

normNLL_nod_aaaaaaaaa <- function(nod_u, nod_sd){
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}

nodfit.aaaaaaaaa <- mle2(normNLL_nod_aaaaaaaaa,start=list(nod_u=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.))
summary(nodfit.aaaaaaaaa)

# All different
normNLL_nod_abcdefghi <- function(nod_u_LL, nod_u_LM, nod_u_LH,
                                  nod_u_ML, nod_u_MM, nod_u_MH,
                                  nod_u_HL, nod_u_HM, nod_u_HH, 
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_LL*lL*lN + nod_u_LM*lL*mN + nod_u_LH*lL*hN + 
    nod_u_ML*mL*lN + nod_u_MM*mL*mN + nod_u_MH*mL*hN + 
    nod_u_HL*hL*lN + nod_u_HM*hL*mN + nod_u_HH*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.abcdefghi <- mle2(normNLL_nod_abcdefghi,start=list(nod_u_LL=1, nod_u_LM=1, nod_u_LH=1,
                                                          nod_u_ML=1, nod_u_MM=1, nod_u_MH=1,
                                                          nod_u_HL=1, nod_u_HM=1, nod_u_HH=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.abcdefghi)

#light treatments different aaa bbb ccc

normNLL_nod_aaabbbccc <- function(nod_u_a, nod_u_b, nod_u_c, 
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_b*mL*lN + nod_u_b*mL*mN + nod_u_b*mL*hN + 
    nod_u_c*hL*lN + nod_u_c*hL*mN + nod_u_c*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaabbbccc <- mle2(normNLL_nod_aaabbbccc,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaabbbccc)

#Nitrogen treatments different abc abc abc

normNLL_nod_abcabcabc <- function(nod_u_a, nod_u_b, nod_u_c, 
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_b*lL*mN + nod_u_c*lL*hN + 
    nod_u_a*mL*lN + nod_u_b*mL*mN + nod_u_c*mL*hN + 
    nod_u_a*hL*lN + nod_u_b*hL*mN + nod_u_c*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.abcabcabc <- mle2(normNLL_nod_abcabcabc,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.abcabcabc)

#High light different aaa aaa bbb

normNLL_nod_aaaaaabbb <- function(nod_u_a, nod_u_b,  
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_a*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_b*hL*lN + nod_u_b*hL*mN + nod_u_b*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaaaaabbb <- mle2(normNLL_nod_aaaaaabbb,start=list(nod_u_a=1, nod_u_b=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaaaaabbb)

# BCC in High light different aaa aaa bcc

normNLL_nod_aaaaaabcc <- function(nod_u_a, nod_u_b, nod_u_c, 
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_a*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_b*hL*lN + nod_u_c*hL*mN + nod_u_c*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaaaaabcc <- mle2(normNLL_nod_aaaaaabcc,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaaaaabcc)

# BCD in High light different aaa aaa bcd

normNLL_nod_aaaaaabcd <- function(nod_u_a, nod_u_b, nod_u_c, nod_u_d,
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_a*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_b*hL*lN + nod_u_c*hL*mN + nod_u_d*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaaaaabcd <- mle2(normNLL_nod_aaaaaabcd,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_u_d=1,nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaaaaabcd)

# BCA in High light different aaa aaa bca

normNLL_nod_aaaaaabca <- function(nod_u_a, nod_u_b, nod_u_c, 
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_a*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_b*hL*lN + nod_u_c*hL*mN + nod_u_a*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaaaaabca <- mle2(normNLL_nod_aaaaaabca,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaaaaabca)

# Med light, low N same as High light, high N: aaa baa cdb

normNLL_nod_aaabaacdb <- function(nod_u_a, nod_u_b, nod_u_c, nod_u_d,
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_b*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_c*hL*lN + nod_u_d*hL*mN + nod_u_b*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaabaacdb <- mle2(normNLL_nod_aaabaacdb,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_u_d=1,nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaabaacdb)

# Med light, low N and all High light different: aaa baa cde

normNLL_nod_aaabaacde <- function(nod_u_a, nod_u_b, nod_u_c, nod_u_d, nod_u_e,
                                  nod_sd,lL,mL,hL,lN,mN,hN){
  nod_u <- nod_u_a*lL*lN + nod_u_a*lL*mN + nod_u_a*lL*hN + 
    nod_u_b*mL*lN + nod_u_a*mL*mN + nod_u_a*mL*hN + 
    nod_u_c*hL*lN + nod_u_d*hL*mN + nod_u_e*hL*hN
  -sum(dnorm(nodule,mean=nod_u,sd=nod_sd,log=TRUE))
}
nodfit.aaabaacde <- mle2(normNLL_nod_aaabaacde,start=list(nod_u_a=1, nod_u_b=1, nod_u_c=1, nod_u_d=1,nod_u_e=1,nod_sd=.05),
                         data=list(nodule=ndfa$nodule.biomass..g.,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(nodfit.aaabaacde)

deltaAICs(c(AIC(nodfit.aaaaaaaaa),AIC(nodfit.aaabbbccc),AIC(nodfit.abcabcabc),AIC(nodfit.aaaaaabbb),AIC(nodfit.aaaaaabcc),AIC(nodfit.aaaaaabcd),
            AIC(nodfit.abcdefghi),AIC(nodfit.aaaaaabca),AIC(nodfit.aaabaacdb),AIC(nodfit.aaabaacde)))
#Looks like the aaa aaa bcd model is the best fit by quite a bit

##################################################################################################################
##### Now MLE Models for Nfixd (g of N fixed per plant)#############################################################################################################
##################################################################################################################

# All the same 

normNLL_fixd_aaaaaaaaa <- function(fixd_u, fixd_sd){
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}

fixdfit.aaaaaaaaa <- mle2(normNLL_fixd_aaaaaaaaa,start=list(fixd_u=1, fixd_sd=.05),
                         data=list(fixd=ndfa$Nfixd))
summary(fixdfit.aaaaaaaaa)

# All different
normNLL_fixd_abcdefghi <- function(fixd_u_LL, fixd_u_LM, fixd_u_LH,
                                   fixd_u_ML, fixd_u_MM, fixd_u_MH,
                                   fixd_u_HL, fixd_u_HM, fixd_u_HH, 
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_LL*lL*lN + fixd_u_LM*lL*mN + fixd_u_LH*lL*hN + 
    fixd_u_ML*mL*lN + fixd_u_MM*mL*mN + fixd_u_MH*mL*hN + 
    fixd_u_HL*hL*lN + fixd_u_HM*hL*mN + fixd_u_HH*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.abcdefghi <- mle2(normNLL_fixd_abcdefghi,start=list(fixd_u_LL=1, fixd_u_LM=1, fixd_u_LH=1,
                                                            fixd_u_ML=1, fixd_u_MM=1, fixd_u_MH=1,
                                                            fixd_u_HL=1, fixd_u_HM=1, fixd_u_HH=1, fixd_sd=.05),
                         data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.abcdefghi)

#light treatments different aaa bbb ccc

normNLL_fixd_aaabbbccc <- function(fixd_u_a, fixd_u_b, fixd_u_c, 
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_a*lL*lN + fixd_u_a*lL*mN + fixd_u_a*lL*hN + 
    fixd_u_b*mL*lN + fixd_u_b*mL*mN + fixd_u_b*mL*hN + 
    fixd_u_c*hL*lN + fixd_u_c*hL*mN + fixd_u_c*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.aaabbbccc <- mle2(normNLL_fixd_aaabbbccc,start=list(fixd_u_a=1, fixd_u_b=1, fixd_u_c=1, fixd_sd=.05),
                          data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                    mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.aaabbbccc)

#Nitrogen treatments different abc abc abc

normNLL_fixd_abcabcabc <- function(fixd_u_a, fixd_u_b, fixd_u_c, 
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_a*lL*lN + fixd_u_b*lL*mN + fixd_u_c*lL*hN + 
    fixd_u_a*mL*lN + fixd_u_b*mL*mN + fixd_u_c*mL*hN + 
    fixd_u_a*hL*lN + fixd_u_b*hL*mN + fixd_u_c*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.abcabcabc <- mle2(normNLL_fixd_abcabcabc,start=list(fixd_u_a=1, fixd_u_b=1, fixd_u_c=1, fixd_sd=.05),
                          data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                    mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.abcabcabc)

#High light different aaa aaa bbb

normNLL_fixd_aaaaaabbb <- function(fixd_u_a, fixd_u_b,  
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_a*lL*lN + fixd_u_a*lL*mN + fixd_u_a*lL*hN + 
    fixd_u_a*mL*lN + fixd_u_a*mL*mN + fixd_u_a*mL*hN + 
    fixd_u_b*hL*lN + fixd_u_b*hL*mN + fixd_u_b*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.aaaaaabbb <- mle2(normNLL_fixd_aaaaaabbb,start=list(fixd_u_a=1, fixd_u_b=1, fixd_sd=.05),
                         data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.aaaaaabbb)

# BCC in High light different aaa aaa bcc

normNLL_fixd_aaaaaabcc <- function(fixd_u_a, fixd_u_b, fixd_u_c, 
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_a*lL*lN + fixd_u_a*lL*mN + fixd_u_a*lL*hN + 
    fixd_u_a*mL*lN + fixd_u_a*mL*mN + fixd_u_a*mL*hN + 
    fixd_u_b*hL*lN + fixd_u_c*hL*mN + fixd_u_c*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.aaaaaabcc <- mle2(normNLL_fixd_aaaaaabcc,start=list(fixd_u_a=1, fixd_u_b=1, fixd_u_c=1, fixd_sd=.05),
                         data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.aaaaaabcc)

# BCD in High light different aaa aaa bcd

normNLL_fixd_aaaaaabcd <- function(fixd_u_a, fixd_u_b, fixd_u_c, fixd_u_d,
                                   fixd_sd,lL,mL,hL,lN,mN,hN){
  fixd_u <- fixd_u_a*lL*lN + fixd_u_a*lL*mN + fixd_u_a*lL*hN + 
    fixd_u_a*mL*lN + fixd_u_a*mL*mN + fixd_u_a*mL*hN + 
    fixd_u_b*hL*lN + fixd_u_c*hL*mN + fixd_u_d*hL*hN
  -sum(dnorm(fixd,mean=fixd_u,sd=fixd_sd,log=TRUE))
}
fixdfit.aaaaaabcd <- mle2(normNLL_fixd_aaaaaabcd,start=list(fixd_u_a=1, fixd_u_b=1, fixd_u_c=1, fixd_u_d=1,fixd_sd=.05),
                         data=list(fixd=ndfa$Nfixd,lL=ndfa$lL,mL=ndfa$mL,hL=ndfa$hL,lN=ndfa$lN,
                                   mN=ndfa$mN,hN=ndfa$hN))
summary(fixdfit.aaaaaabcd)

deltaAICs(c(AIC(fixdfit.aaaaaaaaa),AIC(fixdfit.aaabbbccc),AIC(fixdfit.abcabcabc),AIC(fixdfit.aaaaaabbb),
            AIC(fixdfit.aaaaaabcc),AIC(fixdfit.aaaaaabcd), AIC(fixdfit.abcdefghi)))

##################################################################################################################
##################################################################################################################
##### Now MLE Models FIELD-SAMPLED SEEDLINGS #############################################################################################################
##################################################################################################################
##################################################################################################################
sdl<-read.csv("Seedling Nodule Sampling Data_2017.csv")#read in field-sampled seedling data
nod<-sdl$prop.nod.bio
tot.n<-sdl$tot.n
nh<-sdl$ammonium
no<-sdl$nitrate
lit<-sdl$X..Trans.Tot
rt<-sdl$root.biomass..g.

p <- length(nod[nod==0])/length(nod)
u<-mean(log(nod[nod>0])) # log space mean of nodule biomass (in mg) if there are nodules
s <- sd(log(nod[nod>0])) # log space sd of nodule biomass (in mg) if there are nodules
knownmedian <- (1-p)*exp(u)
#First, plotting the data just to see what they look like
plot(jitter(tot.n),nod,xlab="Total Nitrogen (mg/g)",ylab="Nodulation (%)")
#Writing some functions to calculate and compare AICc's
AICc <- function(fitobj,dat){
  K <- length(coef(fitobj))+1
  n <- length(dat)#make this Nrow of data
  AICc <- AIC(fitobj) + 2*K*(K+1)/(n-K-1)
}

deltaAICs <- function(AICvec){
  for(i in 1:length(AICvec)){
    print(AICvec[i]-min(AICvec))
  }
}

# Starting values
pstart <- 0.5
ustart <- log(5)
sstart <- log(2)

############# DOES NODULATION VARY WITH TOTAL INORGANIC N??? #################
##############################################################################
### Running a Null Model fitting a flat mean to the data #####################
##############################################################################
ZinflnormNLL_nod_null <- function(sdlognod,ulognod,p_no_nod,nod_dat){
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}
fit_nod_null <- mle2(ZinflnormNLL_nod_null,
                     start=list(sdlognod=sstart,ulognod=ustart,p_no_nod=pstart),
                     data=list(nod_dat=nod))
print(summary(fit_nod_null))
predmedian<-(1-coef(fit_nod_null)[3])*exp(coef(fit_nod_null)[2])#for plotting
##################################################################
### Now running a linear model for total soil N #####################
##################################################################
ZinflnormNLL_nod_nit <- function(sdlognod,c,b,e,d,nod_dat,nit){
  ulognod<-c+b*nit
  p_no_nod<-alogitfn(e+d*nit)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_nit <- mle2(ZinflnormNLL_nod_nit,
                    start=list(sdlognod=sstart,c=10,b=-.1,e=.8,d=-.0001),
                    data=list(nod_dat=nod,nit=tot.n))
print(summary(fit_nod_nit))
#################################################################################################
### Now running a linear model for stand age where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_nit <- function(sdlognod,c,b,p_no_nod,nod_dat,nit){
  ulognod<-c+b*nit
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_nit <- mle2(ZinflnormNLL_nod_u_nit,
                      start=list(sdlognod=sstart,c=10,b=-.1,p_no_nod=.5),
                      data=list(nod_dat=nod,nit=tot.n))
print(summary(fit_nod_u_nit))
#################################################################################################
### Now running a linear model for stand age where only the probability of 0 fluctuates#############
#################################################################################################

ZinflnormNLL_nod_p_nit <- function(ulognod,sdlognod,e,d,nod_dat,nit){
  p_no_nod<-alogitfn(e+d*nit)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}

fit_nod_p_nit <- mle2(ZinflnormNLL_nod_p_nit,
                      start=list(ulognod=ustart,sdlognod=sstart,e=.8,d=.01),
                      data=list(nod_dat=nod,nit=tot.n))
print(summary(fit_nod_p_nit))
#################################################################################################
### Now running a hump-shaped model for stand age where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_quadnit <- function(sdlognod,c,b,b2,p_no_nod,nod_dat,nit){
  ulognod<-c+b*nit+b2*nit^2
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_quadnit <- mle2(ZinflnormNLL_nod_u_quadnit,
                          start=list(sdlognod=sstart,c=10,b=-.1,b2=-.1,p_no_nod=.5),
                          data=list(nod_dat=nod,nit=tot.n))
print(summary(fit_nod_u_quadnit))

##### Comparing AIC's ##############################################
deltaAICs(c(AICc(fit_nod_null,nod),AICc(fit_nod_nit,nod),AICc(fit_nod_u_nit,nod),AICc(fit_nod_p_nit,nod),AICc(fit_nod_u_quadnit,nod)))
####################################################################


############# DOES NODULATION VARY WITH LIGHT (% DIRECT TRANSMITTANCE)??? #################
##############################################################################
### Running a Null Model fitting a flat mean to the data #####################
##############################################################################
ZinflnormNLL_nod_null <- function(sdlognod,ulognod,p_no_nod,nod_dat){
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}
fit_nod_null <- mle2(ZinflnormNLL_nod_null,
                     start=list(sdlognod=sstart,ulognod=ustart,p_no_nod=pstart),
                     data=list(nod_dat=nod))
print(summary(fit_nod_null))
##################################################################
### Now running a linear model for light #####################
##################################################################
ZinflnormNLL_nod_lit <- function(sdlognod,c,b,e,d,nod_dat,lit){
  ulognod<-c+b*lit
  p_no_nod<-alogitfn(e+d*lit)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_lit <- mle2(ZinflnormNLL_nod_lit,
                    start=list(sdlognod=sstart,c=10,b=.1,e=.8,d=-.0001),
                    data=list(nod_dat=nod,lit=lit))
print(summary(fit_nod_lit))
#################################################################################################
### Now running a linear model for stand age where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_lit <- function(sdlognod,c,b,p_no_nod,nod_dat,lit){
  ulognod<-c+b*lit
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_lit <- mle2(ZinflnormNLL_nod_u_lit,
                      start=list(sdlognod=sstart,c=10,b=.1,p_no_nod=.5),
                      data=list(nod_dat=nod,lit=lit))
print(summary(fit_nod_u_lit))
#################################################################################################
### Now running a linear model for stand age where only the probability of 0 fluctuates#############
#################################################################################################

ZinflnormNLL_nod_p_lit <- function(ulognod,sdlognod,e,d,nod_dat,lit){
  p_no_nod<-alogitfn(e+d*lit)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}

fit_nod_p_lit <- mle2(ZinflnormNLL_nod_p_lit,
                      start=list(ulognod=ustart,sdlognod=sstart,e=.8,d=.01),
                      data=list(nod_dat=nod,lit=lit))
print(summary(fit_nod_p_lit))
#################################################################################################
### Now running a hump-shaped model for stand age where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_quadlit <- function(sdlognod,c,b,b2,p_no_nod,nod_dat,lit){
  ulognod<-c+b*lit+b2*lit^2
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_quadlit <- mle2(ZinflnormNLL_nod_u_quadlit,
                          start=list(sdlognod=sstart,c=10,b=-.1,b2=-.1,p_no_nod=.5),
                          data=list(nod_dat=nod,lit=lit))
print(summary(fit_nod_u_quadlit))

##### Comparing AIC's ##############################################
deltaAICs(c(AICc(fit_nod_null,nod),AICc(fit_nod_lit,nod),AICc(fit_nod_u_lit,nod),AICc(fit_nod_p_lit,nod),AICc(fit_nod_u_quadlit,nod)))
####################################################################

##############################################################################
############# DOES NODULATION VARY WITH PLANT SIZE??? #################
##############################################################################
##############################################################################
##############################################################################

with(ndfa, plot(nod.prop.bg~Total.biomass, col=Env.Trt,pch=16,cex=1.5))
tst2<-ndfa[ndfa$Total.biomass>4&ndfa$Total.biomass<10,]
with(tst2, summary(aov(nod.prop.bg~Light.Treatment+Total.biomass+Nitrogen.Treatment)))
with(tst2, plot(nod.prop.bg~Total.biomass, col=Light.Treatment,pch=16,cex=1.5))

### Now the effect of Plant size in Field Seedlings ##########################
##############################################################################
### Running a Null Model fitting a flat mean to the data #####################
##############################################################################
ZinflnormNLL_nod_null <- function(sdlognod,ulognod,p_no_nod,nod_dat){
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}
fit_nod_null <- mle2(ZinflnormNLL_nod_null,
                     start=list(sdlognod=sstart,ulognod=ustart,p_no_nod=pstart),
                     data=list(nod_dat=nod))
print(summary(fit_nod_null))
##################################################################
### Now running a linear model for root mass (plant size) #####################
##################################################################
ZinflnormNLL_nod_rt <- function(sdlognod,c,b,e,d,nod_dat,rt){
  ulognod<-c+b*rt
  p_no_nod<-alogitfn(e+d*rt)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_rt <- mle2(ZinflnormNLL_nod_rt,
                    start=list(sdlognod=sstart,c=10,b=.1,e=.8,d=-.0001),
                    data=list(nod_dat=nod,rt=rt))
print(summary(fit_nod_rt))
#################################################################################################
### Now running a linear model for root mass where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_rt <- function(sdlognod,c,b,p_no_nod,nod_dat,rt){
  ulognod<-c+b*rt
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_rt <- mle2(ZinflnormNLL_nod_u_rt,
                      start=list(sdlognod=sstart,c=10,b=.1,p_no_nod=.5),
                      data=list(nod_dat=nod,rt=rt))
print(summary(fit_nod_u_rt))
#################################################################################################
### Now running a linear model for root mass where only the probability of 0 fluctuates#############
#################################################################################################

ZinflnormNLL_nod_p_rt <- function(ulognod,sdlognod,e,d,nod_dat,rt){
  p_no_nod<-alogitfn(e+d*rt)
  -sum(log(p_no_nod[nod==0])) - sum(log(1 - p_no_nod[nod>0])) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod,sd=exp(sdlognod),log=TRUE))
}

fit_nod_p_rt <- mle2(ZinflnormNLL_nod_p_rt,
                      start=list(ulognod=ustart,sdlognod=sstart,e=.8,d=.01),
                      data=list(nod_dat=nod,rt=rt))
print(summary(fit_nod_p_rt))
#################################################################################################
### Now running a hump-shaped model for root mass where only the non-zero mean fluctuates#############
#################################################################################################
ZinflnormNLL_nod_u_quadrt <- function(sdlognod,c,b,b2,p_no_nod,nod_dat,rt){
  ulognod<-c+b*rt+b2*rt^2
  J <- sum(nod==0)
  K <- sum(nod>0)
  -J*log(p_no_nod) - K*log(1-p_no_nod) - 
    sum(dlnorm(nod_dat[nod_dat>0],mean=ulognod[nod_dat>0],sd=exp(sdlognod),log=TRUE))
}

fit_nod_u_quadrt <- mle2(ZinflnormNLL_nod_u_quadrt,
                          start=list(sdlognod=sstart,c=10,b=-.1,b2=-.1,p_no_nod=.5),
                          data=list(nod_dat=nod,rt=rt))
print(summary(fit_nod_u_quadrt))

##### Comparing AIC's ##############################################
deltaAICs(c(AICc(fit_nod_null,nod),AICc(fit_nod_rt,nod),AICc(fit_nod_u_rt,nod),AICc(fit_nod_p_rt,nod),AICc(fit_nod_u_quadrt,nod)))
####################################################################

##################################################################################################################
#### END #########################################################################################################
##################################################################################################################
##################################################################################################################