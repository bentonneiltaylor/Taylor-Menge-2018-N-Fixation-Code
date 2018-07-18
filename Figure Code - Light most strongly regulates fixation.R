#### Light Most Strongly Regulates N Fixation ###################################
### Figure Code - Started 9/20/2017 - Ben Taylor ################################
setwd("C:/Users/Benton/Desktop/Work Files/Research/Costa Rica/Greenhouse Competition Study/Taylor & Menge_Light Regulates Fixation/Data and Analysis Files")
#Read in the data and load packages
library("lattice", lib.loc="~/R/win-library/3.2")
library(ggplot2)
library(multcomp)
library(lsmeans)
library(grid)
library(ggpubr)
source("C:\\RCode\\Homegrown Functions.R")


ndfa.dat<-data.frame("Light.Treatment"=c(rep("L",3),rep("M",3),rep("H",3)),
                     "Nitogen.Treatment"=c(rep(c("L","M","H"),3)),
                     "Ndfa"=c(0,0,0,0,0,0,34.64,20.03, 20.49),
                     "Ndfa.se"=c(4.99,2.02,3.23,5.13,1.86,2.96,4.73,1.90,3.02),
                     "Nfixd"=c(0,0,0,0,0,0,.1084,.0884,.0759),
                     "Nfixd.se"=c(.00085,.00047,.00032,.0015,.0015,.0017,.0208,.0192,.0163))

gh<-read.csv("Ndfa mixing model output.csv")

#First we split things up by light treatments
hl<-gh[gh$Light.Treatment=="H",]
hl$Nitrogen.Treatment<-factor(hl$Nitrogen.Treatment, levels=c("L","M","H"))
ml<-gh[gh$Light.Treatment=="M",]
ml$Nitrogen.Treatment<-factor(ml$Nitrogen.Treatment, levels=c("L","M","H"))
#ml<-ml[ml$Plant.Number!="62.a.a",]
ll<-gh[gh$Light.Treatment=="L",]
ll$Nitrogen.Treatment<-factor(ll$Nitrogen.Treatment, levels=c("L","M","H"))

#Field Seedling Data
sdl<-read.csv("Seedling Nodule Sampling Data_2017.csv")
sdl$nod<-ifelse(sdl$nod.num>0, yes=1, no=0)

############################################################################################
#### Figure 1 - Plant Biomass, Nodulation, and %Ndfa #######################################
############################################################################################

box1<-function(x,y,main,xlab,ylab,ylm,bgcol,pan,let,recup,lh){
  nn<-c("0.51", "20","40")
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,-20,4.2,recup,col=bgcol,lty=1,border=NA)
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F, add=T, notch=F)
  text(x=c(1,2,3),y=lh,labels=let,col=1,cex = 1.7)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

box1.xlb<-function(x,y,main,xlab,ylab,ylm,bgcol,pan,let,recup,lh){
  nn<-c("0.51", "20","40")
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,-20,4.2,recup,col=bgcol,lty=1,border=NA)
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F, add=T, notch=F)
  mtext(side=1,text=expression("Nitrogen Added"~(g~N~m^{-2}~yr^{-1})),cex = 1.3,line=4)
  text(x=c(1,2,3),y=lh,labels=let,col=1,cex = 1.7)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

#png(filename = "Figure 1_Boxplot Biomass-Nodules-Ndfa.png", width = 11, height = 8, units="in",res = 300)
pdf(file = "Figure 1_Boxplot Biomass-Nodules-Ndfa.pdf", width = 11, height = 8)
par(mfrow=c(3,3))
#Row1
par(mar=c(2.1, 5.5, 2.5,1))
xlb<-"Nitrogen Added "
ylb<-"Plant Biomass (g)"
ylm<-c(0,40)
barcol<-c("deepskyblue","dodgerblue","dodgerblue4")
with(ll, box1(Nitrogen.Treatment, Total.biomass, main="8% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50", pan="a",let=c("a","a","a"),recup=50,lh=c(15,15,15)))
with(ml, box1(Nitrogen.Treatment, Total.biomass, main="16% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="b",let=c("b","b","b"),recup=50,lh=c(15,15,15)))
with(hl, box1(Nitrogen.Treatment, Total.biomass, main="40% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="c",let=c("c","d","d"),recup=50,lh=c(15,15,15)))
#Row2
par(mar=c(2.1, 5.5, .75,1))
xlb<-"Nitrogen Added"
ylb<-"Nodule Biomass (%)"
ylm<-c(0,30)
with(ll, box1(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50",pan="d",let=c("a","a","a"),recup=35,lh=c(10,10,10)))
with(ml, box1(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="e",let=c("a","a","a"),recup=35,lh=c(10,10,10)))
with(hl, box1(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="f",let=c("b","c","d"),recup=35,lh=c(11,11,11)))
#Row3
par(mar=c(4.75, 5.5, .75,1))
xlb<-"Nitrogen Added"
ylb<-expression("N"["dfa"]~"(%)")
ylm<-c(-15,65)
with(ll, box1.xlb(Nitrogen.Treatment, Ndfa, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50",pan="g",let=c("a","a","a"),recup=70,lh=c(22,22,22)))
with(ml, box1.xlb(Nitrogen.Treatment, Ndfa, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="h",let=c("a","a","a"),recup=70,lh=c(22,22,22)))
with(hl, box1.xlb(Nitrogen.Treatment, Ndfa, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="i",let=c("b","c","d"),recup=70,lh=c(50,50,50)))
dev.off()
par(mfrow=c(1,1))

########################################################################################
####### Figure 2 - Amount of N Fixed per Plant ############################################
########################################################################################

#png(filename = "Figure 2_Boxplot N fixed per plant.png",width=11, height=5,units = "in",res = 300)
pdf(file = "Figure 2_Boxplot N fixed per plant.pdf",width=11, height=5)
par(mfrow=c(1,3))
par(mar=c(6.1, 6.7, 4.1,3.5))
xlb<-"Nitrogen Added"
ylb<-expression("Nitrogen Fixed"~(g~N~plant^{-1}))
ylm<-c(-.01,.34)
with(ll, box1.xlb(Nitrogen.Treatment, Nfixd, main="8% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50",pan="a",let=c("a","a","a"),recup=70,lh=c(.12,.12,.12)))
with(ml, box1.xlb(Nitrogen.Treatment, Nfixd, main="16% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="b",let=c("a","a","a"),recup=70,lh=c(.12,.12,.12)))
with(hl, box1.xlb(Nitrogen.Treatment, Nfixd, main="40% Total Transmittance", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="c",let=c("b","c","c"),recup=70,lh=c(.12,.12,.12)))
dev.off()
par(mfrow=c(1,1))#sets things back to plot out one thing at a time
par(mar=c(5.1, 4.1, 4.1, 2.1))#sets the margins back to their defaults

########################################################################################
####### Figure 3 - Field Seedling Sampling ############################################
########################################################################################
#This uses the "probability of nodulation" plots from version 1 of this figure
#in the code above
#### Fitting an exponential model to Nodulation vs. Light data ##########
ex.nod.light<-with(sdl, lm(log(prop.nod.bio+1)~X..Trans.Tot))
summary(ex.nod.light)
lightvalues<-seq(6.3,21.1,0.01)
nod.light.pred<-exp(predict(ex.nod.light, list(X..Trans.Tot=lightvalues)))
expred.df<-data.frame("light"=lightvalues,"nod"=nod.light.pred-1)
#### Light plot with Logged y-axis #####################################
litlogsup<-ggplot(sdl, aes(x=X..Trans.Tot, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)",labels=c("0"="0 ", "5"="5 ", "10"="10","15"="15","20"="20"))+
  scale_x_continuous(breaks=c(6,12,18,24.5),labels=c("6","12","18","40"))+
  geom_segment(aes(x=5,xend=21,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=20),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=0,yend=0), colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=5,yend=5),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=10,yend=10),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=15,yend=15),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=20,yend=20),colour="black")+
  geom_segment(aes(x=23,xend=26,y=0,yend=0),colour="black")+
  geom_segment(aes(x=20.75,xend=21.25,y=-.5,yend=.5),colour="black")+
  geom_segment(aes(x=22.75,xend=23.25,y=-.5,yend=.5),colour="black")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5)+
  geom_point(aes(x=7.7,y=0),fill="deepskyblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=8,y=0),fill="dodgerblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=8.3,y=0),fill="dodgerblue4", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=16,y=.592),fill="deepskyblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=15.7,y=0),fill="dodgerblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=16.3,y=0),fill="dodgerblue4", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=24.3,y=21.32),fill="deepskyblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=24.5,y=9.319),fill="dodgerblue", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=24.7,y=1.538),fill="dodgerblue4", size=4.75, shape=22, stroke=1.35)+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

litlogins<-ggplot(sdl, aes(x=X..Trans.Tot, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans = "log", breaks=c(0,.025,.25,2.5,25),labels=c("0","0.025","0.25","2.5","25"))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=25),colour="white")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5,linetype=1)+
  theme(panel.border=element_rect(colour="black",fill=NA))

pnodlt<-1-(alogitfn(3.63220826+(-0.25519186)*lightvalues))
pnl.df<-data.frame("light"=lightvalues, "pnod"=pnodlt)
litprob<-ggplot(pnl.df,aes(x=lightvalues,y=pnod))+
  geom_line(size=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_blank(),axis.title.y=element_text(size=20),axis.ticks.y = element_blank())+
  theme(panel.background=element_blank(),
        plot.margin=unit(c(.2,1,.52,.05),"cm"))+
  ylab("Nodulation \n(probability)")+
  xlab(expression("Light"^{}~"(% Total"[]~"Transmittance)"))+
  geom_jitter(data=sdl, mapping=aes(x=X..Trans.Tot,y=nod), height=.02)+
  geom_segment(aes(x=5,xend=21,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=1),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=0,yend=0), colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.25,yend=.25),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.5,yend=.5),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.75,yend=.75),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=1,yend=1),colour="black")+
  annotate("text",x=3.45,y=0,label="0.00", size=5.6)+
  annotate("text",x=3.45,y=.25,label="0.25", size=5.6)+
  annotate("text",x=3.45,y=.5,label="0.50", size=5.6)+
  annotate("text",x=3.45,y=.75,label="0.75", size=5.6)+
  annotate("text",x=3.45,y=1,label="1.00", size=5.6)

littst<-ggarrange(litlogsup,litprob,heights=c(1.9,0.9),ncol=1,nrow=2,align="v")

#### Nitrogen plot with Logged y-axis ##################################
nod.n.lm<-with(sdl, lm(prop.nod.bio~tot.n+height))
nitlogsup<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab(expression("Soil Inorganic Nitrogen"~(mg~N~kg~soil^{-1})))+
  ylab(" ")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=18, colour="black"))+
  theme(axis.text.y=element_text(hjust=2.5,size=20, colour="black"),axis.ticks.y = element_blank())+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(labels=c("0"="0 ", "5"="5 ", "10"="10","15"="15","20"="20"))+
  scale_x_continuous(breaks=c(20,30,40,52),labels=c("20","30","40","78"))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=20),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=0,yend=0), colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=5,yend=5),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=10,yend=10),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=15,yend=15),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=20,yend=20),colour="black")+
  geom_smooth(method="lm", colour="black", size=1.5, se=F, linetype=2)+
  geom_segment(aes(x=50,xend=54,y=0,yend=0),colour="black")+
  geom_segment(aes(x=47.65,xend=48.35,y=-.5,yend=.5),colour="black")+
  geom_segment(aes(x=49.65,xend=50.35,y=-.5,yend=.5),colour="black")+
  geom_point(aes(x=13.8,y=0),fill="gray50", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=14.2,y=.59199),fill="gray70", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=14,y=21.32),fill="gray95", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=45.75,y=0),fill="gray50", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=46.25,y=0),fill="gray70", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=46,y=9.319),fill="gray95", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=51.75,y=0),fill="gray50", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=52.25,y=0),fill="gray70", size=4.75, shape=22, stroke=1.35)+
  geom_point(aes(x=52,y=1.538),fill="gray95", size=4.75, shape=22, stroke=1.35)+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin=unit(c(5.5,8.5,5.5,5.5),"points"))

nitlogins<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(axis.title=element_blank())+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans="log", breaks=c(0,.025,.25,2.5,25),labels=c("0","0.025","0.25","2.5","25"))+
  geom_segment(aes(x=14,xend=48,y=0,yend=0),colour="white")+
  geom_segment(aes(x=14,xend=14,y=0,yend=25),colour="white")+
  geom_segment(aes(x=15,xend=47,y=0.359,yend=0.359),colour="black",size=1.5,linetype=2)+
  theme(panel.border=element_rect(colour="black",fill=NA))

nitvec<-seq(15.21,47.84,.01)
pnn.df<-data.frame("nit"=nitvec,"pnod"=rep(.6299,3264))
nitprob<-ggplot(pnn.df,aes(x=nit,y=pnod))+
  geom_line(size=1)+
  geom_jitter(data=sdl, mapping=aes(x=tot.n,y=nod), height=.02)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_blank())+
  theme(axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank())+
  xlab(expression("Soil Inorganic Nitrogen"~(mg~N~kg~soil^{-1})))+
  scale_x_continuous(limits=c(10.5,48), breaks=c(20,30,40))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=1),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=0,yend=0), colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.25,yend=.25),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.5,yend=.5),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.75,yend=.75),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=1,yend=1),colour="black")+
  annotate("text",x=10.5,y=0,label="0.00", size=5.6)+
  annotate("text",x=10.5,y=.25,label="0.25", size=5.6)+
  annotate("text",x=10.5,y=.5,label="0.50", size=5.6)+
  annotate("text",x=10.5,y=.75,label="0.75", size=5.6)+
  annotate("text",x=10.5,y=1,label="1.00", size=5.6)

nittst<-ggarrange(nitlogsup,nitprob,heights=c(2,0.9),ncol=1,nrow=2,align="v")

v1<-viewport(width=.21,height=.23,x=.358,y=.8)
v2<-viewport(width=.21,height=.23,x=.89,y=.8)
#png(filename = "Figure 3_Inset version.png",width = 12.2, height = 9, units = "in",res = 300)
pdf(file = "Figure 3_Inset version.pdf",width = 12.2, height=9)
multiplot(littst, nittst, cols=2)
print(litlogins, vp=v1)
print(nitlogins, vp=v2)
dev.off()

############################################################################
########### Figure 4. Conceptual Theory Figure #############################
############################################################################
#png(filename = "Figure 4_Theory Conceptual Figure_V2.png",width = 12,height=6,units = "in",res=300)
pdf(file = "Figure 4_Theory Conceptual Figure_V2.pdf",width = 12,height=6)
par(mar=c(5.1,5.1,4.1,2.1),mfrow=c(1,2))
###Panel A
### Theory For N Limitation ###
plot(1, type="n", bty="l",xlab="N Limitation", ylab=expression(paste("N Fixation (%N"[dfa],")")), 
     xlim=c(0,90),ylim=c(0,100), xaxt="n", yaxt="n", cex.lab=2)
axis(1, at=c(10,45,80), labels=c("N limited","Co-Limited","Not N Limited"),cex.axis=1.2)
#polygon(x=c(30,30,60),y=c(1,79,1),border=NA, col="lightblue")
polygon(x=c(30,60,90,90,60),y=c(79,29.5,29.5,1,1),border=NA, col="pink")
segments(60,0,90,0,col="blue",lwd=4)
#segments(30,0,60,0,col="blue",lwd=4, lty=3)
segments(60,1,90,1,col="forestgreen",lwd=4)
#segments(30,1,60,1,col="forestgreen",lwd=4, lty=3)
segments(1,79,30,79,col="red",lwd=4)
segments(30,79,60,30,col="red",lwd=4)
segments(60,30,90,30,col="red",lwd=4)
segments(1,80,30,80,col="black",lwd=4, lty=1)
segments(30,80,60,2,col="black",lwd=4, lty=1)
segments(60,2,90,2,col="black",lwd=4, lty=1)
legend(x=59, y=100, legend=c("Theory","High Light","Medium Light","Low Light"),
       col=c("black", "red", "forestgreen","blue"), lty=c(1,1,1,1),box.lty=0, lwd=3,cex=1)
abline(v=30, col="black", lty=3, lwd=1.5)
abline(v=60, col="black", lty=3, lwd=1.5)
put.fig.letter(label="a", location="topleft",cex=2, font=2)
points(c(65,75,85,30,70,80),c(0,0,0,79,30,30),col=1,pch=16,cex=2.2)
points(c(65,75,85),c(0,0,0), col="blue", pch=16, cex=2)
points(c(30,70,80),c(79,30,30), col="red", pch=16, cex=2)
points(c(66,76,86),c(1,1,1), col="black", pch=16, cex=2.2)
points(c(66,76,86),c(1,1,1), col="forestgreen", pch=16, cex=2)


###Panel B (Recreation of Fig 3d from Duncan's Nature Plants Paper)
bdat<-data.frame("x"=c(10,20,30,40,50,60,90),"y"=c(22,43,74,92,99,91,0))
par(mar=c(5.1,5.5,4.1,2.1))
plot(1, type="n", bty="l",xlab="SNF Regulation", ylab=expression(paste("Available N Exports (kg N "~ha^-1*~y^-1,")")), 
     xlim=c(0,100),ylim=c(0,100), xaxt="n", cex.lab=1.9, cex.axis=1.5,axes=F)
axis(2)
axis(1, at=c(10,52,90), col=NA,labels=c("Obligate","Incomplete\nDownregulation","Facultative"),cex.axis=1.2)
segments(-5,0,100,0,col="black",lwd=1)
polygon(x=c(15,15,20,30,40,50,60,80,80),y=c(0,32.5,43,74,92,99,91,29.75,0),border=NA, col="pink")
points(bdat$x, bdat$y, lwd=2, type = "l")
#points(bdat$x, bdat$y, cex=2.2, col=1,pch=16)
#points(bdat$x, bdat$y, cex=2, col=c(1,"red","red","red","red","red","blue"),pch=16)
#points(91,1,col=1,cex=2.2,pch=16)
#points(91,1,col="forestgreen",cex=2,pch=16)
put.fig.letter(label="b", location="topleft",cex=2, font=2)
dev.off()
############################################################################################
############################################################################################
############################################################################################
#### Supplemental Figure 1 - Plant Biomass, Nodulation, and %Ndfa with points #######################################
############################################################################################
#### First, Plant Biomass #######################################
seBars1<-function(x,y,main,xlab,ylab,ylm,bgcol,pan,let){
  model<-lm(y~factor(x))
  reps<-length(y)/length(levels(x))
  sem<-summary(model)$sigma/sqrt(reps)
  m<-as.vector(tapply(y,x,mean, na.rm=T))
  #m<-m[c(2,3,1)]
  se<-as.vector(tapply(y,x,sefun))
  #se<-se[c(2,3,1)]
  nn<-c("14", "46","78")
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,0,4.2,24,col=bgcol,lty=1,border=NA)
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F,add=T)
  for(i in 1:length(xs)) {
    arrows(xs[i], m[i]+se[i], xs[i],m[i]-se[i],angle=90,code=3,length=0.1)}
  #mtext(side=1,text=expression("Nitrogen Added"~(g~m^{-2}~yr^{-1})),cex = 1.5,line=4)
  text(x=c(.7,1.9,3.1),y=c(7,7,7),labels=let,col=1,cex = 1.7)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

#### Next, Nodulation ############################################
seBars2<-function(x,y,main,xlab,ylab,ylm,bgcol,pan,let){
  model<-lm(y~factor(x))
  reps<-length(y)/length(levels(x))
  sem<-summary(model)$sigma/sqrt(reps)
  m<-as.vector(tapply(y,x,mean, na.rm=T))
  #m<-m[c(2,3,1)]
  se<-as.vector(tapply(y,x,sefun))
  #se<-se[c(2,3,1)]
  nn<-c("14", "46","78")
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,0,4.2,24,col=bgcol,lty=1,border=NA)
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F, add=T)
  for(i in 1:length(xs)) {
    arrows(xs[i], m[i]+se[i], xs[i],m[i]-se[i],angle=90,code=3,length=0.1)}
  #mtext(side=1,text=expression("Nitrogen Added"~(g~m^{-2}~yr^{-1})),cex = 1.5,line=4)
  text(x=c(.7,1.9,3.1),y=c(14,14,14),labels=let,col=1,cex = 1.7)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

#### Next, %Ndfa ################################################
seBars3<-function(x,main,xlab,ylab,ylm,bgcol,pan,let){
  #model<-lm(y~factor(x))
  #reps<-length(y)/length(levels(x))
  #sem<-summary(model)$sigma/sqrt(reps)
  m<-as.vector(x$Ndfa)
  #m<-m[c(2,3,1)]
  se<-as.vector(x$Ndfa.se)
  #se<-se[c(2,3,1)]
  nn<-c("14", "46","78")
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,0,4.2,50,col=bgcol,lty=1,border=NA)
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.5, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F, add=T)
  for(i in 1:length(xs)) {
    arrows(xs[i], m[i]+se[i], xs[i],m[i]-se[i],angle=90,code=3,length=0.1)}
  mtext(side=1,text=expression("Soil N"~(mg~N~kg~soil^{-1})),cex = 1.4,line=4)
  text(x=c(.7,1.9,3.1),y=c(25.8,25.8,25.8),labels=let,col=1,cex = 1.7)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

#### Now actually making the figure ####################

png(filename = "Supplemental Figure 1_Biomass-Nodules-Ndfa.png", width = 11, height = 8, units="in",res = 300)
par(mfrow=c(3,3))
#Row1
par(mar=c(2.1, 5.5, 2.5,1))
xlb<-"Nitrogen Added "
ylb<-"Plant Biomass (g)"
ylm<-c(0,22)
with(ll, seBars1(Nitrogen.Treatment, Total.biomass, main="8% Full Sunlight", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50", pan="a",let=c("a","a","a")))
with(ml, seBars1(Nitrogen.Treatment, Total.biomass, main="16% Full Sunlight", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="b",let=c("b","b","b")))
with(hl, seBars1(Nitrogen.Treatment, Total.biomass, main="40% Full Sunlight", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="c",let=c("c","d","d")))
#Row2
par(mar=c(2.1, 5.5, .75,1))
xlb<-"Nitrogen Added"
ylb<-"Nodule Biomass (%)"
ylm<-c(0,24)
with(ll, seBars2(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50",pan="d",let=c("a","a","a")))
with(ml, seBars2(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="e",let=c("a","a","a")))
with(hl, seBars2(Nitrogen.Treatment, nod.prop.bg, main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95",pan="f",let=c("b","c","d")))
#Row3
par(mar=c(4.75, 5.5, .75,1))
xlb<-"Nitrogen Added"
ylb<-expression("N"["dfa"]~"(%)")
ylm<-c(0,45)
seBars3(ndfa.dat[c(1:3),], main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50", pan="g",let=c("a","a","a"))
seBars3(ndfa.dat[c(4:6),], main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70", pan="h",let=c("a","a","a"))
seBars3(ndfa.dat[c(7:9),], main="", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray95", pan="i",let=c("b","c","c"))
dev.off()
par(mfrow=c(1,1))

########################################################################################
####### Supplemental Figure 2 - Amount of N Fixed per Plant with Points ############################################
########################################################################################


seBars4<-function(x,main,xlab,ylm,bgcol,pan,let){
  m<-as.vector(x$Nfixd)
  se<-as.vector(x$Nfixd.se)
  nn<-c("14", "46","78")
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.7, cex.names=1.7, bty="l", xlab="", 
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,0,4.2,50,col=bgcol,lty=1,border=NA)
  xs<-barplot(m,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.7, cex.names=1.7, bty="l", xlab="", 
              cex.lab=2, cex.main=2, xpd=F, add=T)
  for(i in 1:length(xs)) {
    arrows(xs[i], m[i]+se[i], xs[i],m[i]-se[i],angle=90,code=3,length=0.1)}
  mtext(side=1,text=expression("Soil N"~(mg~N~kg~soil^{-1})),cex = 1.4,line=4)
  mtext(side=2,text=expression("Nitrogen Fixed"~(g~N~plant^{-1})),cex = 1.4,line=4)
  text(x=c(.7,1.9,3.1),y=c(.02,.02,.02),labels=let,col=1,cex = 1.9)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

png(filename = "Supplemental Figure 2_N fixed per plant.png",width=11, height=5,units = "in",res = 300)
par(mfrow=c(1,3))
par(mar=c(6.1, 6.7, 4.1,3.5))
xlb<-"Nitrogen Added"
ylm<-c(0,.15)
seBars4(ndfa.dat[c(1:3),], main="8% Full Sunlight",xlab=xlb,ylm=ylm, bgcol="gray53", pan="a",let=c("a","a","a"))
seBars4(ndfa.dat[c(4:6),], main="16% Full Sunlight",xlab=xlb,ylm=ylm, bgcol="gray73", pan="b",let=c("a","a","a"))
seBars4(ndfa.dat[c(7:9),], main="40% Full Sunlight",xlab=xlb,ylm=ylm, bgcol="gray93", pan="c",let=c("b","c","c"))
dev.off()
par(mfrow=c(1,1))#sets things back to plot out one thing at a time
par(mar=c(5.1, 4.1, 4.1, 2.1))#sets the margins back to their defaults

########################################################################################
####### Supplemental Figure 3 - Relationship between plant size and nodulation ############################################
########################################################################################
png(filename = "Supplementary Figure 3_Effects of Env Trts on Same-size plants.png",width=12, height=7,units = "in",res = 300)
par(mar=c(5.1, 4.8, 4.1, 2.1))
sub<-gh[gh$Total.biomass>4&gh$Total.biomass<10,]
with(sub, plot(nod.prop.bg~Total.biomass, col=Light.Treatment,pch=c(1,2,3)[Nitrogen.Treatment],cex=2,bty="l",xlab="Total Plant Biomass (g)",
              ylab="Belowground Allocation to Nodules (%)",cex.axis=2, cex.lab=2))
legend(8,15,legend = c("Low", "Medium", "High"),title="Light", col=c(2,3,1),pch=1,box.lty=0, cex=1.5,bg = NA)
legend(9,15,legend = c("Low", "Medium", "High"),title="Nitrogen",pch=c(2,3,1),box.lty=0, cex=1.5,bg = NA)
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#sets the margins back to their defaults
#########################################################################
##### Fig 3 - Alternate version with linear plots and log insets ##########################
########################################################################
#This uses the "probability of nodulation" plots from version 1 of this figure
#in the code above
#### Fitting an exponential model to Nodulation vs. Light data ##########
ex.nod.light<-with(sdl, lm(log(prop.nod.bio+1)~X..Trans.Dir))
summary(ex.nod.light)
lightvalues<-seq(5.3,23.1,0.01)
nod.light.pred<-exp(predict(ex.nod.light, list(X..Trans.Dir=lightvalues)))
expred.df<-data.frame("light"=lightvalues,"nod"=nod.light.pred-1)
#### Light plot with Logged y-axis #####################################
litlogsup<-ggplot(sdl, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=1))+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=20),colour="black")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5)+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

litlogins<-ggplot(sdl, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans = "log", breaks=c(0,.025,.25,2.5,25),labels=c(0,.025,.25,2.5,25))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=25),colour="black")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5,linetype=1)
  
littst<-ggarrange(litlogsup,litprob,heights=c(1.9,0.9),ncol=1,nrow=2,align="v")

#### Nitrogen plot with Logged y-axis ##################################
nod.n.lm<-with(sdl, lm(prop.nod.bio~tot.n+height))
nitlogsup<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab(expression("Soil Inorganic Nitrogen (mgN"~kg^{-1}~")"))+
  ylab(" ")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=18, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=20),colour="black")+
  geom_smooth(method="lm", colour="black", size=1.5, se=F, linetype=2)+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin=unit(c(5.5,8.5,5.5,5.5),"points"))

nitlogins<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(axis.title=element_blank())+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans="log", breaks=c(0,.025,.25,2.5,25),labels=c(0,.025,.25,2.5,25))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=25),colour="black")+
  geom_segment(aes(x=15,xend=47,y=0.359,yend=0.359),colour="black",size=1.5,linetype=2)

nittst<-ggarrange(nitlogsup,nitprob,heights=c(2,0.9),ncol=1,nrow=2,align="v")

v1<-viewport(width=.21,height=.23,x=.39,y=.8)
v2<-viewport(width=.21,height=.23,x=.903,y=.8)
png(filename = "Figure 3_Inset version.png",width = 12, height = 9, units = "in",res = 300)
multiplot(littst, nittst, cols=2)
print(litlogins, vp=v1)
print(nitlogins, vp=v2)
dev.off()


########################################################################################################################
####### Supplementary Figure 4 - 9-panel version of Field Seedling Sampling ############################################
########################################################################################################################
#First we'll break the data into the 9 components
sll<-sdl[sdl$X..Trans.Tot<11&sdl$tot.n<26,]
sll$Nitrogen.Treatment<-"Low"
slm<-sdl[sdl$X..Trans.Tot<11&sdl$tot.n>26&sdl$tot.n<30.5,]
slm$Nitrogen.Treatment<-"Medium"
slh<-sdl[sdl$X..Trans.Tot<11&sdl$tot.n>30.5,]
slh$Nitrogen.Treatment<-"High"
sml<-sdl[sdl$X..Trans.Tot>11&sdl$X..Trans.Tot<13&sdl$tot.n<22.75,]
sml$Nitrogen.Treatment<-"Low"
smm<-sdl[sdl$X..Trans.Tot>11&sdl$X..Trans.Tot<13&sdl$tot.n>22.75&sdl$tot.n<27.5,]
smm$Nitrogen.Treatment<-"Medium"
smh<-sdl[sdl$X..Trans.Tot>11&sdl$X..Trans.Tot<13&sdl$tot.n>27.5,]
smh$Nitrogen.Treatment<-"High"
shl<-sdl[sdl$X..Trans.Tot>13&sdl$tot.n<26,]
shl$Nitrogen.Treatment<-"Low"
shm<-sdl[sdl$X..Trans.Tot>13&sdl$tot.n>26&sdl$tot.n<30.5,]
shm$Nitrogen.Treatment<-"Medium"
shh<-sdl[sdl$X..Trans.Tot>13&sdl$tot.n>30.5,]
shh$Nitrogen.Treatment<-"High"

sl<-rbind(sll,slm,slh)
sm<-rbind(sml,smm,smh)
sh<-rbind(shl,shm,shh)

box1.sup<-function(x,y,main,xlab,ylab,ylm,bgcol,pan,let,recup,lh,nn){
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.3, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F)
  rect(-2,-20,4.2,recup,col=bgcol,lty=1,border=NA)
  xs<-boxplot(y~x,ylim=ylm, names=nn, main=main, col=c("deepskyblue","dodgerblue","dodgerblue4"), 
              axis.lty=1, cex.axis=1.3, cex.names=1.5, bty="l", xlab="", ylab=ylab,
              cex.lab=2, cex.main=2, xpd=F, add=T, notch=F)
  mtext(side=1,text=expression("Soil N Group"~(mg~N~kg~soil^{-1})),cex = 1.3,line=4)
  text(x=c(1,2,3),y=lh,labels=let,col=1,cex = 1.3)
  abline(h=0, lty=3)
  put.fig.letter(label=pan, location="topleft",cex=2, font=2)
}

png(filename = "Supplementary Fig 4_Boxplot Nodulation in Field Seedlings.png",width=12, height=5,units = "in",res = 300)
par(mfrow=c(1,3))
par(mar=c(6.1, 6.7, 4.1,3.5))
xlb<-"Nitrogen Group"
ylb<-"Belowground Allocation to Nodules (%)"
ylm<-c(0,2.5)
with(sl, box1.sup(Nitrogen.Treatment, prop.nod.bio, main="Low Light (6-11%)", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray50",pan="a",
                  let=c(" "," "," "),recup=70,lh=c(.12,.12,.12),nn=c("18-26","26-30.5","30.5-48")))
with(sm, box1.sup(Nitrogen.Treatment, prop.nod.bio, main="Medium Light (11-13%)", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray60",pan="b",
                  let=c(" "," "," "),recup=70,lh=c(.12,.12,.12),nn=c("16.5-22.75","22.75-27.5","27.5-40.5")))
with(sh, box1.sup(Nitrogen.Treatment, prop.nod.bio, main="High Light (13-18%)", xlab=xlb, ylab=ylb,ylm=ylm, bgcol="gray70",pan="c",
                  let=c(" "," "," "),recup=70,lh=c(.12,.12,.12),nn=c("15-26","26-30.5","30.5-46.5")))
dev.off()
par(mfrow=c(1,1))#sets things back to plot out one thing at a time
par(mar=c(5.1, 4.1, 4.1, 2.1))#sets the margins back to their defaults
  
par(mfrow=c(3,3))
ggplot(ll, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=8,xend=13,y=0,yend=0),colour="black")+
  geom_segment(aes(x=8,xend=8,y=0,yend=2),colour="black")+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))
  
ggplot(lm, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=5,xend=12,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=4),colour="black")+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(lh, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=8,xend=13,y=0,yend=0),colour="black")+
  geom_segment(aes(x=8,xend=8,y=0,yend=4.5),colour="black")+
  ggtitle("c")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(ml, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=10,xend=15,y=0,yend=0),colour="black")+
  geom_segment(aes(x=10,xend=10,y=0,yend=4),colour="black")+
  ggtitle("d")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(mm, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=10,xend=16,y=0,yend=0),colour="black")+
  geom_segment(aes(x=10,xend=10,y=0,yend=4),colour="black")+
  ggtitle("e")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(mh, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=11,xend=17,y=0,yend=0),colour="black")+
  geom_segment(aes(x=11,xend=11,y=0,yend=17),colour="black")+
  ggtitle("f")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(hl, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=13,xend=24,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=9),colour="black")+
  ggtitle("g")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(hm, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=13,xend=20,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=5),colour="black")+
  ggtitle("h")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

ggplot(hh, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=9,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)")+
  geom_segment(aes(x=13,xend=22,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=8),colour="black")+
  ggtitle("i")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

#This uses the "probability of nodulation" plots from version 1 of this figure
#in the code above
#### Fitting an exponential model to Nodulation vs. Light data ##########
ex.nod.light<-with(sdl, lm(log(prop.nod.bio+1)~X..Trans.Dir))
summary(ex.nod.light)
lightvalues<-seq(5.3,23.1,0.01)
nod.light.pred<-exp(predict(ex.nod.light, list(X..Trans.Dir=lightvalues)))
expred.df<-data.frame("light"=lightvalues,"nod"=nod.light.pred-1)
#### Light plot with Logged y-axis #####################################
litlogsup<-ggplot(sdl, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab("Light (% Total Transmittance)")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=20, colour="black"),axis.ticks.y = element_blank())+
  scale_y_continuous(name="Belowground Allocation to Nodules (%)",labels=c("0"="0 ", "5"="5 ", "10"="10","15"="15","20"="20"))+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=20),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=0,yend=0), colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=5,yend=5),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=10,yend=10),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=15,yend=15),colour="black")+
  geom_segment(aes(x=4.5,xend=5,y=20,yend=20),colour="black")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5)+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))

litlogins<-ggplot(sdl, aes(x=X..Trans.Dir, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans = "log", breaks=c(0,.025,.25,2.5,25),labels=c(0,.025,.25,2.5,25))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=25),colour="black")+
  geom_line(data=expred.df, aes(light, nod), colour="black", size=1.5,linetype=1)

pnodlt<-1-(alogitfn(3.63220826+(-0.25519186)*lightvalues))
pnl.df<-data.frame("light"=lightvalues, "pnod"=pnodlt)
litprob<-ggplot(pnl.df,aes(x=lightvalues,y=pnod))+
  geom_line(size=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_blank(),axis.title.y=element_text(size=18),axis.ticks.y = element_blank())+
  theme(panel.background=element_blank(),
        plot.margin=unit(c(.2,1,.52,.05),"cm"))+
  ylab("Nodulation \n(probability)")+
  xlab(expression("Light"^{}~"(% Total"[]~"Transmittance)"))+
  geom_segment(aes(x=5,xend=23,y=0,yend=0),colour="black")+
  geom_segment(aes(x=5,xend=5,y=0,yend=1),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=0,yend=0), colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.25,yend=.25),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.5,yend=.5),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=.75,yend=.75),colour="black")+
  geom_segment(aes(x=4.7,xend=5,y=1,yend=1),colour="black")+
  geom_text(aes(x=3.2,y=0),label="0.00", size=5.6)+
  geom_text(aes(x=3.2,y=.25),label="0.25", size=5.6)+
  geom_text(aes(x=3.2,y=.5),label="0.50", size=5.6)+
  geom_text(aes(x=3.2,y=.75),label="0.75", size=5.6)+
  geom_text(aes(x=3.2,y=1),label="1.00", size=5.6)

littst<-ggarrange(litlogsup,litprob,heights=c(1.9,0.9),ncol=1,nrow=2,align="v")

#### Nitrogen plot with Logged y-axis ##################################
nod.n.lm<-with(sdl, lm(prop.nod.bio~tot.n+height))
nitlogsup<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  xlab(expression("Soil Inorganic Nitrogen"~(mg~N~kg~soil^{-1})))+
  ylab(" ")+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=18, colour="black"))+
  theme(axis.text.y=element_text(hjust=2.5,size=20, colour="black"),axis.ticks.y = element_blank())+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(labels=c("0"="0 ", "5"="5 ", "10"="10","15"="15","20"="20"))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=20),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=0,yend=0), colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=5,yend=5),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=10,yend=10),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=15,yend=15),colour="black")+
  geom_segment(aes(x=12.4,xend=13,y=20,yend=20),colour="black")+
  geom_smooth(method="lm", colour="black", size=1.5, se=F, linetype=2)+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin=unit(c(5.5,8.5,5.5,5.5),"points"))

nitlogins<-ggplot(sdl, aes(x=tot.n, y=prop.nod.bio))+
  geom_point(size=3.5, shape=1)+
  theme(axis.title=element_blank())+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(colour="black", size=15))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  scale_y_continuous(trans="log", breaks=c(0,.025,.25,2.5,25),labels=c(0,.025,.25,2.5,25))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=25),colour="black")+
  geom_segment(aes(x=15,xend=47,y=0.359,yend=0.359),colour="black",size=1.5,linetype=2)

nitvec<-seq(15.21,47.84,.01)
pnn.df<-data.frame("nit"=nitvec,"pnod"=rep(.6299,3264))
nitprob<-ggplot(pnn.df,aes(x=nit,y=pnod))+
  geom_line(size=1)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.background=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y = element_blank())+
  xlab(expression("Soil Inorganic Nitrogen"~(mg~N~kg~soil^{-1})))+
  scale_x_continuous(limits=c(10.5,48), breaks=c(20,30,40))+
  geom_segment(aes(x=13,xend=48,y=0,yend=0),colour="black")+
  geom_segment(aes(x=13,xend=13,y=0,yend=1),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=0,yend=0), colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.25,yend=.25),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.5,yend=.5),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=.75,yend=.75),colour="black")+
  geom_segment(aes(x=12.7,xend=13,y=1,yend=1),colour="black")+
  geom_text(aes(x=10.5,y=0),label="0.00", size=5.6)+
  geom_text(aes(x=10.5,y=.25),label="0.25", size=5.6)+
  geom_text(aes(x=10.5,y=.5),label="0.50", size=5.6)+
  geom_text(aes(x=10.5,y=.75),label="0.75", size=5.6)+
  geom_text(aes(x=10.5,y=1),label="1.00", size=5.6)

nittst<-ggarrange(nitlogsup,nitprob,heights=c(2,0.9),ncol=1,nrow=2,align="v")


#################################################################################################
############################# END ###############################################################
#################################################################################################