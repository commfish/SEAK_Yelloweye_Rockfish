###############################################################################
# LIBRARIES
###############################################################################
library(extrafont)
library(ggplot2)
library(plyr)
library(latticeExtra)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(coda)


#------------------------------------------------------------------------------
# Graphics defaults
#------------------------------------------------------------------------------
#font_import()
#loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Calibri"))
theme_set(theme_bw(base_size=18,base_family="Calibri")+
            #theme(legend.position = c(0.25,0.9))+
            theme(legend.direction = "vertical")+
            theme(legend.margin = unit(1,"cm"))+
            theme(axis.title.y=element_text(vjust=1))+
            theme(axis.title.x=element_text(vjust=-0.5))+
            theme(legend.title = element_blank())+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


#------------------------------------------------------------------------------
# Data calls
#------------------------------------------------------------------------------
.REP <- "model" 
source("globals.R")


#------------------------------------------------------------------------------
# Model comparisons
#
# 'model'   = global model with estimated M
# 'model_M' = global model with fixed M
#------------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# ADMB report files - output defined in FINAL_CALCS section of .tpl file
#
# NOTE: you must have all four ADMB files written here: .par, .std, .cor, .rep
#
#------------------------------------------------------------------------------
model    <- read.admb("model_output/model")
modelM   <- read.admb("model_output/modelM")


# ADMB .std files (contain parameter values and standard deviation)
#       - all defined parameters with whatever 'sdreport_vector / matrix'
#         calls are defined in the .tpl file
par   <- read.delim("model_output/model.std",header=TRUE,sep="")$value
var   <- read.delim("model_output/model.std",header=TRUE,sep="")$std
parM  <- read.delim("model_output/modelM.std",header=TRUE,sep="")$value
varM  <- read.delim("model_output/modelM.std",header=TRUE,sep="")$std

# MCMC - primary parameters to check distributions and bounds
theta    <- read.delim("model_output/theta.ps",header=FALSE,sep="")
thetaM   <- read.delim("model_output/thetaM.ps",header=FALSE,sep="")

# MCMC - derived quantities
biomass    <- read.delim("model_output/bm.ps",header=FALSE,sep="")
biomassM   <- read.delim("model_output/bmM.ps",header=FALSE,sep="")
density    <- read.delim("model_output/dens.ps",header=FALSE,sep="")
densityM   <- read.delim("model_output/densM.ps",header=FALSE,sep="")
recruit    <- read.delim("model_output/rec.ps",header=FALSE,sep="")
recruitM   <- read.delim("model_output/recM.ps",header=FALSE,sep="")
spbiomass  <- read.delim("model_output/sbm.ps",header=FALSE,sep="")
spbiomassM <- read.delim("model_output/sbmM.ps",header=FALSE,sep="")



#------------------------------------------------------------------------------
#
# RETROSPECTIVE ANALYSES 
#   You'll have to run these individually and manually
#
#-------------------------------------------------------------------------------
par_R1m   <- read.delim("model_output/model_r1m.std",header=TRUE,sep="")$value
var_R1m   <- read.delim("model_output/model_r1m.std",header=TRUE,sep="")$std
par_R2m   <- read.delim("model_output/model_r2m.std",header=TRUE,sep="")$value
var_R2m   <- read.delim("model_output/model_r2m.std",header=TRUE,sep="")$std
par_R4m   <- read.delim("model_output/model_r4m.std",header=TRUE,sep="")$value
var_R4m   <- read.delim("model_output/model_r4m.std",header=TRUE,sep="")$std
par_R5m   <- read.delim("model_output/model_r5m.std",header=TRUE,sep="")$value
var_R5m   <- read.delim("model_output/model_r5m.std",header=TRUE,sep="")$std
par_R9m   <- read.delim("model_output/model_r9m.std",header=TRUE,sep="")$value
var_R9m   <- read.delim("model_output/model_r9m.std",header=TRUE,sep="")$std
par_R10m  <- read.delim("model_output/model_r10m.std",header=TRUE,sep="")$value
var_R10m  <- read.delim("model_output/model_r10m.std",header=TRUE,sep="")$std

par_R1   <- read.delim("model_output/model_r1.std",header=TRUE,sep="")$value
var_R1   <- read.delim("model_output/model_r1.std",header=TRUE,sep="")$std
par_R2   <- read.delim("model_output/model_r2.std",header=TRUE,sep="")$value
var_R2   <- read.delim("model_output/model_r2.std",header=TRUE,sep="")$std
par_R6   <- read.delim("model_output/model_r6.std",header=TRUE,sep="")$value
var_R6   <- read.delim("model_output/model_r6.std",header=TRUE,sep="")$std
par_R9   <- read.delim("model_output/model_r9.std",header=TRUE,sep="")$value
var_R9   <- read.delim("model_output/model_r9.std",header=TRUE,sep="")$std
par_R10  <- read.delim("model_output/model_r10.std",header=TRUE,sep="")$value
var_R10  <- read.delim("model_output/model_r10.std",header=TRUE,sep="")$std


#------------------------------------------------------------------------------
# Parse and sort MCMC for variances
#   These already have every 100th saved, and if I recall correctly, as ADMB
#   begins the draws from the approximately multivariate normal of the Hessian
#   matrix, a deleted burn-in isn't stricly necessary, but as I am a
#   statistical dunce, I still do it.
#
# MCMC outputs of 1,000,000 draws w/100th saved =  10,000
#  Trim first 20% = 8,000 remaining
#------------------------------------------------------------------------------
theta      <- theta[c(2001:10000),]
thetaM     <- thetaM[c(2001:10000),]
biomass    <- biomass[c(2001:10000),]
biomassM   <- biomassM[c(2001:10000),]
density    <- density[c(2001:10000),]
densityM   <- densityM[c(2001:10000),]
recruit    <- recruit[c(2001:10000),]
recruitM   <- recruitM[c(2001:10000),]
spbiomass  <- spbiomass[c(2001:10000),]
spbiomassM <- spbiomassM[c(2001:10000),]

# Sort each MCMC matrix by column and select 5% and 95%
# ASSUMES ALL MCMC MATRICES ARE OF THE SAME DIMENSIONS
min<-0.05*length(theta[,1])
max<-0.95*length(theta[,1])

for(i in 1:10)
{
  theta[,i]<-sort(theta[,i]) 
  thetaM[,i]<-sort(thetaM[,i]) 
}

for(i in 1:31)
{
  biomass[,i]    <- sort(biomass[,1])
  biomassM[,i]   <- sort(biomassM[,1])
  density[,i]    <- sort(density[,1])
  densityM[,i]   <- sort(densityM[,1])
  recruit[,i]    <- sort(recruit[,1])
  recruitM[,i]   <- sort(recruitM[,1])
  spbiomass[,i]  <- sort(spbiomass[,1])
  spbiomassM[,i] <- sort(spbiomassM[,1])
}



#------------------------------------------------------------------------------
# THETA PARAMETERS
#------------------------------------------------------------------------------
parMcmc <-data.frame(Par  = c(unlist(theta[,1:10]),unlist(thetaM[,1:10])),
                      Name = c(rep(c(rep("M",8000),rep("Mean Recruitment",8000),
                               rep("Mean Year 1 N",8000), rep("Sigma Year 1",8000),
                               rep("Mean F",8000),rep("Mean F sport",8000),
                               rep("Slope (selectivity)",8000),
                               rep("Age 50% sel.",8000),
                               rep("Catchability (commercial)",8000),
                               rep("Catchability (IPHC survey)",8000)),2)),
                      Model =  c(rep("Estimated M",80000),
                                 rep("Fixed M",80000)))

plot.parmc <- function(df) {
  ggplot(subset(df, Name %in% c("M","Mean Recruitment", "Mean Year 1 N",
                                      "Sigma Year 1", "Mean F","Mean F sport",
                                      "Slope (selectivity)", "Age 50% sel.",
                                      "Catchability (commercial)",
                                      "Catchability (IPHC survey)")),
         aes(Par,group=interaction(Name,Model))) + geom_density(aes(Par, y=..scaled.., fill = Model),alpha=.2)+
         facet_wrap(~Name,scales="free",ncol=2)+
         labs(y="")+
         labs(x="")+
         theme(legend.margin=unit(-0.25, "cm"))+
         guides(fill = guide_legend(ncol=3))+
         theme(legend.position = 'bottom')+
         theme(legend.key.size = unit(0.5, "cm"))+
         scale_fill_manual(name='', values=c('Estimated M' = 'black', 
                                             'Fixed M' = 'blue'))
}



#------------------------------------------------------------------------------
# DENSITY GLOBAL
#------------------------------------------------------------------------------
dfD <- data.frame(Year  = rep(seq(1985,2015,by=1),2), 
                  Density = c(par[293:323],parM[292:322]),
                  Stdev   = c(var[293:323],varM[292:322]),
                  Model = c(rep("Estimated M",31), rep("Fixed M",31)),
                  Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),2),
                  Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),2))

plot.density <- function(df) {
  df <- df %>%
    mutate(p.var = 1.96 * Var) %>%
    mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
  df <- df %>% 
    mutate(err= 1.96 * Stdev) %>%
    mutate(lower = Density - err, upper = Density + err)
  
  
  ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    geom_point(aes(x=Year,y= Rov), colour="black",size=4)+
    geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
    scale_colour_manual(name='', values=c('Estimated M' = 'black', 
                                          'Fixed M' = "blue"))+
    scale_fill_manual(name='', values=c('Estimated M' = 'black', 
                                        'Fixed M' = "blue"))+
    labs(x="Year")+
    labs(y="Density (ind. per square kilometer)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))
  
}


#------------------------------------------------------------------------------
# SPAWNING BIOMASS 
#------------------------------------------------------------------------------
dfS <- data.frame(Year  = rep(seq(1985,2015,by=1),2), 
                  Density = c(par[190:220],parM[189:219]),
                  Stdev   = c(var[190:220],varM[189:219]),
                  Model = c(rep("Estimated M",31),rep("Fixed M",31)))

plot.spawners <- function(df) {
  df <- df %>% 
    mutate(err= 1.96 * Stdev) %>%
    mutate(lower = Density - err, upper = Density + err)
  
  
  ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    scale_colour_manual(name='', values=c('Estimated M' = 'black', 
                                          'UnEstimated M' = 'darkorange3',
                                          'Fixed M' = "blue"))+
    scale_fill_manual(name='', values=c('Estimated M' = 'black', 
                                        'UnEstimated M' = 'darkorange3',
                                        'Fixed M' = "blue"))+
    labs(x="Year")+
    labs(y="Spawning biomass (tons)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    theme(legend.key.size = unit(1., "cm"))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))
  
}


#------------------------------------------------------------------------------
# RECRUITMENT 
#------------------------------------------------------------------------------
dfR <- data.frame(Year  = rep(seq(1985,2015,by=1),2), 
                  Density = c(par[159:189],parM[158:188]),
                  Stdev   = c(var[159:189],varM[158:188]),
                  Model = c(rep("Estimated M",31),rep("Fixed M",31)))

plot.recruit <- function(df) {
  df <- df %>% 
    mutate(err= 1.96 * Stdev) %>%
    mutate(lower = Density - err, upper = Density + err)
  
  
  ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    scale_colour_manual(name='', values=c('Estimated M' = 'black', 
                                          'Fixed M' = "blue"))+
    scale_fill_manual(name='', values=c('Estimated M' = 'black', 
                                        'Fixed M' = "blue"))+
    labs(x="Year")+
    labs(y="Recruitment (thousands)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    theme(legend.key.size = unit(1., "cm"))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))
  
}


#------------------------------------------------------------------------------
# Abundance at age
#------------------------------------------------------------------------------
yearT <- seq(1985,2015,by=1)
yearX <- c('1988','1989',seq(1991,2005,by=1),seq(2008,2015,by=1))
age <- seq(8,75,by=1)


SD <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
                 Age = rep(as.numeric(age), each=length(yearT)),
                 p = as.vector(model$natage))

M <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
                Age = rep(as.numeric(age), each=length(yearT)),
                p = as.vector(modelM$natage))

Tr <- data.frame(Year = rep(as.numeric(yearX),length(age)), 
                Age = rep(as.numeric(age), each=length(yearX)),
                p = as.vector(modelM$oac_fish))

png(file='figures/age1.png', res=300, width=7, height=9, units ="in", bg="transparent")
par(mfrow = c(2,1))
symbols(Tr$Year, Tr$Age, Tr$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
        xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition")
symbols(SD$Year, SD$Age, SD$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
        xlab="Year", ylab = "Estimated age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Abundance at age - Estimated M")
  dev.off()
  
  png(file='figures/age2.png', res=300, width=7, height=9, units ="in", bg="transparent")
  par(mfrow = c(2,1))
  symbols(Tr$Year, Tr$Age, Tr$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
          xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition")
  symbols(M$Year, M$Age, M$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
          xlab="Year", ylab = "Estimated age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Abundance at age - fixed M")
  dev.off()
  
  
  #------------------------------------------------------------------------------
  # Catch at age (looking for recruitment indicators)
  #------------------------------------------------------------------------------
  # yearC <- c('1988',seq(1991,2004,by=1),seq(2008,2015,by=1))
  # yearS <- c('1985','1987','1988','1989',	seq(1991,2004,by=1),'2008','2009','2010','2011','2012','2013')
  # yearE <- c(seq(1991,2001,by=1),'2004','2005',seq(2008,2015,by=1))
  # age <- seq(8,75,by=1)
  # 
  # 
  # C <- data.frame(Year = rep(as.numeric(yearC),length(age)), 
  #                  Age = rep(as.numeric(age), each=length(yearC)),
  #                  p = as.vector(model_s$oac_fishC))
  # 
  # E <- data.frame(Year = rep(as.numeric(yearE),length(age)), 
  #                 Age = rep(as.numeric(age), each=length(yearE)),
  #                 p = as.vector(model_s$oac_fishE))
  # 
  # S <- data.frame(Year = rep(as.numeric(yearS),length(age)), 
  #                  Age = rep(as.numeric(age), each=length(yearS)),
  #                  p = as.vector(model_s$oac_fishS))
  # 
  # png(file='../../figures/catch_age_dataC.png', res=300, width=9, height=7, units ="in", bg="transparent")
  # symbols(C$Year, C$Age, C$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
  #         xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition CSEO")
  # dev.off()
  # 
  # png(file='../../figures/catch_age_dataE.png', res=300, width=9, height=7, units ="in", bg="transparent")
  # symbols(E$Year, E$Age, E$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75), 
  #         xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition EYKT")
  # dev.off()
  # 
  # png(file='../../figures/catch_age_dataS.png', res=300, width=9, height=7, units ="in", bg="transparent")
  # symbols(S$Year, S$Age, S$p, inches=0.2, xlim=c(1985,2015), ylim=c(8,75),
  #         xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition SSEO")
  # 
  # dev.off()
  # 




  
  #------------------------------------------------------------------------------
  # Selectivity
  #------------------------------------------------------------------------------
  dfSel <- data.frame(Year =  rep(seq(8,75,by=1),2),
                      Sel  =  c(model$fish_sel,
                                modelM$fish_sel),
                      Model =  c(rep("Estimated M",68),rep("Fixed M",68)))
  
  plot.sel <- function(df) {
    
    ggplot(df,aes(Year,Sel, group=Model, colour=Model)) + geom_line(aes(Year,Sel,group=Model, colour=Model),size=1)+
      scale_colour_manual(name='', values=c('Estimated M' = 'black', 
                                            'UnEstimated M' = 'darkorange3',
                                            'Fixed M' = 'blue'))+
      labs(x="Age")+
      labs(y="Selectivity")+
      theme(legend.position = c(0.75,0.8))+
      theme(legend.direction = "vertical")+
      theme(legend.key = element_rect(fill = "white"))
    
  }
  
  
  #------------------------------------------------------------------------------
  # Fishing mortality
  #------------------------------------------------------------------------------
  dfF <- data.frame(Year =   rep(seq(1985,2015,by=1),2),
                    Mort  =  c(model$Fmort,
                               modelM$Fmort),
                    Model =  c(rep("Estimated M",31),rep("Fixed M",31)))
  
  plot.Fmt <- function(df) {
    
    ggplot(df,aes(Year,Mort, group=Model, colour=Model)) + 
      geom_line(aes(Year,Mort,group=Model, colour=Model))+
      geom_hline(yintercept=0.02)+
      scale_colour_manual(name='', values=c('Estimated M' = 'black', 
                                            'UnEstimated M' = 'darkorange3',
                                            'Fixed M' = 'blue'))+
      labs(x="Year")+
      labs(y="Mortality")+
      theme(legend.position = c(0.75,0.85))+
      theme(legend.direction = "vertical")+
      theme(legend.key = element_rect(fill = "white"))
    
  }


  #------------------------------------------------------------------------------
  # Catch age residuals
  #------------------------------------------------------------------------------
  yearT <- c('1988','1989',seq(1991,2005,by=1),seq(2008,2015,by=1))
  age <- seq(8,75,by=1)
  

  dfs <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
                    Age = rep(as.numeric(age), each=length(yearT)),
                    p = as.vector(model$res_age))
  
  dfm <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
                    Age = rep(as.numeric(age), each=length(yearT)),
                    p = as.vector(modelM$res_age))
  
  
  plot.age <- function(dfs,dfm) {
    
 A <-   ggplot(df,aes(Year,Age,size=p))+
      geom_point(aes(colour=p))+
      scale_size_area(max_size = 5, guide=FALSE)+
      scale_colour_gradient2(limits=c(-0.04,0.08),low="red", high="blue", mid = "grey", midpoint = 0.0)+
      xlab("Age")+
      ylab("Residuals (o - p)")+
      theme(legend.margin=unit(0.2, "cm"))+
      theme(legend.position = "right")
    
B <-    ggplot(dfm,aes(Year,Age,size=p))+
      geom_point(aes(colour=p))+
      scale_size_area(max_size = 5, guide=FALSE)+
      scale_colour_gradient2(limits=c(-0.04,0.08),low="red", high="blue", mid = "grey", midpoint = 0.0)+
      xlab("Age")+
      ylab("Residuals (o - p)")+
      theme(legend.margin=unit(0.2, "cm"))+
      theme(legend.position = "right")

grid.arrange(A,B)
    
  }
  
  
  #------------------------------------------------------------------------------
  # COMMERCIAL CPUE 
  #------------------------------------------------------------------------------
  dfcp <- data.frame(Year  =  as.factor(c(rep(c(seq(1992,2005,by=1),seq(2008,2015,by=1)),3))),
                     CPUE   =  c(model$obs_cpue,model$pred_cpue,modelM$pred_cpue),
                     Var    =  c(model$var_cpue,var[324:345],varM[323:344]),
                     Model  =  c(rep("Observed",22),rep("Estimated M",22),rep("Fixed M",22)))
  
  
  plot.cpue <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Var) %>%
      mutate(lower = CPUE - err, upper = CPUE + err)
    
    ggplot(df,aes(Year,CPUE, group=Model, colour=Model)) + 
      geom_line(data=df,aes(Year,CPUE,group=Model, colour=Model),lwd=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model, colour=NA),alpha=0.2, show.legend=FALSE) + 
      geom_point(data=df, aes(Year,CPUE,group=Model),size=2)+
      labs(y="CPUE (pounds per hook)")+
      scale_x_discrete(breaks=seq(1992,2015,by=4))+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
      theme(legend.position = c(0.4,0.85))+
      theme(legend.direction="vertical")+
      scale_colour_manual(name='', values=c('Estimated M' = 'black', "Observed" = 'green', 'Fixed M' = 'blue'))+
      scale_fill_manual(name='', values=c('Estimated M' = 'black', "Observed" = 'green', 'Fixed M' = 'blue'))
    
  }
  
  
  
  #------------------------------------------------------------------------------
  # IPHC CPUE
  #------------------------------------------------------------------------------
  dfcpi<- data.frame(Year  =  as.factor(c(rep(seq(1998,2015,by=1),3))),
                     CPUE   =  c(model$obs_cpue_iphc,
                                 par[346:363],parM[345:362]),
                     Var    =  c(model$var_cpue_iphc,var[346:363],
                                 varM[345:362]),
                     Model  =  c(rep("Observed",18),rep("Estimated M",18),
                                 rep("Fixed M",18)))
  
  plot.cpuei <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Var) %>%
      mutate(lower = CPUE - err, upper = CPUE + err)
    
    ggplot(df,aes(Year,CPUE, group=Model, colour=Model)) + 
      geom_line(data=df,aes(Year,CPUE,group=Model, colour=Model), lwd=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model, colour=NA),alpha=0.2, show.legend=FALSE) + 
      geom_point(data=df,aes(Year,CPUE,group=Model),size=2)+
      labs(y=("CPUE (individuals per hook)"))+
      scale_x_discrete(breaks=seq(1992,2015,by=4))+
      theme(legend.position = c(0.7,0.85))+
      theme(legend.background = element_blank())+
      theme(legend.margin = unit(0.5,"cm"))+
      theme(legend.direction="vertical")+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
      scale_colour_manual(name='', values=c('Estimated M' = 'black', "Observed" = 'green', "Fixed M" = "blue"))+
      scale_fill_manual(name='', values=c('Estimated M' = 'black', "Observed" = 'green', "Fixed M" = 'blue'))
    
  }
  

  
#------------------------------------------------------------------------------
# Spawning biomass projections 
#------------------------------------------------------------------------------
  dfproj <- data.frame(Year  = rep(seq(2016,2030,by=1),2), 
                       Density = c(par[367:381],parM[366:380]),
                       Error  =   c(var[367:381],varM[366:380]),
                       Model = c(rep("Recommended F = 0.022 (F_65) Estimated M",15),
                                 rep("Recommended F = 0.022 (F_65) Fixed M",15)))
  
  
  plot.proj<- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      scale_colour_manual(name='', values=c('Recommended F = 0.022 (F_65) Estimated M' = 'black', 
                                            'Recommended F = 0.022 (F_65) Fixed M' = 'blue'))+
      scale_fill_manual(name='', values=c('Recommended F = 0.022 (F_65) Estiated M' = 'black',  
                                          'Recommended F = 0.022 (F_65) Fixed M' = 'blue'))+
      labs(x="Year")+
      ylim(0,6000)+
      labs(y="Projected spawning biomass")+
      theme(legend.position = c(0.5,0.2))+
      theme(legend.direction = "vertical")+
      theme(legend.key.size = unit(1., "cm"))+
      theme(legend.background = element_blank())+
      #ggtitle("Projected spawning biomass (mt)")+
      theme(legend.key = element_rect(fill = "white"))
    
  }
  
  
  #------------------------------------------------------------------------------
  # DENSITY GLOBAL RETROSPECTIVE
  #------------------------------------------------------------------------------
  dfDR <- data.frame(Year  = rep(seq(1985,2015,by=1),6), 
                     Density = c(par[283:313],
                                 par_R1[217:246],NA,
                                 par_R2[213:241],NA,NA,
                                 par_R6[197:221],NA,NA,NA,NA,NA,NA,
                                 par_R9[185:206],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                 par_R10[181:201],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                     Error = c(var[283:313],
                               var_R1[217:246],NA,
                               var_R2[213:241],NA,NA,
                               var_R6[197:221],NA,NA,NA,NA,NA,NA,
                               var_R9[185:206],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                               var_R10[181:201],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                     Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),rep("retro 6",31),
                               rep("retro 9",31), rep("retro 10",31)),
                     Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),6),
                     Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),6))
  
  dfDR$Model <- factor(dfDR$Model, as.character(dfDR$Model))
  
  
  plot.retrodensity <- function(df) {
    df <- df %>%
      mutate(p.var = 1.96 * Var) %>%
      mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      geom_point(aes(x=Year,y= Rov), colour="black",size=3)+
      geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                            'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                          'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      labs(x="Year")+
      labs(y="Density (ind. per square kilometer)")+
      theme(legend.position = c(0.87,0.75))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
  }  
  
  
  #------------------------------------------------------------------------------
  # DENSITY GLOBAL RETROSPECTIVE FIXED M
  #------------------------------------------------------------------------------
  dfDRm <- data.frame(Year  = rep(seq(1985,2015,by=1),7), 
                     Density = c(parM[282:312],
                                 par_R1m[216:245],NA,
                                 par_R2m[212:240],NA,NA,
                                 par_R4m[204:230],NA,NA,NA,NA,
                                 par_R5m[200:225],NA,NA,NA,NA,NA,
                                 par_R9m[184:205],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                 par_R10m[180:200],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                     Error = c(varM[282:312],
                               var_R1m[216:245],NA,
                               var_R2m[212:240],NA,NA,
                               var_R4m[204:230],NA,NA,NA,NA,
                               var_R5m[200:225],NA,NA,NA,NA,NA,
                               var_R9m[184:205],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                               var_R10m[180:200],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                     Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),rep("retro 4",31),rep("retro 5",31),
                               rep("retro 9",31), rep("retro 10",31)),
                     Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),7),
                     Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),7))
  
  dfDRm$Model <- factor(dfDRm$Model, as.character(dfDRm$Model))
  
  
  plot.retrodensitym <- function(df) {
    df <- df %>%
      mutate(p.var = 1.96 * Var) %>%
      mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      geom_point(aes(x=Year,y= Rov), colour="black",size=3)+
      geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
                                            'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
                                          'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      labs(x="Year")+
      labs(y="Density (ind. per square kilometer)")+
      theme(legend.position = c(0.87,0.75))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
  }
  
  
  
  #------------------------------------------------------------------------------
  # SPAWNING BIOMASS GLOBAL RETROSPECTIVE
  #------------------------------------------------------------------------------
  dfSR<-data.frame(Year  = rep(seq(1985,2015,by=1),6), 
                   Density = c(par[190:220],
                               par_R1[187:216],NA,
                               par_R2[184:212],NA,NA,
                               par_R6[172:196],NA,NA,NA,NA,NA,NA,
                               par_R9[163:184],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                               par_R10[160:180],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   Error = c(var[190:220],
                             var_R1[187:216],NA,
                             var_R2[184:212],NA,NA,
                             var_R6[172:196],NA,NA,NA,NA,NA,NA,
                             var_R9[163:184],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                             var_R10[160:180],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),
                             rep("retro 6",31),rep("retro 9",31),
                             rep("retro 10",31)))
  
  dfSR$Model <- factor(dfSR$Model, as.character(dfSR$Model))
  
  plot.retrospawn <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                            'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                          'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      #coord_cartesian(xlim=c(1994,2015))+
      #scale_y_continuous(limits=c(500,6000))+
      labs(x="Year")+
      labs(y="Spawning biomass (tons)")+
      theme(legend.position = c(0.85,0.70))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
    
  }
  

  
  #------------------------------------------------------------------------------
  # SPAWNING BIOMASS GLOBAL RETROSPECTIVE FIXED M
  #------------------------------------------------------------------------------
  dfSRm <- data.frame(Year  = rep(seq(1985,2015,by=1),7), 
                      Density = c(parM[189:219],
                                  par_R1m[186:215],NA,
                                  par_R2m[183:211],NA,NA,
                                  par_R4m[177:203],NA,NA,NA,NA,
                                  par_R5m[174:199],NA,NA,NA,NA,NA,
                                  par_R9m[162:183],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                  par_R10m[159:179],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      Error = c(varM[189:219],
                                  var_R1m[186:215],NA,
                                  var_R2m[183:211],NA,NA,
                                  var_R4m[177:203],NA,NA,NA,NA,
                                  var_R5m[174:199],NA,NA,NA,NA,NA,
                                  var_R9m[162:183],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                  var_R10m[159:179],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),rep("retro 4",31),rep("retro 5",31),
                                rep("retro 9",31), rep("retro 10",31)))
  
  dfSRm$Model <- factor(dfSRm$Model, as.character(dfSRm$Model))
  
  
  plot.retrospawnm <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
                                            'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
                                          'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      labs(x="Year")+
      labs(y="Spawning biomass (tons)")+
      theme(legend.position = c(0.87,0.75))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
  }
  
  
  
  
  
  #------------------------------------------------------------------------------
  # Recruitment GLOBAL RETROSPECTIVE
  #------------------------------------------------------------------------------
  dfRR<-data.frame(Year  = rep(seq(1985,2015,by=1),6), 
                   Density = c(par[159:189],
                               par_R1[157:186],NA,
                               par_R2[155:183],NA,NA,
                               par_R6[147:171],NA,NA,NA,NA,NA,NA,
                               par_R9[141:162],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                               par_R10[139:159],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   Error= c(var[159:189],
                            var_R1[157:186],NA,
                            var_R2[155:183],NA,NA,
                            var_R6[147:171],NA,NA,NA,NA,NA,NA,
                            var_R9[141:162],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                            var_R10[139:159],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),
                             rep("retro 6",31),rep("retro 9",31),
                             rep("retro 10",31)))
  
  
  dfRR$Model <- factor(dfRR$Model, as.character(dfRR$Model))
  
  
  plot.retrorec <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                            'retro 6' = 'green',  'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
                                          'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      labs(x="Year")+
      labs(y="Age 8 recruits (thousands)")+
      theme(legend.position = c(0.8,0.7))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
    
  }
  
  
  #------------------------------------------------------------------------------
  # RECRUITMENT GLOBAL RETROSPECTIVE FIXED M
  #------------------------------------------------------------------------------
  dfRRm <- data.frame(Year  = rep(seq(1985,2015,by=1),7), 
                      Density = c(parM[158:188],
                                  par_R1m[156:185],NA,
                                  par_R2m[154:182],NA,NA,
                                  par_R4m[150:176],NA,NA,NA,NA,
                                  par_R5m[148:173],NA,NA,NA,NA,NA,
                                  par_R9m[140:161],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                  par_R10m[138:158],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      Error= c(varM[158:188],
                                  var_R1m[156:185],NA,
                                  var_R2m[154:182],NA,NA,
                                  var_R4m[150:176],NA,NA,NA,NA,
                                  var_R5m[148:173],NA,NA,NA,NA,NA,
                                  var_R9m[140:161],NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                  var_R10m[138:158],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      Model = c(rep("Full",31),rep("retro 1",31),rep("retro 2",31),rep("retro 4",31),rep("retro 5",31),
                                rep("retro 9",31), rep("retro 10",31)))
  
  dfRRm$Model <- factor(dfRRm$Model, as.character(dfRRm$Model))
  
  
  plot.retrorecm <- function(df) {
    df <- df %>% 
      mutate(err= 1.96 * Error) %>%
      mutate(lower = Density - err, upper = Density + err)
    
    
    ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
      scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
                                            'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
                                          'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
      labs(x="Year")+
      labs(y="Age 8 recruits (thousands)")+
      theme(legend.position = c(0.87,0.75))+
      theme(legend.direction = "vertical")+
      guides(fill=guide_legend(ncol=2))
  }
  
  
  
#------------------------------------------------------------------------------
# SDNR - Global
#------------------------------------------------------------------------------
  # dfsdnr <- data.frame(Year  = rep(yearT,4), 
  #                      Result= sdnr$Result,
  #                      Iteration = sdnr$Iteration,
  #                      Metric = sdnr$Metric)
  # 
  # plot.sample <- function(df) {
  #   ggplot(df, aes(Year,Result,group=Iteration)) + geom_line(aes(Year,Result,colour=Iteration),size=1)+
  #     facet_wrap(~Metric, scales = "free", ncol = 1)+
  #     scale_colour_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                           'Third pass' = 'green'))+
  #     scale_fill_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                         'Third pass' = 'green'))+
  #     labs(x="Year")+
  #     labs(y="")+
  #     theme(legend.margin=unit(0.2, "cm"))+
  #     theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 12))+
  #     theme(legend.direction = "vertical")+
  #     guides(fill=guide_legend(ncol=2))
  # }
  # 
  
  #------------------------------------------------------------------------------
  # SDNR - Fixed M
  # #------------------------------------------------------------------------------
  # dfsdnrM <- data.frame(Year  = rep(yearT,4), 
  #                      Result= sdnrm$Result,
  #                      Iteration = sdnrm$Iteration,
  #                      Metric = sdnrm$Metric)
  # 
  # plot.sample <- function(df) {
  #   ggplot(df, aes(Year,Result,group=Iteration)) + geom_line(aes(Year,Result,colour=Iteration),size=1)+
  #     facet_wrap(~Metric, scales = "free", ncol = 1)+
  #     scale_colour_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                           'Third pass' = 'green'))+
  #     scale_fill_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                         'Third pass' = 'green'))+
  #     labs(x="Year")+
  #     labs(y="")+
  #     theme(legend.margin=unit(0.2, "cm"))+
  #     theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 12))+
  #     theme(legend.direction = "vertical")+
  #     guides(fill=guide_legend(ncol=2))
  # }
  # 
  
  
  #------------------------------------------------------------------------------
  # GELMAN
  #------------------------------------------------------------------------------
  # png(file='../figures/FUCKYOU.png', res=300, width=6, height=4, units ="in", bg="transparent")
  # plot(gr_diag[,1],type="l",ylim=c(1,1.8), xlab="Parameter",ylab = "Gelman-Rubin score")
  # lines(gr_diag[,2],col="red")
  # dev.off()
  # 
  
  
#------------------------------------------------------------------------------
# FUNCTION CALLS
#------------------------------------------------------------------------------
d1 <- plot.parmc(df = parMcmc) 
d2 <- plot.density(df = dfD) 
d3 <- plot.spawners(df = dfS) 
d4 <- plot.recruit(df = dfR) 
d5 <- plot.sel(df = dfSel) 
d6 <- plot.Fmt(df = dfF) 
d7 <- plot.age(dfs = dfs, dfm = dfm)
d8 <- plot.cpue(df = dfcp)
d9 <- plot.cpuei(df = dfcpi)
d10 <- plot.proj(df = dfproj)
d11 <- plot.retrodensity(df = dfDR)
d12 <- plot.retrospawn(df = dfSR)
d13 <- plot.retrorec(df = dfRR)
#d14 <- plot.sample(df = dfsdnr)
#d15 <- plot.sample(df = dfsdnrM)
d16 <- plot.retrodensitym(df = dfDRm)
d17 <- plot.retrospawnm(df = dfSRm)
d18 <- plot.retrorecm(df = dfRRm)
#d19 <- plot.age(df = dfx)

#------------------------------------------------------------------------------
# EXPORT
#------------------------------------------------------------------------------
ggsave(d1,file='figures/theta.png', dpi=300, width=7, height=8, units ="in", bg="transparent") 
ggsave(d2,file='figures/density.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d3,file='figures/spawners.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d4,file='figures/recruit.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d5,file='figures/selectivity.png', dpi=300, width=7, height=4., units ="in", bg="transparent") 
ggsave(d6,file='figures/Fmort.png', dpi=300, width=7, height=4., units ="in", bg="transparent") 
ggsave(d7,file='figures/age_res.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d8,file='figures/cpue.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d9,file='figures/cpuei.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d10,file='figures/spbm_proj.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d11,file='figures/retro_den.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d12,file='figures/retro_spawn.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d13,file='figures/retro_rec.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
#ggsave(d14,file='../../figures/sample.png', dpi=300, width=7, height=4, units ="in", bg="transparent") 
#ggsave(d15,file='../../figures/sample2.png', dpi=300, width=7, height=4, units ="in", bg="transparent") 
ggsave(d16,file='figures/retro_denm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d17,file='figures/retro_spawnm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d18,file='figures/retro_recm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
#ggsave(d19,file='../../figures/age_res_uncorrected.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 










# Ds <- read.delim("../model_output/dens_sdnr.ps",header=FALSE,sep="")
# Dg <- read.delim("../model_output/dens_g.ps",header=FALSE,sep="")
# Do <- read.delim("../model_output/dens_4.ps",header=FALSE,sep="")
# 
# 
# Ds<-Ds[c(10001:25000),]
# Dg<-Dg[c(430358:445357),]
# Do<-Do[c(5001:20000),]
# 
# for(i in 1:31)
# {
#   Ds[,i]<-sort(Ds[,i]) 
#   Dg[,i]<-sort(Dg[,i]) 
#   Do[,i]<-sort(Do[,i]) 
# }
# 
# MenD<-as.matrix(c(colMeans(Ds),colMeans(Dg),colMeans(Do)))
# MinD<-as.matrix(c(Ds[min,],Dg[min,],Do[min,]))
# MaxD<-as.matrix(c(Ds[max,],Dg[max,],Do[max,]))
# 
# dfD <- data.frame(Year  = rep(seq(1985,2015,by=1),3), 
#                   Density = MenD,
#                   MinD  = unlist(MinD),
#                   MaxD  = unlist(MaxD),
#                   Model = c(rep("S",31),rep("G",31),rep("O",31)),
#                   Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),3),
#                   Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),3))
# 
# plot.density <- function(dfD) {
#   dfD <- dfD %>%
#     mutate(p.var = 1.96 * Var) %>%
#     mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
#   
#   
#   ggplot(dfD,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#     geom_ribbon(aes(ymin = MinD, ymax = MaxD, fill=Model),alpha=0.2, show.legend=FALSE) + 
#     geom_point(aes(x=Year,y= Rov), colour="black",size=4)+
#     geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
#     scale_colour_manual(name='', values=c('S' = '#000000', 'G' = '#009E73', 'O' = '#e79f00'))+
#     scale_fill_manual(name='', values=c('S' = '#000000', 'G' = '#009E73', 'O' = '#e79f00'))+
#     labs(x="Year")+
#     labs(y="Density (ind. per square kilometer)")+
#     theme(legend.position = c(0.8,0.8))+
#     theme(legend.direction = "vertical")+
#     #theme(legend.position = c(0.8,0.8))+
#     theme(legend.key.size = unit(1., "cm"))+
#     #ggtitle("Total density (individuals per square kilometer)")+
#     theme(legend.key = element_rect(fill = "white"))
#   
#}

# dev_mcmc <-data.frame(Par  = c(unlist(mcmc_Gx[,11:41]),unlist(mcmc_2x[,11:41])),
#                       Time =   c(rep(c(rep(seq(1,31,by=1),each=15000)),2)),
#                       Model =  c(rep("Global",465000),rep("Uncorrected global",465000)))
# 
# plot.devmc <- function(df) {
#   ggplot(dev_mcmc,
#          aes(Par,group = interaction(Time,Model))) + geom_density(aes(Par, y=..scaled.., fill = Model),alpha=.2)+
#     facet_wrap(~Time,ncol=5)+
#     labs(y="")+
#     scale_fill_manual(name='', values=c('Global' = '#009E73',
#                                         'Uncorrected global' = '#e79f00'))
# }
