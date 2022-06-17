#Analysis of Global SR Data for PCR inclusion in MMRT models using this website 
#Jim Den Uyl
#Created 5-24-22

#Setup####

library(sjPlot)
library(dplyr)
library(lme4)
library(data.table)
library(ggpubr)
library(mice)
library(devtools)
library(bbmle)
library(factoextra)
library(tibble)
library(scales)

setwd("/Users/jdenuyl/Downloads/RFiles")


#PCA####

dat1 <- read.csv("Universal_av(6-1-22).csv")

dat2 <- dat1[complete.cases(dat1$mean.temp), ]

imputed_Data <- mice(dat2)
dat3 <- complete(imputed_Data, 2)

comball <- subset(dat3, select = c(ug_NH4_week_cm2 ,ug_NO3_week_cm2, Inorganic_N,
                                   ug_TP_week_cm2, NP, pH, Soil_Percent_C, Soil_Percent_N, Soil_CN, 
                                   Species_Richness, Total_Plant_Cover_All, VWC, mean.temp, Total_Plant_Cover_All
))

comball$ug_NH4_week_cm2 <- rescale(comball$ug_NH4_week_cm2, to = c(0,100))
comball$ug_NO3_week_cm2 <- rescale(comball$ug_NO3_week_cm2, to = c(0,100))
comball$Inorganic_N <- rescale(comball$Inorganic_N, to = c(0,100))
comball$ug_TP_week_cm2 <- rescale(comball$ug_TP_week_cm2, to = c(0,100))
comball$NP <- rescale(comball$NP, to = c(0,100))
comball$pH <- rescale(comball$pH, to = c(0,100))
comball$Soil_Percent_C <- rescale(comball$Soil_Percent_C, to = c(0,100))
comball$Soil_Percent_N <- rescale(comball$Soil_Percent_N, to = c(0,100))
comball$Soil_CN <- rescale(comball$Soil_CN, to = c(0,100))
comball$Species_Richness <- rescale(comball$Species_Richness, to = c(0,100))
comball$VWC <- rescale(comball$VWC, to = c(0,100))
comball$mean.temp <- rescale(comball$mean.temp, to = c(0,100))
comball$Total_Plant_Cover_All <- rescale(comball$Total_Plant_Cover_All, to = c(0,100))

covar.matrix <- cov(comball)

pca <- prcomp(scale(comball),center = TRUE)
var <- get_pca_var(pca)
cos2.var <- var$cos2 #Imprortance of factors in PCs
ind <- get_pca_ind(pca)
coord.ind <- ind$coord #Coordinates for individuals

#Merge
coord.ind.tab1 <- enframe(coord.ind[,"Dim.1"])
coord.ind.tab2 <- enframe(coord.ind[,"Dim.2"])
coord.ind.tab3 <- enframe(coord.ind[,"Dim.3"])
coord.ind.tab4 <- enframe(coord.ind[,"Dim.4"])
coord.ind.tab5 <- enframe(coord.ind[,"Dim.5"])

dat3$Dim.1 <- coord.ind.tab1$value
dat3$Dim.2 <- coord.ind.tab2$value
dat3$Dim.3 <- coord.ind.tab3$value
dat3$Dim.4 <- coord.ind.tab4$value
dat3$Dim.5 <- coord.ind.tab5$value

#Preparing Data####
comball <- dat3

comball$value<-comball$EFFLUX

comball$Dim.1 <- rescale(comball$Dim.1, to = c(0,100))
comball$Dim.2 <- rescale(comball$Dim.2, to = c(0,100))
comball$Dim.3 <- rescale(comball$Dim.3, to = c(0,100))
comball$Dim.4 <- rescale(comball$Dim.4, to = c(0,100))
comball$Dim.5 <- rescale(comball$Dim.5, to = c(0,100))


comballHtest <- subset(comball, Elevation=='High') # 
comballLtest <- subset(comball, Elevation=='Low') # 
comballwarmie <- subset(comball, Warming=='W')
comballcoolie <- subset(comball, Warming=='A')
highwarm <- subset(comballwarmie, Elevation=='High')
highambient <- subset(comballcoolie, Elevation=='High')
lowwarm <- subset(comballwarmie, Elevation=='Low')
lowambient <- subset(comballcoolie, Elevation=='Low')

mean(comballHtest$mean.EFFLUX)
mean(comballLtest$mean.EFFLUX)
mean(comballwarmie$mean.EFFLUX)
mean(comballcoolie$mean.EFFLUX)
mean(highwarm$mean.EFFLUX)
mean(highambient$mean.EFFLUX)
mean(lowwarm$mean.EFFLUX)
mean(lowambient$mean.EFFLUX)

#Find common temp range

comballe <- subset(comball, mean.temp>5.1)
comballe <- subset(comballe, mean.temp<20)
comballe <- subset(comballe, mean.EFFLUX>0)
comballe$logvalue <- log(comballe$mean.EFFLUX)
comballe$K <- comballe$mean.temp+273.15

comballe$ug_NH4_week_cm2 <- rescale(comballe$ug_NH4_week_cm2, to = c(0,100))
comballe$ug_NO3_week_cm2 <- rescale(comballe$ug_NO3_week_cm2, to = c(0,100))
comballe$Inorganic_N <- rescale(comballe$Inorganic_N, to = c(0,100))
comballe$ug_TP_week_cm2 <- rescale(comballe$ug_TP_week_cm2, to = c(0,100))
comballe$NP <- rescale(comballe$NP, to = c(0,100))
comballe$pH <- rescale(comballe$pH, to = c(0,100))
comballe$Soil_Percent_C <- rescale(comballe$Soil_Percent_C, to = c(0,100))
comballe$Soil_Percent_N <- rescale(comballe$Soil_Percent_N, to = c(0,100))
comballe$Soil_CN <- rescale(comballe$Soil_CN, to = c(0,100))
comballe$Species_Richness <- rescale(comballe$Species_Richness, to = c(0,100))
comballe$VWC <- rescale(comballe$VWC, to = c(0,100))
comballe$mean.temp <- rescale(comballe$mean.temp, to = c(0,100))
comballe$Total_Plant_Cover_All <- rescale(comballe$Total_Plant_Cover_All, to = c(0,100))


comballt <- subset(comball, mean.temp>1.49)
comballt <- subset(comballt, mean.temp<22)
comballt <- subset(comballt, mean.EFFLUX>0)
comballt$logvalue <- log(comballt$mean.EFFLUX)
comballt$K <- comballt$mean.temp+273.15


comballhigh <- subset(comballe, Elevation=='High')
comballlow <- subset(comballe, Elevation=='Low')
comballhighw <- subset(comballhigh, Warming=='W')
comballhigha <- subset(comballhigh, Warming=='A')
comballloww <- subset(comballlow, Warming=='W')
comballlowa <- subset(comballlow, Warming=='A')
comballr <- subset(comballe, Removal=='R')
comballc <- subset(comballe, Removal=='C')

#Create Plain MMRT models####


fun <- function(K,H,C,S){log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
    ((H+C*(K-298))/(8.31446261815324*K))+
    ((S+C*(log(K)-log(298)))/8.31446261815324)}

model.E <- nls(logvalue ~ fun(K,H,C,S), 
               data=comballe,
               start = list(H = 70182.62, C = -12.37, S = 1))

model.E.high <- nls(logvalue ~ fun(K,H,C,S), 
                    data=comballhigh,
                    start = list(H = 70182.62, C = -12.37, S = 1))

model.E.low <- nls(logvalue ~ fun(K,H,C,S), 
                   data=comballlow,
                   start = list(H = 70182.62, C = -12.37, S = 1))

model.E.r <- nls(logvalue ~ fun(K,H,C,S), 
                 data=comballr,
                 start = list(H = 70182.62, C = -12.37, S = 1))

model.E.c <- nls(logvalue ~ fun(K,H,C,S), 
                 data=comballc,
                 start = list(H = 70182.62, C = -12.37, S = 1))

model.T <- nls(logvalue ~ fun(K,H,C,S), 
               data=comballt,
               start = list(H = 70182.62, C = -12.37, S = 1))

model.E.highw <- nls(logvalue ~ fun(K,H,C,S), 
                     data=comballhighw,
                     start = list(H = 70182.62, C = -12.37, S = 1))

model.E.loww <- nls(logvalue ~ fun(K,H,C,S), 
                    data=comballloww,
                    start = list(H = 70182.62, C = -12.37, S = 1))

model.E.higha <- nls(logvalue ~ fun(K,H,C,S), 
                     data=comballhigha,
                     start = list(H = 70182.62, C = -12.37, S = 1))

model.E.lowa <- nls(logvalue ~ fun(K,H,C,S), 
                    data=comballlowa,
                    start = list(H = 352097.7, C = 33551.1, S = 959.5),
                    control = nls.control(minFactor = 1/6999999999,
                                          maxiter = 199999))

summary(model.E)
summary(model.E.high)
summary(model.E.low)
summary(model.E.r)
summary(model.E.c)
summary(model.E.highw)
summary(model.E.loww)
summary(model.E.higha)
summary(model.E.lowa)
summary(model.T) #Use H,C,S estimates

fun.gg <- function(K,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                           ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                           ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+b}

model.E <- nls(logvalue ~ fun.gg(K,b), 
               data=comballe,
               start = list(b = .5))

comballe$pred <- predict(model.E, newdata = comballe)

ggplot() + 
  geom_point(data = comballe,
             aes(x = K,
                 y = logvalue),
             color = "blue") +
  geom_line(data= comballe,
            aes(x = K, 
                y = pred), 
            color = "blue")  +
  ggtitle("Global Soil Respiration Sensitivity")+
  xlab("Soil Temperature (K)")+
  ylab("Natural Log of Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##

comballe$pred.r <- predict(model.E.r, newdata = comballe)
comballe$pred.c <- predict(model.E.c, newdata = comballe)

ggplot() + 
  geom_point(data = comballe,
             aes(x = K,
                 y = logvalue)) +
  geom_line(data= comballe,
            aes(x = K, 
                y = pred.r),
            color = "green") + 
  geom_line(data= comballe,
            aes(x = K, 
                y = pred.c),
            color = "orange")  +
  ggtitle("Global Soil Respiration Sensitivity")+
  xlab("Soil Temperature (K)")+
  ylab("Natural Log of Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#Create MMRT * Dim.1 models####

fun1 <- function(K,a,Dim.1,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Dim.1)+b}

model.1.E <- nls(logvalue ~ fun1(K,a,Dim.1,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))



summary(model.1.E)

comballe$pred1 <- predict(model.1.E, newdata = comballe)

ggplot() + 
  geom_point(data = comballe,
             aes(x = K,
                 y = logvalue),
             color = "blue") +
  geom_line(data= comballe,
            aes(x = K, 
                y = pred1), 
            color = "blue")  +
  ggtitle("Global Soil Respiration Sensitivity")+
  xlab("Soil Temperature (K)")+
  ylab("Natural Log of Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Create MMRT * Dim.2 models####

fun2 <- function(K,a,Dim.2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Dim.2)+b}

model.2.E <- nls(logvalue ~ fun2(K,a,Dim.2,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))

summary(model.2.E)


#Create MMRT * Dim.3 models####

fun3 <- function(K,a,Dim.3,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Dim.3)+b}

model.3.E <- nls(logvalue ~ fun3(K,a,Dim.3,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))



summary(model.3.E)

#Create MMRT * Dim.4 models####

fun4 <- function(K,a,Dim.4,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Dim.4)+b}

model.4.E <- nls(logvalue ~ fun4(K,a,Dim.4,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))


summary(model.4.E)

#Create MMRT * Dim.5 models####

fun5 <- function(K,a,Dim.5,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Dim.5)+b}

model.5.E <- nls(logvalue ~ fun5(K,a,Dim.5,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))


summary(model.5.E)


#Create MMRT * pH models####

fun6 <- function(K,a,pH,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                              ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                              ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*pH)+b}

model.p.E <- nls(logvalue ~ fun6(K,a,pH,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))


summary(model.p.E)


#Create MMRT * CN models####

fun7 <- function(K,a,Soil_CN,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                   ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                   ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Soil_CN)+b}

model.cn.E <- nls(logvalue ~ fun7(K,a,Soil_CN,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))


summary(model.cn.E)


#Create MMRT * Total_Cover_Scaled models####

fun8 <- function(K,a,Total_Cover_Scaled,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                              ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                              ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Total_Cover_Scaled)+b}

model.t.E <- nls(logvalue ~ fun8(K,a,Total_Cover_Scaled,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))



summary(model.t.E)


#Create MMRT * Soil Nitrogen models####

fun9 <- function(K,a,Soil_Percent_N,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                          ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                          ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Soil_Percent_N)+b}

model.soiln.E <- nls(logvalue ~ fun9(K,a,Soil_Percent_N,b), 
                     data=comballe,
                     start = list(a=1, b=.5),
                     control = nls.control(minFactor = 1/6999999999,
                                           maxiter = 199999))

summary(model.soiln.E)

fun.N <- function(K,H,C,S,Soil_Percent_N){log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
    ((H+C*(K-298))/(8.31446261815324*K))+
    ((S+C*(log(K)-log(298)))/8.31446261815324)*(1.95071*Soil_Percent_N)+0.05377}

model.soiln.all <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                       data=comballe,
                       start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.high <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                        data=comballhigh,
                        start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.low <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                       data=comballlow,
                       start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.highw <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                         data=comballhighw,
                         start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.loww <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                        data=comballloww,
                        start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.higha <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                         data=comballhigha,
                         start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.lowa <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                        data=comballlowa,
                        start = list(H = 70182.62, C = -12.37, S = 1))

summary(model.soiln.all)
summary(model.soiln.high)
summary(model.soiln.low)
summary(model.soiln.highw)
summary(model.soiln.loww)
summary(model.soiln.higha)
summary(model.soiln.lowa)


#Create MMRT * VWC models####

fun10 <- function(K,a,VWC,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*VWC)+b}

model.m.E <- nls(logvalue ~ fun10(K,a,VWC,b), 
                 data=comballe,
                 start = list(a=1, b=.5),
                 control = nls.control(minFactor = 1/6999999999,
                                       maxiter = 199999))




summary(model.m.E)


#Create MMRT * Soil Carbon models####

fun11 <- function(K,a,Soil_Percent_C,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                           ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                           ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Soil_Percent_C)+b}

model.soilc.E <- nls(logvalue ~ fun11(K,a,Soil_Percent_C,b), 
                     data=comballe,
                     start = list(a=1, b=.5),
                     control = nls.control(minFactor = 1/6999999999,
                                           maxiter = 199999))



summary(model.soilc.E)



#Create MMRT * NH4 models####

fun12 <- function(K,a,ug_NH4_week_cm2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                            ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                            ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*ug_NH4_week_cm2)+b}

model.nh4.E <- nls(logvalue ~ fun12(K,a,ug_NH4_week_cm2,b), 
                   data=comballe,
                   start = list(a=1, b=.5),
                   control = nls.control(minFactor = 1/6999999999,
                                         maxiter = 199999))


summary(model.nh4.E)



#Create MMRT * NO3 models####

fun13 <- function(K,a,ug_NO3_week_cm2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                            ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                            ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*ug_NO3_week_cm2)+b}

model.no3.E <- nls(logvalue ~ fun13(K,a,ug_NO3_week_cm2,b), 
                   data=comballe,
                   start = list(a=1, b=.5),
                   control = nls.control(minFactor = 1/6999999999,
                                         maxiter = 199999))


summary(model.no3.E)




#Create MMRT * Species Richness models####

fun14 <- function(K,a,Species_Richness,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                             ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                             ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*Species_Richness)+b}

model.spr.E <- nls(logvalue ~ fun14(K,a,Species_Richness,b), 
                   data=comballe,
                   start = list(a=1, b=.5),
                   control = nls.control(minFactor = 1/6999999999,
                                         maxiter = 199999))


summary(model.spr.E)



#Create MMRT * NP models####

fun15 <- function(K,a,NP,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                               ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                               ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(a*NP)+b}

model.np.E <- nls(logvalue ~ fun15(K,a,NP,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))


summary(model.np.E)




#Create MMRT + Dim.1 models####

fun1 <- function(K,a,Dim.1,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Dim.1)+b}

model.1.Eb <- nls(logvalue ~ fun1(K,a,Dim.1,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))



summary(model.1.Eb)

comballe$pred1 <- predict(model.1.E, newdata = comballe)

ggplot() + 
  geom_point(data = comballe,
             aes(x = K,
                 y = logvalue),
             color = "blue") +
  geom_line(data= comballe,
            aes(x = K, 
                y = pred1), 
            color = "blue")  +
  ggtitle("Global Soil Respiration Sensitivity")+
  xlab("Soil Temperature (K)")+
  ylab("Natural Log of Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Create MMRT + Dim.2 models####

fun2 <- function(K,a,Dim.2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Dim.2)+b}

model.2.Eb <- nls(logvalue ~ fun2(K,a,Dim.2,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))

summary(model.2.Eb)


#Create MMRT + Dim.3 models####

fun3 <- function(K,a,Dim.3,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Dim.3)+b}

model.3.Eb <- nls(logvalue ~ fun3(K,a,Dim.3,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))



summary(model.3.Eb)

#Create MMRT + Dim.4 models####

fun4 <- function(K,a,Dim.4,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Dim.4)+b}

model.4.Eb <- nls(logvalue ~ fun4(K,a,Dim.4,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))


summary(model.4.Eb)

#Create MMRT + Dim.5 models####

fun5 <- function(K,a,Dim.5,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                 ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                 ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Dim.5)+b}

model.5.Eb <- nls(logvalue ~ fun5(K,a,Dim.5,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))


summary(model.5.Eb)


#Create MMRT + pH models####

fun6 <- function(K,a,pH,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                              ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                              ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*pH)+b}

model.p.Eb <- nls(logvalue ~ fun6(K,a,pH,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))


summary(model.p.Eb)


#Create MMRT + CN models####

fun7 <- function(K,a,Soil_CN,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                   ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                   ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Soil_CN)+b}

model.cn.Eb <- nls(logvalue ~ fun7(K,a,Soil_CN,b), 
                   data=comballe,
                   start = list(a=1, b=.5),
                   control = nls.control(minFactor = 1/6999999999,
                                         maxiter = 199999))


summary(model.cn.Eb)


#Create MMRT + Total_Cover_Scaled models####

fun8 <- function(K,a,Total_Cover_Scaled,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                              ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                              ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Total_Cover_Scaled)+b}

model.t.Eb <- nls(logvalue ~ fun8(K,a,Total_Cover_Scaled,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))



summary(model.t.Eb)


#Create MMRT + Soil Nitrogen models####

fun9 <- function(K,a,Soil_Percent_N,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                          ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                          ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Soil_Percent_N)+b}

model.soiln.Eb <- nls(logvalue ~ fun9(K,a,Soil_Percent_N,b), 
                      data=comballe,
                      start = list(a=1, b=.5),
                      control = nls.control(minFactor = 1/6999999999,
                                            maxiter = 199999))

summary(model.soiln.Eb)

fun.N <- function(K,H,C,S,Soil_Percent_N){log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
    ((H+C*(K-298))/(8.31446261815324*K))+
    ((S+C*(log(K)-log(298)))/8.31446261815324)+(1.84918*Soil_Percent_N)+0.05746}

model.soiln.highb <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                         data=comballhigh,
                         start = list(H = 70182.62, C = -12.37, S = 1))

model.soiln.lowb <- nls(logvalue ~ fun.N(K,H,C,S,Soil_Percent_N), 
                        data=comballlow,
                        start = list(H = 70182.62, C = -12.37, S = 1))

summary(model.soiln.highb)
summary(model.soiln.lowb)



#Create MMRT + VWC models####

fun10 <- function(K,a,VWC,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*VWC)+b}

model.m.Eb <- nls(logvalue ~ fun10(K,a,VWC,b), 
                  data=comballe,
                  start = list(a=1, b=.5),
                  control = nls.control(minFactor = 1/6999999999,
                                        maxiter = 199999))




summary(model.m.Eb)


#Create MMRT + Soil Carbon models####

fun11 <- function(K,a,Soil_Percent_C,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                           ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                           ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Soil_Percent_C)+b}

model.soilc.Eb <- nls(logvalue ~ fun11(K,a,Soil_Percent_C,b), 
                      data=comballe,
                      start = list(a=1, b=.5),
                      control = nls.control(minFactor = 1/6999999999,
                                            maxiter = 199999))



summary(model.soilc.Eb)



#Create MMRT + NH4 models####

fun12 <- function(K,a,ug_NH4_week_cm2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                            ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                            ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*ug_NH4_week_cm2)+b}

model.nh4.Eb <- nls(logvalue ~ fun12(K,a,ug_NH4_week_cm2,b), 
                    data=comballe,
                    start = list(a=1, b=.5),
                    control = nls.control(minFactor = 1/6999999999,
                                          maxiter = 199999))


summary(model.nh4.Eb)



#Create MMRT + NO3 models####

fun13 <- function(K,a,ug_NO3_week_cm2,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                            ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                            ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*ug_NO3_week_cm2)+b}

model.no3.Eb <- nls(logvalue ~ fun13(K,a,ug_NO3_week_cm2,b), 
                    data=comballe,
                    start = list(a=1, b=.5),
                    control = nls.control(minFactor = 1/6999999999,
                                          maxiter = 199999))


summary(model.no3.Eb)




#Create MMRT + Species Richness models####

fun14 <- function(K,a,Species_Richness,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                             ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                             ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*Species_Richness)+b}

model.spr.Eb <- nls(logvalue ~ fun14(K,a,Species_Richness,b), 
                    data=comballe,
                    start = list(a=1, b=.5),
                    control = nls.control(minFactor = 1/6999999999,
                                          maxiter = 199999))


summary(model.spr.Eb)



#Create MMRT + NP models####

fun15 <- function(K,a,NP,b){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                               ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                               ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+(a*NP)+b}

model.np.Eb <- nls(logvalue ~ fun15(K,a,NP,b), 
                   data=comballe,
                   start = list(a=1, b=.5),
                   control = nls.control(minFactor = 1/6999999999,
                                         maxiter = 199999))


summary(model.np.Eb)




#Create Non-MMRT models for elevation data####

mod.temp <- lm(mean.EFFLUX ~ K, data = comballe)
mod.pH <- lm(mean.EFFLUX ~ pH, data = comballe)
mod.CN <- lm(mean.EFFLUX ~ Soil_CN, data = comballe)
mod.TCS <- lm(mean.EFFLUX ~ Total_Cover_Scaled, data = comballe)
mod.SN <- lm(mean.EFFLUX ~ Soil_Percent_N, data = comballe)
mod.VWC <- lm(mean.EFFLUX ~ VWC, data = comballe)
mod.SC <- lm(mean.EFFLUX ~ Soil_Percent_C, data = comballe)
mod.NH4 <- lm(mean.EFFLUX ~ ug_NH4_week_cm2, data = comballe)
mod.NO3 <- lm(mean.EFFLUX ~ ug_NO3_week_cm2, data = comballe)
mod.SPR <- lm(mean.EFFLUX ~ Species_Richness, data = comballe)
mod.NP <- lm(mean.EFFLUX ~ NP, data = comballe)

summary(mod.temp)
summary(mod.pH)
summary(mod.CN)
summary(mod.TCS)
summary(mod.SN)
summary(mod.VWC)
summary(mod.SC)
summary(mod.NH4)
summary(mod.NO3)
summary(mod.SPR)
summary(mod.NP)

#Elevation check
mod.temp.high <- lm(mean.EFFLUX ~ K, data = comballhigh)
mod.pH.high <- lm(mean.EFFLUX ~ pH, data = comballhigh)
mod.CN.high <- lm(mean.EFFLUX ~ Soil_CN, data = comballhigh)
mod.TCS.high <- lm(mean.EFFLUX ~ Total_Cover_Scaled, data = comballhigh)
mod.SN.high <- lm(mean.EFFLUX ~ Soil_Percent_N, data = comballhigh)
mod.VWC.high <- lm(mean.EFFLUX ~ VWC, data = comballhigh)
mod.SC.high <- lm(mean.EFFLUX ~ Soil_Percent_C, data = comballhigh)
mod.NH4.high <- lm(mean.EFFLUX ~ ug_NH4_week_cm2, data = comballhigh)
mod.NO3.high <- lm(mean.EFFLUX ~ ug_NO3_week_cm2, data = comballhigh)
mod.SPR.high <- lm(mean.EFFLUX ~ Species_Richness, data = comballhigh)
mod.NP.high <- lm(mean.EFFLUX ~ NP, data = comballhigh)

mod.temp.low <- lm(mean.EFFLUX ~ K, data = comballlow)
mod.pH.low <- lm(mean.EFFLUX ~ pH, data = comballlow)
mod.CN.low <- lm(mean.EFFLUX ~ Soil_CN, data = comballlow)
mod.TCS.low <- lm(mean.EFFLUX ~ Total_Cover_Scaled, data = comballlow)
mod.SN.low <- lm(mean.EFFLUX ~ Soil_Percent_N, data = comballlow)
mod.VWC.low <- lm(mean.EFFLUX ~ VWC, data = comballlow)
mod.SC.low <- lm(mean.EFFLUX ~ Soil_Percent_C, data = comballlow)
mod.NH4.low <- lm(mean.EFFLUX ~ ug_NH4_week_cm2, data = comballlow)
mod.NO3.low <- lm(mean.EFFLUX ~ ug_NO3_week_cm2, data = comballlow)
mod.SPR.low <- lm(mean.EFFLUX ~ Species_Richness, data = comballlow)
mod.NP.low <- lm(mean.EFFLUX ~ NP, data = comballlow)


#Model Selection####

Cand.models <- list("MMRT" = model.E,
                    "MMRT * pH" = model.p.E,
                    "MMRT * CN" = model.cn.E,
                    "MMRT * Total Veg Cover" = model.t.E,
                    "MMRT * Soil N" = model.soiln.E,
                    "MMRT * Soil C" = model.soilc.E,
                    "MMRT * VWC" = model.m.E,
                    "MMRT * Species Richness" = model.spr.E,
                    "MMRT * NP" = model.np.E,
                    "MMRT * NH4" = model.nh4.E,
                    "MMRT * NO3" = model.no3.E,
                    "MMRT + pH" = model.p.Eb,
                    "MMRT + CN" = model.cn.Eb,
                    "MMRT + Total Veg Cover" = model.t.Eb,
                    "MMRT + Soil N" = model.soiln.Eb,
                    "MMRT + Soil C" = model.soilc.Eb,
                    "MMRT + VWC" = model.m.Eb,
                    "MMRT + Species Richness" = model.spr.Eb,
                    "MMRT + NP" = model.np.Eb,
                    "MMRT + NH4" = model.nh4.Eb,
                    "MMRT + NO3" = model.no3.Eb,
                    "Temperature" = mod.temp,
                    "pH" = mod.pH,
                    "Total Veg Cover" = mod.TCS,
                    "CN" = mod.CN,
                    "Soil N" = mod.SN,
                    "Soil C" = mod.SC,
                    "Soil VWC" = mod.VWC,
                    "Species Richness" = mod.SPR,
                    "NP" = mod.NP,
                    "NO3" = mod.NO3,
                    "NH4" = mod.NH4)

Cand.models <- list("MMRT" = model.E,
                    "MMRT * Dim.1" = model.1.E,
                    "MMRT * Dim.2" = model.2.E,
                    "MMRT * Dim.3" = model.3.E,
                    "MMRT * Dim.4" = model.4.E,
                    "MMRT * Dim.5" = model.5.E,
                    "MMRT * pH" = model.p.E,
                    "MMRT * CN" = model.cn.E,
                    "MMRT * Total Veg Cover" = model.t.E,
                    "MMRT * Soil N" = model.soiln.E,
                    "MMRT * Soil C" = model.soilc.E,
                    "MMRT * VWC" = model.m.E,
                    "MMRT * Species Richness" = model.spr.E,
                    "MMRT * NP" = model.np.E,
                    "MMRT * NH4" = model.nh4.E,
                    "MMRT * NO3" = model.no3.E,
                    "MMRT + Dim.1" = model.1.Eb,
                    "MMRT + Dim.2" = model.2.Eb,
                    "MMRT + Dim.3" = model.3.Eb,
                    "MMRT + Dim.4" = model.4.Eb,
                    "MMRT + Dim.5" = model.5.Eb,
                    "MMRT + pH" = model.p.Eb,
                    "MMRT + CN" = model.cn.Eb,
                    "MMRT + Total Veg Cover" = model.t.Eb,
                    "MMRT + Soil N" = model.soiln.Eb,
                    "MMRT + Soil C" = model.soilc.Eb,
                    "MMRT + VWC" = model.m.Eb,
                    "MMRT + Species Richness" = model.spr.Eb,
                    "MMRT + NP" = model.np.Eb,
                    "MMRT + NH4" = model.nh4.Eb,
                    "MMRT + NO3" = model.no3.Eb,
                    "Temperature" = mod.temp,
                    "pH" = mod.pH,
                    "Total Veg Cover" = mod.TCS,
                    "CN" = mod.CN,
                    "Soil N" = mod.SN,
                    "Soil C" = mod.SC,
                    "Soil VWC" = mod.VWC,
                    "Species Richness" = mod.SPR,
                    "NP" = mod.NP,
                    "NO3" = mod.NO3,
                    "NH4" = mod.NH4)

Cand.models <- list("MMRT" = model.E,
                    "MMRT + pH" = model.p.E,
                    "MMRT + CN" = model.cn.E,
                    "MMRT + Total Veg Cover" = model.t.E,
                    "MMRT + Soil N" = model.soiln.E,
                    "MMRT + Soil C" = model.soilc.E,
                    "MMRT + VWC" = model.m.E,
                    "MMRT + Species Richness" = model.spr.E,
                    "MMRT + NP" = model.np.E,
                    "MMRT + NH4" = model.nh4.E,
                    "MMRT + NO3" = model.no3.E,
                    "pH" = mod.pH,
                    "Total Veg Cover" = mod.TCS,
                    "CN" = mod.CN,
                    "Soil N" = mod.SN,
                    "Soil C" = mod.SC,
                    "Soil VWC" = mod.VWC,
                    "Species Richness" = mod.SPR,
                    "NP" = mod.NP,
                    "NO3" = mod.NO3,
                    "NH4" = mod.NH4)

Cand.models <- list("MMRT" = model.E,
                    "Temperature" = mod.temp,
                    "pH" = mod.pH,
                    "Total Veg Cover" = mod.TCS,
                    "CN" = mod.CN,
                    "Soil N" = mod.SN,
                    "Soil C" = mod.SC,
                    "Soil VWC" = mod.VWC,
                    "Species Richness" = mod.SPR,
                    "NP" = mod.NP,
                    "NO3" = mod.NO3,
                    "NH4" = mod.NH4)


Cand.models <- list("MMRT" = model.E.high,
                    "Temperature" = mod.temp.high,
                    "pH" = mod.pH.high,
                    "Total Veg Cover" = mod.TCS.high,
                    "CN" = mod.CN.high,
                    "Soil N" = mod.SN.high,
                    "Soil C" = mod.SC.high,
                    "Soil VWC" = mod.VWC.high,
                    "Species Richness" = mod.SPR.high,
                    "NP" = mod.NP.high,
                    "NO3" = mod.NO3.high,
                    "NH4" = mod.NH4.high)

Cand.models <- list("MMRT" = model.E.low,
                    "Temperature" = mod.temp.low,
                    "pH" = mod.pH.high,
                    "Total Veg Cover" = mod.TCS.low,
                    "CN" = mod.CN.low,
                    "Soil N" = mod.SN.low,
                    "Soil C" = mod.SC.low,
                    "Soil VWC" = mod.VWC.low,
                    "Species Richness" = mod.SPR.low,
                    "NP" = mod.NP.low,
                    "NO3" = mod.NO3.low,
                    "NH4" = mod.NH4.low)

AICctab(Cand.models,
        base=T, weights=T, logLik=T)

AIC(model.1.E,model.2.E,model.3.E,model.4.E,model.5.E, model.p.E, model.cn.E, model.t.E, model.soilc.E, model.soiln.E,
    model.m.E, model.spr.E, model.nh4.E, model.no3.E, model.np.E, mod.pH, mod.TCS, mod.CN, mod.SN, mod.SC, mod.VWC, mod.SPR, mod.NP, mod.NO3, mod.NH4)

AIC(model.soiln.higha, model.soiln.lowa,model.soiln.highw, model.soiln.loww)

mean(comballe$Soil_Percent_N)
mean(comballe$Soil_Percent_C)
mean(comballe$ug_NH4_week_cm2)
mean(comballe$ug_NO3_week_cm2)
mean(comballe$NP)

#Graph####

#Extras
summary(model.soiln.E) #Selected as best model
summary(model.E) #Selected as best model
high <- subset(comballe, Elevation=="High") # 
low <- subset(comballe, Elevation=="Low") # 
mean(high$Soil_Percent_N)
mean(low$Soil_Percent_N)

#Functions to graph
fun.dim2.veryhigh <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                    ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                    ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(1.5)+0.11534}
fun.dim2.high <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(1.2)+0.11534}
fun.dim2.verylow <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                                   ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                                   ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(.25)+0.11534}
fun.dim2.low <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                               ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                               ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(.48)+0.11534}
funMMRT <- function(K){exp((log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                              ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                              ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))*(1)+0.11534)}
funn <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                       ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                       ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))}
funp1 <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                        ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                        ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+.2}
funp2 <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                        ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                        ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))+.4}
funm1 <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                        ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                        ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))-.2}
funm2 <- function(K){(log(((1.38064852*(10^-23))*K)/(6.62607004*(10^-34)))-
                        ((-194865+-15455.9*(K-298))/(8.31446261815324*K))+
                        ((-908.5+-15455.9*(log(K)-log(298)))/8.31446261815324))-.4}


#Not Natural Log - interact
curve(exp(fun(x)), from=278, to=300, col=2,ylim=c(-0,2.5), xlab = "Soil Temperature (K)", ylab = "Soil Respiration (µmol CO2 m-2 s-1)", lwd=3)
curve(exp(fun.dim2.veryhigh(x)), from=278, to=300, add=T)
curve(exp(fun.dim2.high(x)), from=278, to=300, add=T)
curve(exp(fun.dim2.verylow(x)), from=278, to=300, add=T)
curve(exp(fun.dim2.low(x)), from=278, to=300, add=T)

#Not Natural Log - addition
curve(exp(funn(x)), from=278, to=300, col=2,ylim=c(-0,2.5), xlab = "Soil Temperature (K)", ylab = "Soil Respiration (µmol CO2 m-2 s-1)", lwd=3)
curve(exp(funp1(x)), from=278, to=300, add=T)
curve(exp(funp2(x)), from=278, to=300, add=T)
curve(exp(funm1(x)), from=278, to=300, add=T)
curve(exp(funm2(x)), from=278, to=300, add=T)

#Natural log
curve(fun(x), from=278, to=300, main="SoilN=0", ylim=c(-0.5,1.2), col=2, lwd=3, xlab = "Soil Temperature (K)", ylab = "Natural Log of Soil Respiration (µmol CO2 m-2 s-1)")
curve(fun.dim2.80(x), add=T)
curve(fun.dim2.50(x), add=T)
curve(fun.dim2.30(x), add=T)
curve(fun.dim2.15(x), add=T)

#ggplot

mod.1 <- nls(logvalue ~ fun.dim2.veryhigh(K), 
             data=comballe,
             start = list(K=290),
             control = nls.control(minFactor = 1/699999999999999999999,
                                   maxiter = 199999))

comballe$pred1 <- exp(predict(mod.1, newdata = comballe))

mod.2 <- nls(logvalue ~ fun.dim2.high(K), 
             data=comballe,
             start = list(K=1),
             control = nls.control(minFactor = 1/6999999999,
                                   maxiter = 199999))

comballe$pred2 <- exp(predict(mod.2, newdata = comballe))

#mod.3 <- nls(logvalue ~ fun.dim2.verylow(K), 
data=comballe,
start = list(K=285),
control = nls.control(minFactor = 1/699999999999999999,
                      maxiter = 199999))

#comballe$pred3 <- exp(predict(mod.3, newdata = comballe))

mod.4 <- nls(logvalue ~ fun.dim2.low(K), 
             data=comballe,
             start = list(K=290),
             control = nls.control(minFactor = 1/999999999999999999999,
                                   maxiter = 199999))

comballe$pred4 <- exp(predict(mod.4, newdata = comballe))

mod.high <- nls(logvalue ~ fun.dim2.high(K), 
                data=comballe,
                start = list(K=1),
                control = nls.control(minFactor = 1/6999999999,
                                      maxiter = 199999))

comballe$predhigh <- exp(predict(mod.high, newdata = comballe))

mod.low <- nls(logvalue ~ fun.dim2.low(K), 
               data=comballe,
               start = list(K=1),
               control = nls.control(minFactor = 1/6999999999,
                                     maxiter = 199999))

comballe$predlow <- exp(predict(mod.low, newdata = comballe))

#Split data between categories for SoilN

soilNhigh <- subset(comballe, Soil_Percent_N>40) # 
soilNmid <- subset(comballe, Soil_Percent_N>30) # 
soilNmid <- subset(comballe, Soil_Percent_N<40) # 
soilNlow <- subset(comballe, Soil_Percent_N<30) # 


#Plot - select SoilN values
ggplot() + 
  geom_point(data = soilNhigh,
             aes(x = K,
                 y = exp(logvalue)),
             color = "orange") +
  geom_point(data = soilNmid,
             aes(x = K,
                 y = exp(logvalue)),
             color = "green") +
  geom_point(data = soilNlow,
             aes(x = K,
                 y = exp(logvalue)),
             color = "blue") +
  geom_function(fun = funMMRT, colour = "red")  +
  geom_line(data = comballe,
            aes(x = K, 
                y = pred1))+
  geom_line(data = comballe,
            aes(x = K, 
                y = pred2))+
  geom_line(data = comballe,
            aes(x = K, 
                y = pred4))+
  ggtitle("Global Soil Respiration Sensitivity to Warming and Soil Nitrogen")+
  xlab("Soil Temperature (K)")+
  ylab("Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#Plot - High and low SoilN average values
ggplot() + 
  geom_point(data = comballe,
             aes(x = K,
                 y = exp(logvalue)),
             color = "blue") +
  geom_function(fun = fun, colour = "red")  +
  geom_line(data = comballe,
            aes(x = K, 
                y = predhigh))+
  geom_line(data = comballe,
            aes(x = K, 
                y = predlow))+
  ggtitle("Global Soil Respiration Sensitivity to Warming and Soil Nitrogen by Elevation")+
  xlab("Soil Temperature (K)")+
  ylab("Soil Respiration (µmol CO2 m-2 s-1)")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




