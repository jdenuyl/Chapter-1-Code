#Create Publishable Correlation Matrix
#Date: 6-16-22

install.packages("Hmisc")

#library(Hmisc)

setwd("/Users/jdenuyl/Downloads/RFiles")

#Setup####

dat1 <- read.csv("Universal_av(6-1-22).csv")
 
dat2 <- dat1[complete.cases(dat1$mean.temp), ]

imputed_Data <- mice(dat2)
dat3 <- complete(imputed_Data, 2)

comball <- subset(dat3, select = c(ug_NH4_week_cm2 ,ug_NO3_week_cm2, Inorganic_N,
                                   ug_TP_week_cm2, NP, pH, Soil_Percent_C, Soil_Percent_N, Soil_CN, 
                                   Species_Richness, Total_Plant_Cover_All, VWC, mean.temp
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


#Exploring Data####

mcor<-round(cor(comball),3)
mcor

install.packages('htmltools')

library(xtable)
library("corrplot")                             
corrplot(cor(comball), method = "circle")

#
corrplot(cor(comball), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#

library(tidyverse)  
library(corrr)

res.cor <- correlate(comball)
res.cor

res.cor %>%
  shave() %>% 
  stretch(na.rm = TRUE) %>% 
  ggplot(aes(r)) +
  geom_histogram(bins = 10)

res.cor %>% rplot()

plot <- res.cor %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>% 
  rplot(shape = 15, colours = c("red", "green"))

plot + theme(axis.text.x = element_text(angle = 60, hjust = 1))


res.cor %>% network_plot(min_cor = .4)







