# Statistical validation analyses on the data sets and results of "Remote sensing analysis on primary productivity and forest cover dynamics: a Western Ghats case study"

## Multiple Time Series Analysis of Dry Matter Productivity in the Indian Western Ghats region
## Periods are always Jan 10th for every available year of the data products (Western Ghats dry season)
## Multiple Time Series Analysis on DMP, FAPAR, NDVI V2 1Km spatial scale from VITO catalog 1999-2020 of Central Western Ghats region

# Used R required packages

library(ncdf4)
library(raster)
library(RStoolbox)

#________________________________________________________________________________________________
# DMP import the dataset

setwd("D:/toy/uppangala/cop_dmp_1km/")               

rlist <- list.files (pattern ="DMP")
rlist                                                
import <- lapply(rlist,raster)
dmp.multi <- stack(import) 

# central Western Ghats Region
ext <- c(75.25, 75.55, 12.15, 12.45)                 
dmp.wg <- crop(dmp.multi,ext)

names(dmp.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009",
                   "Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")   


#________________________________________________________________________________________________
# FAPAR import of the dataset

setwd("D:/toy/uppangala/cop_fapar_1km/")                   

rlist2 <- list.files (pattern ="FAPAR")
rlist2                                                     
import2 <- lapply(rlist2,raster)
fapar.multi <- stack(import2)

# Crop the products for my AOI
fapar.wg <- crop(fapar.multi,ext)

names(fapar.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009",
                     "Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")


#________________________________________________________________________________________________
# NDVI 

setwd("D:/toy/uppangala/cop_ndvi_1km/")                    

rlist3 <- list.files (pattern ="NDVI")
rlist3                                                 
import3 <- lapply(rlist3,raster)
ndvi.multi <- stack(import3)

# ext <- c(75.25, 75.55, 12.15, 12.45)                 
ndvi.wg <- crop(ndvi.multi,ext)

names(ndvi.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009",
                    "Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")

#________________________________________________________________________________________________
### to see the general statistics of the cropped for the AOI satellite producs and in particular the medians values

# DMP
summary(dmp.wg)

# FAPAR
summary(fapar.wg)

# NDVI
summary(ndvi.wg)


#________________________________________________________________________________________________
# Pearson's Correlation Analysis amoung our satellite cropped for the Central Western Ghat (CWG) products

# Used R required packages
library(GGally)

# Pearson's Correlation Analyses amoung DMP dataset CWG products
ggpairs(dmp.wg)

# Pearson's Correlation Analyses amoung FAPAR dataset CWG products
ggpairs(fapar.wg)

# Pearson's Correlation Analyses amoung NDVI dataset CWG products
ggpairs(ndvi.wg)


#________________________________________________________________________________________________
# Principal Component Analysis (PCA) for our raster objects datasets (dataset CWG products)


# Raster PCA analssys all the DMP images products to modelling raster objects in a single one evaluating a single component (DMP index in this case)
dmpPCA <- rasterPCA(dmp.wg)
summary(dmpPCA$model)                                                          
plot(dmpPCA$map)

# Raster PCA analysis all the FAPAR images products to modelling raster objects in a single one evaluating a single component (FAPAR index in this case)
faparPCA <- rasterPCA(fapar.wg)
summary(faparPCA$model)                                                         
plot(faparPCA$map)

# Raster PCA analysis all the NDVI images products to modelling raster objects in a single one evaluating a single component (NDVI index in this case)
ndviPCA <- rasterPCA(ndvi.wg)
summary(ndviPCA$model)                                                          
plot(ndviPCA$map)

#________________________________________________________________________________________________________________________________________
## Builds the three data frame used for the statistical validation
# Used required R packages
library(tidyverse)

#___________________DMP data frame_____________________

fn <- system.file("external/test.grd", package="raster")

stc <- stack(
  fn,
  fn
)

stc_df <- fortify(dmp.wg, maxpixels = 5689958400) %>% 
  pivot_longer(
    .,
    cols = -(1:2),
    names_to = "layer",
  )


#___________________Data frame structure_________________
stc_df$value
str(stc_df)


#___________________FAPAR data frame_____________________

fn <- system.file("external/test.grd", package="raster")

stc <- stack(
  fn,
  fn
)

stc_df2 <- fortify(fapar.wg, maxpixels = 5689958400) %>% 
  pivot_longer(
    .,
    cols = -(1:2),
    names_to = "layer",
  )

#___________________Data frame structure_________________
stc_df2

str(stc_df2)

stc_df$value


#___________________NDVI data frame_____________________

fn <- system.file("external/test.grd", package="raster")

stc <- stack(
  fn,
  fn
)

stc_df3 <- fortify(ndvi.wg, maxpixels = 5689958400) %>% 
  pivot_longer(
    .,
    cols = -(1:2),
    names_to = "layer",
  )

#___________________Data frame structure_________________
stc_df3$value

str(stc_df3)



#_____________Linear regression models___________________

#_______________________DMP______________________________

plot(stc_df)
mod1 <- lm(stc_df)
mod1
summary(mod1)
plot(mod1)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)

#_______________________FAPAR___________________________

plot(stc_df2)
mod2 <- lm(stc_df2)
mod2
summary(mod2)
plot(mod2)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)

#________________________NDVI___________________________

plot(stc_df3)
mod3 <- lm(stc_df3)
mod3
summary(mod3)
plot(mod3)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)

# Computing the medians of each products of the datasets
# DMP
summary(dmp.wg)

# FAPAR
summary(fapar.wg)

# NDVI
summary(ndvi.wg)


#________________________________________________________________________________________________
## Linear regression models

#______________________DMP_________________________
plot(stc_df)
mod1 <- lm(stc_df)
mod1
summary(mod1)

#Residual standard error: 0.05669 on 25213 degrees of freedom
#(195 observations deleted due to missingness)
#Multiple R-squared:  0.5818,	Adjusted R-squared:  0.5814 
#F-statistic:  1525 on 23 and 25213 DF,  p-value: < 2.2e-16

plot(mod1)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)

#______________________FAPAR__________________________

plot(stc_df2)
mod2 <- lm(stc_df2)
mod2
summary(mod2)

# Residual standard error: 0.06684 on 25408 degrees of freedom
# Multiple R-squared:  0.4183,	Adjusted R-squared:  0.4178 
# F-statistic: 794.5 on 23 and 25408 DF,  p-value: < 2.2e-16

plot(mod2)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)

#________________________NDVI__________________________
plot(stc_df3)
mod3 <- lm(stc_df3)
mod3
summary(mod3)

# Residual standard error: 0.06957 on 25408 degrees of freedom
# Multiple R-squared:  0.3699,	Adjusted R-squared:  0.3693 
# F-statistic: 648.5 on 23 and 25408 DF,  p-value: < 2.2e-16

plot(mod3)
#abline(a = mod1$coefficients[1], b = mod1$coefficients[2], col = "red", lwd = 3)


#________________________________________________________________________________________________
# Create a matrix for the dmp median and for the time
summary(dmp.wg)
DMP_Medians <- c(93.19,79.89,79.18,81.87,77.41,79.43,83.84,83.86,94.74,88.61,92.42,94.85,81.11,87.39,75.16,80.39,
                 100.36,88.80,96.69,84.68,109.67,85.39)

summary(fapar.wg)
FAPAR_Medians <- c(0.736, 0.708, 0.68, 0.724, 0.692, 0.68, 0.728, 0.728, 0.728, 0.74, 0.756, 0.784, 0.74, 0.752, 
                   0.736, 0.692, 0.804, 0.788, 0.736, 0.76, 0.82, 0.83)

summary(ndvi.wg)
NDVI_Medians <- c(0.752, 0.788, 0.76, 0.772, 0.704, 0.736, 0.752, 0.776, 0.74, 0.78, 0.776, 0.82, 0.756, 0.776, 0.788, 
                  0.756, 0.816, 0.808, 0.74, 0.728, 0.816, 0.836)


Time <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020)


# Annual Precipitazioni in Uppangala

Precipitation <- c(5987, 4675, 5492, 5522, 5473, 5371, 5043, 5415, 3907, 4271, 4945, 5510, 5167, 5766, 4674, 5556, 5087, 5115, 3993, 6479, 
       5042, 4320, 4463, 5026, 7220, 5375)
time<- c(1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 
         2014, 2015, 2016, 2017, 2018, 2019)

#________________________________________________________________________________________________

plot(Time, DMP_Medians)

plot(Time, DMP_Medians, type = "b", pch = 19, 
     col = "red", xlab = "x", ylab = "y")
plot(Time,DMP_Medians, type="l", col="green", lwd=5, xlab="time", ylab="Kg / ha / day", main="January Western Ghats Dry Matter Productivity")

### line fitting the DMP_Medians values

line(DMP_Medians)

#   Call:
#     line(DMP_Medians)

#   Coefficients:
#      intercept      Angolar coeff.
#     [1]  80.5205   0.4439

### line fitting the FAPAR_Medians

line(FAPAR_Medians)

#    Call:
#      line(FAPAR_Medians)

#    Coefficients:
#      intercept      Angolar coeff.
#      [1]  0.696714  0.004143

### line fitting the FAPAR_Medians

line(NDVI_Medians)

#    Call:
#      line(NDVI_Medians)

#    Coefficients:
#      intercept   Angolar coeff.
#      [1]  0.746  0.003


#___________________________________________________________________________________
# Data Frames creation for Time and Variables

df1 <- data.frame(Time,DMP_Medians)
df2 <- data.frame(Time,FAPAR_Medians)
df3 <- data.frame(Time,NDVI_Medians)
df4 <- data.frame(Time,Precipitation)

# Data Frames creation for Variables and Variables

df5 <- data.frame(Precipitation,DMP_Medians)
df6 <- data.frame(Precipitation,FAPAR_Medians)
df7 <- data.frame(Precipitation,NDVI_Medians)
df8 <- data.frame(FAPAR_Medians,DMP_Medians)
df9 <- data.frame(NDVI_Medians,DMP_Medians)
df10 <- data.frame(NDVI_Medians,FAPAR_Medians)

#_____________________________________________________________________________________________________________
## Linear regression models computation and visualization between Time (dipendent) and Variables (indipendent)

# Time ~ DMP
lm_dmp <- lm(df1)
summary(lm_dmp)

# Time ~ FAPAR
lm_fapar <- lm(df2)
summary(lm_fapar)

# time ~ NDVI
lm_ndvi <- lm(df3)
summary(lm_ndvi)

# Time ~ Annual Precipitation
lm_prec <- lm(df4)
summary(lm_prec)

#_____________________________________________________________________________________________________________
## Linear regression models computation and visualization between Variable (dipendent) and Variable (indipendent)

# Annual Precipitation ~ DMP
lm_dmp_prec <- lm(df5)
summary(lm_dmp_prec)

# Annual Precipitation ~ FAPAR
lm_fapar_prec <- lm(df6)
summary(lm_fapar_prec)

# Annual Precipitation ~ NDVI
lm_ndvi_prec <- lm(df7)
summary(lm_ndvi_prec)

# FAPAR ~ DMP 
lm_fapar_dmp <- lm(df8)
summary(lm_fapar_dmp)

# NDVI ~ DMP
lm_ndvi_dmp <- lm(df9)
summary(lm_ndvi_dmp)

# NDVI ~ FAPAR
lm_ndvi_fapar <- lm(df10)
summary(lm_ndvi_fapar)

#_______________________________________________________________________________________________________________
# Once we've identified this model as the best, we can proceed to fit the model and analyze the results including the R-squared value and the beta coefficients to determine the exact relationship between the set of predictor variables and the response variable.

# Used required R packages 
library(ggplot2)
library(gridExtra)

# Method Linear Regression Models between Time (dependent) and Variables (independent)

a1 <- ggplot(df1, aes(Time,DMP), xlab = ("Time (year)"), ylab = ("DMP (kg/ha/day)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a2 <- ggplot(df2, aes(Time,FAPAR), xlab = ("Time (year)"), ylab = ("FAPAR (n.index)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a3 <- ggplot(df3, aes(Time,NDVI), xlab = ("Time (year)"), ylab = ("NDVI (m.index)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a4 <- ggplot(df4, aes(Time,Precipitation), xlab = ("Time (year)"), ylab = ("Annual Precipitation (mm)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)

# Multi-visualization Scatter plot with LM of Time ~ Variables 

grid.arrange(a1,a2,a3,a4, ncol=2)

#________________________________________________________________________________________________
## Linear regression models computation and visualization between Time (dipendent) and Variables (indipendent)

# Method Linear Regression Models between Variables (dependent) and Variables (independent)

a5 <- ggplot(df5, aes(Precipitation,DMP), xlab = ("Annual Precipitation (mm)"), ylab = ("DMP (kg/ha/day)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a6 <- ggplot(df6, aes(Precipitation,FAPAR), xlab = ("Annual Precipitation (mm)"), ylab = ("FAPAR (n.index)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a7 <- ggplot(df7, aes(Precipitation,NDVI), xlab = ("Annual Precipitation (mm)"), ylab = ("NDVI (m.index)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a8 <- ggplot(df8, aes(FAPAR,DMP), xlab = ("FAPAR (n.index)"), ylab = ("DMP (kg/ha/day)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a9 <- ggplot(df9, aes(NDVI,DMP), xlab = ("NDVI (m.index)"), ylab = ("DMP (kg/ha/day)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a10 <- ggplot(df10, aes(NDVI,FAPAR), xlab = ("NDVI (m.index)"), ylab = ("FAPAR (n.index)")) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)

# Multi-visualization Scatter plot and LM between Time and Variables

grid.arrange(a5,a6,a7,a8,a8,a10, ncol=2)

#_______________________________________________________________________________________________
# Method Loess (locally-weighted polynomial regression)

a1 <- ggplot(df1, aes(DMP, Time)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a2 <- ggplot(df2, aes(FAPAR, Time)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a3 <- ggplot(df3, aes(NDVI, Time)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a4 <- ggplot(df4, aes(Precipitation, Time)) + geom_point() + stat_smooth(method = "loess", size = 1.5)


# Multi-visualization Scatter plot and LOESS model between Time and Variables

grid.arrange(a1,a2,a3,a4, ncol=2)

#________________________________________________________________________________________________
# Used R required packages

# Applied the T test - Student's Test at our data frames

# Time ~ Variables
t.test(Time,DMP_Medians,paired=TRUE)
t.test(Time,FAPAR_Medians,paired=TRUE)
t.test(Time,NDVI_Medians,paired=TRUE)
t.test(Time,Precipitation,paired=TRUE)

# Annual Precipitation ~ Vegetation Indices (VI)
t.test(Precipitation,DMP_Medians,paired=TRUE)
t.test(Precipitation,FAPAR_Medians,paired=TRUE)
t.test(Precipitation,NDVI_Medians,paired=TRUE)

# VI ~ VI
t.test(FAPAR_Medians,DMP_Medians,paired=TRUE)
t.test(NDVI_medians,DMP_Medians,paired=TRUE)
t.test(NDVI_medians,FAPAR_Medians,paired=TRUE)

#________________________________________________________________________________________________
# Akaike information criterion (AIC) test

## AIC test part 
library(AICcmodavg)

models <- list(mod1, mod2, mod3)
mod.names <- c("DMP","FAPAR","NDVI")

aictab(cand.set = models,modnames = mod.names)

#         Model selection based on AICc:
  
#       K   AICc       Delta_AICc  AICcWt Cum.Wt  LL
# DMP   25 -73219.80       0.00      1      1     36634.93
# FAPAR 25 -65412.91    7806.89      0      1     32731.48
# NDVI  25 -63378.27    9841.53      0      1     31714.16

# Here's how to interpret the output:
  
# K: The number of parameters in the model.
# AICc: The AIC value of the model. The lowercase 'c' indicates that the AIC has been calculated from the AIC corrected for small sample sizes.
# Delta_AICc: The difference between the AIC of the best model compared to the current model being compared.
# AICcWt: The proportion of the total predictive power that can be found in the model.
# Cum.Wt: The cumulative sum of the AIC weights.
# LL: The log-likelihood of the model. This tells us how likely the model is, given the data we used.
# The model with the lowest AIC value is always listed first. From the output we can see that the following model has the lowest AIC value and is thus the best fitting model:

# e.g. mpg = β0 + β1(disp) + β2(hp) + β3(wt) + β4(qsec)





