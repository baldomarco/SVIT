## Multiple Time Series Analysis of Dry Matter Productivity in Indian Western Gaths region
## Periods are always Jan 10th for every available year of the data products (Western Ghats dry season)
## Multiple Time Series Analysis on DMP, FAPAR, NDVI V2 1Km spatial scale from VITO catalogues 1999-2020 of South Western Ghats Region

library(ncdf4)
library(raster)
library(RStoolbox)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(rasterVis)
library(GGally)
library(scico)

# DMP

setwd("D:/uppangala/cop_dmp_1km/")                   # if I want to use the dataset at 1km spatial resolution # 22years time resolution

rlist <- list.files (pattern ="DMP")
rlist                                                # 22 years
import <- lapply(rlist,raster)
dmp.multi <- stack(import) 

ext <- c(75.25, 75.55, 12.15, 12.45)                 # South-Center Western Gaths Region
dmp.wg <- crop(dmp.multi,ext)

names(dmp.wg) <- c("Jan 99"," Jan 00"," Jan 01"," Jan 02"," Jan 03","Jan 04","Jan 05","Jan 06","Jan 07","Jan 08","Jan 09","Jan 10","Jan 11","Jan 12","Jan 13","Jan 14","Jan 15","Jan 16","Jan 17","Jan 18","Jan 19","Jan 20")
ggpairs(dmp.wg)

levelplot(dmp.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(dmp.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")

#names(dmp.wg) <- c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020")
#ggpairs(dmp.wg, names= c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))

# FAPAR

setwd("D:/uppangala/cop_fapar_1km/")                   # if I want to use the dataset at 1km spatial resolution # 22years time resolution

rlist2 <- list.files (pattern ="FAPAR")
rlist2                                                 # 22 years
import2 <- lapply(rlist2,raster)
fapar.multi <- stack(import2)

# ext <- c(75.25, 75.55, 12.15, 12.45)                 # South-Central Western Gaths Region
fapar.wg <- crop(fapar.multi,ext)

names(fapar.wg) <- c("Jan 99"," Jan 00"," Jan 01"," Jan 02"," Jan 03","Jan 04","Jan 05","Jan 06","Jan 07","Jan 08","Jan 09","Jan 10","Jan 11","Jan 12","Jan 13","Jan 14","Jan 15","Jan 16","Jan 17","Jan 18","Jan 19","Jan 20")
ggpairs(fapar.wg)

levelplot(fapar.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(fapar.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")


# NDVI
setwd("D:/uppangala/cop_ndvi_1km/")                   # if I want to use the dataset at 1km spatial resolution # 22years time resolution
x2 <- raster("vc_gls_NDVI_202001110000_GLOBE_PROBAV_V2.2.1.nc")
plot(x2)
x2 <- reclassify(x2, cbind(253:255, NA)) 


rlist3 <- list.files (pattern ="NDVI")
rlist3                                                 # 22 years
import3 <- lapply(rlist3,raster)
ndvi.multi <- stack(import3)

# ext <- c(75.25, 75.55, 12.15, 12.45)                 # South-Center Western Gaths Region
ndvi.wg <- crop(ndvi.multi,ext)

names(ndvi.wg) <- c("Jan 99"," Jan 00"," Jan 01"," Jan 02"," Jan 03","Jan 04","Jan 05","Jan 06","Jan 07","Jan 08","Jan 09","Jan 10","Jan 11","Jan 12","Jan 13","Jan 14","Jan 15","Jan 16","Jan 17","Jan 18","Jan 19","Jan 20")
ggpairs(ndvi.wg)

levelplot(ndvi.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(ndvi.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")


## costruisco i 3 DF

library(tidyverse)
library(raster)
library(RStoolbox)


# DF dry matter productivity SW ghats 

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

stc_df$value

str(stc_df)

# DF fapar sw ghats

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

stc_df2

str(stc_df2)

stc_df$value

# DF ndvi sw ghats

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

stc_df3$value

str(stc_df3)

plot(stc_df3)

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

### to see the medians

summary(dmp.wg)

summary(fapar.wg)

summary(ndvi.wg)

# trying the principal component analysis
dmp_pca <- rasterPCA(dmp.wg)
plot(dmp_pca$map) 
summary(dmp_pca$model)
#dmp_pca$model

# Create a matrix for the dmp median and for the time
summary(dmp.wg)
DMP_Medians <- c(93.19,79.89,79.18,81.87,77.41,79.43,83.84,83.86,94.74,88.61,92.42,94.85,81.11,87.39,75.16,80.39,100.36,88.80,96.69,84.68,109.67,85.39)

summary(fapar.wg)
FAPAR_Medians <- c(0.736, 0.708, 0.68, 0.724, 0.692, 0.68, 0.728, 0.728, 0.728, 0.74, 0.756, 0.784, 0.74, 0.752, 0.736, 0.692, 0.804, 0.788, 0.736, 0.76, 0.82, 0.83)

summary(ndvi.wg)
NDVI_Medians <- c(0.752, 0.788, 0.76, 0.772, 0.704, 0.736, 0.752, 0.776, 0.74, 0.78, 0.776, 0.82, 0.756, 0.776, 0.788, 0.756, 0.816, 0.808, 0.74, 0.728, 0.816, 0.836)


Time <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020)

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




###########################__________________________________________###########################################

Tilde operator is used to define the relationship between dependent variable and independent variables in a statistical model formula. The variable on the left-hand side of tilde operator is the dependent variable and the variable(s) on the right-hand side of tilde operator is/are called the independent variable(s). So, tilde operator helps to define that dependent variable depends on the independent variable(s) that are on the right-hand side of tilde operator.
Example

> Regression_Model <- lm(y~ x1 + x2 + x3)

Here, the object Regression_Model stores the formula for linear regression model created by using function lm and y is the dependent variable and x1, x2, and x3 are independent variables.

This model can be created by using a dot (.) if we want to include all the independent variables but for this purpose, we should have all the variables stored in a data frame.
Example

> Regression_Data <- data.frame(x1, x2, x3, y)
> Regression_Model_New < - lm(y~ . , data = Regression_Data)

This will have the same output as the previous model, but we cannot use tilde with dot if we want to create a model with few variables.

Suppose you want to create a new model with x1 and x3 only then it can be done as follows −

> Regression_Model_New1 <- lm(y~ x1 + x3, data = Regression_Data)

But we cannot do it using dot with tilde as −

> Regression_Model_New2_Incorrect <- lm(y~ . + x3, data = Regression_Data)

https://www.tutorialspoint.com/how-to-write-text-and-output-it-as-a-text-file-using-r


# COULD BE PRESENT SOME ERROR STARTING FROM THE DIPENDENT (VEG INDICES) AND INDIPENDENT VARIABLE (TIME) OF THE SISTEM!!! 
# IN MORE WE SHOULD EVALUATE THE STATISTICS COMPERING THE THREE INDECES.


#     https://stackoverflow.com/questions/27539033/r-apply-lm-on-each-data-frame-row



## Linear model time-dry matter productivity medians (same result)
m <- lm(Time ~ DMP_Medians)
m
plot(Time, DMP_Medians)
summary(m)

## LINEAR MODELS BETWEEN INDICES
m1 <- lm(DMP_Medians ~ FAPAR_Medians)
m1
plot(DMP_Medians, FAPAR_Medians)
summary(m1)


#Call:
#  lm(formula = DMP_Medians ~ FAPAR_Medians)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-13.5209  -2.9593  -0.5434   3.0980  12.0996 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -12.35      25.18   -0.49 0.629164    
#FAPAR_Medians   134.05      33.85    3.96 0.000773 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 6.514 on 20 degrees of freedom
#Multiple R-squared:  0.4395,	Adjusted R-squared:  0.4114 
#F-statistic: 15.68 on 1 and 20 DF,  p-value: 0.0007727



m2 <- lm(FAPAR_Medians ~ NDVI_Medians)
m2
plot(FAPAR_Medians, NDVI_Medians)
summary(m2)

#Call:
#  lm(formula = FAPAR_Medians ~ NDVI_Medians)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.051794 -0.019110  0.007048  0.014681  0.058523 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.01177    0.14187   0.083    0.935    
#NDVI_Medians  0.94741    0.18369   5.158 4.79e-05 ***
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02819 on 20 degrees of freedom
#Multiple R-squared:  0.5708,	Adjusted R-squared:  0.5494 
#F-statistic:  26.6 on 1 and 20 DF,  p-value: 4.789e-05

m3 <- lm(DMP_Medians ~ NDVI_Medians)
m3
plot(DMP_Medians, NDVI_Medians)
summary(m3)

#Call:
#  lm(formula = t3 ~ ndvi_ed)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-10.2367  -3.5376   0.6828   3.2230   8.8819 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1884.19      39.52  47.677 5.55e-16 ***
#  ndvi_ed       159.96      51.06   3.133  0.00793 ** 
#  ---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 5.141 on 13 degrees of freedom
#Multiple R-squared:  0.4302,	Adjusted R-squared:  0.3864 
#F-statistic: 9.816 on 1 and 13 DF,  p-value: 0.007927

## Edit the matrix for thanks the loess function, cutting the out-layers 

dmp_ed <- c(93.19,79.89,79.18,81.87,77.41,79.43,83.84,83.86,94.74,88.61,92.42,94.85,81.11,87.39,80.39,88.80,96.69,85.39)
t1 <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2014,2016,2017,2020)
mod4 <- lm(t1 ~ dmp_ed)
summary(mod4)

#Call:
  #lm(formula = t1 ~ dmp_ed)

#Residuals:
 # Min      1Q  Median      3Q     Max 
#-11.379  -3.487  -1.194   4.377  12.223 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 1979.2936    19.9613  99.157   <2e-16 ***
#  dmp_ed         0.3336     0.2314   1.442    0.169    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 5.95 on 16 degrees of freedom
#Multiple R-squared:  0.115,	Adjusted R-squared:  0.05965 
#F-statistic: 2.078 on 1 and 16 DF,  p-value: 0.1687

fapar_ed <- c(0.736, 0.708, 0.724, 0.692, 0.728, 0.728, 0.728, 0.74, 0.756, 0.74, 0.752, 0.736, 0.82, 0.83)
t2 <- c(1999,2000,2002,2003,2005,2006,2007,2008,2009,2011,2012,2013,2019,2020)
mod5 <- lm(t2 ~ fapar_ed)
summary(mod5)

#lm(formula = t2 ~ fapar_ed)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-7.9695 -0.8428  0.0548  2.0739  6.0305 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1900.91      19.08  99.615  < 2e-16 ***
#  fapar_ed      144.10      25.61   5.626 0.000111 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 3.503 on 12 degrees of freedom
#Multiple R-squared:  0.7251,	Adjusted R-squared:  0.7022 
#F-statistic: 31.65 on 1 and 12 DF,  p-value: 0.0001115

ndvi_ed <- c(0.752, 0.788, 0.76, 0.772, 0.736, 0.752, 0.776, 0.74, 0.78, 0.776, 0.776, 0.788, 0.756, 0.816, 0.836)
t3 <- c(1999,2000,2001,2002,2004,2005,2006,2007,2008,2009,2012,2013,2014,2019,2020)
mod6 <- lm(t3 ~ ndvi_ed)
summary(mod6)

#Call:
#  lm(formula = t3 ~ ndvi_ed)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-10.2367  -3.5376   0.6828   3.2230   8.8819 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1884.19      39.52  47.677 5.55e-16 ***
#  ndvi_ed       159.96      51.06   3.133  0.00793 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 5.141 on 13 degrees of freedom
#Multiple R-squared:  0.4302,	Adjusted R-squared:  0.3864 
#F-statistic: 9.816 on 1 and 13 DF,  p-value: 0.007927

#### QQplot  linear regression model applied to quartiles

qqnorm(DMP_Medians, pch = 1, frame = FALSE)


## Create a DF and use ggplot visualization

# DMP

DMP_df <- data.frame(Time, DMP_Medians)
DMP_df

ggplot(DMP_df, aes(Time, DMP_Medians)) + geom_point() + stat_smooth()

# FAPAR

FAPAR_df <- data.frame(Time, FAPAR_Medians)
FAPAR_df

ggplot(DMP_df, aes(Time, FAPAR_Medians)) + geom_point() + stat_smooth()

# NDVI

NDVI_df <- data.frame(Time, NDVI_Medians)
NDVI_Medians

ggplot(DMP_df, aes(Time, NDVI_Medians)) + geom_point() + stat_smooth()


## AIC test part 
install.packages("Rtools")
library(Rtools)
install.packages("AICcmodavg")
library(AICcmodavg)

models <- list(mod1, mod2, mod3)
mod.names <- c("DMP","FAPAR","NDVI")

aictab(cand.set = models,modnames = mod.names)

#Model selection based on AICc:
  
#       K   AICc       Delta_AICc  AICcWt Cum.Wt  LL
# DMP   25 -73219.80       0.00      1      1     36634.93
# FAPAR 25 -65412.91    7806.89      0      1     32731.48
# NDVI  25 -63378.27    9841.53      0      1     31714.16

#Here’s how to interpret the output:
  
# K: The number of parameters in the model.
# AICc: The AIC value of the model. The lowercase ‘c’ indicates that the AIC has been calculated from the AIC corrected for small sample sizes.
# Delta_AICc: The difference between the AIC of the best model compared to the current model being compared.
# AICcWt: The proportion of the total predictive power that can be found in the model.
# Cum.Wt: The cumulative sum of the AIC weights.
# LL: The log-likelihood of the model. This tells us how likely the model is, given the data we used.
# The model with the lowest AIC value is always listed first. From the output we can see that the following model has the lowest AIC value and is thus the best fitting model:
  
#  mpg = β0 + β1(disp) + β2(hp) + β3(wt) + β4(qsec)

# Once we’ve identified this model as the best, we can proceed to fit the model and analyze the results including the R-squared value and the beta coefficients to determine the exact relationship between the set of predictor variables and the response variable.



