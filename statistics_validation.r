# Statistical validation analyses on the data sets and results of "Remote sensing analysis on primary productivity and forest cover dynamics: a Western Ghats case study"

## Multiple Time Series Analysis of Dry Matter Productivity in the Indian Western Ghats region
## Periods are always Jan 10th for every available year of the data products (Western Ghats dry season)
## Multiple Time Series Analysis on DMP, FAPAR, NDVI V2 1Km spatial scale from VITO catalog 1999-2020 of central Western Ghats region

# Used R packages

library(ncdf4)
library(raster)
library(RStoolbox)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(rasterVis)
library(GGally)
library(scico)

#________________________________________________________________________________________________
# DMP

setwd("D:/toy/uppangala/cop_dmp_1km/")               

rlist <- list.files (pattern ="DMP")
rlist                                                
import <- lapply(rlist,raster)
dmp.multi <- stack(import) 

ext <- c(75.25, 75.55, 12.15, 12.45)                 # central Western Ghats Region
dmp.wg <- crop(dmp.multi,ext)

names(dmp.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009","Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")   
ggpairs(dmp.wg)

levelplot(dmp.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(dmp.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")

#names(dmp.wg) <- c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020")
#ggpairs(dmp.wg, names= c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))

#________________________________________________________________________________________________
# FAPAR

setwd("D:/toy/uppangala/cop_fapar_1km/")                   

rlist2 <- list.files (pattern ="FAPAR")
rlist2                                                     
import2 <- lapply(rlist2,raster)
fapar.multi <- stack(import2)

# Crop the products for my AOI
fapar.wg <- crop(fapar.multi,ext)

names(fapar.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009","Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")
ggpairs(fapar.wg)

levelplot(fapar.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(fapar.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")

#________________________________________________________________________________________________
# NDVI
setwd("D:/toy/uppangala/cop_ndvi_1km/")                   
x2 <- raster("vc_gls_NDVI_202001110000_GLOBE_PROBAV_V2.2.1.nc")
plot(x2)
x2 <- reclassify(x2, cbind(253:255, NA)) 


rlist3 <- list.files (pattern ="NDVI")
rlist3                                                 
import3 <- lapply(rlist3,raster)
ndvi.multi <- stack(import3)

# ext <- c(75.25, 75.55, 12.15, 12.45)                 
ndvi.wg <- crop(ndvi.multi,ext)

names(ndvi.wg) <- c("Jan 1999"," Jan 2000"," Jan 2001"," Jan 2002"," Jan 2003","Jan 2004","Jan 2005","Jan 2006","Jan 2007","Jan 2008","Jan 2009","Jan 2010","Jan 2011","Jan 2012","Jan 2013","Jan 2014","Jan 2015","Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020")
ggpairs(ndvi.wg)

levelplot(ndvi.wg, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"))
boxplot(ndvi.wg, outline=F, horizontal=F, axes=T, names=c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"), col="gold")


#________________________________________________________________________________________________
# Principal Component Analysis (PCA) for raster objects datasets
# PCA DMP
dmp_pca <- rasterPCA(dmp.wg)
plot(dmp_pca$map) 
summary(dmp_pca$model)

# PCA FAPAR
fapar_pca <- rasterPCA(fapar.wg)
plot(fapar_pca$map) 
summary(fapar_pca$model)

# PCA DMP
ndvi_pca <- rasterPCA(ndvi.wg)
plot(ndvi_pca$map) 
summary(ndvi_pca$model)


#________________________________________________________________________________________________
### to see the medians

summary(dmp.wg)

summary(fapar.wg)

summary(ndvi.wg)


#________________________________________________________________________________________________
## Data frame extraction

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
DMP_Medians <- c(93.19,79.89,79.18,81.87,77.41,79.43,83.84,83.86,94.74,88.61,92.42,94.85,81.11,87.39,75.16,80.39,100.36,88.80,96.69,84.68,109.67,85.39)

summary(fapar.wg)
FAPAR_Medians <- c(0.736, 0.708, 0.68, 0.724, 0.692, 0.68, 0.728, 0.728, 0.728, 0.74, 0.756, 0.784, 0.74, 0.752, 0.736, 0.692, 0.804, 0.788, 0.736, 0.76, 0.82, 0.83)

summary(ndvi.wg)
NDVI_Medians <- c(0.752, 0.788, 0.76, 0.772, 0.704, 0.736, 0.752, 0.776, 0.74, 0.78, 0.776, 0.82, 0.756, 0.776, 0.788, 0.756, 0.816, 0.808, 0.74, 0.728, 0.816, 0.836)


Time <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020)


# Precipitazioni

P <- c(5987, 4675, 5492, 5522, 5473, 5371, 5043, 5415, 3907, 4271, 4945, 5510, 5167, 5766, 4674, 5556, 5087, 5115, 3993, 6479, 5042, 4320, 4463, 5026, 7220, 5375)
time<- c(1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)
P
time

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

#________________________________________________________________________________________________
## Linear model time-dry matter productivity medians 
m <- lm(Time ~ DMP_Medians)
m
ggplot(df1, aes(y1, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m)

# FAPAR

m1 <- lm(Time ~ FAPAR_Medians)
m1
ggplot(df1, aes(y1, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m1)

# NDVI

m2 <- lm(Time ~ NDVI_Medians )
m2
ggplot(df1, aes(y1, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m2)

#________________________________________________________________________________________________
## Create a DF and use ggplot visualization

# DMP

DMP_df <- data.frame(Time, DMP_Medians)
DMP_df

lm_dmp <- lm(DMP_df)
summary(lm_dmp)

ggplot(DMP_df, aes(Time, DMP_Medians)) + geom_point() + stat_smooth()

# FAPAR

FAPAR_df <- data.frame(Time, FAPAR_Medians)
FAPAR_df

ggplot(DMP_df, aes(Time, FAPAR_Medians)) + geom_point() + stat_smooth()

# NDVI

NDVI_df <- data.frame(Time, NDVI_Medians)
NDVI_Medians

ggplot(DMP_df, aes(Time, NDVI_Medians)) + geom_point() + stat_smooth()

#________________________________________________________________________________________________
## LINEAR MODELS BETWEEN INDICES

m3 <- lm(DMP_Medians ~ FAPAR_Medians)
m3
df1 <- data.frame(DMP_Medians, FAPAR_Medians)
p1 <- ggplot(df1, aes(DMP_Medians, FAPAR_Medians)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m3)

m4 <- lm(DMP_Medians ~ NDVI_Medians)
m4
df2 <- data.frame(DMP_Medians, NDVI_Medians)
p2 <- ggplot(df2, aes(DMP_Medians, NDVI_Medians)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m4)

m5 <- lm(NDVI_Medians ~ FAPAR_Medians)
m5
df3 <- data.frame(NDVI_Medians, FAPAR_Medians)
p3 <- ggplot(df3, aes(NDVI_Medians, FAPAR_Medians)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
summary(m5)


m6 <- lm(FAPAR_Medians ~ NDVI_Medians)
m6
plot(FAPAR_Medians, NDVI_Medians)
summary(m6)


m7 <- lm (DMP_Medians ~ NDVI_Medians)
m7
plot(DMP_Medians, NDVI_Medians)
summary(m7)


m8 <- lm (time ~ P)
m8
plot(time, P)
summary(m8)


#Call:
#  lm(formula = t3 ~ ndvi_ed)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-10.2367  -3.5376   0.6828   3.2230   8.8819 

#  Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  1884.19      39.52  47.677 5.55e-16 ***
#  ndvi_ed       159.96      51.06   3.133  0.00793 ** 
#  ------
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#  Residual standard error: 5.141 on 13 degrees of freedom
#  Multiple R-squared:  0.4302,	Adjusted R-squared:  0.3864 
#  F-statistic: 9.816 on 1 and 13 DF,  p-value: 0.007927


#________________________________________________________________________________________________
# Precipitazioni

P_df <- data.frame(time, P)
P_df

l<- lm(P_df)
summary(l)

ggplot(P_df, aes(time, P)) + geom_point() + stat_smooth()

#________________________________________________________________________________________________
# Akaike information criterion (AIC) test part

## AIC test part 
library(Rtools)
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

# Once we've identified this model as the best, we can proceed to fit the model and analyze the results including the R-squared value and the beta coefficients to determine the exact relationship between the set of predictor variables and the response variable.

#________________________________________________________________________________________________
library(gridExtra)

df1 <- data.frame(y1, y5)
df2 <- data.frame(y2, y5)
df3 <- data.frame(y3, y5)
df4 <- data.frame(y4, y5)

#________________________________________________________________________________________________
# Method LM
a1 <- ggplot(df1, aes(y1, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a2 <- ggplot(df2, aes(y2, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a3 <- ggplot(df3, aes(y3, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
a4 <- ggplot(df4, aes(y4, y5)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x)
      
grid.arrange(a1,a2,a3,a4, ncol=2)

#________________________________________________________________________________________________
# Method Loess
a1 <- ggplot(df1, aes(y1, y5)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a2 <- ggplot(df2, aes(y2, y5)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a3 <- ggplot(df3, aes(y3, y5)) + geom_point() + stat_smooth(method = "loess", size = 1.5)
a4 <- ggplot(df4, aes(y4, y5)) + geom_point() + stat_smooth(method = "loess", size = 1.5)

grid.arrange(a1,a2,a3,a4, ncol=2)

#________________________________________________________________________________________________
# LM with time as indipendent variable (at the growth of the time, the indices 've grown)

m_dmp <- lm (y3 ~ y1)
summary(m_dmp)

m_fapar <- lm  (y4 ~ y3)   
summary(m_fapar)

m_ndvi <- lm (y5 ~ y3)
summary(m_ndvi)

m_p <- lm (y5 ~ y4)
summary(m_p)

#________________________________________________________________________________________________
# T test - Student test

t.test(y1,y2,paired=TRUE)
t.test(y1,y3,paired=TRUE)
t.test(y2,y3,paired=TRUE)
t.test(y1,y4,paired=TRUE)
t.test(y2,y4,paired=TRUE)
t.test(y3,y4,paired=TRUE)

#________________________________________________________________________________________________
# LM with info in the graphics 
# https://rpkgs.datanovia.com/ggpubr/reference/stat_regline_equation.html

install.packages("ggpmisc")
library(ggpmisc)


# Simple scatter plot with correlation coefficient and
# regression line
#::::::::::::::::::::::::::::::::::::::::::::::::::::


ggscatter(df1, x = "y5", y = "y1", add = "reg.line") +
  geom_smooth(method = "lm", formula = y ~ x)


CORRELATIONP3 <-CORRELATIONP2[product=='',]

x<-CORRELATIONP3$b
y<-CORRELATIONP3$p


df <- data.frame(x = x)
m <- lm(y ~ x, data = df)
p <- ggplot(data = df, aes(x = x, y = y)) +
  scale_x_continuous("b (%)") +
  scale_y_continuous("p (%)")+
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_point()
p

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(        a = format(coef(m)[1], digits = 4),
                              b = format(coef(m)[2], digits = 4),
                              r2 = format(summary(m)$r.squared, digits = 3)))

dftext <- data.frame(x = 3, y = 0.2, eq = as.character(as.expression(eq)))

p + geom_text(aes(label = eq), data = dftext, parse = TRUE)

#________________________________________________________________________________________________
# LM with p and r2 values in the graphics

CORRELATIONP3 <-CORRELATIONP2[product=='a',]

x<-CORRELATIONP3$b
y<-CORRELATIONP3$p


df <- data.frame(x = x)
m <- lm(y ~ x, data = df)
p <- ggplot(data = df, aes(x = x, y = y)) +
  scale_x_continuous("b (%)") +
  scale_y_continuous("p (%)")+
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_point()
p

eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), pvalue = format(summary(m)$coefficients[2,4], digits = 3)))


dftext <- data.frame(x = 30, y = 0.4, eq = as.character(as.expression(eq)))
p + geom_text(aes(label = eq), data = dftext, parse = TRUE)

#________________________________________________________________________________________________
# Put in a same figure the LM scatter plots
grid.arrange(a1,a2,a3,a4,layout_matrix=rbind(c(2,2)))

grid.arrange(gg1,gg2,gg3,layout_matrix=rbind(c(1,1),c(2,3)))

grid.arrange(arrangeGrob(gg1,gg2,ncol=2), gg3)


#________________________________________________________________________________________________
# Compare models
fit1 <- lm(y5 ~  + y1)
fit2 <- lm(y5 ~ y2)
anova(fit1)

# LM with time as indipendent variable (at the growth of the time, the indices have grown)

m_dmp <- lm (y3 ~ y1)
summary(m_dmp)

m_fapar <- lm  (y3 ~ y2)
summary(m_fapar)

m_ndvi <- lm (y5 ~ y3)
summary(m_ndvi)

m_p <- lm (y5 ~ y4)
summary(m_p)

#________________________________________________________________________________________________
# T test

t.test(y1,y2,paired=TRUE)
t.test(y1,y3,paired=TRUE)
t.test(y2,y3,paired=TRUE)
t.test(y1,y4,paired=TRUE)
t.test(y2,y4,paired=TRUE)
t.test(y3,y4,paired=TRUE)

matrix(data, nrow, ncol, byrow = F)
#________________________________________________________________________________________________

# LM and multiple LM 

set.seed(1)

DF <- data.frame(A=rnorm(50, 100, 3),
                 B=rnorm(50, 100, 3))

resultlist   <- apply(DF, 1, function(y) lm(y1 ~ y5))
resultcoeffs <- apply(DF, 1, function(y) lm(y1 ~ y5)$coefficients)

# It is just one observation per row. Note that you get NA estimates as there are not enough degrees of freedom.

mapply(function(x,y) lm(y~x)$coefficients, DF[,1], DF[,2])

# or

apply(DF1, 1, function(x) lm(x[2]~x[1])$coefficients)

# Suppose, you have many observations per row i.e. x and y variables span over many columns

mapply(function(x,y) lm(y~x)$coefficients, as.data.frame(t(DFNew[1:3])),
       as.data.frame(t(DFNew[4:6])))
# or

apply(DFNew, 1, function(x) lm(x[4:6]~x[1:3])$coefficients)

# data

set.seed(25)
DFNew <- as.data.frame(matrix(sample(1:50,10*6, replace=TRUE), ncol=6))

# END___________________________________________________________________________________________





