# initial work conducted by S. Decker
# imported to rep by K. Davidson June 2020
# minor changes made to run code on this wd 

# the issue with this code is that it isn't transparent where the 'esc_data_expanded.txt' file comes from, and if any processing occurred between
# exporting it from the original parent excel file before it was imported to R. As well, a lot of columns appear to be products of calculations,
# which were done in excel (or elsewhere?) which is also not transparent. 
# this analysis will attempt to be re-created in another script () to allow for reproducability and transparency in analysis. 

setwd("~/Documents/ANALYSIS/Data/calibration")


rm(list=ls(all=TRUE))

#Read in the dataframe; select appropriate .txt file for year of interest 
data <- read.table("esc_data_expanded.txt", header=T)
names(data)
str(data)
  
  
#cube-root transformation does the best job of normalizing the data
data$CubeFenceEsc <- (data$FenceEsc)^(1/3)
data$CubeMaxLiveDead2 <- (data$MaxLiveDead2)^(1/3)
data$CubeTotalCarc <- (data$TotalCarc)^(1/3)
data$CubeDeadvsLive <- (data$DeadvsLive)^(1/3)
data$CubeMaxLive <- (data$MaxLive)^(1/3)
data$CubeMaxLiveDead1 <- (data$MaxLiveDead1)^(1/3)
data$ExpandMaxLiveDead1 <- data$MaxLiveDead1*1.66
data$ResidExpandMaxLiveDead1 <- data$ExpandMaxLiveDead1-data$FenceEsc
data$StandResidExpandMaxLiveDead1 <- data$ResidExpandMaxLiveDead1/data$FenceEsc

# There are trends over time (not linear) in the ratio of total carcs to live + cum. dead, whcih means there will be trends in the degree of bias 
# in escaspement estimates
  
dataForfar <- data[data$Stream == "Forfar",]
dataKynock <- data[data$Stream == "Kynock",]
dataGluske <- data[data$Stream == "Gluske",]
  
par(mfrow=c(2,2))
scatter.smooth(x=dataKynock$Year, y=dataKynock$DeadvsLive, main="Kynock", xlab="Year", ylab="TotalCarc/PeakLive+CumDead")
scatter.smooth(x=dataForfar$Year, y=dataForfar$DeadvsLive, main="Forfar", xlab="Year", ylab="TotalCarc/PeakLive+CumDead")
scatter.smooth(x=dataGluske$Year, y=dataGluske$DeadvsLive, main="Gluske", xlab="Year", ylab="TotalCarc/PeakLive+CumDead")


data=na.omit(data)
#look at histograms of input variables to assess normal distribution assumption (rawa parameters are skewed but FenceEsc/MaxLive1 is normal)
# FenceEsc / MAxLive2 is also somewhat skewed
par(mfrow=c(2,2))
hist(data$FenceEsc, main= "Fence Count")
hist(data$TotalCarc, main= "Total carcasses recovered")
hist(data$MaxLiveDead1, main= "max(live + cumulative dead)")
hist(data$FenceEsc/data$MaxLiveDead1, main= "Fence Count/(peak live + cum. dead)")


#Best candidate models of the bunch; ModelDeadTotal has slightly higher r2 but the 2 dependant var. are highly correlated; they
#are not correlated in ModelDeadRatio

ModelDeadRatio <- lm(CubeFenceEsc ~ CubeDeadvsLive + CubeMaxLiveDead2, data=data)
ModelDeadTotal <- lm(CubeFenceEsc ~ CubeTotalCarc + CubeMaxLiveDead2, data=data)
summary(ModelDeadRatio)
summary(ModelDeadTotal)

par(mfrow=c(2,2)) # 4 charts per frame
plot(ModelDeadRatio) # generate deafualt R diagnostic graphs to evaluate model fit
plot(ModelDeadTotal) # generate deafualt R diagnostic graphs to evaluate model fit


#generate models predictions Of Fence Count (back-transformed from cube root-space) and residuals
PredictModelDeadRatio <- (predict(ModelDeadRatio,list(CubeDeadvsLive=data$CubeDeadvsLive,CubeMaxLiveDead2=data$CubeMaxLiveDead2)))^3
PredictModelDeadTotal <- (predict(ModelDeadTotal,list(CubeTotalCarc=data$CubeTotalCarc,CubeMaxLiveDead2=data$CubeMaxLiveDead2)))^3
ResidModelDeadRatio <- residuals(ModelDeadRatio)
ResidModelDeadTotal <- residuals(ModelDeadTotal)
RealResidModelDeadTotal <- (PredictModelDeadTotal-data$FenceEsc)
RealResidModelDeadRatio <- (PredictModelDeadRatio-data$FenceEsc)
RealStandResidModelDeadTotal <- RealResidModelDeadTotal/data$FenceEsc
RealStandResidModelDeadRatio <- RealResidModelDeadRatio/data$FenceEsc

# this chunk is just to see if adding total carcs actually results in a meaningful improvement in model predictions (doesnt improve r2 sig., but does shrink the range of residuals and slightly decease bias) 
ModelMaxLiveDead2 <- lm(CubeFenceEsc ~ CubeMaxLiveDead2, data=data)
summary(ModelMaxLiveDead2)
par(mfrow=c(2,2)) # 4 charts per frame
plot(ModelMaxLiveDead2)
PredictModelMaxLiveDead2 <- (predict(ModelMaxLiveDead2,list(CubeMaxLiveDead2=data$CubeMaxLiveDead2)))^3
ResidModelMaxLiveDead2 <- residuals(ModelMaxLiveDead2)
RealResidModelMaxLiveDead2 <- (PredictModelMaxLiveDead2-data$FenceEsc)
RealStandResidModelMaxLiveDead2 <- RealResidModelMaxLiveDead2/data$FenceEsc

#model predictions vs fence counts (with 1:1 line) (real space)
par(mfrow=c(2,2))
plot(PredictModelDeadTotal, data$FenceEsc, main="Dead Total Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelDeadTotal)
abline(fit, col = "red")
plot(PredictModelDeadRatio, data$FenceEsc, main="Dead Ratio Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelDeadRatio)
abline(fit, col = "red")
plot(data$ExpandMaxLiveDead1, data$FenceEsc, main="Current 1.65 Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ data$ExpandMaxLiveDead1)
abline(fit, col = "red")
plot(PredictModelMaxLiveDead2, data$FenceEsc, main="Dead Max Live 2 Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelMaxLiveDead2)
abline(fit, col = "red")

# model residuals (real space)
par(mfrow=c(2,2))
plot(data$FenceEsc, RealResidModelDeadTotal, main="Dead Total Model residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000))
plot(data$FenceEsc, data$ResidExpandMaxLiveDead1, main="Current 1.65 method", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000))
plot(data$FenceEsc, RealResidModelDeadRatio, main="Dead Ratio Model", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000))
plot(data$FenceEsc, RealResidModelMaxLiveDead2, main="MaxLiveDead2 Model", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000))

# model standardized residuals with trend line for bias vs abundance (real space)
par(mfrow=c(2,2))
plot(data$FenceEsc, RealStandResidModelDeadTotal, main="Dead Total Model", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadTotal ~ data$FenceEsc)
abline(fit, col = "red")

plot(data$FenceEsc, data$StandResidExpandMaxLiveDead1, main="Current 1.65 method", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ data$FenceEsc)
abline(fit, col = "red")

plot(data$FenceEsc, RealStandResidModelDeadRatio, main="Dead Ratio Model", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadRatio ~ data$FenceEsc)
abline(fit, col = "red")

plot(data$FenceEsc, RealStandResidModelMaxLiveDead2, main="MaxLiveDead2  Model", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelMaxLiveDead2 ~ data$FenceEsc)
abline(fit, col = "red")


# % bias in estimate for each model as a function of (total recovery / pop estimate (derived by current peak live + cum. dead method))
RatioDeadMaxLiveDead1 <-data$TotalCarc/(data$MaxLiveDead1)
par(mfrow=c(2,2))

plot(RatioDeadMaxLiveDead1, RealStandResidModelDeadTotal, main="max_count + total dead regression", xlab="total carcass / peak live + cum. dead", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadTotal ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")

plot(RatioDeadMaxLiveDead1, RealStandResidModelDeadRatio, main="Dead Total Ratio residuals", xlab="total carcass / peak live + cum. dead", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadRatio ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")
  plot(RatioDeadMaxLiveDead1, data$StandResidExpandMaxLiveDead1, main="Current 1.65 method", xlab="total carcass / peak live + cum. dead", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")

plot(RatioDeadMaxLiveDead1, RealStandResidModelMaxLiveDead2, main="Max count only regression", xlab="total carcass / peak live + cum. dead", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelMaxLiveDead2 ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")


# test whether Cube transformation is important to model fit (it is beacause of skew in TotalCarc, MaxLive2 and FenceEsc)
# results in this chunk show that regression without the cube root transformation results in poor fit
ModelDeadTotalUT <- lm(FenceEsc ~ TotalCarc + MaxLiveDead2, data=data)
ModelDeadTotal <- lm(CubeFenceEsc ~ CubeTotalCarc + CubeMaxLiveDead2, data=data)
summary(ModelDeadTotalUT)
summary(ModelDeadTotal)

par(mfrow=c(2,2)) # 4 charts per frame
plot(ModelDeadTotalUT) # generate deafualt R diagnostic graphs to evaluate model fit
plot(ModelDeadTotal) # generate deafualt R diagnostic graphs to evaluate model fit

PredictModelDeadTotal <- (predict(ModelDeadTotal,list(CubeTotalCarc=data$CubeTotalCarc,CubeMaxLiveDead2=data$CubeMaxLiveDead2)))^3
ResidModelDeadTotal <- residuals(ModelDeadTotal)
RealResidModelDeadTotal <- (PredictModelDeadTotal-data$FenceEsc)
RealStandResidModelDeadTotal <- RealResidModelDeadTotal/data$FenceEsc

PredictModelDeadTotalUT <- predict(ModelDeadTotalUT,list(TotalCarc=data$TotalCarc,MaxLiveDead2=data$MaxLiveDead2))
RealResidModelDeadTotalUT <- (PredictModelDeadTotalUT-data$FenceEsc)
RealStandResidModelDeadTotalUT <- RealResidModelDeadTotalUT/data$FenceEsc


par(mfrow=c(1,2))
plot(data$FenceEsc, RealStandResidModelDeadTotal, main="Dead Total Model", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadTotal ~ data$FenceEsc)
abline(fit, col = "red")

plot(data$FenceEsc, RealStandResidModelDeadTotalUT, main="Dead Total Model un-transformed", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000))
fit <- lm(RealStandResidModelDeadTotalUT ~ data$FenceEsc)
abline(fit, col = "red")  

plot(RatioDeadMaxLiveDead1, RealStandResidModelDeadTotalUT, main="Dead Total Model residuals", xlab="ratio total carcasses to peak live + dead * 1.66", ylab="std. residuals", abline(0,0), ylim=c(-.5,.5))
fit <- lm(RealStandResidModelDeadTotalUT ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")


#test effect of cube root transformation on estimation of simple peak+cum. dead expansion factor
# result here is that performing regression thru origin on transformed MaxLive1 data doesn't really improve the fit compared to the
#conventional (untransformed) 1.66 expansion.  This is because FenceEsc/MaxLive1 is not skewed to begin with (unlike TotalCarc, MaxLive1,2, etc)
(mean(data$CubeFenceEsc/data$CubeMaxLiveDead1))^3
# 1.64 expansion (doesn't change much from non-transformed estimate)
ModelSimRegExpansion <- lm(CubeFenceEsc ~ CubeMaxLiveDead1 -1, data=data)
summary(ModelSimRegExpansion)
PredictModelSimRegExpansion <- (predict(ModelSimRegExpansion,list(CubeMaxLiveDead1=data$CubeMaxLiveDead1)))^3
RealResidModelSimRegExpansion <- (PredictModelSimRegExpansion-data$FenceEsc)
RealStandResidModelSimRegExpansion <- (RealResidModelSimRegExpansion/data$FenceEsc)


par(mfrow=c(2,2))
plot(PredictModelSimRegExpansion, data$FenceEsc, main="Simple Regression thru origin Model trans. data", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelSimRegExpansion)
abline(fit, col = "red")

plot(data$FenceEsc, RealStandResidModelSimRegExpansion, main="Simple Regression thru origin Model trans. data", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelSimRegExpansion ~ data$FenceEsc)
abline(fit, col = "red")

#current expannsion model for comparison
plot(data$ExpandMaxLiveDead1, data$FenceEsc, main="Current 1.65 Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ data$ExpandMaxLiveDead1)
abline(fit, col = "red")

plot(data$FenceEsc, data$StandResidExpandMaxLiveDead1, main="Current 1.65 method", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ data$FenceEsc)
abline(fit, col = "red")


plot(RatioDeadMaxLiveDead1, RealStandResidModelSimRegExpansion, main="Simple Regression thru origin Model trans. data", xlab="ratio total carcasses to conventional 1.66 expansion estimate", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelSimRegExpansion ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")

#current expansion model for comparison
plot(RatioDeadMaxLiveDead1, data$StandResidExpandMaxLiveDead1, main="Expanded peak + cum. dead model residuals", xlab="ratio total carcasses to conventional 1.66 expansion estimate ", ylab="std. residuals", abline(0,0), ylim=c(-.6,.6))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ RatioDeadMaxLiveDead1)
abline(fit, col = "red")



rm(list=ls(all=TRUE))


#extra redundant bits

par(mfrow=c(2,2))
plot(PredictModelDeadRatio, data$FenceEsc, main="Dead Ratio Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
plot(PredictModelDeadTotal, data$FenceEsc, main="Dead Total Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
plot(data$FenceEsc, ResidModelDeadRatio, main="Dead Ratio Model residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-3,3))
plot(data$FenceEsc, ResidModelDeadTotal, main="Dead Total Model residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-3,3))
par(mfrow=c(2,2))
plot(data$FenceEsc, RealResidModelMaxLiveDead2, main="Dead Max Live 2 Model residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000))
plot(data$FenceEsc, RealResidModelDeadTotal, main="Dead Total Model residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-6000,6000)) 


plot(data$FenceEsc, RealStandResidModelMaxLiveDead2, main="Dead Max Live 2 Model std. residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.5,.5))
fit <- lm(RealStandResidModelMaxLiveDead2 ~ data$FenceEsc)
abline(fit, col = "red")
plot(data$FenceEsc, RealStandResidModelDeadTotal, main="Dead Total Model std. residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.5,.5)) 
fit <- lm(RealStandResidModelDeadTotal ~ data$FenceEsc)
abline(fit, col = "red")



#look at how the DeadTotal and MaxLiveDead2 models are superior to current method (ie, more precise and unbiased)
plot(data$FenceEsc, data$StandResidExpandMaxLiveDead1, main="Current 1.65 Model std. residuals", xlab="fence estimate", ylab="residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.5,.5))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ data$FenceEsc)
abline(fit, col = "red")

par(mfrow=c(1,1))
plot(data$Year[data$S], data$DeadvsLiveAndDead, pch=21, 
bg=c("red","green3","blue")[unclass(data$Stream)], 
main="Trend in TotalCarc/MaxLive1 over time")


par(mfrow=c(2,2))
plot(PredictModelDeadTotal, data$FenceEsc, main="max_count + total dead regression", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelDeadTotal)
abline(fit, col = "red")
plot(data$ExpandMaxLiveDead1, data$FenceEsc, main="Current 1.65 Model", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ data$ExpandMaxLiveDead1)
abline(fit, col = "red")
plot(PredictModelMaxLiveDead2, data$FenceEsc, main="max count only regression", xlab="prediction", ylab="fence estimate", abline(0,1), xlim=c(0,30000), ylim=c(0,30000))
fit <- lm(data$FenceEsc ~ PredictModelMaxLiveDead2)
abline(fit, col = "red")



par(mfrow=c(2,2))
plot(data$FenceEsc, RealStandResidModelDeadTotal, main="max_count + total dead regression", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelDeadTotal ~ data$FenceEsc)
abline(fit, col = "red")

plot(data$FenceEsc, data$StandResidExpandMaxLiveDead1, main="Current 1.65 method", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(data$StandResidExpandMaxLiveDead1 ~ data$FenceEsc)
abline(fit, col = "red")


plot(data$FenceEsc, RealStandResidModelMaxLiveDead2, main="max count only regression", xlab="fence estimate", ylab="std. residuals", abline(0,0), xlim=c(0,30000), ylim=c(-.6,.6))
fit <- lm(RealStandResidModelMaxLiveDead2 ~ data$FenceEsc)
abline(fit, col = "red")


