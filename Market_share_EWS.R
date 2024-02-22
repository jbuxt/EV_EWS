library(stats)

#Import market share time series
#Data is from Jan 2009 to Dec 2022
market_share_data <- read.csv("ICEV_market_shares.csv")
date <- as.Date(unlist(market_share_data['date']))
national_markets <- c("Germany_UK_France", "China", "USA")

#The code below implements the following steps:
#1. Call in each market share time series and cut it at the start of 2020
#2. Remove the long term trend using a Kernel Regression Smoother with a bandwidth of 8
#3. Calculate the Lag-1 Autocorrelation (AR(1)) and variance of this residual time series
#4. Calculate the AR(1) and variance trend with the Mann-Kendall Tau
mkt_share <- matrix(array(NA,dim=3*168),nrow=168,ncol=3)
mkt_share_trend <- matrix(array(NA,dim=3*132),nrow=132,ncol=3)
mkt_share_res <- matrix(array(NA,dim=3*132),nrow=132,ncol=3)
mkt_share_ar1 <- matrix(array(NA,dim=3*66),nrow=66,ncol=3)
mkt_share_vari <- matrix(array(NA,dim=3*66),nrow=66,ncol=3)
mkt_share_ar1_tau <- array(NA, dim=3)
mkt_share_vari_tau <- array(NA, dim=3)

for (i in 1:length(national_markets)){
  market_share = unlist(market_share_data[national_markets[i]])
  mkt_share[,i] <- market_share
  #Cut market share at start of 2020
  ts = market_share[1:132]
  #Remove long term trend using a Kernel Regression Smoother 
  #Bandwidth is set at 8, but this can be altered to change the fit
  ts_res <- ts - ksmooth(1:length(ts), ts, bandwidth=8)$y
  mkt_share_trend[,i] <- ksmooth(1:length(ts), ts, bandwidth=8)$y
  mkt_share_res[,i] <- ts_res
  
  #Calculate the AR(1) and variance
  #Here we use a moving window with length equal to half the time series length
  l <- length(ts_res)
  wl <- l/2
  ar1 <- array(NA,dim=(l-wl))
  vari <- array(NA,dim=(l-wl))
  for (j in 1:(l-wl)){
    ar1[j] <- cor(ts_res[j:(j+wl-1)],ts_res[(j+1):(j+wl)])
    vari[j] <- var(ts_res[j:(j+wl-1)])
  }
  mkt_share_ar1[,i] <- ar1
  mkt_share_vari[,i] <- vari
  
  #We can now calculate the Mann-Kenall Tau Value for each of these
  mkt_share_ar1_tau[i] <- cor.test(1:length(ar1),ar1,test='kendall')$estimate
  mkt_share_vari_tau[i] <- cor.test(1:length(vari),vari,test='kendall')$estimate
}

#We now can plot these results

par(mar=c(4,4,2, 2)) # c(bottom, left, top, right)

#First, we plot the market share plots for each country
par(mfrow=c(4,3))
plot(date,mkt_share[,1],xlab='Date',ylab='ICV Market Share',main='Germany, UK and France',type='l',ylim=c(0,100),bty='l')
lines(date[(1:132)],mkt_share_trend[,1],col='red')
abline(v=date[132],col='blue',lty=3)

plot(date,mkt_share[,2],xlab='Date',ylab='ICV Market Share',main='China',type='l',ylim=c(0,100),bty='l')
lines(date[(1:132)],mkt_share_trend[,2],col='red')
abline(v=date[132],col='blue',lty=3)

plot(date,mkt_share[,3],xlab='Date',ylab='ICV Market Share',main='USA',type='l',ylim=c(0,100),bty='l')
lines(date[(1:132)],mkt_share_trend[,3],col='red')
abline(v=date[132],col='blue',lty=3)

#Here we plot the residual time series as a barplot
barplot(c(mkt_share_res[,1],rep(NA,36)),ylab='Detrended residual')
barplot(c(mkt_share_res[,2],rep(NA,36)),ylab='Detrended residual')
barplot(c(mkt_share_res[,3],rep(NA,36)),ylab='Detrended residual')


#Then we plot the AR(1) time series with the tau value
plot(date,c(rep(NA,33),mkt_share_ar1[,1],rep(NA,69)),type='l',xlab='Date',ylab='AR(1)',main=paste('AR(1) τ =', round(mkt_share_ar1_tau[1],2)),bty='l',font.main = 1)
plot(date,c(rep(NA,33),mkt_share_ar1[,2],rep(NA,69)),type='l',xlab='Date',ylab='AR(1)',main=paste('AR(1) τ =', round(mkt_share_ar1_tau[2],2)),bty='l',font.main = 1)
plot(date,c(rep(NA,33),mkt_share_ar1[,3],rep(NA,69)),type='l',xlab='Date',ylab='AR(1)',main=paste('AR(1) τ =', round(mkt_share_ar1_tau[3],2)),bty='l',font.main = 1)

#Finally we plot the variance time series with the tau value
plot(date,c(rep(NA,33),mkt_share_vari[,1],rep(NA,69)),type='l',xlab='Date',ylab='Variance',main=paste('Variance τ =', round(mkt_share_vari_tau[1],2)),bty='l',font.main = 1)
plot(date,c(rep(NA,33),mkt_share_vari[,2],rep(NA,69)),type='l',xlab='Date',ylab='Variance',main=paste('Variance τ =', round(mkt_share_vari_tau[2],2)),bty='l',font.main = 1)
plot(date,c(rep(NA,33),mkt_share_vari[,3],rep(NA,69)),type='l',xlab='Date',ylab='Variance',main=paste('Variance τ =', round(mkt_share_vari_tau[3],2)),bty='l',font.main = 1)





