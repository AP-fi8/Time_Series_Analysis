library(data.table)
library(forecast)


data =  fread("AmtrakBig_raw.csv")
data<-data[,Month:=as.Date(Month,format="%m/%d/%Y")]


dmy <- caret::dummyVars(" ~ Season", data = data, fullRank = T) 
dat_transformed <- data.frame(predict(dmy, newdata = data)) 

my_data<-data.table(cbind.data.frame(data,dat_transformed))


train<-my_data[Month>="1991-01-01" & Month<="2001-07-01"]
test<-my_data[Month>"2001-07-01"]



# trial model ARIMA
First_differeneced_Ridership_data<-diff(train[,Ridership],1)
acf(First_differeneced_Ridership_data)
pacf(First_differeneced_Ridership_data,36)

my_arima<-Arima(train[,Ridership],c(6,1,6),
                include.mean=T,
                seasonal=list(order=c(1,0,1),period=3))
my_arima
# trail model ARIMAX
my_arimax<-Arima(train[,Ridership],c(6,1,6),
                include.mean=T,
                seasonal=list(order=c(1,0,1),period=3),
                xreg=as.matrix(train[,5:ncol(my_data)]))
my_arimax
############################################
xreg=as.matrix(train[,5:ncol(my_data)])
write.csv(xreg,file="dummies.csv")
############################################
frcst.arima<-forecast(my_arima,h=51)
frcst.arimax.test<-forecast(my_arimax,h=39,
                            xreg=as.matrix(test[,5:ncol(my_data)]))
frcst.arimax<-forecast(my_arimax,h=12,
                       xreg=as.matrix(test[,5:ncol(my_data)]))

data_SARIMA<-c(data[,Ridership],rep(NA,12))
plot(data_SARIMA,type="l",main="SARIMA")
lines(my_arima$fitted,col="blue")
lines(frcst.arima$mean,col="red")

data_SARIMAX<-c(data[,Ridership],rep(NA,12))
plot(data_SARIMAX,type="l",main="SARIMAX")
lines(my_arimax$fitted,col="blue")
lines(frcst.arimax.test$mean,col="red")
lines(frcst.arimax$mean,col="green")

data_lr<-c(data[,Ridership],rep(NA,12))
plot(data_lr,type="l",main="LINEAR REGRESSION")
lines(my_lm$fitted.values,col="blue")
lines(c(rep(NA,120),frcst.lm),col="red")

data<-data[,`:=`(arima_projections=c(my_arima$fitted,frcst.arima$mean),
                 arimax_projections=c(my_arimax$fitted,frcst.arimax$mean),
                 lr_projections=c(my_lm$fitted.values,frcst.lm))]
#write.csv(data, file = "tsdata.csv")

val.df<-data[(Month>"2001-07-01") & (!is.na(Ridership))]
train.df<-data[Month>="1991-01-01" & Month<="2001-07-01"]


apply(val.df[,5:ncol(val.df)],2,
      function(x) MLmetrics::RMSE(x,val.df[,Ridership]))
apply(train.df[,5:ncol(val.df)],2,
      function(x) MLmetrics::RMSE(x,train.df[,Ridership]))
###############################################################################
Dm <- ts
plot.ts(mdata)
Rider <- ts(mdata, frequency = 12, start = c(1991,1))
Rider
decmp<-decompose(Rider)
plot(decmp)
seasonal<-decmp$seasonal

plot(seasonal[0:12],type="l")


plot(Rider)
scatter.smooth(mdata)

############################ checking for randomness ###########################
n=127
u=0
for(i in 2:(n-1))
{
  if(((mdata$Ridership[i]>mdata$Ridership[i+1])
      &&(mdata$Ridership[i]>mdata$Ridership[i-1]))
     ||((mdata$Ridership[i]<mdata$Ridership[i+1])
       &&(mdata$Ridership[i]<mdata$Ridership[i-1]))) {
    u=u+1
  }
}
u
e_u=2*(n-2)/3
v_u=(16*n-29)/90
z=(u-e_u)/sqrt(v_u)
z
qnorm(0.975)

################ checking presence of trend, relative ordering test ############
rltv_oderng_test=function(y){
  cnt=0
  for(i in 1:(length(y)-1)){
    for(j in (i+1):length(y)){
      if(y[i]>y[j]){
        cnt=cnt+1
      }
    }
  }
  exp_q=length(y)*(length(y)-1)/4
  tau=1-4*(cnt/(length(y)*(length(y)-1)))
  V_tau=2*(2*length(y)+5)/(9*length(y)*(length(y)-1))
  tst_stat=tau/sqrt(V_tau)
  z_alpha=qnorm(0.025)
  if(abs(tst_stat)>abs(z_alpha)){
    cat("Reject Null Hypothesis")
  }
  else {
    cat("Accept Null Hypothesis")
  }
}
rltv_oderng_test(mdata$Ridership)

##################### checking for seasonality #################################
calc_chi=function(m,c,r){
  x=0
  den=c*(r+1)
  num=c*(r+1)*0.5
  for(i in 1:r){
    x=x+((m[i]-num)*(m[i]-num))
  }
  return(x/den)
}
sesnlty_check=function(y){
  dt_mtrx= matrix(y,nrow=12)
  dt_rankedmtrx=matrix(0,nrow=12,ncol=ncol(dt_mtrx))
  for(i in 1:ncol(dt_mtrx)){
    dt_rankedmtrx[,i]=rank(dt_mtrx[,i])
  }
  m_i=rowSums(dt_rankedmtrx)
  chi=calc_chi(m_i,ncol(dt_rankedmtrx),nrow(dt_rankedmtrx))
  tab_Chi=qchisq(.95,nrow(dt_mtrx)-1)
  if(chi>tab_Chi){
    cat("data shows the presence of seansonality",chi,tab_Chi)
  }
  else{
    cat("Based on the data we may say that there is no presence of sea-
sonality in the data")
  }
}
sesnlty_check(mdata$Ridership)
acf(mdata$Ridership,lag.max=40,main="Auto Correlation
Function Plot")

################# Estimation and Elimination of Trend #########################
trnd_est=function(x){
  trnd=c(rep(0,(length(x)-12)))
  for(i in 7:(length(x)-6)){
    trnd[i-6]=(1/12)*((0.5)*x[i-6]+x[i-5]+x[i-4]+x[i-3]+x[i-2]
                      +x[i-1]+x[i]+x[i+1]+x[i+2]+x[i+3]+x[i+4]
                      +x[i+5]+(0.5)*x[i+6])
  }
  return(trnd)
}
estmd_trnd=trnd_est(mdata$Ridership)
estmd_trnd
length(estmd_trnd)
plot(estmd_trnd,col="black",main = "Plot of Trend component",ylab
     = "estimated trend component")
remv_trnd=function(x,y){
  trnd_rmvd_dt=x[7:121]-y
  return(trnd_rmvd_dt)
}
trnd_rmvd_dt=remv_trnd(mdata$Ridership,estmd_trnd)
plot(trnd_rmvd_dt,type = "o",col= "black",ylab= "Data points",
     main= "plot after removing trend") 

################## estimation and elimination of seasonality ###################
differencing=function(x){
  diff=c(rep(0,length(x)-12))
  for(i in 1:(length(x)-12))
  {
    diff[i]=x[12+i]-x[i]
  }
  return(diff)
}
desesnlizd_detrndd_dt=differencing(trnd_rmvd_dt)
sesnlty_check(desesnlizd_detrndd_dt)
acf(desesnlizd_detrndd_dt,lag.max=30,main="Auto Correlation Func-
tion Plot")

plot(differencing(mdata$Ridership),type = "o",col="black",
     main = "Plot of Seasonal component",
     ylab = "estimated seasonality")

############################## white noise test ################################
rho_cal=function(x,h,n){
  gammah=0
  xbarn=mean(x)
  c=n-h
  for(i in 1:c){
    gammah=gammah+((x[i]-xbarn)*(x[i+h]-xbarn))
  }
  return(gammah)
}
rho = rho_cal(desesnlizd_detrndd_dt,1,
              length(desesnlizd_detrndd_dt))/rho_cal(desesnlizd_detrndd_dt,0,
                                                length(desesnlizd_detrndd_dt))
rho

test_stat_for_whitenoise = sqrt(length(desesnlizd_detrndd_dt))*rho
test_stat_for_whitenoise

######################### plot acf and pacf ####################################
acf(desesnlizd_detrndd_dt,lag.max=30,main="Auto Correlation Function Plot")
pacf(desesnlizd_detrndd_dt,lag.max=30,main="Partial Auto Correlation Function
     Plot")

##################################### AIC ######################################
#len <- 13
#aic=array(0,dim=c(len,len))
#for(i in 0:(len-1))
#{
  #for(j in 0:(len-1))
  #{
    #aic[(i+1),(j+1)]=0
    #model<-arima(x=desesnlizd_detrndd_dt,order=c(i,0,j),method = "ML",
                 #optim.control=list(maxit=1000),optim.method
                 #= "BFGS")
    #aic[(i+1),(j+1)]=model$aic
  #}
#}
################################################

############################# model identification #############################
desesnlizd_detrndd_dt_df<-data.frame(desesnlizd_detrndd_dt)

train=desesnlizd_detrndd_dt_df[(1:(.80*nrow(desesnlizd_detrndd_dt_df))),]
test=desesnlizd_detrndd_dt_df[-(1:(.80*nrow(desesnlizd_detrndd_dt_df))),]

auto_1 = auto.arima(train,seasonal=TRUE)
summary(auto_1)

########### Estimates of model and noise parameters of ARMA ####################
model<-arima(x=train,order=c(1,0,0),method = "ML",seasonal=TRUE)
summary(model)
############# (checking for stationarity) ######################################
polyroot(c(1,-1*model$coef[1],- 1*model$coef[2],-1*model$coef[3],
           - 1*model$coef[4]))

############## forcasting model ################################################

fcast<-forecast(model, h=12)
autoplot(fcast)
fcast
fcast<-data.frame(fcast)
fcast
#Point_fcast<-data.frame(fcast$Point.Forecast)
#Point_fcast
#test$Ridership
#desesnlizd_detrndd_dt
estmd_trnd
r=lm(estmd_trnd~seq(1,115))
summary(r)

auto_2=auto.arima(mdata,seasonal=TRUE)
fcast1<-forecast(auto_2,h=32)
autoplot(fcast1)


library(datasets); library(xts); library(forecast)

fit <- ets()
season <- fit$states[,"s1"]

plot(as.xts(season), major.format = "%Y-%m", auto.grid=F)
################# trend #################
tre
BIC(trend5)nd<-lm(estmd_trnd~seq(7,121))
BIC(trend)
trend2<-lm(estmd_trnd~seq(7,121)+I(seq(7,121)^2))
BIC(trend2)
trend3<-lm(estmd_trnd~seq(7,121)+I(seq(7,121)^2)+I(seq(7,121)^3))
BIC(trend3)
trend4<-lm(estmd_trnd~seq(7,121)+I(seq(7,121)^2)+I(seq(7,121)^3)+I(seq(7,121)^4))
BIC(trend4)
summary(trend4)##############4
trend5<-lm(estmd_trnd~seq(7,121)+I(seq(7,121)^2)+I(seq(7,121)^3)+I(seq(7,121)^4)
           +I(seq(7,121)^5))
############ seasonality ########################
model1<-lm(seasonal[1:12]~seq(1,12))
BIC(model1)
model2<-lm(seasonal[1:12]~seq(1,12)+I(seq(1,12)^2))
BIC(model2)
model3<-lm(seasonal[1:12]~seq(1,12)+I(seq(1,12)^2)+I(seq(1,12)^3))
BIC(model3)
summary(model3)################3
model4<-lm(seasonal[1:12]~seq(1,12)+I(seq(1,12)^2)+I(seq(1,12)^3)+I(seq(1,12)^4))
BIC(model4)
model5<-lm(seasonal[1:12]~seq(1,12)+I(seq(1,12)^2)+I(seq(1,12)^3)+I(seq(1,12)^4)
           +I(seq(1,12)^5))
BIC(model5)

predval=c(rep(1659,12))+c(18.82*seq(128,139))+
  c(-0.6708*seq(128,139)^2)+c(0.007655*seq(128,139)^3)+
  c(-0.00002679*seq(128,139)^4)+c(rep(-503.7277,12))+
+c(seq(1,12)*257.0195)+c((seq(1,12)^2)*-34.4947)+c((seq(1,12)^3)*1.3838)

predval1=c(rep(1659,12))+c(18.82*seq(140,151))+
  c(-0.6708*seq(140,151)^2)+c(0.007655*seq(140,151)^3)+
  c(-0.00002679*seq(140,151)^4)+c(rep(-503.7277,12))+
  c(seq(1,12)*257.0195)+c(-34.4947*(seq(1,12)^2))+c((seq(1,12)^3)*1.3838)

predval2=c(rep(1659,8))+c(18.82*seq(152,159))+
  c(-0.6708*seq(152,159)^2)+c(0.007655*seq(152,159)^3)+
  c(-0.00002679*seq(152,159)^4)+c(rep(-503.7277,8))+
  c(seq(1,8)*257.0195)+c(-34.4947*(seq(1,8)^2))+c((seq(1,8)^3)*1.3838)

predval
predval1
predval2
results<-c(predval,predval1,predval2)

auto.arima(mdata,seasonal="TRUE")
plot(results,type="l")
scatter = ggplot(test, aes(x=seq(1,32), y= Ridership))+
  labs(title = "Scatter Plot of Monthly Amtrak mdata Ridership",
       x="Time points",
       y="Number of Riders travelled")+geom_line()
scatter
p = ggplot() + 
  geom_line(data = prescription1, aes(x = dates, y = Difference), color = "blue")
+
  geom_line(data = prescription2, aes(x = dates, y = Difference), color = "red")
+
  xlab('Dates') +
  ylab('percent.change')
prediction_f=data.frame(prediction)

test_df<-data.frame(test)
test_df
fcast<-data.frame(fcast)
colnames(fcast)[colnames(fcast)=="Point Forecast"] <- "Point_Forecast"

p = ggplot() + 
  geom_line(data = test_df, aes(x = seq(1,27), y = test), color = "blue") +ylim(c(0,5000))
  geom_line(data = fcast, aes(x = seq(1,27), y = Point.Forecast), color = "red") +
  xlab('Dates') +
  ylab('percent.change')
p

fcast
test_df
