############## COPULA
library(dplyr)

datos <- read.csv("https://raw.githubusercontent.com/karen0c/Invest/main/Codigo%202.0/Datos/Prediccion_rf.csv")## RF

# Fijar semilla para reproducibilidad
set.seed(1234)

library(psych)
corr.test(datos$TCI, datos$WMA) # coef corr es 0.69

library(fitdistrplus)
y_v = as.vector(datos$WMA)
tci_v = as.vector(datos$TCI)
par(mfrow=c(1,1))


plotdist(y_v, histo = TRUE, demp = TRUE)
descdist(y_v, discrete=FALSE, obs.col = "red") ## 

fb <- fitdist(y_v, "beta")
summary(fb)
fg <- fitdist(y_v, "gamma")
fe <- fitdist(y_v, "exp")

par(mfrow = c(2, 2))
plot.legend <- c("beta", "gamma", "expo")
denscomp(list(fb, fg, fe), legendtext = plot.legend, cex = 0.8)
qqcomp(list(fb, fg, fe), legendtext = plot.legend, cex = 0.8)
cdfcomp(list(fb, fg, fe), legendtext = plot.legend, cex = 0.8)
ppcomp(list(fb, fg, fe), legendtext = plot.legend, cex = 0.8)


## pruebas
library(fitdistrplus)

# Realizar la prueba de Kolmogorov-Smirnov para una distribución beta
ks_beta <- ks.test(y_v,  "pbeta", shape1 = fb$estimate[1], shape2 = fb$estimate[2])
ks_beta #pvalor maypr al nivel de significancia, los datos siguen la distribución beta. H0:los datos siguen la distribución especificada.

# Realizar la prueba de Kolmogorov-Smirnov para una distribución gamma
ks_gamma <- ks.test(y_v, "pgamma", shape =fg$estimate[1] , rate = fg$estimate[2])
ks_gamma 

ks_exponencial <- ks.test(y_v, "pexp", rate = fe$estimate[1])
ks_exponencial

############################ TCI
par(mfrow = c(1, 1))
plotdist(tci_v, histo = TRUE, demp = TRUE)
descdist(tci_v, discrete=FALSE,  obs.col = "red")  ##  es gamma

fw3 <- fitdist(tci_v, "weibull")
fg3 <- fitdist(tci_v, "gamma")
fl3 <- fitdist(tci_v, "lnorm")
fn3 <- fitdist(tci_v, "beta")
summary(fl3)
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "lognormal", "normal")
denscomp(list(fw3, fg3, fl3, fn3), legendtext = plot.legend, cex = 0.8)
qqcomp(list(fw3, fg3, fl3, fn3), legendtext = plot.legend, cex = 0.8)
cdfcomp(list(fw3, fg3, fl3, fn3), legendtext = plot.legend, cex = 0.8)
ppcomp(list(fw3, fg3, fl3, fn3), legendtext = plot.legend, cex = 0.8)


par(mfrow = c(1, 1))

tciecdf <- ecdf(tci_v)(tci_v)
yecdf<-ecdf(y_v)(y_v)   #Empirical Cumulative Distribution Functio 
plot(tciecdf,yecdf,xlab="tci",ylab="y_pred",main="tci and Y_pred CDFs")


ks_weibull <- ks.test(tci_v, "pweibull", shape = fw3$estimate[1], scale =fw3$estimate[2] )
ks_weibull # no se rechaza

# Realizar la prueba de Kolmogorov-Smirnov para una distribución gamma
ks_gamma3 <- ks.test(tci_v, "pgamma", shape =fg3$estimate[1] , rate = fg3$estimate[2])
ks_gamma3 #no se rechaza

ks_lognormal <- ks.test(tci_v, "plnorm", meanlog = fl3$estimate[1], sdlog = fl3$estimate[2])
ks_lognormal # no se rechaza



############################ TCI y Y_pred ########################################3

##funcisn BiCopSelect() 
library(VineCopula)
library(copula)
var_b=pobs(y_v)
var_a2=pobs(tci_v)

library(VineCopula)
CopulaFit_tci<-BiCopSelect(var_a2,var_b, method="mle")
summary(CopulaFit_tci)

par(mfrow = c(1,1))
plot(tciecdf,yecdf,xlab="tci",ylab="y_pred",main="TCI and Y_pred CDFs")
CopulaSim2<-BiCopSim(1000, CopulaFit_tci$family, CopulaFit_tci$par)
points(CopulaSim2, col = rgb(0,0,1,alpha=0.3), pch=16)

#install.packages("E:/1/VC2copula_0.1.2.tar.gz")

#install.packages('VC2copula')
library(VC2copula)
library(copula)
#tawnCopula()

#cop_model_spei <-rotCopula(gumbelCopula(), flip = c(TRUE,FALSE)) ###########################################
#cop_model <- r90GumbelCopula(param = -1.34) ## qué parametro exactamente debo poner aquí?
# Crear un objeto de copula t-Student
cop_model2 <- normalCopula(param = 0.81, dim = 2, dispstr = "ex") ##estarán correctos los parametros,mas adelante me dice qeu df debe ser entero por lo que reemplazo 2.99 por 3


##Simulations:
my_dist_pr2 <- mvdc(cop_model2, margins = c("lnorm","beta"), 
                    paramMargins = list(list(meanlog=fl3$estimate[1], sdlog=fl3$estimate[2]), list(shape1=fb$estimate[1], shape2=fb$estimate[2])))
library(ggplot2)
t_cop=rCopula(1000,cop_model2)
colnames(t_cop)=c("TCI", "Y")
p_tci <- qplot(t_cop[,1], t_cop[,2], colour = t_cop[,2], main="Copula Gauss", xlab = "TCI", ylab = "Y")

library(grid)
print(p_tci, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
# More graphics.Generate random sample observations from the multivariate distribution

#######################################COPULA
mycopula2 <- cop_model2
u2 <- rCopula(2000,cop_model2)
# Compute the density
pdf2 <- dCopula(u2, mycopula2)
# Compute the CDF
cdf2 <- pCopula(u2, mycopula2)
# Generate random sample observations from the multivariate distribution
v2 <- rMvdc(2000, my_dist_pr2)
# Compute the density
pdf_mvd2 <- dMvdc(v2, my_dist_pr2)
# Compute the CDF
cdf_mvd2 <- pMvdc(v2, my_dist_pr2)

par(mfrow = c(1, 2))
#scatterplot3d(u[,1], u[,2], pdf_, color="blue", main="Density SPEI6 and Yield", xlab ="SPEI6", ylab="Yiedl", zlab="dCopula", pch=".")
persp(mycopula2, dCopula, main ="Density TCI and Y",xlab ="TCI", ylab="Y") 
contour(mycopula2, dCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot", xlab ="TCI", ylab="Y")

persp(mycopula2, pCopula, main = "CDF",  xlab = "TCI", ylab="Y")
contour(mycopula2, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot" , xlab = "TCI", ylab="Y")

library(scatterplot3d) 
scatterplot3d(v2[,1],v2[,2],pdf_mvd2,highlight.3d = T)
scatterplot3d(v2[,1],v2[,2],cdf_mvd2,highlight.3d = T)


################################### TCI - Y, QR ######################################################3
ns <- 10000 #sample size

xysTCI <- rMvdc(ns, my_dist_pr2)
par(mfrow = c(1, 1))
#Inspect data
plot(xysTCI[,1], xysTCI[,2], xlab="TCI", ylab="Y", main="Artificial Data TCI")

library(MASS)


##Fit a bivariate copula to the artificial data
betaMLE_2 <- fitdistr(xysTCI[,2], densfun="beta", list(shape1=fb$estimate[1], shape2=fb$estimate[2]))
betaShape2_1 <- betaMLE_2$estimate[1]
betaShape2_2 <- betaMLE_2$estimate[2]


#Obtain start values for fitting the copula
lnormMLE <- fitdistr(xysTCI[,1], densfun="lognormal")# no se pudo añadir valores iniciales
lnormMean<- lnormMLE$estimate[1]
lnormSd <- lnormMLE$estimate[2]


start2 <- c(lnormMean, lnormSd, betaShape2_1, betaShape2_2,0.81 )

fitCop2 <- fitMvdc(xysTCI,my_dist_pr2, start2)

#copula2<-fitCop2


library(dplyr)
xys11=as.data.frame(xysTCI)
plot(xys11)

colnames(xys11)<- c("tci","yield")
xpoints01=xys11%>%group_by(tci)%>%summarise(promedio=mean(yield))
head(xpoints01)
xpoints001 <- xpoints01[,1]
xpoints1 <- as.data.frame(xpoints001)
head(xpoints1)
#x=xpoints1$tci

#compute the median
ymed2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.5,
                                                    marginal=2, x=k))
y05_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.05,
                                                    marginal=2, x=k))

#compute the .1th-quantile
y15_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.15,
                                                   marginal=2, x=k))

#compute the .2th-quantile
y2_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.2,
                                                   marginal=2, x=k))
#compute the .3 th-quantile
y3_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.3,
                                                   marginal=2, x=k))
y7_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.7,
                                                   marginal=2, x=k))
y85_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.85,
                                                   marginal=2, x=k))
y9_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.9,
                                                   marginal=2, x=k))
y95_2 <- sapply(xpoints1$tci, function (k) qrCopula(copula=fitCop2, p=.95,
                                                    marginal=2, x=k))


#Plot the results
par(mfrow=c(1,1))
plot(xysTCI[,1], xysTCI[,2], xlab="TCI", ylab="Y", col=gray(.2,.2), main="Critical Quatiles on Marginals From Gaussian Copula")
#add the median to the plot
#lines(xpoints, ymed, lty=1, col="black", lwd=1.5)
#add the 90% prediction interval: approximately 20 (out of 200) points
#should lie outside the boundaries of the 90% prediction interval
lines( xpoints1$tci,y05_2, lty=2, col="firebrick3", lwd=2)
#lines(xpoints1$tci, y15_2, lty=1, col="pink", lwd=2)
lines(xpoints1$tci, y3_2, lty=2, col="dodgerblue4", lwd=2)
lines( xpoints1$tci,y7_2, lty=2, col="#CD6600", lwd=2)
#lines(xpoints1$tci, y85_2, lty=1, col="pink", lwd=2)
lines( xpoints1$tci,y95_2, lty=2, col="#9A32CD", lwd=2)
lines(xpoints1$tci, ymed2, lty=1, col="forestgreen", lwd=2)
x2=xpoints1$tci
str(x2)

#add a legend
legend("topright", lty=c( 2,1, 2,1, 2,1,2),
#       col=c( "firebrick3","pink","dodgerblue4", "forestgreen", "#CD6600","pink", "#9A32CD"),
#       legend=c( "5% quantile","15% quantile", "30% quantile", "50% quantile", "70% quantile","85% quantile", "95% quantile"),
        col=c( "firebrick3","dodgerblue4", "forestgreen", "#CD6600", "#9A32CD"),
        legend=c( "5% quantile", "30% quantile", "50% quantile", "70% quantile", "95% quantile"),

       bty="n", y.intersp=1.2, cex=.9)
# Convert data into data frame to manipulate subsets:
xys_df2=as.data.frame(xysTCI)
##
#ymed=as.data.frame(ymed)
x2<- as.data.frame(x2)
xymed_2 =cbind(x2, ymed2)
head(xymed_2)
colnames(xymed_2)=c("TCI", "Y")
#plot(xymed_2)
#

xy05_2 =cbind(x2, y05_2)
head(xy05_2)
colnames(xy05_2)=c("TCI", "Y")
#

xy15_2 =cbind(x2, y15_2)
colnames(xy15_2)=c("TCI", "Y")
head(xy15_2)

xy3_2 =cbind(x2, y3_2)
colnames(xy3_2)=c("TCI", "Y")
head(xy3_2)
#

xy7_2 =cbind(x2, y7_2)
colnames(xy7_2)=c("TCI", "Y")
tail(xy7_2)

xy85_2 =cbind(x2, y85_2)
colnames(xy85_2)=c("TCI", "Y")
tail(xy85_2)

xy95_2 =cbind(x2, y95_2)
colnames(xy95_2)=c("TCI", "Y")
tail(xy95_2)


par(mfrow=c(1,5))
#Datos ordenados con base en TCI
datos2 <- arrange(xys_df2,V1)
colnames(datos2)<-c("TCI", "Y")
tail(datos2)
QA_tci <- arrange(xy05_2,TCI)
tail(QA_tci)
plot(QA_tci)
QB_tci <- arrange(xy3_2,TCI)
tail(QB_tci)
plot(QB_tci)
#QMED_tci<-arrange(xymed_2,TCI)
#tail(QMED_tci)
#plot(QMED_tci)
QMED1_tci<-arrange(xy15_2,TCI)
tail(QMED1_tci)
plot(QMED1_tci)


QC_tci <- arrange(xy7_2,TCI)
tail(QC_tci)
plot(QC_tci)
QD_tci <- arrange(xy95_2,TCI)
tail(QD_tci)
plot(QD_tci)

QMED2_tci<-arrange(xy85_2,TCI)
tail(QMED2_tci)
plot(QMED2_tci)


##########

lower_limit_TCI<- subset(datos2, datos2$Y <= QA_tci$Y)
lower_threshold_TCI= subset(datos2, datos2$Y < QB_tci$Y & datos2$Y>QA_tci$Y)
lower_TCI <- subset(datos2, datos2$Y < QB_tci$Y )
mean_values_TCI<-subset(datos2, datos2$Y >= QB_tci$Y & datos2$Y<=QC_tci$Y)
upper_treshold_TCI=subset(datos2, datos2$Y < QD_tci$Y & datos2$Y>QC_tci$Y)
upper_limit_TCI <- subset(datos2, datos2$Y>= QD_tci$Y)
upper_TCI <- subset(datos2,  datos2$Y>QC_tci$Y)

par(mfrow=c(1,5))
#plot(datos2, main="TCI and Yield")
plot(lower_limit_TCI[,1:2], main="Zone A")
plot(lower_threshold_TCI[,1:2], main="Zone B")
plot(mean_values_TCI[,1:2], main ="Zone C")
plot(upper_treshold_TCI[,1:2], main="Zone D")
plot(upper_limit_TCI[,1:2], main="Zone E")


par(mfrow=c(1,1))
# Gráfico original (sin dividir en zonas)
plot(datos2, main="TCI and Y", col="black")

# Añade puntos de cada zona con colores diferentes
points(lower_limit_TCI[,1:2], col="firebrick3")
points(lower_threshold_TCI[,1:2], col="dodgerblue4")
points(mean_values_TCI[,1:2], col="forestgreen")
points(upper_treshold_TCI[,1:2], col="#CD6600")
points(upper_limit_TCI[,1:2], col="#9A32CD")

lines(xpoints1$tci, y15_2, lty=1, col="pink", lwd=2)
lines(xpoints1$tci, y85_2, lty=1, col="pink", lwd=2)

#add a legend
legend("topleft", pch=c(19, 19, 19, 19, 19),
       col=c( "firebrick3","dodgerblue4", "forestgreen", "#CD6600", "#9A32CD"),
       legend=c( "5% quantile", "30% quantile", "50% quantile", "70% quantile", "95% quantile"),
       bty="n", y.intersp=1.2, cex=.9)

describe(lower_limit_TCI)
describe(lower_threshold_TCI)
describe(mean_values_TCI)
describe(upper_treshold_TCI)
describe(upper_limit_TCI)

plot(datos2, main="TCI and Y", col="#838B8B")
points(lower_TCI[,1:2], col="#528B8B")
points(upper_TCI[,1:2], col="darkorange4")
lines(xpoints1$tci, y15_2, lty=1, col="pink", lwd=2)
lines(xpoints1$tci, y85_2, lty=1, col="pink", lwd=2)

####### ####################------------- Indemnities:

#C?lculo de indemnizaciones: Due to the lack of information of monthly yields, we derived a yield measure from production, by averaging the historic harvest area in Samana (2007-2020)
PrecioInternoCop= 2730000/125 #precio COP interno promedio 2023 café por kg
exchange=4500    #(COL/USD) 2023
CostoProduccion =PrecioInternoCop*0.65/exchange #Costo produccion entre promedio TRM 2023 USD/kg
Sum_insured=3750*CostoProduccion #kg/ha * usd/kg=USD/ha
max_payment= 3.75*0.24*CostoProduccion*1000   #ton/ha * *USD/kg
#HA_Samana=read.csv("stats.samana.csv", sep = ";", header=TRUE)
#prom_HA=3522.161429  # Hectareas (ha) dato 2007-2020
price= PrecioInternoCop*1000/exchange          #(US/ton) 
max_indemnity=2493*(1+0.1312)  #(US/ha) ###Ajustar
#k=mean(mean_values$V2)  #historical average production (ton)
#k           #305 ton   vs 1.091 ton/ha
####### ####################------------- Indemnities TCI:
#mean_values_TCI=mean_values_TCI*1000

k_tci  = mean_values_TCI%>%group_by(TCI)%>%summarise(Y=mean(Y))
k_y_tci = (k_tci)*price
k_y_tci          ## Strike: historical average yield in US/ha   4654.103 USD/ha
str(k_y_tci)
# Expected Shortfalls:
library(dplyr)
ES_A_tci<- lower_limit_TCI%>%group_by(TCI)%>%summarise(Y=mean(Y))
ES_B_tci<- lower_threshold_TCI%>%group_by(TCI)%>%summarise(Y=mean(Y))
ES_D_tci<- upper_treshold_TCI%>%group_by(TCI)%>%summarise(Y=mean(Y))
ES_E_tci<- upper_limit_TCI%>%group_by(TCI)%>%summarise(Y=mean(Y))

k_qr_tci = k_tci
#plot(xymed)
k_y_qr_tci=(k_qr_tci)*price
tail(k_y_qr_tci)            ## Strike: historical yield in US/ha

######################                 Indemnities

#km_tci=QMED_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))
km1_tci=QMED1_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))
km2_tci=QMED2_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))



ES_A1_tci<-ES_A_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))

ES_B1_tci<-ES_B_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))

ES_D1_tci<-ES_D_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))

ES_E1_tci<-ES_E_tci%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))


###### ZONA A   (Upper Limit)
ZonaA_TCI=merge(km1_tci,ES_A1_tci, by="TCI")
ZonaA_TCI
Indemnity_tci_lower_limit_qr= ifelse(abs(ZonaA_TCI$Y.x - ZonaA_TCI$Y.y) > 0,
                                     abs(ZonaA_TCI$Y.x - ZonaA_TCI$Y.y) , 0)
Premium_tci_lower_limit_qr=mean(Indemnity_tci_lower_limit_qr)
Premium_tci_lower_limit_qr
Premium_tci_lower_limit_qr*price            # 

indemnities_A_tci=cbind(ZonaA_TCI, Indemnity_tci_lower_limit_qr)
colnames(indemnities_A_tci) <- c("TCI", "k", "Miu*", "I_A")
head(indemnities_A_tci)


###### ZONA B:   (Lower Threshold)
ZonaB_TCI=merge(km1_tci,ES_B1_tci, by="TCI")
ZonaB_TCI
Indemnity_tci_lower_thre_qr= ifelse(abs(ZonaB_TCI$Y.x - ZonaB_TCI$Y.y) > 0,
                                    abs(ZonaB_TCI$Y.x - ZonaB_TCI$Y.y) , 0)
Premium_tci_lower_thre_qr=mean(Indemnity_tci_lower_thre_qr)
Premium_tci_lower_thre_qr
Premium_tci_lower_thre_qr*price            # 


indemnities_B_TCI=cbind(ZonaB_TCI, Indemnity_tci_lower_thre_qr)
colnames(indemnities_B_TCI) <- c("TCI", "k", "Miu*", "I_B")
head(indemnities_B_TCI)

##
###### ZONA D:(Upper Threshold)
ZonaD_TCI=merge(km2_tci,ES_D1_tci, by="TCI")
ZonaD_TCI
Indemnity_tci_upper_thre_qr= ifelse(abs(ZonaD_TCI$Y.x - ZonaD_TCI$Y.y) > 0,
                                    abs(ZonaD_TCI$Y.x - ZonaD_TCI$Y.y) , 0)
Premium_tci_upper_thre_qr=mean(Indemnity_tci_upper_thre_qr)
Premium_tci_upper_thre_qr
Premium_tci_upper_thre_qr*price            # 

indemnities_D_tci=cbind(ZonaD_TCI, Indemnity_tci_upper_thre_qr)
colnames(indemnities_D_tci) <- c("TCI", "k", "Miu*", "I_D")
head(indemnities_D_tci)
str(Indemnity_tci_upper_thre_qr)
###
###### ZONA E:(Upper limit)
ZonaE_TCI=merge(km2_tci,ES_E1_tci, by="TCI")
str(ZonaE_TCI)
Indemnity_tci_upper_lim_qr= ifelse(abs(ZonaE_TCI$Y.x - ZonaE_TCI$Y.y)> 0,
                                   abs(ZonaE_TCI$Y.x - ZonaE_TCI$Y.y) , 0)

Premium_tci_upper_limit_qr=mean(Indemnity_tci_upper_lim_qr)
Premium_tci_upper_limit_qr
Premium_tci_upper_limit_qr*price            # 

indemnities_E_tci=cbind(ZonaE_TCI, Indemnity_tci_upper_lim_qr)
colnames(indemnities_E_tci) <- c("TCI", "k", "Miu*", "I_E")
head(indemnities_E_tci)


##PREMIUMS LOWER: (ZONA A +ZONA B)
Premium_qr_all_lower_TCI= Premium_tci_lower_limit_qr + Premium_tci_lower_thre_qr
Premium_qr_all_lower_TCI*price    #532.78 usd/HA     vs gumbel 90: 349.4 USD/ha  vs GR): 323.94
Premium_qr_all_lower_TCI

###PREMIUMS UPPER: (ZONA D +ZONA E)
Premium_qr_all_upper_TCI= Premium_tci_upper_thre_qr +  Premium_tci_upper_limit_qr
Premium_qr_all_upper_TCI*price         #567.06 usd/ha          vs gumbel 90: 412.9063/ha vs: 384.76 (RG)
Premium_qr_all_upper_TCI

#Fair_insurance_premium=sum(Premium_spei_upper_thre_qr, Premium_spei_upper_limit_qr)*price
Fair_insurance_premium1_TCI = append(Indemnity_tci_upper_thre_qr, Indemnity_tci_upper_lim_qr)
Fair_insurance_premium2_TCI = append(Indemnity_tci_lower_limit_qr, Indemnity_tci_lower_thre_qr)
Fair_insurance_premium3_TCI = append(Fair_insurance_premium1_TCI, Fair_insurance_premium2_TCI) 
mean(Fair_insurance_premium3_TCI)
mean(Fair_insurance_premium3_TCI)*price/max_indemnity           #Tasa IR correcta, con todas las zonas

str(Indemnity_tci_lower_limit_qr)
str(Indemnity_tci_lower_thre_qr)
str(Indemnity_tci_upper_thre_qr)
str(Indemnity_tci_upper_lim_qr)

#----promediando todas las zonas
#mean(Indemnity_spei_upper_lim_qr)*price

###############Uninsured representa la zona de no disbuursals:
uninsured_TCI = mean_values_TCI
uninsuredAB_TCI = lower_TCI%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))
uninsuredDE_TCI = upper_TCI%>% mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))

str(uninsured_TCI)
## ----------------------------Insured Yields:
#mean_values_redon = mean_values%>%mutate(across(c('PR', 'Y'), round, 5))
#mean_values_redon_TCI = QMED_tci%>%mutate(across(c('TCI', 'Y'), round, 4))%>%group_by(TCI)%>%summarise(Y=mean(Y))
mean_values_redon1_TCI = QMED1_tci%>%mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))
mean_values_redon2_TCI = QMED2_tci%>%mutate(across(c('TCI', 'Y'), round, 5))%>%group_by(TCI)%>%summarise(Y=mean(Y))


#str(mean_values_redon_TCI)

#Indemnities organizados:
#Ind_D_TCI=cbind(Indemnity_tci_upper_thre_qr,ES_D1_tci)
Indemnity_tci_upper_lim_qr
##
head(indemnities_D_tci)
head(indemnities_E_tci)
indemnities_DE_tci=merge(indemnities_D_tci, indemnities_E_tci, by="TCI")[,-c(2,3,5,6)]
str(indemnities_DE_tci)
head(indemnities_DE_tci)

indemnities_AB_tci=merge(indemnities_A_tci, indemnities_B_TCI, by="TCI")[,-c(2,3,5,6)]
str(indemnities_AB_tci)
head(indemnities_AB_tci)


base_indenmities_TCI = merge(indemnities_DE_tci, mean_values_redon2_TCI, by="TCI")
head(base_indenmities_TCI)
str(base_indenmities_TCI)         ############Datos con yield insured y no insured (que es la zona media o C)

###########--------------#################333

base_indenmities2_TCI = merge(indemnities_AB_tci, mean_values_redon1_TCI, by="TCI")
head(base_indenmities2_TCI)
###########--------------#####################
##indemnities_E 
## ----------------------------Insured & Uninsured premiums:
Premium_tci_lower_limit_qr
Premium_tci_lower_thre_qr
Premium_tci_upper_thre_qr*price
Premium_tci_upper_limit_qr*price

#yield_insured:
Y_Insured_TCI_DE= base_indenmities_TCI%>% mutate(Y_ins = Y+I_D+I_E-Premium_qr_all_upper_TCI)
head(Y_Insured_TCI_DE)
str(Y_Insured_TCI_DE)

#yield_insured AB:
Y_Insured_TCI_AB= base_indenmities2_TCI%>% mutate(Y_ins = Y+I_A+I_B-Premium_qr_all_lower_TCI)
head(Y_Insured_TCI_AB)
str(Y_Insured_TCI_AB)

##Uninsured Average
Y_uninsured_gorro_TCI = mean(uninsured_TCI$Y)
Y_uninsured_gorroAB_TCI = mean(uninsuredAB_TCI$Y)
Y_uninsured_gorroDE_TCI = mean(uninsuredDE_TCI$Y)
##  -------------------------------------------------            Uninsured SV:
##Operacion diferencia cuadr?tica Yield_uninsured menos Y_uninsured_gorro
SV_Uninsured_0_AB_tci=ifelse(uninsuredAB_TCI$Y - Y_uninsured_gorro_TCI != 0,1,0)
SV_Uninsured_0_AB_tci
SV_uninsured_00_AB_tci=cbind(uninsuredAB_TCI, SV_Uninsured_0_AB_tci)
tail(SV_uninsured_00_AB_tci)
SV_uninsured_1_AB_tci=subset(SV_uninsured_00_AB_tci, SV_Uninsured_0_AB_tci == 1)


SV_unsured2_AB_tci= SV_uninsured_1_AB_tci %>% mutate(SV_un_dif = (Y- Y_uninsured_gorro_TCI)^2)
SV_uninsured_AB_tci=mean(SV_unsured2_AB_tci$SV_un_dif)
SV_uninsured_AB_tci                                  ######SV uninsured
SV_uninsured_AB_tci*price



SV_Uninsured_0_DE_tci=ifelse(uninsuredDE_TCI$Y - Y_uninsured_gorro_TCI !=0,1,0)
SV_Uninsured_0_DE_tci
SV_uninsured_00_DE_tci=cbind(uninsuredDE_TCI, SV_Uninsured_0_DE_tci)
tail(SV_uninsured_00_DE_tci)
SV_uninsured_1_DE_tci=subset(SV_uninsured_00_DE_tci, SV_Uninsured_0_DE_tci == 1)


SV_unsured2_DE_tci= SV_uninsured_1_DE_tci %>% mutate(SV_un_dif = (Y- Y_uninsured_gorro_TCI)^2)
SV_uninsured_DE_tci=mean(SV_unsured2_DE_tci$SV_un_dif)
SV_uninsured_DE_tci                                  ######SV uninsured
SV_uninsured_DE_tci*price


###############################################
##   -----------------------------------------------          insured SV UPPER:
##Operacion diferencia cuadr?tica Yield_insured menos Y_uninsured_gorro
SV_insured_0_tci_DE=ifelse(Y_Insured_TCI_DE$Y_ins - Y_uninsured_gorroDE_TCI !=0,1,0)
SV_insured_0_tci_DE
SV_insured_00_tci_DE=cbind(Y_Insured_TCI_DE, SV_insured_0_tci_DE)
tail(SV_insured_00_tci_DE)
SV_insured_1_tci_DE=subset(SV_insured_00_tci_DE, SV_insured_0_tci_DE == 1)
head(SV_insured_1_tci_DE) 

SV_insured2_tci_DE= SV_insured_1_tci_DE %>% mutate(SV_ins_dif = (Y_ins- Y_uninsured_gorroDE_TCI)^2)
SV_insured_tci_DE=mean(SV_insured2_tci_DE$SV_ins_dif)
SV_insured_tci_DE                                  ######SV uninsured
SV_insured_tci_DE*price

##   -----------------------------------------------          insured SV LOWER:
##Operacion diferencia cuadr?tica Yield_insured menos Y_uninsured_gorro
SV_insured_0_tci_AB=ifelse(Y_Insured_TCI_AB$Y_ins - Y_uninsured_gorroAB_TCI !=0,1,0)
SV_insured_0_tci_AB
SV_insured_00_tci_AB=cbind(Y_Insured_TCI_AB, SV_insured_0_tci_AB)
tail(SV_insured_00_tci_AB)
SV_insured_1_tci_AB=subset(SV_insured_00_tci_AB, SV_insured_0_tci_AB == 1)
head(SV_insured_1_tci_AB) 

SV_insured2_tci_AB= SV_insured_1_tci_AB %>% mutate(SV_ins_dif = (Y_ins- Y_uninsured_gorroAB_TCI)^2)
SV_insured_tci_AB=mean(SV_insured2_tci_AB$SV_ins_dif)
SV_insured_tci_AB                                  ######SV uninsured
SV_insured_tci_AB*price

#Hedging Effectiveness:
HE_qr_TCI_DE = (SV_insured_tci_DE- SV_uninsured_DE_tci)/SV_uninsured_DE_tci
HE_qr_TCI_DE*100
#HE_qr_TCI = 1- (SV_insured_2/SV_uninsured) ###es esta opcion correcta?######################33


HE_qr_TCI_AB = (SV_insured_tci_AB- SV_uninsured_AB_tci)/SV_uninsured_AB_tci
HE_qr_TCI_AB*100




