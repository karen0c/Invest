############## COPULA
library(dplyr)

datos <- read.csv("https://raw.githubusercontent.com/karen0c/Invest/main/Codigo%202.0/Datos/Prediccion_rf.csv")## RF

# Fijar semilla para reproducibilidad
set.seed(1234)


 #datos$pr[13]<-0.0000000000001 #PARA LSTM
 #datos$pr[14]<-0.0000000000001 #PARA LSTM

library(psych)
corr.test(datos$pr, datos$WMA, method = "spearman") # coef corr es 0.98, #spearman no lineal

library(fitdistrplus)
pr_v=as.vector(datos$pr)%>%na.omit()
y_v = as.vector(datos$WMA)
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

######################### PR
par(mfrow = c(1, 1))
plotdist(pr_v, histo = TRUE, demp = TRUE)
descdist(pr_v, discrete=FALSE,  obs.col = "red")  #

fb2 <- fitdist(pr_v, "beta")
summary(fb2)
fg2 <- fitdist(pr_v, "gamma")

par(mfrow = c(2, 2))
plot.legend <- c("beta", "gamma")
denscomp(list(fb2, fg2), legendtext = plot.legend, cex = 0.8)
qqcomp(list(fb2, fg2), legendtext = plot.legend, cex = 0.8)
cdfcomp(list(fb2, fg2), legendtext = plot.legend, cex = 0.8)
ppcomp(list(fb2, fg2), legendtext = plot.legend, cex = 0.8)


# Realizar la prueba de Kolmogorov-Smirnov para una distribución beta
ks_beta2 <- ks.test(pr_v,  "pbeta", shape1 = fb2$estimate[1], shape2 = fb2$estimate[2])
ks_beta2 #pvalor maypr al nivel de significancia, los datos siguen la distribución beta. H0:los datos siguen la distribución especificada.

# Realizar la prueba de Kolmogorov-Smirnov para una distribución gamma
ks_gamma2 <- ks.test(pr_v, "pgamma", shape =fg2$estimate[1] , rate = fg2$estimate[2])
ks_gamma2


par(mfrow = c(1, 1))

precdf<-ecdf(pr_v)(pr_v)  #Empirical Cumulative Distribution 
yecdf<-ecdf(y_v)(y_v)   #Empirical Cumulative Distribution Functio 
plot(precdf,yecdf,xlab="pr",ylab="y_pred",main="pr and Y_pred CDFs")
#plot(tciecdf,yecdf,xlab="tci",ylab="y_pred",main="tci and Y_pred CDFs")





##funcisn BiCopSelect() 
library(VineCopula)
library(copula)
var_a=pobs(pr_v)
var_b=pobs(y_v)

CopulaFit_pr<-BiCopSelect(var_a,var_b, method="mle")
summary(CopulaFit_pr)


CopulaSim<-BiCopSim(5000, CopulaFit_pr$family, CopulaFit_pr$par)
points(CopulaSim, col = rgb(0,0,1,alpha=0.3), pch=16)

#install.packages("E:/1/VC2copula_0.1.2.tar.gz")

#install.packages('VC2copula')
library(VC2copula)
#install.packages('VC2copula')
#tawnCopula()

#cop_model_spei <-rotCopula(gumbelCopula(), flip = c(TRUE,FALSE)) ###########################################
#cop_model <- r90GumbelCopula(param = -1.34) ## qué parametro exactamente debo poner aquí?
cop_model <- frankCopula(param = 20.46 , dim = 2)

##Simulations:                               "pr", "y"
my_dist_pr3 <- mvdc(cop_model, margins = c("beta","beta"), 
                    paramMargins = list(list(shape1=fb2$estimate[1], shape2=fb2$estimate[2]), list(shape1=fb$estimate[1], shape2=fb$estimate[2])))

library(grid)
library(ggplot2)
frankC=rCopula(1000,cop_model)
colnames(frankC)=c("PR", "Y")
p_pr <- qplot(frankC[,1], frankC[,2], colour = frankC[,2], main="Frank Copula Random Samples", xlab = "PR", ylab = "Y")
print(p_pr, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) 
# More graphics.Generate random sample observations from the multivariate distribution



mycopula <- cop_model
u <- rCopula(2000,cop_model)
# Compute the density
pdf_ <- dCopula(u, mycopula)
# Compute the CDF
cdf <- pCopula(u, mycopula)
# Generate random sample observations from the multivariate distribution
v <- rMvdc(2000, my_dist_pr3)
# Compute the density
pdf_mvd <- dMvdc(v, my_dist_pr3)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist_pr3)

par(mfrow = c(1,2))
# 3D plain scatterplot of the density, plot of the density and contour plot
library(scatterplot3d)
#scatterplot3d(u[,1], u[,2], pdf_, color="blue", main="Density SPEI6 and Yield", xlab ="SPEI6", ylab="Yiedl", zlab="dCopula", pch=".")
persp(mycopula, dCopula, main ="Density PR and Y",xlab ="PR", ylab="Y") 
contour(mycopula, dCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot", xlab ="PR", ylab="Y")

persp(mycopula, pCopula, main = "CDF",  xlab = "PR", ylab="Y")
contour(mycopula, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot" , xlab = "PR", ylab="Y")


scatterplot3d(v[,1],v[,2],pdf_mvd,highlight.3d = T)
scatterplot3d(v[,1],v[,2],cdf_mvd,highlight.3d = T)


##### Frank

ns <- 10000 #sample size
xys <- rMvdc(ns, my_dist_pr3)
par(mfrow = c(1, 1))
plot(xys)
#Inspect data
plot(xys[,1], xys[,2], xlab="PR", ylab="Y", main="Artificial Data")

library(MASS)


##Fit a bivariate copula to the artificial data

#Obtain start values for fitting the copula
betaMLE_1 <- fitdistr(xys[,1], densfun="beta", list(shape1=fb2$estimate[1], shape2=fb2$estimate[2]))
betaShape1 <- betaMLE_1$estimate[1]
betaShape2 <- betaMLE_1$estimate[2]

betaMLE_2 <- fitdistr(xys[,2], densfun="beta", list(shape1=fb$estimate[1], shape2=fb$estimate[2]))
betaShape2_1 <- betaMLE_2$estimate[1]
betaShape2_2 <- betaMLE_2$estimate[2]


start <- c(betaShape1, betaShape2, betaShape2_1, betaShape2_2, 20.46)

fitCop <- fitMvdc(xys,my_dist_pr3, start)


###################------------------
##Quantile regression: regress y(Yield) on x(PR)


library(dplyr)
xys1=as.data.frame(xys)
plot(xys1)

colnames(xys1)<- c("pr","y")
xpoints0=xys1%>%group_by(pr)%>%summarise(promedio=mean(y))
head(xpoints0)
xpoints00 <- xpoints0[,1]
xpoints <- as.data.frame(xpoints00)
head(xpoints)
#x=xpoints$pr


#compute the median
ymed <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.5,
                                                    marginal=2, x=k))
y05 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.05,
                                               marginal=2, x=k))

#compute the .15 th-quantile
y15 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.15,
                                                    marginal=2, x=k))

#compute the .2th-quantile
y2 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.2,
                                                    marginal=2, x=k))
#compute the .3 th-quantile
y3 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.3,
                                                   marginal=2, x=k))
y7 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.7,
                                                    marginal=2, x=k))
y85 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.85,
                                               marginal=2, x=k))
y9 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.9,
                                               marginal=2, x=k))
y95 <- sapply(xpoints$pr, function (k) qrCopula(copula=fitCop, p=.95,
                                               marginal=2, x=k))


#mean(y30)
#Plot the results
par(mfrow=c(1,1))
plot(xys[,1], xys[,2], xlab="PR", ylab="Y", col=gray(.2,.2), main="Critical Quatiles on Marginals From Frank Copula")
#add the median to the plot
#lines(xpoints, ymed, lty=1, col="black", lwd=1.5)
#add the 90% prediction interval: approximately 20 (out of 200) points
#should lie outside the boundaries of the 90% prediction interval
lines( xpoints$pr,y05, lty=2, col="firebrick3", lwd=2)
#lines(xpoints$pr, y15, lty=1, col="pink", lwd=2)
lines(xpoints$pr, y3, lty=2, col="dodgerblue4", lwd=2)
lines( xpoints$pr,y7, lty=2, col="#CD6600", lwd=2)
#lines(xpoints$pr, y85, lty=1, col="pink", lwd=2)
lines( xpoints$pr,y95, lty=2, col="#9A32CD", lwd=2)
lines(xpoints$pr, ymed, lty=1, col="forestgreen", lwd=2)
x=xpoints$pr
str(x)

#add a legend
legend("topright", lty=c( 2,1, 2,1, 2,1,2),
#       col=c( "firebrick3","pink","dodgerblue4", "forestgreen", "#CD6600","pink", "#9A32CD"),
#       legend=c( "5% quantile","15% quantile", "30% quantile", "50% quantile", "70% quantile","85% quantile", "95% quantile"),
        col=c( "firebrick3","dodgerblue4", "forestgreen", "#CD6600", "#9A32CD"),
        legend=c( "5% quantile", "30% quantile", "50% quantile", "70% quantile", "95% quantile"),
       bty="n", y.intersp=1.2, cex=.9)
# Convert data into data frame to manipulate subsets:
xys_df=as.data.frame(xys)
##
#ymed=as.data.frame(ymed)
x<- as.data.frame(x)
xymed =cbind(x, ymed)
head(xymed)
colnames(xymed)=c("PR", "Y")
plot(xymed)
#

xy05 =cbind(x, y05)
head(xy05)
colnames(xy05)=c("PR", "Y")
#

xy15 =cbind(x, y15)
colnames(xy15)=c("PR", "Y")
head(xy15)
#

xy3 =cbind(x, y3)
colnames(xy3)=c("PR", "Y")
head(xy3)
#

xy7 =cbind(x, y7)
colnames(xy7)=c("PR", "Y")
tail(xy7)

xy85 =cbind(x, y85)
colnames(xy85)=c("PR", "Y")
tail(xy85)

xy95 =cbind(x, y95)
colnames(xy95)=c("PR", "Y")
tail(xy95)



#Datos ordenados con base en SPEI
datos <- arrange(xys_df,V1)
colnames(datos)<-c("PR", "Y")
tail(datos)
QA <- arrange(xy05,PR)
tail(QA)
plot(QA)
QB <- arrange(xy3,PR)
tail(QB)
#QMED<-arrange(xymed,PR)
#tail(QMED)
#plot(QMED)
QMED1<-arrange(xy15,PR)
tail(QMED1)
plot(QMED1)


QC <- arrange(xy7,PR)
tail(QC)
QD <- arrange(xy95,PR)
tail(QD)
plot(QD)

QMED2<-arrange(xy85,PR)
tail(QMED2)
plot(QMED2)

##########

lower_limit<- subset(datos, datos$Y <= QA$Y)
lower_threshold= subset(datos, datos$Y < QB$Y & datos$Y>QA$Y)
lower <- subset(datos, datos$Y < QB$Y)
mean_values<-subset(datos, datos$Y >= QB$Y & datos$Y<=QC$Y)
upper_treshold=subset(datos, datos$Y < QD$Y & datos$Y>QC$Y)
upper_limit <- subset(datos, datos$Y>= QD$Y)
upper<-subset(datos, datos$Y > QC$Y)

par(mfrow=c(1,5))
#plot(datos, main="PR and Y")
plot(lower_limit[,1:2], main="Zone A")
plot(lower_threshold[,1:2], main="Zone B")
plot(mean_values[,1:2], main ="Zone C")
plot(upper_treshold[,1:2], main="Zone D")
plot(upper_limit[,1:2], main="Zone E")

par(mfrow=c(1,1))
# Gráfico original (sin dividir en zonas)
plot(datos, main="PR and Y", col="black")

# Añade puntos de cada zona con colores diferentes
points(lower_limit[,1:2], col="firebrick3")
points(lower_threshold[,1:2], col="dodgerblue4")
points(mean_values[,1:2], col="forestgreen")
points(upper_treshold[,1:2], col="#CD6600")
points(upper_limit[,1:2], col="#9A32CD")

lines(xpoints$pr, y15, lty=1, col="pink", lwd=2)
lines(xpoints$pr, y85, lty=1, col="pink", lwd=2)

#add a legend
legend("topright", pch=c(19, 19, 19, 19, 19),
       col=c( "firebrick3","dodgerblue4", "forestgreen", "#CD6600", "#9A32CD"),
       legend=c( "5% quantile", "30% quantile", "50% quantile", "70% quantile", "95% quantile"),
       bty="n", y.intersp=1.2, cex=.9)

describe(lower_limit)
describe(lower_threshold)
describe(mean_values)
describe(upper_treshold)
describe(upper_limit)

plot(datos, main="PR and Y", col="#838B8B")
points(lower[,1:2], col="#528B8B")
points(upper[,1:2], col="darkorange4")
lines(xpoints$pr, y15, lty=1, col="pink", lwd=2)
lines(xpoints$pr, y85, lty=1, col="pink", lwd=2)
####### ####################------------- Indemnities PR:

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
k  = mean_values%>%group_by(PR)%>%summarise(Y_datos=mean(Y))
k_y = (k)*price
k_y          ## Strike: historical average yield in US/ha   4654.103 USD/ha
str(k_y)
# Expected Shortfalls:
library(dplyr)
ES_A<- lower_limit%>%group_by(PR)%>%summarise(Y=mean(Y))
ES_B<- lower_threshold%>%group_by(PR)%>%summarise(Y=mean(Y))
ES_D<- upper_treshold%>%group_by(PR)%>%summarise(Y=mean(Y))
ES_E<- upper_limit%>%group_by(PR)%>%summarise(Y=mean(Y))

k_qr = k
#plot(xymed)
k_y_qr=(k_qr)*price
tail(k_y_qr)            ## Strike: historical yield in US/ha

######################                 Indemnities
#km=as.data.frame(k)%>% mutate(across(c('PR', 'Y'), round, 4))
#km=QMED%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))
km1=QMED1%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))
km2=QMED2%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))


ES_A1<-ES_A%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))

ES_B1<-ES_B%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y = mean(Y))

ES_D1<-ES_D%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))

ES_E1<-ES_E%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))

###### ZONA A   (Upper Limit)
ZonaA=merge(km1,ES_A1, by="PR")


Indemnity_pr_lower_limit_qr= ifelse(abs(ZonaA$Y.x- ZonaA$Y.y) > 0,
                                     (ZonaA$Y.x - ZonaA$Y.y ), 0)
Premium_pr_lower_limit_qr=mean(Indemnity_pr_lower_limit_qr)
Premium_pr_lower_limit_qr
Premium_pr_lower_limit_qr*price   

########---------------################3
indemnities_A=cbind(ZonaA, Indemnity_pr_lower_limit_qr)
colnames(indemnities_A) <- c("PR", "k", "Miu*", "I_A")
head(indemnities_A)
########----------------################

#

###### ZONA B:   (Lower Threshold)
ZonaB=merge(km1,ES_B1, by="PR")

Indemnity_pr_lower_thre_qr= ifelse(abs(ZonaB$Y.x - ZonaB$Y.y) > 0,
                                   abs(ZonaB$Y.x - ZonaB$Y.y) , 0)
Premium_pr_lower_thre_qr=mean(Indemnity_pr_lower_thre_qr)
Premium_pr_lower_thre_qr
Premium_pr_lower_thre_qr*price            # 

########----------------################
indemnities_B=cbind(ZonaB, Indemnity_pr_lower_thre_qr)
colnames(indemnities_B) <- c("PR", "k", "Miu*", "I_B")
head(indemnities_B)
########----------------################

##
###### ZONA D:(Upper Threshold)
ZonaD=merge(km2,ES_D1, by="PR")

Indemnity_pr_upper_thre_qr= ifelse(abs(ZonaD$Y.x - ZonaD$Y.y) > 0,
                                   abs(ZonaD$Y.x - ZonaD$Y.y) , 0)
Premium_pr_upper_thre_qr=mean(Indemnity_pr_upper_thre_qr)
Premium_pr_upper_thre_qr
Premium_pr_upper_thre_qr*price            # 

indemnities_D=cbind(ZonaD, Indemnity_pr_upper_thre_qr)
colnames(indemnities_D) <- c("PR", "k", "Miu*", "I_D")
head(indemnities_D)
str(Indemnity_pr_upper_thre_qr)
###
###### ZONA E:(Upper limit)
ZonaE=merge(km2,ES_E1, by="PR")

Indemnity_pr_upper_lim_qr= ifelse(abs(ZonaE$Y.x - ZonaE$Y.y) > 0,
                                  abs(ZonaE$Y.x - ZonaE$Y.y) , 0)

Premium_pr_upper_limit_qr=mean(Indemnity_pr_upper_lim_qr)
Premium_pr_upper_limit_qr
Premium_pr_upper_limit_qr*price            # 

indemnities_E=cbind(ZonaE, Indemnity_pr_upper_lim_qr)
colnames(indemnities_E) <- c("PR", "k", "Miu*", "I_E")
head(indemnities_E)


##PREMIUMS LOWER: (ZONA A +ZONA B)
Premium_qr_all_lower= Premium_pr_lower_limit_qr + Premium_pr_lower_thre_qr
Premium_qr_all_lower*price    #532.78 usd/HA     vs gumbel 90: 349.4 USD/ha  vs GR): 323.94
Premium_qr_all_lower

##PREMIUMS UPPER: (ZONA D +ZONA E)
Premium_qr_all_upper= Premium_pr_upper_thre_qr +  Premium_pr_upper_limit_qr
Premium_qr_all_upper*price         #567.06 usd/ha          vs gumbel 90: 412.9063/ha vs: 384.76 (RG)
Premium_qr_all_upper

#Fair_insurance_premium=sum(Premium_spei_upper_thre_qr, Premium_spei_upper_limit_qr)*price
Fair_insurance_premium1 = append(Indemnity_pr_upper_thre_qr, Indemnity_pr_upper_lim_qr)
Fair_insurance_premium2 = append(Indemnity_pr_lower_limit_qr, Indemnity_pr_lower_thre_qr)
Fair_insurance_premium3 = append(Fair_insurance_premium1, Fair_insurance_premium2) 
mean(Fair_insurance_premium3)
mean(Fair_insurance_premium3)*price/max_indemnity           #Tasa IR correcta, con todas las zonas

str(Indemnity_pr_lower_limit_qr)
str(Indemnity_pr_lower_thre_qr)
str(Indemnity_pr_upper_thre_qr)
str(Indemnity_pr_upper_lim_qr)

#----promediando todas las zonas
#mean(Indemnity_spei_upper_lim_qr)*price

###############Uninsured representa la zona de no disbuursals:
uninsured = mean_values
uninsuredAB = lower%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))
uninsuredDE = upper%>% mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))



## ----------------------------Insured Yields:

#mean_values_redon = mean_values%>%mutate(across(c('PR', 'Y'), round, 5))
mean_values_redon1 = QMED1%>%mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))
mean_values_redon2 = QMED2%>%mutate(across(c('PR', 'Y'), round, 5))%>%group_by(PR)%>%summarise(Y=mean(Y))



#Indemnities organizados:
#Ind_D=cbind(Indemnity_pr_upper_thre_qr,ES_D1)
Indemnity_pr_upper_lim_qr
##
head(indemnities_D)
head(indemnities_E)
indemnities_DE=merge(indemnities_D, indemnities_E, by="PR")[,-c(2,3,5,6)]
str(indemnities_DE)
head(indemnities_DE)
base_indenmities = merge(indemnities_DE, mean_values_redon2, by="PR")
head(mean_values_redon1)
      ############Datos con yield insured y no insured (que es la zona media o C)

#####----------------###############3
indemnities_AB=merge(indemnities_A, indemnities_B, by="PR")[,-c(2,3,5,6)]
str(indemnities_AB)
head(indemnities_AB)

base_indenmities2 = merge(indemnities_AB, mean_values_redon1, by="PR")
head(base_indenmities2)


##########-----------#################3

## ----------------------------Insured & Uninsured premiums:
Premium_pr_lower_limit_qr
Premium_pr_lower_thre_qr
Premium_pr_upper_thre_qr*price
Premium_pr_upper_limit_qr*price

#yield_insured:
Y_Insured_DE= base_indenmities%>% mutate(Y_ins = Y+I_D+I_E-Premium_qr_all_upper)
head(Y_Insured_DE)
str(Y_Insured_DE)


#yield_insured AB:
Y_Insured_AB= base_indenmities2%>% mutate(Y_ins = Y+I_A+I_B-Premium_qr_all_lower)
head(Y_Insured_AB)


##Uninsured Average
Y_uninsured_gorroAB = mean(uninsuredAB$Y)
Y_uninsured_gorroDE = mean(uninsuredDE$Y)
Y_uninsured_gorro = mean(uninsured$Y)

##  -------------------------------------------------            Uninsured SV AB:
##Operacion diferencia cuadr?tica Yield_uninsured menos Y_uninsured_gorro
SV_Uninsured_0_AB=ifelse(uninsuredAB$Y - Y_uninsured_gorro!=0,1,0)
SV_Uninsured_0_AB
SV_uninsured_00_AB=cbind(uninsuredAB, SV_Uninsured_0_AB)
tail(SV_uninsured_00_AB)
SV_uninsured_1_AB=subset(SV_uninsured_00_AB, SV_Uninsured_0_AB == 1)
head(SV_uninsured_1_AB) 

SV_unsured2_AB= SV_uninsured_1_AB %>% mutate(SV_un_dif = (Y- Y_uninsured_gorro)^2)
SV_uninsuredAB=mean(SV_unsured2_AB$SV_un_dif)
SV_uninsuredAB                                  ######SV uninsured
SV_uninsuredAB*price
###############################################

##  -------------------------------------------------            Uninsured SV:
##Operacion diferencia cuadr?tica Yield_uninsured menos Y_uninsured_gorro
SV_Uninsured_0_DE=ifelse(uninsuredDE$Y - Y_uninsured_gorro!=0,1,0)
SV_Uninsured_0_DE
SV_uninsured_00_DE=cbind(uninsuredDE, SV_Uninsured_0_DE)
tail(SV_uninsured_00_DE)
SV_uninsured_1_DE=subset(SV_uninsured_00_DE, SV_Uninsured_0_DE == 1)
head(SV_uninsured_1_DE) 

SV_unsured2_DE= SV_uninsured_1_DE %>% mutate(SV_un_dif = (Y- Y_uninsured_gorro)^2)
SV_uninsured_DE=mean(SV_unsured2_DE$SV_un_dif)
SV_uninsured_DE                                ######SV uninsured
SV_uninsured_DE*price
###############################################


##   -----------------------------------------------          insured SV Upper
##Operacion diferencia cuadr?tica Yield_insured menos Y_uninsured_gorro
SV_insured_0_DE=ifelse(Y_Insured_DE$Y_ins - Y_uninsured_gorroDE !=0,1,0)
SV_insured_0_DE
SV_insured_00_DE=cbind(Y_Insured_DE, SV_insured_0_DE)
tail(SV_insured_00_DE)
SV_insured_1_DE=subset(SV_insured_00_DE, SV_insured_0_DE == 1)
head(SV_insured_1_DE) 

SV_insured2_DE= SV_insured_1_DE %>% mutate(SV_ins_dif = (Y_ins- Y_uninsured_gorroDE)^2)
SV_insured_DE=mean(SV_insured2_DE$SV_ins_dif)
SV_insured_DE                                  ######SV uninsured upper
SV_insured_DE*price
## --------------------------------------------------       insured SV Lower
##Operacion diferencia cuadr?tica Yield_insured menos Y_uninsured_gorro
SV_insured_0_AB=ifelse(Y_Insured_AB$Y_ins - Y_uninsured_gorroAB !=0,1,0)
SV_insured_0_AB
SV_insured_00_AB=cbind(Y_Insured_AB, SV_insured_0_AB)
tail(SV_insured_00_AB)
SV_insured_1_AB=subset(SV_insured_00_AB, SV_insured_0_AB == 1)
head(SV_insured_1_AB) 

SV_insured2_AB= SV_insured_1_AB %>% mutate(SV_ins_dif = (Y_ins- Y_uninsured_gorroAB)^2)
SV_insured_AB=mean(SV_insured2_AB$SV_ins_dif)
SV_insured_AB                                  ######SV uninsured upper
SV_insured_AB*price

#--------------
#Hedging Effectiveness:
HE_qr_DE = (SV_insured_DE- SV_uninsured_DE)/SV_uninsured_DE
#HE_qr_DE = 1- (SV_insured_DE/SV_uninsured) ###es esta opcion correcta?######################33
HE_qr_DE*100


HE_qr_AB = (SV_insured_AB- SV_uninsuredAB)/SV_uninsuredAB
#HE_qr = 1- (SV_insured/SV_uninsured) ###es esta opcion correcta?######################33
HE_qr_AB*100

