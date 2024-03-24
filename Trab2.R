# Leitura dos dados -------------------------------------------------------
dados <- read.table("http://sobrevida.fiocruz.br/dados/ctinca.dat", header = TRUE)
head(dados)


# Pacotes -----------------------------------------------------------------
require(survival)
require(muhaz)
require(pander)
require(flexsurv)
require(hnp)
require(ggplot2)
require(survminer)
require(parmsurvfit)
require(pander)
require(hnp)
library(MASS)

#Funções -----------------------------------------------------------------
# Função para Teste da Razão de Verossimilhanças (TRV)
TRV <- function(m0, m1) {
  # m0: Modelo nulo
  # m1: Modelo alternativo (saturado)
  TRV = 2 * (logLik(m1) - logLik(m0))
  gl = length(coef(m1)) - length(coef(m0))
  p = 1 - pchisq(TRV, gl)
  res = cbind(TRV, gl, p)
  return(res)
}



# Traducao dos dados ------------------------------------------------------
dados$sexo <- ifelse(dados$sexo == "Male","Mas","Fem")
dados$desnut <- ifelse(dados$desnut == "n", "Não", "Sim")
dados$comorbi <- ifelse(dados$comorbi == "n", "Não", "Sim")
dados$leucopenia <- ifelse(dados$leucopenia == "n", "Não", "Sim")

dados$sexo <- factor(dados$sexo, levels = c("Mas","Fem"))
dados$desnut <- factor(dados$desnut, levels = c("Sim","Não"))
dados$comorbi <- factor(dados$comorbi, levels = c("Sim","Não"))
dados$leucopenia <- factor(dados$leucopenia, levels = c("Sim","Não"))
dados$gptumor <- factor(dados$gptumor, levels = c("Loco", "Mtx", "Hemato"))

summary(dados)

str(dados)

# Descritiva --------------------------------------------------------------
km1 <- survfit(Surv(tempo, status)~gptumor, dados)
plot(km1, xlab='Tempo', ylab='Prob de sobrevivencia', lwd=2, col=1:3)
legend('bottomleft', title='Tumores', c("Loco", "Mtx", "Hemato"), col=1:3, lwd=2, bty='n')

km2 <- survfit(Surv(tempo, status)~sexo, dados)
plot(km2, xlab='Tempo', ylab='Prob de sobrevivencia', lwd=2, col=1:2)
legend('bottomleft', title='Sexo', c("Mas","Fem"), col=1:2, lwd=2, bty='n')

km3 <- survfit(Surv(tempo, status)~desnut, dados)
plot(km3, xlab='Tempo', ylab='Prob de sobrevivencia', lwd=2, col=1:3)
legend('bottomleft', title='Desnutricao', c("Sim","Não"), col=1:3, lwd=2, bty='n')

km4 <- survfit(Surv(tempo, status)~comorbi, dados)
plot(km4, xlab='Tempo', ylab='Prob de sobrevivencia', lwd=2, col=1:3)
legend('bottomleft', title='Comorbidade', c("Sim","Não"), col=1:3, lwd=2, bty='n')

km5 <- survfit(Surv(tempo, status)~leucopenia, dados)
plot(km5, xlab='Tempo', ylab='Prob de sobrevivencia', lwd=2, col=1:3)
legend('bottomleft', title='Leucopenia', c("Sim","Não"), col=1:3, lwd=2, bty='n')

km1; km2; km3; km4; km5

# Regressão ---------------------------------------------------------------
fit1 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='exponential')

fit2 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='weibull')

fit3 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='gamma')

fit4 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='lnorm')

fit5 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='llogis')

fit6 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='gompertz')

fit7 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='gengamma')

fit8 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='genf')

fit9 <- flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi+leucopenia, data=dados, dist='exp')

fit1$AIC; fit2$AIC; fit3$AIC; fit4$AIC; fit5$AIC; fit6$AIC; fit7$AIC; fit8$AIC; fit9$AIC

# Retirar sexo?


aic_bic_valores <- data.frame(
  Modelo = c("Exponential", "Weibull", "Gamma", "Log-Normal", "Log-Logistic", "Gompertz", "Generalized Gamma", "Generalized F", "Exponential 2"),
  AIC = c(AIC(fit1), AIC(fit2), AIC(fit3), AIC(fit4), AIC(fit5), AIC(fit6), AIC(fit7), AIC(fit8), AIC(fit9)),
  BIC = c(BIC(fit1), BIC(fit2), BIC(fit3), BIC(fit4), BIC(fit5), BIC(fit6), BIC(fit7), BIC(fit8), BIC(fit9))
)

pander(aic_bic_valores)




s1 <- c() ; s2 <- c() 

for(i in 1:nrow(dados)) {
  s2[i] <- summary(fit4, type='survival', t=dados$tempo[i], newdata=dados[i,], tidy=T)$est
  if(i%%10==0) print(i)
}



a2 <- runif(length(s2), min=0, max=s2)
r2 <- ifelse(dados$status==1, qnorm(s2), qnorm(a2))


qqnorm(r2, main='Weibull')     ; abline(0,1)



hnp(r2, main='Weibull', half=F,print.on = T)    

########################

for(i in 1:nrow(dados)) {
  s2[i] <- summary(fit6, type='survival', t=dados$tempo[i], newdata=dados[i,], tidy=T)$est
  if(i%%10==0) print(i)
}



a2 <- runif(length(s2), min=0, max=s2)
r2 <- ifelse(dados$status==1, qnorm(s2), qnorm(a2))


qqnorm(r2, main='Weibull')     ; abline(0,1)



hnp(r2, main='Weibull', half=F,print.on = T)  
 



### teste variaveis
fit_nulo = flexsurvreg(Surv(tempo, status)~ 1, data=dados, dist='lnorm')

fit_SEM_sexo = flexsurvreg(Surv(tempo, status)~gptumor+desnut+comorbi+leucopenia, data=dados, dist='lnorm')
fit_SEM_gptumor = flexsurvreg(Surv(tempo, status)~gptumor+sexo+comorbi+leucopenia, data=dados, dist='lnorm')
fit_SEM_comorbi = flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+leucopenia, data=dados, dist='lnorm')
fit_SEM_desnut = flexsurvreg(Surv(tempo, status)~gptumor+sexo+comorbi+leucopenia, data=dados, dist='lnorm')
fit_SEM_leucopenia = flexsurvreg(Surv(tempo, status)~gptumor+sexo+desnut+comorbi, data=dados, dist='lnorm')

TRV(fit_SEM_sexo, fit4) 
# tirar sexo
# TRV gl         p
# [1,] 1.964037  1 0.1610822

TRV(fit_SEM_gptumor, fit4) # gptumor IMPORTANTE

TRV(fit_SEM_comorbi, fit4) # comorbi IMPORTANTE

TRV(fit_SEM_desnut, fit4) # desnut IMPORTANTE

TRV(fit_SEM_leucopenia, fit4) # leucopenia IMPORTANTE


####fit_SEM_sexo é o melhor
fit_SEM_gptumor.2 = flexsurvreg(Surv(tempo, status)~desnut+comorbi+leucopenia, data=dados, dist='lnorm')
fit_SEM_comorbi.2 = flexsurvreg(Surv(tempo, status)~gptumor+desnut+leucopenia, data=dados, dist='lnorm')
fit_SEM_desnut.2 = flexsurvreg(Surv(tempo, status)~gptumor+comorbi+leucopenia, data=dados, dist='lnorm')
fit_SEM_leucopenia.2 = flexsurvreg(Surv(tempo, status)~gptumor+desnut+comorbi, data=dados, dist='lnorm')

        
TRV(fit_SEM_gptumor.2, fit_SEM_sexo)
TRV(fit_SEM_comorbi.2, fit_SEM_sexo)
TRV(fit_SEM_desnut.2, fit_SEM_sexo)
TRV(fit_SEM_leucopenia.2, fit_SEM_sexo)



##### fit_SEM_sexo é o melhor
fit_APENAS_gptumor = flexsurvreg(Surv(tempo, status)~gptumor, data=dados, dist='lnorm')
fit_APENAS_comorbi = flexsurvreg(Surv(tempo, status)~comorbi, data=dados, dist='lnorm')
fit_APENAS_desnut = flexsurvreg(Surv(tempo, status)~desnut, data=dados, dist='lnorm')
fit_APENAS_leucopenia = flexsurvreg(Surv(tempo, status)~leucopenia, data=dados, dist='lnorm')


TRV(fit_APENAS_gptumor, fit_SEM_sexo)
TRV(fit_APENAS_comorbi, fit_SEM_sexo)
TRV(fit_APENAS_desnut, fit_SEM_sexo)
TRV(fit_APENAS_leucopenia, fit_SEM_sexo)


##### fit_COM_idade é o melhor, sem sexo

fit_SEM_sexo = flexsurvreg(Surv(tempo, status)~gptumor+desnut+comorbi+leucopenia, data=dados, dist='lnorm')
fit_COM_idade = flexsurvreg(Surv(tempo, status)~gptumor+desnut+comorbi+leucopenia+idade, data=dados, dist='lnorm')

TRV(fit_SEM_sexo, fit_COM_idade)



AIC(fit4)
AIC(fit_COM_idade)

