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