dados <- read.table("http://sobrevida.fiocruz.br/dados/ctinca.dat", header = TRUE)
head(dados)

# Pacotes
require(survival)
require(muhaz)
require(pander)
require(flexsurv)
require(hnp)
require(ggplot2)
require(survminer)
require(parmsurvfit)


teste 2