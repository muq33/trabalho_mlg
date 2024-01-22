if(!require(statmod)){install.packages("statmod")}
if(!require(fitdistrplus)){install.packages("fitdistrplus")}
if(!require(doBy)){install.packages("doBy")}
if(!require(Rcmdr)){install.packages("Rcmdr")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(lmtest)){install.packages("lmtest")}
if(!require(xtable)){install.packages("xtable")}
if(!require(car)){install.packages("car")}

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

dados <- read.table("dados.txt", header = T, sep = ",", dec = ".")
dados$grupo <- factor(dados$grupo)

#Descritivo

plot(fitdist(dados$tempo_duracao, distr = "gamma"))
plot(fitdist(dados$tempo_duracao, distr = "exp"))
summaryBy(tempo_duracao~grupo, data=dados, FUN=c(mean, sd))
plotMeans(dados$tempo_duracao,  dados$grupo,  error.bars="se")
with(dados, plot(densidade_máxima, tempo_duracao))

ggplot(dados, aes(x=densidade_máxima, y=tempo_duracao, color=grupo)) + 
  geom_point(size=6) + ggtitle("Gráfico de dispersão dos dados")

#Modelo

ajust <- vector("list", 3)
formulas <- c(tempo_duracao~densidade_máxima, tempo_duracao~densidade_máxima + grupo, tempo_duracao~densidade_máxima*grupo)

ajust[[1]] <- glm(formulas[[1]], data = dados, family = Gamma(link = "identity"))
ajust[[2]] <- glm(formulas[[2]], data = dados, family = Gamma(link = "identity"))
ajust[[3]] <- glm(formulas[[3]], data = dados, family = Gamma(link = "identity"))

#Testes

testes <- vector("list", 2)
testes[[1]] <- lrtest(ajust[[1]],ajust[[2]])
testes[[2]] <- lrtest(ajust[[2]],ajust[[3]])

formulas_str <- as.character(formulas)
saida_teste <- data.frame(
  formula = c(paste0(formulas_str[[1]], " : ", formulas_str[[2]]),
              paste0(formulas_str[[2]], " : ", formulas_str[[3]])),
  Df = c(testes[[1]]$Df[[2]],
         testes[[2]]$Df[[2]]),
  Chisq = c(testes[[1]]$Chisq[[2]],
            testes[[2]]$Chisq[[2]]),
  p = c(testes[[1]]$`Pr(>Chisq)`[[2]],
        testes[[2]]$`Pr(>Chisq)`[[2]])
); xtable(saida_teste)


linearHypothesis(ajust[[3]],hypothesis.matrix=c(0,0,0,0,0,1,1,1))

  
#Diagnóstico
fit.model <- ajust[[3]]
source("envelope_gama_ident.R")
