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
plot(fitdist(dados$tempo_duracao, distr = "invgauss", start = list(mean = 1, shape = 1)))
 
summaryBy(tempo_duracao~grupo, data=dados, FUN=c(mean, sd))
plotMeans(dados$tempo_duracao,  dados$grupo,  error.bars="se")
with(dados, plot(densidade_máxima, tempo_duracao))

ggplot(dados, aes(x=densidade_máxima, y=tempo_duracao, color=grupo)) + 
  geom_point(size=6) + ggtitle("Gráfico de dispersão dos dados")

#Modelo - Gamma

ajust <- vector("list", 3)
formulas <- c(tempo_duracao~densidade_máxima, tempo_duracao~densidade_máxima + grupo, tempo_duracao~densidade_máxima*grupo)

ajust[[1]] <- glm(formulas[[1]], data = dados, family = Gamma(link = "identity"))
ajust[[2]] <- glm(formulas[[2]], data = dados, family = Gamma(link = "identity"))
ajust[[3]] <- glm(formulas[[3]], data = dados, family = Gamma(link = "identity"))

#Testes - Gamma

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


linearHypothesis(ajust[[3]],hypothesis.matrix=c(0,0,0,0,1,1))
linearHypothesis(ajust[[3]],hypothesis.matrix=c(0,0,0,0,1,1))
  
#Diagnóstico - Gamma
fit.model <- ajust[[3]]
source("envelope_gama_ident.R")

#--grafico influencia
yajust <- fitted.values(ajust[[3]])
X <- model.matrix(ajust[[3]])
w <- ajust[[3]]$weights
W <- diag(w)
H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W) #Matriz de projeção
h <- diag(H) #Elementos da matriz de projeção H


plot(h, xlab = "Indices", ylab = "Elementos da Matriz de Projeção H",
     pch = 19, cex = 1.3)
title("Pontos de Alavanca (Leverage)")
which(h > 0.02)

par(mfrow = c(1,2))
influenceIndexPlot(ajust[[3]], vars=c("Cook"), pch=19, cex=1.2, main="Distância de Cook")
influenceIndexPlot(ajust[[3]], vars=c("Studentized"), pch=19, cex=1.2, main = "Resíduos Padronizados")

#--residuos quantilicos
par(mfrow=c(2,1))
residuos <- qresiduals(ajust[[3]])
plot(residuos, ylab = "Resíduos Quantilicos", pch=19, cex=1.2)
qqnorm(residuos, pch=19, cex=1.2)
qqline(residuos, col = 2)
shapiro.test(residuos)
par(mfrow=c(1,1))
shapiro.test(residuos)

#--residuos vs valores ajustados
par(mfrow = c(1,2))
residuos <- qresid(ajust[[3]])
ajustados <- predict(ajust[[3]])
plot(residuos ~ ajustados, pch = 20, cex = 1.4, col = 'blue', main = "Resíduos vs Valores Ajustados")


#--quantil quantil normal
par(mfrow = c(1,1))
qqnorm(residuos, pch = 20, cex = 1.4, col = 'blue')
qqline(residuos)

#--componentes da deviance
phi = 1
res = residuals(ajust[[3]], type = 'response') # Resíduos Ordinários

rd = residuals(ajust[[3]], type = "deviance") # Componente da Deviance
rdl = rd/sqrt(phi*(1-h)) # Componente da Dev. padronizado

rp = residuals(ajust[[3]], type = "pearson") # Residuo de Pearson
rpl = rp/sqrt(phi*(1-h)) # Residuo Pearson Padronizado

LD = (h/(1-h))*(rpl)^2 # Afastamento da Vero. para Dist. de Cook

# Componentes da Deviance
par(mfrow = c(1,2))
plot(rdl, xlab = "Indices", ylab = "Componentes da Deviance",
     pch = 1, cex = 1)
title("Residuo - Componente da Deviance Padronizada")

# Distancia de Cook
plot(LD, xlab = "Indices", ylab = "Distância de Cook",
     pch = 19, cex = 1.3)
title("Distância de Cook vs Indices")
par(mfrow = c(1,1))

#Modelo - Normal Inversa

ajust_ninv <- vector("list", 3)

ajust_ninv[[1]] <- glm(formulas[[1]], data = dados, family = inverse.gaussian(link = "identity"))
ajust_ninv[[2]] <- glm(formulas[[2]], data = dados, family = inverse.gaussian(link = "identity"))
ajust_ninv[[3]] <- glm(formulas[[3]], data = dados, family = inverse.gaussian(link = "identity"))

#Testes - Normal Inversa

testes_ninv <- vector("list", 2)
testes_ninv[[1]] <- lrtest(ajust_ninv[[1]],ajust_ninv[[2]])
testes_ninv[[2]] <- lrtest(ajust_ninv[[2]],ajust_ninv[[3]])

formulas_str <- as.character(formulas)
saida_teste_ninv <- data.frame(
  formula = c(paste0(formulas_str[[1]], " : ", formulas_str[[2]]),
              paste0(formulas_str[[2]], " : ", formulas_str[[3]])),
  Df = c(testes_ninv[[1]]$Df[[2]],
         testes_ninv[[2]]$Df[[2]]),
  Chisq = c(testes_ninv[[1]]$Chisq[[2]],
            testes_ninv[[2]]$Chisq[[2]]),
  p = c(testes_ninv[[1]]$`Pr(>Chisq)`[[2]],
        testes_ninv[[2]]$`Pr(>Chisq)`[[2]])
); xtable(saida_teste_ninv)


linearHypothesis(ajust_ninv[[3]],hypothesis.matrix=c(0,0,0,0,1,1))
linearHypothesis(ajust_ninv[[3]],hypothesis.matrix=c(0,0,0,0,1,1))

#Diagnóstico - Normal Inversa
fit.model <- ajust_ninv[[3]]
source("envelope_ginv_ident.R")


