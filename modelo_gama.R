if(!require(statmod)){install.packages("statmod")}
if(!require(fitdistrplus)){install.packages("fitdistrplus")}
if(!require(doBy)){install.packages("doBy")}
if(!require(Rcmdr)){install.packages("Rcmdr")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(lmtest)){install.packages("lmtest")}
if(!require(xtable)){install.packages("xtable")}
if(!require(car)){install.packages("car")}

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("utils.R")

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
formulas <- c(tempo_duracao~1,tempo_duracao~densidade_máxima, tempo_duracao~densidade_máxima + grupo, tempo_duracao~densidade_máxima*grupo)


ajust[[1]] <- glm(formulas[[1]], data = dados, family = Gamma(link = "log"))
ajust[[2]] <- glm(formulas[[2]], data = dados, family = Gamma(link = "log"))
ajust[[3]] <- glm(formulas[[3]], data = dados, family = Gamma(link = "log"))
ajust[[4]] <- glm(formulas[[4]], data = dados, family = Gamma(link = "log"))

# 0 é referencia nas 2

#Testes

testes <- vector("list", 2)
testes[[1]] <- lrtest(ajust[[1]],ajust[[2]])
testes[[2]] <- lrtest(ajust[[2]],ajust[[3]])
testes[[3]] <- lrtest(ajust[[3]],ajust[[4]])

formulas_str <- c("Densidade", "Grupo | Densidade", "Grupo:Densidade | Grupo, Densidade")
saida_teste <- data.frame(
  formula = formulas_str,
  Df = c(testes[[1]]$Df[[2]],
         testes[[2]]$Df[[2]],
         testes[[3]]$Df[[2]]),
  Chisq = c(testes[[1]]$Chisq[[2]],
            testes[[2]]$Chisq[[2]],
            testes[[3]]$Chisq[[2]]),
  p = c(testes[[1]]$`Pr(>Chisq)`[[2]],
        testes[[2]]$`Pr(>Chisq)`[[2]],
        testes[[3]]$`Pr(>Chisq)`[[2]])
); xtable(saida_teste)


#Comparação de modelos
ajust[[5]] <- glm(formulas[[4]], data = dados, family = Gamma(link = "identity"))
saida_comparacao <- data.frame(
  modelo = c("Gamma", "Gamma"),
  ligacao = c("Logaritmo", "Identidade"),
  AIC = c(AIC(ajust[[4]]), AIC(ajust[[5]])),
  BIC = c(BIC(ajust[[4]]), BIC(ajust[[5]])),
  TRV = c(lrtest(ajust[[4]], ajust[[5]])$`Pr(>Chisq)`)
) ; xtable(saida_comparacao)

#Contrastes (precisamos de hipóteses melhores)


#Diagnóstico - Gamma (Ident)
fit.model <- ajust[[5]]
source("envelope_gama_ident.R")


#Gráfico influencia
par(mfrow = c(1,2))
grafico_inf <- influence_graph(ajust[[4]], plot =T)
grafico_inf <- influence_graph(ajust[[5]], plot =T)

which(grafico_inf > 0.02)

#Distancia Cook
influenceIndexPlot(ajust[[5]], vars=c("Cook"), pch=19, cex=1.2, main="Distância de Cook") #Modelo ident
#Residuos estudentizados
influenceIndexPlot(ajust[[5]], vars=c("Studentized"), pch=19, cex=1.2, main = "Resíduos Estudentizados")#Modelo ident

#Residuos quantilicos
residuos <- qresiduals(ajust[[5]])
plot(residuos, ylab = "Resíduos Quantilicos", pch=19, cex=1.2)
qqnorm(residuos, pch=19, cex=1.2)
qqline(residuos, col = 2)
shapiro.test(residuos)


#Residuos vs valores ajustados
resid_pred_ident <- res_ajust(ajust[[5]])
plot(resid_pred_ident$resid ~ resid_pred_ident$ajust, pch = 20, cex = 1.4, col = 'blue', main = "Resíduos vs Valores Ajustados")

#Quantil quantil normal
qqnorm(resid_pred_ident$resid, pch = 20, cex = 1.4, col = 'blue')
qqline(resid_pred_ident$resid)
