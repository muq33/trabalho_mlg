#------------------------------------------------------------#
# Para rodar este programa deixe no objeto fit.model a saída      
# do ajuste da regressão com erros normal inversa e ligação
# log.  Deixe também os dados disponíveis através do
# comando attach(...).Depois use o comando source(...) no R ou
# S-Plus para executar o programa. A sequência de comandos é 
# a seguinte:               
#
#     fit.model <- ajuste
#     attach(dados)
#     source("envel_ninv")                                      
#
# A saída será o gráfico de envelope para o resíduo componente 
# do desvio padronizado. Para colocar um título no gráfico use
# o comando title("..."). Para  usar  outras  ligações
# mudar no programa abaixo o termo family=inverse.gaussian
# (link=log)  por inverse.gaussian no caso de ligação
# recíproca ao quadrado  ou por family=inverse.gaussian
# (link=identidade) no caso de ligação identidade.
#------------------------------------------------------------#
# Função para gerar observações de uma distribuição
# normal inversa

rig <- function(n, mu = stop("no shape arg"), lambda = 1)
{
  #  Random variates from inverse Gaussian distribution
  #  Reference:
  #      Chhikara and Folks, The Inverse Gaussian Distribution,
  #      Marcel Dekker, 1989, page 53.
  #  GKS  15 Jan 98
  #
  if(any(mu<=0)) stop("mu must be positive")
  if(any(lambda<=0)) stop("lambda must be positive")
  if(length(n)>1) n <- length(n)
  if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
  if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
  y2 <- rchisq(n,1)
  u <- runif(n)
  r1 <- mu/(2*lambda) * (2*lambda + mu*y2 - sqrt(4*lambda*mu*y2 + mu^2*y2^2))
  r2 <- mu^2/r1
  ifelse(u < mu/(mu+r1), r1, r2)
}
#------------------------------------------------------------#                     
#
par(mfrow=c(1,1))
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
w <- fit.model$weights
mu <-predict(fit.model,type="response")
W <- diag(w)
H <- solve(t(X)%*%W%*%X)
H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
h <- diag(H)
ro <- resid(fit.model,type="response")
fi <- (n-p)/sum((ro^2)/(fitted(fit.model)^3))
td <- resid(fit.model,type="deviance")*sqrt(fi/(1-h))
#
e <- matrix(0,n,100)
#
for(i in 1:100){
  resp <- rig(n,mu,fi)
  fit <- glm(resp ~ X, family=inverse.gaussian(link=identity))
  w <- fit$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  ro <- resid(fit,type="response")
  phi <- (n-p)/sum((ro^2)/(fitted(fit)^3))
  e[,i] <- sort(resid(fit,type="deviance")*sqrt(phi/(1-h)))}
#
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:n){
  eo <- sort(e[i,])
  e1[i] <- (eo[2]+eo[3])/2
  e2[i] <- (eo[97]+eo[98])/2}
#
med <- apply(e,1,mean)
faixa <- range(td,e1,e2)
par(pty="s")
qqnorm(td,xlab="Percentil da N(0,1)",
       ylab="Componente do Desvio", ylim=faixa, pch=16, main="")
par(new=T)
#
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
par(new=T)
qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2, main="")
#------------------------------------------------------------#  

