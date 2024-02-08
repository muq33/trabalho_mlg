influence_graph <- function(model, plot = F){
  yajust <- fitted.values(model)
  X <- model.matrix(model)
  w <- model$weights
  W <- diag(w)
  H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W) #Matriz de projeção
  h <- diag(H) #Elementos da matriz de projeção H
  
  if(plot == T){
    plot(h, xlab = "Indices", ylab = "Elementos da Matriz de Projeção H",
         pch = 19, cex = 1.3)
  }
  return(h)
}

res_ajust <- function(model){
  residuos <- qresid(ajust[[5]])
  ajustados <- predict(ajust[[5]])
  return(list(
    resid = residuos,
    ajust = ajustados
  ))
}