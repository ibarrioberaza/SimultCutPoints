# Functions for the optimal location of the cut-off points based on GAM ---

# CortesGAM ---------------------------------------------------------------
  # Function that calculates the associated error for each cut-off point
  # (comparison between GAM and piece-wise functions)

  # Arguments:
  # x: data of the variable to be categorized.
  # y: GAM model prediction with terms
  # w: se.fit^(-2) (inverse of variance)
  # pred: prediction of the piece-wise function.
  # ymed: mean of "y" at each level of the factor.

  CortesGAM <- function(x, y, w = rep(1,length(y)), cortes, nmin = 0.02*length(x)){
    cortes = unique(sort(c(min(x), cortes, max(x))))
    pred = numeric(length(y))
    
    x2 = cut(x, include.lowest = TRUE, breaks = cortes)
    niveles = levels(x2)
    nn = length(niveles)
    minimo = min(table(x2)) # Minimum frequency with the new cut-off points
    
    error = Inf
    if (minimo >= nmin) {
      ymed = tapply(y, x2, mean); 
      for (j in 1:nn) pred[x2 == niveles[j]] = ymed[j]
      error = sum(w*(y - pred)**2)}
    
    list(x = x, y = y, w = w, cortes = cortes, pred = pred, err = error)
}

# Opt.cpoints.GAM ---------------------------------------------------------
  # optimal cut-off points with the BackAddFor search. 
  # The model is fitted outside the function.
  # nc: Number of cut-off points to search for.
  
  Opt.cpoints.GAM <- function(x, y, w, nfino=100, nc, eps = 0.001, repmax=5 ){
    cfino = seq(min(x), max(x), length = nfino) # Grid to search
    mat.aux <- NULL
    
    # Initial nodes
    c0 = seq(quantile(x, 0.025), quantile(x, 0.975), length = nc + 2)[-c(1, nc + 2)] # Removing the first and the last ones.
    err0 = CortesGAM(x = x , y = y , w = w , c = c0)$err
    
    # Sequential search
    cfino = seq(min(x), max(x), length = nfino + 2)[-c(1, nfino + 2)]
    
    seguir = TRUE
    nrep = 0
    err_old = err0
    
    while (seguir) {
      nrep = nrep + 1 
      err_old = err0
      for (j in 1:nc) { 
        for (ifino in 1:nfino) {
          c1 = c0
          c1[j] = cfino[ifino]
          err1 = CortesGAM(x = x , y = y , w = w , c = c1)$err
          mat.aux = rbind(mat.aux, c(c1, err1))
          if (err1 < err0) {err0 = err1; c0 = c1}
        }
      }
      if (nrep > 1 & abs((err0 - err_old) / err_old) < eps) {seguir = FALSE}
      if (nrep > repmax) {seguir = FALSE} 
      err_old = err0  
    }
    
    list(res = CortesGAM(x, y, w, c0), errores = mat.aux)
  }

# Pseudo_BIC function -----------------------------------------------------
  # Created from the logLik.gam{mgcv}
  Pseudo_BIC2 <- function(object, nc){ 
    b <- BIC(object) + nc*log(length(object$y))
    b
  }

# SimLOG.multi ------------------------------------------------------------
  # Function to simulate the theoretical logistic model                   
  # Double multivariate search
  SimLOG.multi <- function(cortes1 = c(0.2,0.5), cortes2 = c(0.2,0.5),
                        fx1 = c(0.5,1.5,2), fx2 = c(0.5,1.5,2),
                        x1 = runif(200), x2 = runif(200), z = runif(200)){
    n = length(x1); K1 = length(cortes1) ; K2 = length(cortes2) ; 
    pred1 = numeric(n); pred2 = numeric(n); pred = numeric(n);
    pred1[x1 <= cortes1[1]] = fx1[1]
    pred1[x1 > cortes1[K1]] = fx1[K1+1]
    pred2[x2 <= cortes2[1]] = fx2[1]
    pred2[x2 > cortes2[K2]] = fx2[K2+1]
    for (k in 1:(K1-1)) pred1[cortes1[k] < x1 & x1 <= cortes1[k+1]] = fx1[k+1]
    for (j in 1:(K2-1)) pred2[cortes2[j] < x2 & x2 <= cortes2[j+1]] = fx2[j+1]
    
    pred = pred1 + pred2 + 0.1*z
    
    pred = exp(pred)/(1+exp(pred))
    
    list(x1 = x1, x2 = x2, z = z, pred = pred)
    }
