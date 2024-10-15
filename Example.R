
# Sources required --------------------------------------------------------
source("SimultCutPoints.R")

# Data --------------------------------------------------------------------
  # Urine analysis data from boot package
  data(urine, package = "boot")
  # Converting r (Indicator of the presence of calcium oxalate crystals, 
  # with values 0,1) into factor
  urine$r <- as.factor(urine$r) 
  GGally::ggpairs(urine[,c(2:7,1)], ggplot2::aes(colour = r))
  
# Example 1 ---------------------------------------------------------------
  # Binomial predictor variable (r)
  # Simultaneously categorise three continuous covariates (ph, calc and urea)
  Ex1 <- SimultCutPoints(formula = r ~ 1, data = urine,
                         cat.vars = c("ph","calc","urea"), 
                         family = "binomial", k.min = 1, k.max = 4, k.gam = 10)
  
  summary(Ex1$model.smooth)
  
  par(mfrow = c(1,3))
  plot(Ex1$model.smooth, select = 1)
  abline(v = Ex1$Cuts$ph[1,], col = "gray", lwd = 2) # k = 1 cut-off points
  plot(Ex1$model.smooth, select = 2)
  abline(v = Ex1$Cuts$calc[2,], col = "gray", lwd = 2) # k = 2 cut-off points
  plot(Ex1$model.smooth, select = 3)
  abline(v = Ex1$Cuts$urea[1,], col = "gray", lwd = 2) # k = 1 cut-off points
  
  # Best five BICps
  head(Ex1$BICps, n = 5)
  # Worst five BICps
  tail(Ex1$BICps, n = 5)
  
# Example 2 ---------------------------------------------------------------
  # Gaussian predictor variable (cond)
  # Simultaneously categorise four continuous covariates 
  # (gravity, ph, osmo and urea)
  Ex2 <- SimultCutPoints(formula = cond ~ 1, data = urine,
                         cat.vars = c("gravity","ph","osmo","urea"), 
                         family = "gaussian", k.min = 1, k.max = 4, k.gam = 10)
  
  summary(Ex2$model.smooth)
  
  par(mfrow = c(2,2))
  plot(Ex2$model.smooth, select = 1)
  abline(v = Ex2$Cuts$gravity[1,], col = "gray", lwd = 2) # k = 1 cut-off points
  plot(Ex2$model.smooth, select = 2)
  abline(v = Ex2$Cuts$ph[1,], col = "gray", lwd = 2) # k = 1 cut-off points
  plot(Ex2$model.smooth, select = 3)
  abline(v = Ex2$Cuts$osmo[4,], col = "gray", lwd = 2) # k = 4 cut-off points
  plot(Ex2$model.smooth, select = 4)
  abline(v = Ex2$Cuts$urea[1,], col = "gray", lwd = 2) # k = 1 cut-off points

  # Best five BICps
  head(Ex2$BICps, n = 5)
  # Worst five BICps
  tail(Ex2$BICps, n = 5)


  
  


  
  

  