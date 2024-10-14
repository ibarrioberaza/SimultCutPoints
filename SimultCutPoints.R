# This method allows to find cut-off points for p variables simultaneously.
# The location of the cut-off points works well, but the selection of the number of points is still an issue. 
# At the moment it is based on the minimum BIC_pseudo obtained (an adaptation of the BIC).

# To navigate more easily through the code when using RStudio, we recommend to collapse it: Alt+o (Cmd+Alt+o on the Mac)

# Sources required --------------------------------------------------------
  source("Functions.R")

# Arguments ---------------------------------------------------------------
  # formula: formula with the response variable and other covariates in the data set to be included in the model.
  # cat.vars: a vector with the names of the variables in the data set to be categorised.
  # data: a data frame containing all needed variables.
  # family: family used in the GAM model.
  # k.min: a numeric value which indicates the minimum cut-off points to be considered for each cat.var.
  # k.max: a numeric value which indicates the maximum cut-off points to be considered for each cat.var.

# Values ------------------------------------------------------------------
  # model.smooth: Results for the GAM model
  # data: a data frame containing the original data frame plus the categorised cat.vars.
  # BICps: a data frame with the BICps for all possible models when combine the cat.vars with different k in (k.min,k.max)
  # Cuts: a list with arrays for each cat.var with the cut-off points for k in (k.min,k.max)
  # formula.gam: the formula used for the GAM model
  # formula.cat: the formula used for the model with the optimal (minimum BICps) cut-off points for each cat.var
  
# Function ----------------------------------------------------------------
  SimultCutPoints <- function(formula, cat.vars, data, family, k.min = 1, k.max = 2){ 
    data <- na.omit(data[,c(all.vars(formula), cat.vars)])
    formula.new <- update(formula, as.formula(paste0("~ . + ", paste0("s(", cat.vars, ", k = 30, m = 5, bs = 'ad')", collapse = "+"))))
    # formula.new.sop <- update(formula, as.formula(paste0("~ . -1 + ", paste0("ad(", cat.vars, ", nseg = 30, nseg.sp = 5)", collapse = "+"))))
    
    # 1. GAM model for the smooth relationship
      mod.gam <- mgcv::gam(formula = formula.new, data = data, family = family, method = "REML")
      # summary(mod.gam)
      # plot(mod.gam)
      
      # The following is the same model but using SOP
      # This is much faster than GAM 
      # mod.sop <- SOP::sop(formula.new.sop, family = binomial(link = "logit"), data = data)
      
    # 2. Predictions with the GAM model
      pred.gam  <-  predict(mod.gam, type = "iterms", se.fit = TRUE) 
      # The problem with SOP is the predictions are not the same for the smooth terms :(
      # Check what is the "problem"
      # pred.sop <- predict(mod.sop, type = "terms", se.fit = TRUE)
      # plot(pred.gam$fit[,1], pred.sop$fit[,1])
      
    # 3. Categorisation of p cat.vars with k cut-off points each
      # k in (k.min, k.max)
      # We can also think that each covariate has its own (k.min, k.max)
      # For simplicity we assume the same (k.min, k.max)
      # Optimal location based on the GAM-model proposal.
      cuts <- lapply(cat.vars, function(x){
        Cuts=array(dim = c(k.max, k.max))
        data.cut <- as.data.frame(array(dim = c(nrow(data), k.max)))
        for(k in k.min:k.max){
          Ax = Opt.cpoints.GAM(data[,x], 
                               pred.gam$fit[,paste0("s(",x,")")], 
                               w = 1/(pred.gam$se.fit[,paste0("s(",x,")")]**2), 
                               nc = k, 
                               nfino = 200)
          data.cut[,k] = cut(data[,x], breaks = Ax$res$cortes, include.lowest = TRUE) 
          names(data.cut)[k] <- paste0("k",k)
          Cuts[k,1:k] = Ax$res$cortes[-c(1,k+2)]
        }
        res <- list(data.cut = data.cut, Cuts = Cuts)
        return(res)
        })
      names(cuts) <- cat.vars
      
      # Joining the data.cut for all cat.vars
      data <- cbind(data, do.call("cbind", lapply(cuts, function(x) x$data.cut)))
  
    # 4. Selecting the optimal number of cut-off points
      # It is not taken into account whether the category is significant
      # eval <- expand.grid(lapply(cat.vars, function(x){paste0(x, ".", paste0("k",1:k.max))}), stringsAsFactors = FALSE)
      eval <- expand.grid(lapply(cat.vars, function(x){1:k.max}), stringsAsFactors = FALSE)
      names(eval) <- paste0("k_", cat.vars)
      eval$BICps <- NA
      for(i in 1:nrow(eval)){
        mod.cat <- mgcv::gam(update(formula, as.formula(paste0("~ . +", paste0(cat.vars,".k",eval[i,1:length(cat.vars)],collapse = "+")))), 
                             method = "REML", family = family, data = data)
        # nc <- sum(as.numeric(substring(eval[i,1:length(cat.vars)], first = 5)))
        nc <- sum(eval[i,1:length(cat.vars)])
        eval[i,"BICps"] = Pseudo_BIC2(mod.cat, nc = nc)
      }
      
      best <- (eval[order(eval$BICps),])[1,]
      cat("\n*************************************************\n")
      cat("Formula: \n")
      print(update(formula, as.formula(paste0("~ . +", paste0(cat.vars,".k", best[,1:length(cat.vars)], collapse = "+")))), quote = FALSE)
      cat(paste0("BICps = ", round(best[,length(cat.vars)+1],4)),"\n")
      cat("\n*************************************************\n")
      cat("Number of optimal cut-off points for each cat.var: \n \n")
      lapply(1:length(cat.vars), function(i){
        cat(paste0("For cat.var = ", cat.vars[i], ", take k = ", best[,i], ", with cut-off points at:"), "\n")
        cat(cuts[[cat.vars[i]]]$Cuts[best[,i],1:best[,i]], "\n \n")
      })
      
    # 5. Values to be returned 
      res <- list()
      res$model.smooth <- mod.gam
      res$data <- data
      res$BICps <- eval[order(eval$BICps),]
      res$Cuts <- lapply(cuts, function(x){x$Cuts})
      res$formula.gam <- formula.new
      res$formula.cat <- update(formula, as.formula(paste0("~ . +", paste0(cat.vars,".k", best[,1:length(cat.vars)], collapse = "+"))))
      res
  }
