median_impute <- function(x) {
  # Impute missing samples in vector as median
  x[is.na(x)] <- median(x, na.rm=TRUE)
  return(x)
}

abind3 <- function(...) abind::abind(..., along=3)

irf_fit <- function(x, y, train.id, test.id=NULL, k=2:20, null.cells=FALSE) {
  # Fit random forest model constrained to a select number of cell lines
  if (is.null(test.id)) test.id <- setdiff(1:nrow(x), train.id)
  
  # Fit initial model to assess variable importance
  fit <- iRF(x=x[train.id,], 
             y=as.factor(y)[train.id], 
             type='ranger', 
             n.iter=2)
  
  importance <- fit$rf.list$variable.importance
  
  # Fit models with constrained number of cell lines
  out <- lapply(k, function(kk) {
    # Select optimal cell lines
    if (!null.cells) id.select <- tail(sort(importance), kk)
    if (null.cells) id.select <- sample(importance, kk)
    xselect <- x[,names(id.select)]
    
    # Fit model constraining # of cell lines
    fitk <- iRF(x=xselect[train.id,], 
                y=as.factor(y)[train.id], 
                type='ranger', 
                n.iter=1)
    
    # Evaluate model performance
    ypred <- predict(fitk$rf.list, xselect[test.id,])
    accuracy <- mean(ypred$predictions == y[test.id])
    
    # Return selected cell lines
    selected <- colnames(x) %in% names(id.select)
    names(selected) <- colnames(x)
    return(list(accuracy=accuracy, selected=selected))
  })
  
  accuracy <- sapply(out, function(z) z$accuracy)
  selected <- sapply(out, function(z) z$selected)
  return(list(selected=selected, accuracy=accuracy))
}





################################################################################
# Old utility functions
################################################################################
model_accuracy <- function(fit, x, y, smax=5) {
  # Generate model predictions given a maximum number of cell lines
  
  # Evaluate number of cell lines used for each model (lambda value)
  
  betas <- lapply(coef(fit), function(b) as.matrix(b[-1,]))
  betas <- do.call(abind3, betas)
  
  # Collapse betas and identify maximal model under given constraint
  betas.collapse <- apply(betas != 0, MAR=1:2, any)
  s <- colSums(betas.collapse != 0)
  idmax <- max(which(s <= smax))
  
  # Get model predictions
  ypred <- predict(fit, x, type='response')
  ypred <- ypred[,,idmax]
  
  # Normalize predictions for class imblanace
  ypred.true <- mapply(function(ii, iy) ypred[ii, iy + 1], 1:length(y), y)
  ypred <- apply(ypred, MAR=2, function(z) z / mean(z))
  ypred.class <- apply(ypred, MAR=1, which.max) - 1
  
  # Return beta as average absolute coefficient
  betas <- rowMeans(abs(betas[, idmax, ]))
  
  
  return(list(accuracy=mean(ypred.class == y), betas=betas, ypred=ypred.true))
}

