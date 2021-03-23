median_impute <- function(x) {
  # Impute missing samples in vector as median
  x[is.na(x)] <- median(x, na.rm=TRUE)
  return(x)
}

abind3 <- function(...) abind::abind(..., along=3)

model_accuracy <- function(fit, x, y, smax=5) {
  # Generate model predictions given a maximum number of cell lines
  
  # Evaluate number of cell lines used for each model (lambda value)
  betas <- lapply(coef(fit$glmnet.fit), function(b) as.matrix(b[-1,]))
  betas <- do.call(abind3, betas)
  
  # Collapse betas and identify maximal model under given constraint
  betas.collapse <- apply(betas != 0, MAR=1:2, any)
  s <- colSums(betas.collapse != 0)
  idmax <- max(which(s <= smax))
  
  # Get model predictions
  ypred <- predict(fit$glmnet.fit, x, type='response')
  ypred <- ypred[,,idmax]
  
  # Normalize predictions for class imblanace
  ypred.true <- mapply(function(ii, iy) ypred[ii, iy + 1], 1:length(y), y)
  ypred.class <- apply(ypred, MAR=1, which.max) - 1
  
  # Return beta as average absolute coefficient
  betas <- rowMeans(abs(betas[, idmax, ]))
  
  
  return(list(accuracy=mean(ypred.class == y), betas=betas, ypred=ypred.true))
}

