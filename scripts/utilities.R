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

glm_fit <- function(x, y, train.id, test.id=NULL, k=2:20, null.cells=FALSE) {
  # Fit random forest model constrained to a select number of cell lines
  if (is.null(test.id)) test.id <- setdiff(1:nrow(x), train.id)
  
  # Fit initial model to assess variable importance
  x <- as.matrix(x)
  fit <- cv.glmnet(x=x[train.id,], y=y[train.id], family='multinomial', nfolds=5)

  # Read coefficients at each level of lambda
  betas <- lapply(fit$glmnet.fit$beta, function(b) as.matrix(b != 0))
  betas <- do.call(abind3, betas)
  betas <- apply(betas, MAR=1:2, any)
  counts <- colSums(betas)
  
  # Read predictions at each level of lambda
  ypred <- predict(fit$glmnet.fit, x[test.id,])
  
  # Fit models with constrained number of cell lines
  out <- lapply(k, function(kk) {
    if (!null.cells) id.select <- max(which(counts <= kk))
    if (null.cells) id.select <- sample(1:ncol(x), kk)
    
    ypred.select <- ypred[,,id.select]
    ypred.class <- apply(ypred.select, MAR=1, which.max) - 1
    
    # Evaluate model performance
    accuracy <- mean(ypred.class == y[test.id])
    
    # Return selected cell lines
    selected <- betas[,id.select]
    return(list(accuracy=accuracy, selected=selected))
  })
  
  accuracy <- sapply(out, function(z) z$accuracy)
  selected <- sapply(out, function(z) z$selected)
  return(list(selected=selected, accuracy=accuracy))
}