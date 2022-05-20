### These functions are based on the same function in "lasso.proj" of package "hdi".


calculate.Z<-function (x, wantTheta=F, parallel=F, ncores, verbose, Z, do.ZnZ = FALSE, debug.verbose = FALSE) 
{
  if (is.null(Z)) {
    message("Nodewise regressions will be computed as no argument Z was provided.")
    message("You can store Z to avoid the majority of the computation next time around.")
    message("Z only depends on the design matrix x.")
    nodewiselasso.out <- score.nodewiselasso(x = x, wantTheta = F, verbose=verbose,
                                             parallel = parallel, ncores = ncores, cv.verbose = verbose || 
                                               debug.verbose, do.ZnZ = do.ZnZ)
    Z <- nodewiselasso.out$out$Z
    scaleZ <- nodewiselasso.out$out$scaleZ
  }
  else {
    scaleZ <- rep(1, ncol(Z))
    if (!isTRUE(all.equal(rep(1, ncol(x)), colSums(Z * x)/nrow(x), 
                          tolerance = 10^-8))) {
      rescale.out <- score.rescale(Z = Z, x = x)
      Z <- rescale.out$Z
      scaleZ <- rescale.out$scaleZ
    }
  }
  list(Z = Z, scaleZ = scaleZ)
}

score.nodewiselasso<-function (x, wantTheta = FALSE, verbose = FALSE, lambdaseq = "quantile", 
                               parallel = FALSE, ncores = 8, oldschool = FALSE, lambdatuningfactor = 1, 
                               cv.verbose = FALSE, do.ZnZ = TRUE) 
{
  lambdas <- switch(lambdaseq, quantile = nodewise.getlambdasequence(x), 
                    linear = nodewise.getlambdasequence.old(x, verbose), 
                    stop("invalid 'lambdaseq': ", lambdaseq))
  if (verbose) {
    cat("Using the following lambda values:", lambdas, "\n")
  }
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x, 
                                      parallel = parallel, ncores = ncores, oldschool = oldschool, 
                                      verbose = cv.verbose)
  if (verbose) {
    cat(paste("lambda.min is", cvlambdas$lambda.min), "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), "\n")
  }
  if (do.ZnZ) {
    bestlambda <- improve.lambda.pick(x = x, parallel = parallel, 
                                      ncores = ncores, lambdas = lambdas, bestlambda = cvlambdas$lambda.min, 
                                      verbose = verbose)
    if (verbose) {
      cat("Doing Z&Z technique for picking lambda\n")
      cat("The new lambda is", bestlambda, "\n")
      cat("In comparison to the cross validation lambda, lambda = c * lambda_cv\n")
      cat("c=", bestlambda/cvlambdas$lambda.min, "\n")
    }
  }
  else {
    if (lambdatuningfactor == "lambda.1se") {
      if (verbose) 
        cat("lambda.1se used for nodewise tuning\n")
      bestlambda <- cvlambdas$lambda.1se
    }
    else {
      if (verbose) 
        cat("lambdatuningfactor used is", lambdatuningfactor, 
            "\n")
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  }
  if (verbose) {
    cat("Picked the best lambda:", bestlambda, "\n")
  }
  if (wantTheta) {
    out <- score.getThetaforlambda(x = x, lambda = bestlambda, 
                                   parallel = parallel, ncores = ncores, oldschool = TRUE, 
                                   verbose = verbose)
  }
  else {
    Z <- score.getZforlambda(x = x, lambda = bestlambda, 
                             parallel = parallel, ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out, bestlambda = bestlambda)
  return(return.out)
}

est.stderr.despars.lasso<-function (x, y, Z, betalasso, sigmahat, robust, robust.div.fixed = FALSE) 
{
  if (robust) {
    stderr <- sandwich.var.est.stderr(x = x, y = y, Z = Z, 
                                      betainit = betalasso)
    n <- nrow(x)
    if (robust.div.fixed) 
      stderr <- stderr * n/(n - sum(betalasso != 0))
  }
  else {
    stderr <- (sigmahat * sqrt(diag(crossprod(Z))))/nrow(x)
  }
  stderr
}

sandwich.var.est.stderr<-function (x, y, betainit, Z) 
{
  n <- nrow(x)
  p <- ncol(x)
  if (!isTRUE(all.equal(rep(1, p), colSums(Z * x)/n, tolerance = 10^-8))) {
    rescale.out <- score.rescale(Z = Z, x = x)
    Z <- rescale.out$Z
  }
  if (length(betainit) > ncol(x)) {
    x <- cbind(rep(1, nrow(x)), x)
  }
  else {
    if (all(x[, 1] == 1)) {
    }
    else {
    }
  }
  eps.tmp <- as.vector(y - x %*% betainit)
  eps.tmp <- eps.tmp - mean(eps.tmp)
  sigmahatZ.direct <- sqrt(colSums(sweep(eps.tmp * Z, MARGIN = 2, 
                                         STATS = crossprod(eps.tmp, Z)/n, FUN = `-`)^2))
  sigmahatZ.direct/n
}


score.rescale<-function (Z, x) 
{
  scaleZ <- diag(crossprod(Z, x))/nrow(x)
  Z <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z = Z, scaleZ = scaleZ))
}

nodewise.getlambdasequence<-function (x) 
{
  nlambda <- 100
  p <- ncol(x)
  lambdas <- c()
  for (c in 1:p) {
    lambdas <- c(lambdas, glmnet(x[, -c], x[, c])$lambda)
  }
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}


cv.nodewise.bestlambda<-function (lambdas, x, K = 10, parallel = FALSE, ncores = 8, oldschool = FALSE, 
                                  verbose = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)
  dataselects <- sample(rep(1:K, length = n))
  if (oldschool) {
    message("doing cv.nodewise.error oldschool")
    totalerr <- numeric(l)
    for (c in 1:p) {
      for (i in 1:K) {
        whichj <- dataselects == i
        glmnetfit <- glmnet(x[!whichj, -c, drop = FALSE], 
                            x[!whichj, c, drop = FALSE], lambda = lambdas)
        predictions <- predict(glmnetfit, x[whichj, -c, 
                                            drop = FALSE], s = lambdas)
        totalerr <- totalerr + apply((x[whichj, c] - 
                                        predictions)^2, 2, mean)
      }
    }
    totalerr <- totalerr/(K * p)
    stop("lambda.1se not implemented for oldschool cv.nodewise.bestlamba")
  }
  else {
    if (parallel) {
      totalerr <- mcmapply(cv.nodewise.err.unitfunction, 
                           c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                           x = list(x = x), lambdas = list(lambdas = lambdas), 
                           mc.cores = ncores, SIMPLIFY = FALSE, verbose = verbose, 
                           p = p)
    }
    else {
      totalerr <- mapply(cv.nodewise.err.unitfunction, 
                         c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                         x = list(x = x), lambdas = list(lambdas = lambdas), 
                         SIMPLIFY = FALSE, verbose = verbose, p = p)
    }
    err.array <- array(unlist(totalerr), dim = c(length(lambdas), 
                                                 K, p))
    err.mean <- apply(err.array, 1, mean)
    err.se <- apply(apply(err.array, c(1, 2), mean), 1, sd)/sqrt(K)
  }
  pos.min <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]
  stderr.lambda.min <- err.se[pos.min]
  list(lambda.min = lambda.min, lambda.1se = max(lambdas[err.mean < 
                                                           (min(err.mean) + stderr.lambda.min)]))
}

cv.nodewise.err.unitfunction<-function (c, K, dataselects, x, lambdas, verbose, p) 
{
  if (verbose) {
    interesting.points <- round(c(1/4, 2/4, 3/4, 4/4) * p)
    names(interesting.points) <- c("25%", "50%", "75%", "100%")
    if (c %in% interesting.points) {
      message("The expensive computation is now ", names(interesting.points)[c == 
                                                                               interesting.points], " done")
    }
  }
  cv.nodewise.totalerr(c = c, K = K, dataselects = dataselects, 
                       x = x, lambdas = lambdas)
}

cv.nodewise.totalerr<-function (c, K, dataselects, x, lambdas) 
{
  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  for (i in 1:K) {
    whichj <- dataselects == i
    glmnetfit <- glmnet(x = x[!whichj, -c, drop = FALSE], 
                        y = x[!whichj, c, drop = FALSE], lambda = lambdas)
    predictions <- predict(glmnetfit, newx = x[whichj, -c, 
                                               drop = FALSE], s = lambdas)
    totalerr[, i] <- apply((x[whichj, c] - predictions)^2, 
                           2, mean)
  }
  totalerr
}
score.getZforlambda<-function (x, lambda, parallel = FALSE, ncores = 8, oldschool = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n * p), n)
  if (oldschool) {
    message("doing getZforlambda oldschool")
    for (i in 1:p) {
      glmnetfit <- glmnet(x[, -i], x[, i])
      prediction <- predict(glmnetfit, x[, -i], s = lambda)
      Z[, i] <- x[, i] - prediction
    }
  }
  else {
    if (parallel) {
      Z <- mcmapply(score.getZforlambda.unitfunction, i = 1:p, 
                    x = list(x = x), lambda = lambda, mc.cores = ncores)
    }
    else {
      Z <- mapply(score.getZforlambda.unitfunction, i = 1:p, 
                  x = list(x = x), lambda = lambda)
    }
  }
  Z <- score.rescale(Z, x)
  return(Z)
}

score.getZforlambda.unitfunction<-function (i, x, lambda) 
{
  glmnetfit <- glmnet(x[, -i], x[, i])
  prediction <- predict(glmnetfit, x[, -i], s = lambda)
  return(x[, i] - prediction)
}


score.getThetaforlambda<-function (x, lambda, parallel = FALSE, ncores = 8, oldschool = FALSE, 
                                   verbose = FALSE, oldtausq = TRUE) 
{
  message("Calculating Thetahat by doing nodewise regressions and dropping the unpenalized intercept")
  n <- nrow(x)
  p <- ncol(x)
  C <- diag(rep(1, p))
  T2 <- numeric(p)
  if (oldschool) {
    message("doing getThetaforlambda oldschool")
    for (i in 1:p) {
      glmnetfit <- glmnet(x[, -i], x[, i])
      coeffs <- as.vector(predict(glmnetfit, x[, -i], type = "coefficients", 
                                  s = lambda))[-1]
      C[-i, i] <- -as.vector(coeffs)
      if (oldtausq) {
        T2[i] <- as.numeric(crossprod(x[, i])/n - x[, 
                                                    i] %*% (x[, -i] %*% coeffs)/n)
      }
      else {
        T2[i] <- as.numeric((x[, i] %*% (x[, i] - predict(glmnetfit, 
                                                          x[, -i], s = lambda)))/n)
      }
    }
  }
  else {
    stop("not implemented yet!")
  }
  thetahat <- C %*% solve(diag(T2))
  if (verbose) {
    cat("1/tau_j^2:", solve(diag(T2)), "\n")
  }
  thetahat <- t(thetahat)
  if (all(thetahat[lower.tri(thetahat)] == 0, thetahat[upper.tri(thetahat)] == 
          0) && verbose) 
    cat("Thetahat is a diagonal matrix!\n")
  return(thetahat)
}

prepare.data<-function (x, y, standardize, family) 
{
  x <- scale(x, center = TRUE, scale = standardize)
  dataset <- switch(family, gaussian = {
    list(x = x, y = y)
  }, {
    switch.family(x = x, y = y, family = family)
  })
  x <- scale(dataset$x, center = TRUE, scale = FALSE)
  y <- scale(dataset$y, scale = FALSE)
  y <- as.numeric(y)
  list(x = x, y = y)
}
