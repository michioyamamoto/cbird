## Supplemental functions
OptimCBIRD_C <- function(X, F, A, g, U, mu, Q, info, lambda, pll.record){
  .Call("OptimCBIRD_C",
        as.double(X),
        as.double(F),
        as.double(A),
        as.double(g),
        as.double(U),
        as.double(mu),
        as.double(Q),
        as.double(info),
        as.double(lambda),
        as.double(pll.record),
        package = "OptimCBIRD_C")
}


## Main function
cbird <- function (Y, N.comp, N.clust, lambda=0, N.ite=10000, N.random=1, show.random.ite=FALSE, eps=0.0001, mc.cores=1) {
  N.ite.clust <- 100 ##the maximum number of iterations for k-means algorithm
  N.ite.GP <- 100
  N.ite.alpha <- 30

  ##indices
  N.sub <- dim(Y)[1]
  N.var <- dim(Y)[2]
  vec.1.sub <- rep(1, N.sub)
  Q <- 2 * Y - 1

  ret <- NULL
  ret.val.pllfunc <- Inf
  ptm = proc.time()
  best.pll <- NULL

  ##For checking local optimum and uniqueness of a solution
  pll.record.mat <- matrix(0, N.ite, N.random)
  U.all.mat <- NULL


  ## Use parallel package?
  if (mc.cores > 1) {
    if (requireNamespace("parallel", quietly=TRUE)) {

      temp.solution <- parallel::mclapply(1:N.random, FUN = function(n.random) {
        if (show.random.ite)
          if (n.random %% (N.random / 10) == 0)
            cat(paste(n.random, "; ", sep=""))

        ##Initialization of parameters (use sparse logistic PCA)
        set.seed(n.random)

        mu <- rnorm(N.var)
        F <- qr.Q(qr(matrix(rnorm(N.clust * N.comp), N.clust, N.comp))) ##Orthogonalization
        A <- matrix(rnorm(N.var * N.comp), N.var, N.comp)

        h <- hclust(dist(as.matrix(Y)), method = "average")
        initial <- tapply(as.matrix(Y), list(rep(cutree(h, N.clust), ncol(Y)), col(as.matrix(Y))), mean)
        cluster <- kmeans(Y, initial, nstart=N.ite.clust)$cluster
        g <- numeric(N.clust)
        for (n.clust in 1:N.clust)
          g[n.clust] <- length(which(cluster == n.clust)) / N.sub
        U <- matrix(0, N.sub, N.clust)

        pll.record <- numeric(N.ite)


        ##Optimization by EM + Majorization algorithm
        info <- c(N.sub, N.var, N.comp, N.ite, N.ite.GP, N.ite.alpha, N.clust, eps)
        solution <- OptimCBIRD_C(Y, F, A, g, U, mu, Q, info, lambda, pll.record)

        ans <- list()
        ans$F <- solution[[1]]
        ans$A <- solution[[2]]
        ans$g <- solution[[3]]
        ans$U <- solution[[4]]
        ans$mu <- solution[[5]]
        ans$n.ite <- solution[[6]][1]
        ans$pll <- solution[[6]][2] ##penalized log likelihood
        ans$pll.record <- solution[[7]] ##record of penalized log likelihoods

        return(ans)
      }, mc.cores=mc.cores)

      ##find the solution which has maximum value of penalized log likelihood
      best.pll <- sapply(temp.solution,
                         FUN = function(x) x$pll
                         )
      if (all(is.na(best.pll)))
        stop("\nmywarning-->Could not find a feasible starting point...exiting\n",
             call. = FALSE)
      nb <- which(best.pll == max(best.pll, na.rm = TRUE))[1]

      ##record of penalized log likelihoods
      pll.record.mat <- sapply(temp.solution, FUN = function (x) x$pll.record)

      solution <- temp.solution[[nb]]
      F <- solution$F
      A <- solution$A
      mu <- solution$mu
      U <- solution$U
      g <- solution$g
      n.ite <- solution$n.ite
      val.pll <- solution$pll

    } else { ##if parallel package has not been installed
      if(nchar(system.file(package="parallel")) == 0)
        stop("Package 'parallel' is required to calculate using multicores.")
    }

  } else {
    ## Not use parallel package
    for (n.random in 1:N.random) {
      if (show.random.ite)
        if (n.random %% (N.random / 10) == 0)
          cat(paste(n.random, "; ", sep=""))

      ##Initialization of parameters (use sparse logistic PCA)
      set.seed(n.random)

      mu <- numeric(N.var)
      F <- qr.Q(qr(matrix(rnorm(N.clust * N.comp), N.clust, N.comp))) ## Orthogonalization
      A <- matrix(rnorm(N.var * N.comp), N.var, N.comp)

      h <- hclust(dist(as.matrix(Y)), method = "average")
      initial <- tapply(as.matrix(Y), list(rep(cutree(h, N.clust), ncol(Y)), col(as.matrix(Y))), mean)
      cluster <- kmeans(Y, initial, nstart=N.ite.clust)$cluster
      g <- numeric(N.clust)
      for (n.clust in 1:N.clust)
        g[n.clust] <- length(which(cluster == n.clust)) / N.sub
      U <- matrix(0, N.sub, N.clust)

      pll.record <- numeric(N.ite)

      ##Optimization by EM + Majorization algorithm
      info <- c(N.sub, N.var, N.comp, N.ite, N.ite.GP, N.ite.alpha, N.clust, eps)
      solution <- OptimCBIRD_C(Y, F, A, g, U, mu, Q, info, lambda, pll.record)

      pll.record.mat[, n.random] <- solution[[7]] ##record of penalized log likelihoods

      ##Compare the values of loss functions
      if (solution[[6]][2] < ret.val.pllfunc) {
        ret.val.pllfunc <- solution[length(solution)][[1]][2] ##why use solution[[6]][2]? (130916)

        F.ret <- solution[[1]]
        A.ret <- solution[[2]]
        g.ret <- solution[[3]]
        U.ret <- solution[[4]]
        mu.ret <- solution[[5]]
        n.ite.ret <- solution[[6]][1]
        val.pll.ret <- solution[[6]][2]
      }
    }
    F <- F.ret; A <- A.ret; g <- g.ret; U <- U.ret; mu <- mu.ret; n.ite <- n.ite.ret; val.pll <- val.pll.ret; pll.record <- pll.record.mat
  }


  ##Log likelihood (not with penalty)
  LL <- val.pll + N.sub * lambda * sum(abs(A))

  ##BIC
  bic <- c(-2 * LL + log(N.sub) * (N.var + N.clust * N.comp + length(which(!(abs(A) < 1e-7))) + N.clust))

  if (show.random.ite)
    cat("\n")

  class(A) <- "loadings"
  class(U) <- "loadings"
  cluster <- numeric(N.sub)
  for (n.sub in 1:N.sub)
    cluster[n.sub] <- which(U[n.sub, ] == max(U[n.sub, ]))[1]

  ##proc time
  ptm.ret <- proc.time() - ptm

  ##estimate clusters of subjects for each initial value
  if (!is.null(U.all.mat)) {
    cluster.all <- matrix(0, N.random, N.sub)
    for (n.random in 1:N.random) {
      U.temp <- matrix(U.all.mat[, n.random], N.sub, N.clust)
      for (n.sub in 1:N.sub)
        cluster.all[n.random, n.sub] <- which(U.temp[n.sub, ] == max(U.temp[n.sub, ]))[1]
    }
  }

##  return(list("F"=F, "A"=A, "mu"=mu, "U"=U, "g"=g, "n.ite"=n.ite, "pll"=val.pll, "bic"=bic, "LL"=LL, "cluster"=cluster, "ptime"=ptm.ret, "pll.seq"=best.pll, "pll.record"=pll.record.mat, "cluster.all"=cluster.all, "mu.all"=mu.all.mat, "A.all"=A.all.mat))
  return(list("F"=F, "A"=A, "mu"=mu, "U"=U, "g"=g, "n.ite"=n.ite, "pll"=val.pll, "bic"=bic, "LL"=LL, "cluster"=cluster, "ptime"=ptm.ret, "pll.seq"=best.pll, "pll.record"=pll.record.mat))
}

