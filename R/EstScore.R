## Supplemental functions
EstScore_Orthogonal_C<- function(X, F, A, mu, Q, info){
  .Call("EstScore_Orthogonal_C",
        as.double(X),
        as.double(F),
        as.double(A),
        as.double(mu),
        as.double(Q),
        as.double(info),
        package = "EstScore_C")
}

EstScore_Oblique_C<- function(X, F, A, mu, Q, info){
  .Call("EstScore_Oblique_C",
        as.double(X),
        as.double(F),
        as.double(A),
        as.double(mu),
        as.double(Q),
        as.double(info),
        package = "EstScore_C")
}


## Main function
EstScore <- function (X, A, mu, N.ite=10000, N.random=1, show.random.ite=FALSE, oblique=FALSE, mc.cores=1) {

  N.ite.GP <- 100 ## the number of iterations for GP algorithm
  N.ite.alpha <- 30 ## the number of iterations for searching alpha in GP algorithm

  ##indices
  N.sub <- dim(X)[1]
  N.var <- dim(X)[2]
  N.comp <- dim(A)[2]
  Q <- 2 * X - 1

  ##----- several random starts
  ret <- NULL
  ret.val.lossfunc <- Inf

  ## Use parallel package?
  if (mc.cores > 1) {
    if (requireNamespace("parallel", quietly=TRUE)) {

      temp.solution <- parallel::mclapply(1:N.random, FUN = function(n.random) {
        if (show.random.ite)
          if (n.random %% (N.random / 10) == 0)
            cat(paste(n.random, "; ", sep=""))

        ##Initialization of parameters (use sparse logistic PCA)
        set.seed(n.random)
        F <- qr.Q(qr(matrix(rnorm(N.sub * N.comp), N.sub, N.comp))) ##Orthogonalization

        ##Optimization by ALS algorithm
        info <- c(N.sub, N.var, N.comp, N.ite, N.ite.GP, N.ite.alpha)

        if (oblique) {
          solution <- EstScore_Oblique_C(X, F, A, mu, Q, info)
        } else {
          solution <- EstScore_Orthogonal_C(X, F, A, mu, Q, info)
        }

        ##Solutions
        ans <- list()
        ans$F <- matrix(solution[[1]], N.sub, N.comp)
        ans$n.ite <- solution[[2]][1]
        ans$loss <- solution[[2]][2]


        ##Convergence?
        if (length(which(is.na(solution))) > 1 || length(which(is.nan(ans$F))) > 1) {
          ans$conv <- 0
        } else {
          ans$conv <- 1
        }

        return(ans)
      }, mc.cores=mc.cores)

      best.loss <- sapply(temp.solution,
                          FUN = function(x) if (x$conv != 1) NA
                        else x$loss
                          )
      if (all(is.na(best.loss)))
        stop("\nmywarning-->Could not find a feasible starting point...exiting\n",
             call. = FALSE)
      nb <- which(best.loss == min(best.loss, na.rm = TRUE))[1]
      solution <- temp.solution[[nb]]
      F <- solution$F
      loss <- solution$loss
      n.ite <- solution$n.ite

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
      F <- qr.Q(qr(matrix(rnorm(N.sub * N.comp), N.sub, N.comp))) ##Orthogonalization

      ##Optimization by ALS algorithm
      info <- c(N.sub, N.var, N.comp, N.ite, N.ite.GP, N.ite.alpha)
      if (oblique) {
        solution <- EstScore_Oblique_C(X, F, A, mu, Q, info)
      } else {
        solution <- EstScore_Orthogonal_C(X, F, A, mu, Q, info)
      }

      ##Compare the values of loss function
      if (solution[[2]][2] < ret.val.lossfunc) {
        ret.val.lossfunc <- solution[[2]][2]

        F <- matrix(solution[[1]], N.sub, N.comp)
        n.ite <- solution[[2]][1]
        loss <- solution[[2]][2]
      }
    }
  }

  if (show.random.ite)
    cat("\n")

  class(A) <- "loadings"

  return(list("F"=F, "n.ite"=n.ite, "loss"=loss))
}

