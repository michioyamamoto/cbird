\name{cbird}
\alias{cbird}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Clustering of multivariate binary data with dimension reduction via
	L1-regularized likelihood maximization.
}
\description{
	This function conducts the clustering of binary data with reducing the
	dimensionality (CLUSBIRD) proposed by Yamamoto and Hayashi (2015).
}
\usage{
cbird(Y, N.comp, N.clust, lambda=0, N.ite=10000, N.random=1,
     show.random.ite=FALSE, eps=0.0001, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
		Binary data matrix (N * D), where N denotes sample size and D
		denotes the number of binary variables (0 or 1).
	}
  \item{N.comp}{
		The number of component to be extracted.
	}
  \item{N.clust}{
		The number of mixture components, which corresponds to the number of clusters.
	}
  \item{lambda}{
		A tuning parameter of an L1 penalty for loadings. A non-negative
		real value should be used as the value of lambda.
	}
  \item{N.ite}{
		The number of maximum of iterations for the EM algorithm.
	}
  \item{N.random}{
		The number of random sets of parameters for initial random starts.
	}
  \item{show.random.ite}{
		If \code{"TRUE"}, the number of each iteration is shown on the R console.
	}
  \item{eps}{
		The criterion for the convergence of the alternating least-squares
		algorithm, which should be specified as a positive real value. If
		the difference between the values of penalized log likelihood
		functions of successive iteration is smaller than \code{eps}, then
		\code{cbird} makes a decision about the convergence of the
		algorithm.
	}
  \item{mc.cores}{
		If \code{"parallel"} package has been installed, \code{"cbird"}
		adopts a multithread process for multiple initial random starts. If
		\code{"mc.cores"=1}, \code{"parallel"} package is not needed, and a
		single core process is conducted.
	}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
	\item{F}{An estimated component score matrix for cluster centroids.}
	\item{A}{An estimated loading matrix.}
	\item{mu}{Estimated mean values in the subspace.}
	\item{U}{The cluster assignment matrix (N * \code{N.clust}).}
	\item{g}{The estimated mixture probability.}
	\item{n.ite}{The number of iteration needed for convergence.}
	\item{loss}{The value of log likelihood with L1 penalty.}
	\item{bic}{The value of BIC.}
	\item{LL}{The value of log likelihood.}
	\item{cluster}{Estimated clusters where subjects were assigned to}
	\item{ptime}{Time for calculation}
}
\references{
	Yamamoto, M. and Hayashi, K. (2015). Clustering of multivariate binary
	data with dimension reduction via L1-regularized maximization. Pattern
	Recognition, 48, 3959-3968.
}
\author{
	Michio Yamamoto\cr
	\email{michio.koko@gmail.com}
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
##Random Binary Data (unmeaningful example)
##100 subjects and 20 variables
##Consider three mixture components in the data
set.seed(1)
Y <- matrix(rbinom(100 * 20, 1, 0.5), 100, 20)
out <- cbird(Y, 2, 3)

est <- EstScore(Y, out$A, out$mu)
}

