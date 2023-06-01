\name{EstScore}
\alias{EstScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimate compontent scores for each subject using the result of cbird.
}
\description{
	This function estimates components scores for each subject using the
	result of CLUSBIRD.
}
\usage{
EstScore(X, A, mu, N.ite=10000, N.random=1, show.random.ite=FALSE,
oblique=FALSE, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
		Binary data matrix (N * D).
	}
  \item{A}{
		Loading matrix (D * L) estimated by \code{cbird}.
	}
  \item{mu}{
		A D-length mean vector estimated by \code{cbird}.
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
  \item{oblique}{
		If \code{"TRUE"}, the oblique component scores F are estimated. The
		default is \code{"FALSE"}.
	}
  \item{mc.cores}{
		If \code{"parallel"} package has been installed, \code{"EstScore"}
		adopts a multithread process for multiple initial random starts. If
		\code{"mc.cores"=1}, \code{"parallel"} package is not needed, and a
		single core process is conducted.
	}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
	\item{F}{An estimated component score matrix (N * D) containing scores
				for subjects.}
	\item{n.ite}{The number of iteration needed for convergence.}
	\item{loss}{The value of loss function used in ALS algorithm}
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
##See the example of the function "cbird".
}