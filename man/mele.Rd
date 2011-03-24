\name{mele}
\alias{mele}
\title{
  Maximum Estimated Likelihood Estimator (MELE)
}
\description{
  This is the function to calculate the maximum estimated likelihood estimator
  (MELE) for the treatment-biomarker interaction in two-phase sampling
  nested within randomized clinical trials, incorporating the independence
  between a randomized treatment and the baseline markers.
}
\usage{
mele(data, response, treatment, BaselineMarker, extra = NULL, phase, ind = TRUE, difffactor = 1e-07)
}
\arguments{
  \item{data}{
	A data.frame containing all the variables for the model
  }
  \item{response}{
	Quoted string to indicate the variable of response.  The variable must
	be binary.
  }
  \item{treatment}{
	Quoted string to indicate the variable of treatment in the first
	phase.  The variable must be binary.
  }
  \item{BaselineMarker}{
	Quoted string to indicate the baseline marker in the second phase.  The variable must
	be continuous variable.
  }
  \item{extra}{
	A vector of quoted column names to be used as extra variables in the model.
  }
  \item{phase}{
	Quoted string to indicate the variable of phase.
  }
  \item{ind}{
	TRUE for independent assumption; FALSE for non-independent assumption.
  }
  \item{difffactor}{
	A small number to control the iteration.
  }
%%  \item{maxit}{
%%	integer for the number of iteration.
%%  }
%%  \item{verbose}{
%%	TRUE to turn on the log file.
%%  }
}
\details{
  Please see Dai et al. (2009) for more details.
}
\value{
  It returns a data.frame object with the estimates, standard errors,
  and p-values for all parameters in the model.
}
\references{
J. Y. Dai, M. LeBlanc, and C. Kooperberg. Semiparametric estimation
exploiting co-variate independence in two-phase randomized
trials. Biometrics, 65(1):178-187, Mar 2009.
}
\author{
  James Y. Dai and Ting-Yuan Liu
}
\seealso{
  \code{\link{spmle}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
