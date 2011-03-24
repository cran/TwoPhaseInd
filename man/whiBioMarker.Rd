\name{whiBioMarker}
\alias{whiBioMarker}
\docType{data}
\title{
  Women's Health Initiative (WHI) Biomarker Data
}
\description{
  This is the data used in Kooperberg et al. (2007) and in Dai et al. (2009).
}
\usage{data(whiBioMarker)}
\format{
  A data frame with 16608 observations on the following 10 variables.
  \describe{
    \item{\code{stroke}}{An indicator of stroke (1="Yes", 0="No")}
    \item{\code{hrtdisp}}{An indicator of treatment in the Estrogen Plus
	  PRogestin Trial (1="Estrogen Plus Progestin", 0="placebo")}
    \item{\code{papbl}}{Biomarker PAP (plasmin-antiplasmin complex) in logarithmic scale (base 10)}
    \item{\code{age}}{Age}
    \item{\code{dias}}{Diastolic BP (1="Yes", 0="No")}
    \item{\code{syst}}{Systolic BP}
    \item{\code{diabtrt}}{Diabetes with levels \code{Missing} \code{No} \code{Yes}}
    \item{\code{phase}}{Phase indicator: 1 for phase 1, and 2 for pahse 2.}
  }
}
\details{
  An example dataset adapted from Kooperberg et al. (2007) to demostrate
    the usage of MELE and SPMLE algorithms in Dai et al. (2009).
}
\source{
C. Kooperberg, M. Cushman, J. Hsia, J. G. Robinson, A. K. Aragaki,
J. K. Lynch, A. E. Baird, K. C. Johnson, L. H. Kuller, S. A. Beresford,
and B. Rodriguez. Can biomarkers identify women at increased stroke
risk? the women's health initiative hormone trials. PLoS clinical
trials, 2(6):e28, Jun 15 2007.   
}
\references{
J. Y. Dai, M. LeBlanc, and C. Kooperberg. Semiparametric estimation
exploiting co-variate independence in two-phase randomized
trials. Biometrics, 65(1):178-187, Mar 2009.
}
\examples{
data(whiBioMarker)
colnames(whiBioMarker)
dim(whiBioMarker)
## maybe str(whiBioMarker) ; plot(whiBioMarker) ...
}
\keyword{datasets}
