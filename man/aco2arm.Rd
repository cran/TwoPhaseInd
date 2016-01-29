\name{aco2arm}
\alias{aco2arm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate parameters in augmented case-only designs, the genotype is ascertained
for a random subcohort from both the active treatment arm and the placebo arm
}
\description{
This function estimates parameters of proportional hazards model with gene-treatment
interaction. It employs case-cohort estimation incorporating the case-only estimators. The method was published in Dai et al. (2015) Biometrics.
}
\usage{
aco2arm(data, svtime, event, treatment, BaselineMarker, id, 
subcohort, esttype = 1, extra)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame used to access the following data.
}
  \item{svtime}{
A character string of column name, corresponds to one column of the data frame, which is used to store the failure time variable (numeric).
}
  \item{event}{
A character string of column name, corresponds to one column of the data frame,which is used to store the indicator of failure event (1: failure, 0: not failure).
}
  \item{treatment}{
A character string of column name, corresponds to one column of the data frame, which is used to store the binary vector of treatment variable (1: treatment, 0: placebo).
}
  \item{BaselineMarker}{
A character string of column name, corresponds to one column of the data frame, which is used to store a vector of biomarker.
}
  \item{id}{
A character string of column name, corresponds to one column of the data frame, which is used to store the sample identifier.
}
  \item{subcohort}{
A character string of column name, corresponds to one column of the data frame, which is used to store the indicator of sub-cohort (1: sample belong to the sub-cohort, 0: not belong to the sub-cohort)
}
  \item{esttype}{
The option of estimation methods (1: Self-Prentice estimator, 0: Lin-Ying estimator).
}
  \item{extra}{
A string vector of column name(s), corresponds to more or more column(s) of the data frame, which is/are used to store the extra baseline covariate(s) to be adjusted for in addition to treatment and biomarker.
}
}
\details{
The function returns estimates of the proportional hazards model, and variance of the estimates. The method was published in Dai et al. (2015) Biometrics.
}
\value{
\item{beta}{Estimated parameter}
\item{stder}{Standard error}
\item{pVal}{p value}
}
\references{
J. Y. Dai, XC Zhang,CY Wang, and C. Kooperberg. Augmented Case-Only Designs for Randomized Clinical Trials with Failure Time Endpoints. Biometrics, DOI: 10.1111/biom.12392, 2015.
}
\author{
James Y. Dai
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{aco1arm}}
}
\examples{
## Load the example data
data(acodata)
## Case-cohort + case-only estimators
rfit1 <- aco2arm(data=acodata,
                 svtime="vacc1_evinf",
                 event="f_evinf",
                 treatment="f_treat",
                 BaselineMarker="fcgr2a.3",
                 id="ptid",
                 subcohort="subcoh",
                 esttype=1,
                 extra=c("f_agele30","f_hsv_2","f_ad5gt18","f_crcm",
                 "any_drug","num_male_part_cat","uias","uras")) 
rfit1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ case-only designs }
\keyword{ case-cohort designs }% __ONLY ONE__ keyword per line
