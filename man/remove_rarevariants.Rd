\name{remove_rarevariants}
\alias{remove_rarevariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function used in spmle and acoarm to remove rare-variant covariates
}
\description{
It is used to remove a rare-variant covariates, which can cause divergence problem.
}
\usage{
remove_rarevariants(data, cutoff = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A dataframe composed of covariates.
}
  \item{cutoff}{
Proportion cutoff. If data composed of more than (1-cutoff) proportion of a constant value, we call it rare-variant.
}
}
