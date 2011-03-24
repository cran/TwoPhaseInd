spmle <- function(data,					## dataset (data.frame)
                  response,				## response variable (binary)
                  treatment,			## treatment variable (binary)
                  BaselineMarker,		## environment variable (continuous)
                  extra=NULL,			## extra variable(s)
                  phase,				## variable for phase indicator
                  ind=TRUE,				## independent or non-indepentent
                  difffactor=0.000001,	## 
                  maxit=1000,			## max iteration
                  verbose=FALSE			## TRUE to turn on debug log
                  ){
  ## argument validation
  if(!is.data.frame(data)){
    stop("Argument data must be a data.frame object.")
  }else{
    colNames <- colnames(data)
    if(!(response %in% colNames)){
      stop("Response variable not found in the data.")
    }
    if(!(treatment %in% colNames)){
      stop("Treatment variable not found in the data.")
    }
    if(!(BaselineMarker %in% colNames)){
      stop("BaselineMarker variable not found in the data.")
    }
    if(!is.null(extra)){
      if(!any(extra %in% colNames)){
        extraNotFound <- paste(extra[!(extra %in% colNames)], sep="", collapse=", ")
        stop(paste("Extra variable(s) not found in the data:", extraNotFound))
      }
    }

    if(!(phase %in% colNames)){
      stop("Phase variable not found in the data.")
    }else{
      if(any(levels(factor(data[, phase])) != c("1", "2"))){
        stop("Phase variable must be either 1 or 2 only.")
      }
    }
  }
 

  ## phase 1 
  y <- data[, response]
  x <- data[, treatment]
  z <- data[, BaselineMarker]
  N1 <- sum(y==1)
  N0 <- sum(y==0)

  NXYcount <- c(sum(x==0 & y==0),
                sum(x==1 & y==0),
                sum(x==0 & y==1),
                sum(x==1 & y==1)
                )

  ## phase 2
  phase2 <- data[, phase]==2
  phase2Data <- data[phase2, ]
  
  y <- data[phase2, response]
  x <- data[phase2, treatment]
  z <- data[phase2, BaselineMarker]
  n <- length(y)
  n0 <- sum(y==0)
  n1 <- sum(y==1)

  wgt0 <- ifelse(y==1,N1/n1,N0/n0)

  ##if(!is.null(extra)){
  ##  glmData <- data.frame(cbind(y, x, z, x*z, data[phase2, extra]))
  ##}else{
  ##  glmData <- data.frame(cbind(y, x, z, x*z))
  ##}
  ##
  ##sfit <- glm(y~., data=glmData, family=binomial, weights=wgt0)

  if(!is.null(extra)){
    glmFormula <- as.formula(paste(response, "~",
                                   paste(## main model
                                         paste(treatment,
                                               BaselineMarker,
                                               paste(treatment, BaselineMarker, sep=" * "), sep=" + "
                                               ),
                                         ## extra variables
                                         paste(extra, collapse=" + "),
                                         sep=" + "
                                         )
                                   )
                             )
  }else{
    glmFormula <- as.formula(paste(response, "~",
                                   paste(## main model
                                         paste(treatment,
                                               BaselineMarker,
                                               paste(treatment, BaselineMarker, sep=" * "), sep=" + "
                                               ),
                                         sep=" + "
                                         )
                                   )
                             )
  }

  sfit <- glm(glmFormula, data=data[phase2, ], family=binomial, weights=wgt0)  
  
  beta <- sfit$coef
  varmat <- rep(0,length(beta)^2)

  if(!is.null(extra)){
    extracov <- 1
    ## the following 2 lines don't work if some variables are factoral
    ##nextracov <- length(extra)	
    ##covvec <- c(t(as.matrix(data[, extra])))
    nextracov <- length(beta) - 4 ## 4 factors

    ## re-recognize levels of factors
    phase2Extra <- data[phase2, extra]
    for(var in extra){
      if(is.factor(phase2Extra[, var])){
        phase2Extra[, var] <- factor(phase2Extra[, var])
      }
    }
    ## dim(phase2Extra)

    #### FIXME:  should we constrand on phase2, i.e., data=data[phase2, extra] or data=data[, extra]?
    extraDesMat <- model.matrix(as.formula(paste("~", paste(extra, collapse=" + "))), data=phase2Extra)
    extraDesMat <- extraDesMat[, -1] ## remove the "(Intercept)" column
    covvec <- c(t(extraDesMat))

    ##setdiff(colnames(extraDesMat), names(beta))
    ##setdiff(names(beta), colnames(extraDesMat))

    betaName <- c("(Intercept)", treatment, BaselineMarker, paste(treatment, BaselineMarker, sep=":"),
                   colnames(extraDesMat))


  }else{
    extracov <- 0
    covvec <- 0
    nextracov <- 0

    betaName <- c("(Intercept)", treatment, BaselineMarker, paste(treatment, BaselineMarker, sep=":"))

  }
  
  ## re-order beta to fit extraDesMat
  beta <- beta[betaName]
  

  ## decide which C function to call by independent or not
  if(ind){
    nrCode <- "profile_NR_ind"
  }else{
    nrCode <- "profile_NR_noind"
  }

  nrResult <- .C(nrCode,
                 as.integer(n),		## Number of subject
                 as.integer(1:n),	## Subject ID
                 as.numeric(y),		## Y = Response variable (binary)
                 as.numeric(x),		## X = gene? (binary)
                 as.numeric(z),		## Z = environment? (continuous)
                 as.integer(NXYcount),	## 
                 as.integer(extracov),
                 as.numeric(covvec),
                 as.integer(nextracov),
                 beta=as.numeric(beta),
                 varmat=as.numeric(varmat),
                 as.numeric(difffactor),
                 as.integer(maxit),
                 as.integer(sum(verbose))
                 )

  beta <- nrResult$beta
  ## FIXME: CovMat is not symmetric
  CovMat <- matrix(nrResult$varmat, length(beta), length(beta), byrow=TRUE)
  stder <- sqrt(diag(CovMat))
  pVal <- 2*(1-pnorm(abs(beta/stder)))

  tmpResult <- data.frame(beta=round(beta, 4), stder=round(stder, 4), pVal=round(pVal, 4))

  ## add more information to the row names
  betaName[2] <- paste(treatment, "(Treatment)")
  betaName[3] <- paste(BaselineMarker, "(BaselineMarker)")
  betaName[4] <- paste(treatment, BaselineMarker, sep=":")
  rownames(tmpResult) <- betaName
  
  ##FIXME
  return(tmpResult)

}
