mxtrmodLL <-
function(params,obsY,xVars,zVars,Tvals,includeDiscrete){
  #################################################################################
  # Specifies the mixture distribution of logistic and log normal variables. 
  #
  # Args:
  #   params: A vector of parameter estimates for all paramters in the mixture model.  
  #   obsY:   A vector containing the response variable, which must be normally or log 
  #           normally distributed. If lognormal, log transformed Y's should be input 
  #           as the response variable. Missing Y values should be indicated by NA. 
  #   xVars:  The design matrix for the coavariates included in the discrete portion 
  #           of the model.   
  #   zVars:  The design matrix for the covariates inlcuded in the continuous portion 
  #           of the model. 
  #   Tvals:  A vector of thresholds below which continuous variables are not 
  #           observable.               
  # includeDiscrete: A logical indicator for whether or not to include the discrete 
  #                  portion of the model.  
  # Returns:
  #   The negative log-likelihood of the specified mixture model. 
  ################################################################################# 

  ##########################################################################################
  #If includeDiscrete is false, exclude any missing values and use only the continuous
  #portion of the model likelihood. 
  ##########################################################################################

  if(!includeDiscrete){

  #Determine the number of continuous covariates 
  nContinuous <- as.numeric(dim(zVars)[2])
  
  #Specify the parameter estimates for the model 
  betaVec <- params[1:nContinuous]
  sigma <- params[(nContinuous+1)]

  #Write the log likelihood 
  obsInd <- 1*(!is.na(obsY))
  ZprimeBeta <- as.vector(zVars %*% betaVec)

  theseObs <- which(obsInd==1)
  
  #suppressing warnings for the dnorm and pnorm functions 
  #because these give warnings about NANs produced, which 
  #are expected and can be ignored 
  ans1 <- suppressWarnings(sum(log(dnorm(obsY[theseObs],mean=ZprimeBeta[theseObs],
                         sd=sigma))))

  #Take the negative of the log likelihood because the optimization 
  #function (optimx) finds a minimum. 
  ans <- -ans1
  }
  
  #######################################################################
  #If includeDiscrete is true, include both discrete and continuous
  #portions of the mdoel likelihood. 
  ####################################################################### 

  if(includeDiscrete){

  #Determine number of covariates for the discrete portion of the model
  nDiscrete <- as.numeric(dim(xVars)[2])
  
  #Determine the number of covariates for the continuous portion of the model
  nContinuous <- as.numeric(dim(zVars)[2])
  
  #Specify parameters for the discrete portion
  alphaVec <- params[1:nDiscrete]

  #Specify parameters for the continuous portion
  betaVec <- params[(nDiscrete+1):(nDiscrete+nContinuous)]

  #Specify variance paramter
  sigma <- params[(nDiscrete+nContinuous+1)]
 
  #Write the log likelihood
  #obsInd is an indicator variable, equal to 1 if the response variable
  #was observed, and 0 if censored  
  obsInd <- 1*(!is.na(obsY))
  XprimeAlpha <- as.vector(xVars %*% alphaVec)
  ZprimeBeta <- as.vector(zVars %*% betaVec)

  #part 1
  part1 <- sum(-log(1+exp(XprimeAlpha)))

  #part 2
  cdfPart <- suppressWarnings(pnorm(Tvals,mean=ZprimeBeta,sd=sigma))
  part2 <- suppressWarnings(sum((1-obsInd) * log(1+exp(XprimeAlpha)*cdfPart)))
  
  #part 3
  part3 <- sum(obsInd*XprimeAlpha)

  #part 4
  theseObs <- which(obsInd==1)
  part4 <- suppressWarnings(sum(log(dnorm(obsY[theseObs],mean=ZprimeBeta[theseObs],
                         sd=sigma))))
  
  ans1 <- part1+part2+part3+part4

  #Take the negative of the log likelihood because the optimization 
  #function (optimx) finds a minimum.
  ans <- -ans1
 }

  ans
}
