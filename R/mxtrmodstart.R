mxtrmodstart <-
function(obsY,xVars,zVars,includeDiscrete){
  #######################################################################################
  # Specifies the starting values to use for optimization of the mxtrmodLL function.
  #
  # Args:   
  #   obsY:   A vector containing the response variable, which must be normally or log 
  #           normally distributed. If lognormal, log transformed Y's should be input 
  #           as the response variable. Missing Y values should be indicated by NA. 
  #   xVars:  The design matrix for the coavariates included in the discrete portion 
  #           of the model.   
  #   zVars:  The design matrix for the covariates inlcuded in the continuous portion 
  #           of the model. 
  # includeDiscrete: A logical indicator for whether or not to include the discrete 
  #                  portion of the model.   
  # Returns:
  #   A vector containing the starting values for each parameter in the 
  #   mixture model function, to be used as starting points when optimizing 
  #   the parameter estimates.   
  #########################################################################################
 
  ####################################################################################
  # If includeDiscrete is false, specify starting values for the model excluding the
  # discrete component
  ####################################################################################  

  if(!includeDiscrete){

  #Determine the number of covariates for the continuous portion of the model 
  nContinuous <- dim(zVars)[2]

  #Use the mean of the observed responses as the starting value for the intercept term 
  beta0 <- mean(obsY,na.rm=TRUE)

  #Use 0 as the starting value for the remaining beta parameters. 
  if (nContinuous>1){
    betaStart <- c(beta0,rep(0,nContinuous-1))
  } else {
  betaStart <- beta0
  }

  #Use the standard deviation of observed responses as the starting value for the variance 
  #parameter
  sigmaStart <- sd(obsY,na.rm=TRUE)

  #Combine starting values into one vector
  startVals <- c(betaStart,sigmaStart)
  }

  ####################################################################################
  # If includeDiscrete is true, specify starting values for the model including both
  # discrete and continuous portions
  ####################################################################################

  if(includeDiscrete){

  #Determine number of covariates for the categorical portion of the model
  nDiscrete <- dim(xVars)[2]

  #Determine number of covariates for the continuous portion of the model 
  nContinuous <- dim(zVars)[2]
  
  #Use the log odds of having observed a response as the starting value for 
  #the intercept term of the discrete component of the model. 
  p0 <- sum(!is.na(obsY))/length(obsY)
  alpha0 <- log(p0/(1-p0))

  #Use 0 as the starting value for the remaining alpha parameters. 
  if (nDiscrete>1){
    alphaStart <- c(alpha0,rep(0,nDiscrete-1))
  } else {
  alphaStart <- alpha0
  }  
  #Use the mean of the observed responses as the starting value for the intercept term
  #of the continuous portion 
  beta0 <- mean(obsY,na.rm=TRUE)

  #Use 0 as the starting value for the remaining beta parameters. 
  if (nContinuous>1){
    betaStart <- c(beta0,rep(0,nContinuous-1))
  } else {
  betaStart <- beta0
  }

  #Use the standard deviation of observed responses as the starting value for the variance 
  #parameter
  sigmaStart <- sd(obsY,na.rm=TRUE)

  #Combine starting values into one vector
  startVals <- c(alphaStart,betaStart,sigmaStart)
  }
  startVals
}
