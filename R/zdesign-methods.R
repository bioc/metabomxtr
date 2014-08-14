  #################################################################################
  # zdesign defines the design matrix for the continuous portion of the model 
  #
  # Args:
  #   x: A matrix, data frame, or expression set containing the covariates for 
  #      the continuous portion of the model. 
  #   m: A Formula object of the form ~x1+x2...|z1+z2.... where z's
  #      represent the names of the covariates for the continuous portion of the model. 
  #
  # Returns:
  #   The design matrix for the continuous portion of the model. 
  #################################################################################    
   
  #Set method for expression set inputs 
  setMethod("zdesign",signature("ExpressionSet","ANY"),function(x,m){
  			options(na.action=na.pass)
            zVars<-model.matrix(formula(m,lhs=0,rhs=2),data=pData(x))
            zVars})
            
  #Set method for matrix inputs  
  setMethod("zdesign",signature("matrix","ANY"),function(x,m){
  			options(na.action=na.pass)
            zVars<-model.matrix(formula(m,lhs=0,rhs=2),data=x)
            zVars}) 
            
  #Set method for data frame inputs                 
  setMethod("zdesign",signature("data.frame","ANY"),function(x,m){
  			options(na.action=na.pass)
            zVars<-model.matrix(formula(m,lhs=0,rhs=2),data=x)
            zVars})