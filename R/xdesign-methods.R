  #################################################################################
  # xdesign defines the design matrix for the discrete portion of the model 
  #
  # Args:
  #   x: A matrix, data frame, or expression set containing the covariates for 
  #      the discrete portion of the model. 
  #   m: A Formula object of the form ~x1+x2...|z1+z2.... where x's
  #      represent the names of the covariates for the discrete portion of the model. 
  #
  # Returns:
  #   The design matrix for the discrete portion of the model. 
  ################################################################################# 
 
 #Set method for expression set inputs 
 setMethod("xdesign",signature("ExpressionSet","ANY"),function(x,m){
 			options(na.action=na.pass)
            xVars<-model.matrix(formula(m,lhs=0,rhs=1),data=pData(x))
            xVars})
            
  #Set method for matrix inputs          
  setMethod("xdesign",signature("matrix","ANY"),function(x,m){
  			options(na.action=na.pass)
            xVars<-model.matrix(formula(m,lhs=0,rhs=1),data=x)
            xVars})
            
  #Set method for data frame inputs            
  setMethod("xdesign",signature("data.frame","ANY"),function(x,m){
  			options(na.action=na.pass)
            xVars<-model.matrix(formula(m,lhs=0,rhs=1),data=x)
            xVars})
