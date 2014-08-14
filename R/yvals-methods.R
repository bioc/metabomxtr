  #################################################################################
  # yvals defines a data frame of response variables 
  #
  # Args:
  #   y: A matrix, data frame, or expression set containing the response variables
  #   n: A character vector of the names of the response variables 
  #
  # Returns:
  #   A data frame of response variables, with columns representing different
  #  variables. 
  ################################################################################# 
  
  #Set method for Expression Sets
  setMethod("yvals",signature("ExpressionSet","character"),function(y,n){
            yvals<-as.data.frame(t(exprs(y)[n,,drop=FALSE]))
            yvals })
  #Set method for matrices          
  setMethod("yvals",signature("matrix","character"),function(y,n){
            yvals<-as.data.frame(y[,n,drop=FALSE])})
            
  #Set method for data frames 
  setMethod("yvals",signature("data.frame","character"),function(y,n){
            yvals<-as.data.frame(y[,n,drop=FALSE])})
