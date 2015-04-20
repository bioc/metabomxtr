
  #Generic function for creating data frame of response variables 
  setGeneric("yvals",function(y,n){standardGeneric("yvals")})
  
  #Generic function for creating the design matrix for the discrete portion of the model
  setGeneric("xdesign",function(x,m){standardGeneric("xdesign")} )

  #Generic function for creating the design matrix for the continuous portion of the model
  setGeneric("zdesign",function(x,m){standardGeneric("zdesign")} )
  
  