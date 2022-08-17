removeAllMissingCatVar<-function(cat.varname, mxtrModel){

	################################################################################################
	#This function removes categorical variables with insufficient data to estimate effects 
	#from mixture models. This function must take mixture model of class formula, and will return 
	#a formula of class formula. 
	#
	#Args:
	#	cat.varname: The character string name of a categorical variable with insufficient data to 
	#				  estimate effects, which will be removed from the model 
	#	mxtrModel: A mixture model of class formula containing the variable specified in 
	#			   argument cat.varname
	#
	#Returns:
	#	A mixture model of class formula that no longer includes the variable in cat.varname
	####################################################################################################
	
	#convert mixture model from class formula to class Formula 
	mxtrModelF<-Formula(mxtrModel)

	#need to remove variables from both sides of the mixture model formula 
	#so begin by checking whether the variable in both sides 
	side1.vars<-attr(terms(mxtrModelF, lhs=0, rhs=1), "term.labels")
	side2.vars<-attr(terms(mxtrModelF, lhs=0, rhs=2), "term.labels")
	in.side1<-cat.varname %in% side1.vars
	in.side2<-cat.varname %in% side2.vars
	
	#now update the mixture model formula
	if (in.side1 & in.side2){
	
		new.modelF<-eval(substitute(update(mxtrModelF, ~.-cat.varname|.-cat.varname), list(cat.varname=as.name(cat.varname))))
	
	} else if (in.side1 & ! in.side2){
	
		new.modelF<-eval(substitute(update(mxtrModelF, ~.-cat.varname|.), list(cat.varname=as.name(cat.varname))))
	
	} else if (! in.side1 & in.side2){
	
		new.modelF<-eval(substitute(update(mxtrModelF, ~.|.-cat.varname), list(cat.varname=as.name(cat.varname))))
	
	} else {
	
		new.modelF<-mxtrModelF
	
	}
	
	#return new model of class formula 
	new.model<-formula(new.modelF)
	return(new.model)

}