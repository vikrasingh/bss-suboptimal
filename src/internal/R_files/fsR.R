fsR <- function(x,y,k){
  include_packages <- c("ggplot2", "glmnet", "Matrix","bestsubset");
  lapply(include_packages, library, character.only = TRUE)
  start.time=proc.time()
  fs.obj = fs(x,y,k,intercept=FALSE,verbose=FALSE)  # maxsteps=tmax 
  end.time=proc.time()
  fs.obj$cputime=as.double(end.time[2]-start.time[2]) 
  return(fs.obj)
}

 
