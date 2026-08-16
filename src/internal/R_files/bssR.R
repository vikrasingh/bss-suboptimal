bssR <- function(x,y,k,cpulimit){
  Packages <- c("ggplot2", "glmnet", "Matrix","bestsubset")  # packages to be loaded
  lapply(Packages, library, character.only = TRUE)  # load all the packages
  start.time=proc.time()    # to find the cputime
  result_bss = bs(x,y,intercept=FALSE,k,time.limit=cpulimit,verbose=TRUE)
  end.time=proc.time()
  result_bss$cputime=as.double(end.time[2]-start.time[2])  # add cputime as a new field
  return(result_bss)
}
