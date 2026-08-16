abessUbxR <- function(xMatrix, yArray, tmax1){
    
    Packages <- c("abess")  # packages to be loaded
    lapply(Packages, library, character.only = TRUE)  # load all the packages
    start.time=proc.time()    # to find the cputime
    abess_fit =abess(xMatrix,yArray,family=c("gaussian"),tune.path="sequence",weight=NULL,normalize=0,support.size=tmax1,fit.intercept=FALSE)
    end.time=proc.time()
    ctime=as.double(end.time[2]-start.time[2])  # add cputime as a new field
    beta=coef(abess_fit)
    list(out1 = beta, out2 = ctime)
  }