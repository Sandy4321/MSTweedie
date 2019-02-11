MSTweedie<-function(x,y,w,source,rho=1.5,
                    nlambda=100,lambda.min,lambda,
                    x.normalize=T,eps, sr=T, kktstop=F,
                    reg = c('L2','Linf'),alpha = 0,
                    dfmax=nvars+1,pmax=min(dfmax*1.2,nvars),
                    pf=rep(1,nvars),maxit=1e6)
{
   this.call=match.call()
   t0 <- proc.time()
   #################################################################################
   ### INPUT SETUP #################################################################
   #################################################################################

   ##################
   ## transform data into lists for convenience

   # first check if already in right format (e.g. called from cv.MSTweedie)
   if(is.list(x) && is.list(y) && is.list(w) &&
      !is.data.frame(x) && !is.data.frame(y) && !is.data.frame(w)){
      # all lists, check formats
      nvars <- ncol(x[[1]])
      nobs <- sapply(x, nrow)
      ntasks <- length(x)
      if(!ntasks == length(y)) stop("x and y must have the same number of tasks (list check)")
      if(!ntasks == length(w)) stop("x and w must have the same number of tasks (list check)")
   }else{
   # not in target format, do transformations
      x <- as.data.frame(x)
      nvars <- ncol(x)
      if(missing(w)){
         w<-nvars+1
         x[,w]<-1
      }
      # case 1 : X = (y,src,w,x,...,x) not simultaneous
      if(!missing(source)){
         #should be only one index
         if(length(source)>1)stop("only one source index is accepted, use multiple y instead")
         #should have only one y
         idy <- y
         if(length(idy)>1)stop("source is supplied, only one y index is accepted")
         SRC <- as.factor(x[,source])
         sources <- levels(as.factor(SRC))
         # extract y
         y <- lapply(sources, function(src) subset(x, SRC == src , idy ))
         nobs <- sapply(y, nrow)
         ntasks <- length(y)
         for(k in seq(ntasks)) names(y[[k]]) <- paste(names(y[[k]]),sources[k], sep='.')
         # extract w
            if(length(w)>1)stop("only one w index is accepted")
            idw<-w
            w <- lapply(sources, function(src) subset(x, SRC == src , select = w ))
            x <- subset(x, select =-idw )
         # extract x, drop y nidex, w index and source index
         x <- lapply(sources, function(src) subset(x, SRC == src , -c(idy,source) ))
         nvars <- ncol(x[[1]])
      }
      # case 2 : X = (y,...,y,w,x,...,x)  simultaneous
      if(missing(source)){
         #source is disregarded in this case
         idy <- y
         # extract y
         y <- lapply(idy, function(src) subset(x, T, select = src ))
         nobs <- sapply(y, nrow)
         ntasks <- length(y)
         # extract w
            if(length(w)>1)stop("only one w index is accepted")
            idw <- w
            w <- lapply(idy, function(src) subset(x, T, select = w ))
            x <- subset(x, select =-idw )
         # extract x, drop y nidex, w index and source index
         x <- lapply(idy, function(src) subset(x, T, select =-c(idy) ))
         nvars <- ncol(x[[1]])
      }
   }

   ##treatment of x
   x=lapply(x, function(xm) apply(as.matrix.noquote(xm),2,as.numeric))
   nobs=sapply(x,nrow)
   nobsmax=max(nobs)
   xnames=colnames(x[[1]])
   if(is.null(xnames))xnames=paste("X",seq(nvars),sep="")

   ##treatment of y
   y=lapply(y, as.matrix)
   ynames <- sapply(y,names)
   if(is.null(ynames))ynames=paste("Y",seq(ntasks),sep="")

   ##treatment of w
   w=lapply(w, as.matrix)

   ##coerce x, y and w to arrays of doubles
   xmat <- array(0, dim = c(ntasks, nobsmax, nvars))
   ymat <- array(0, dim = c(ntasks, nobsmax))
   wmat <- array(0, dim = c(ntasks, nobsmax))
   for(k in seq(ntasks)){
      xmat[k,1:nobs[k],]<-as.double(x[[k]])
      ymat[k,1:nobs[k]]<-as.double(y[[k]])
      wmat[k,1:nobs[k]]<-as.double(w[[k]])
   }

   ##treatment of rho
   if(!is.numeric(rho)) stop("Parameter rho must be numeric")
   if(length(rho)>1) stop('cannot specify more than one rho value')
   if(rho <=1 | rho>=2) stop("Parameter rho must be between 1 and 2")

   ##treatment of lambda
   nlambda=as.integer(nlambda)
   if(missing(lambda))
   {
      if(missing(lambda.min)){
         lambda.min <- 1e-3 / switch(reg, 'L2' = sqrt(ntasks), 'Linf' = ntasks)
      }
      if(lambda.min>=1) stop("lambda.min should be less than 1")
      if(lambda.min<1.0E-6) stop("lambda.min is too small")
      lambda.min=as.double(lambda.min)
      userlam=double(1) #ulam=0 if lambda is missing
   }
   else
   {
      #flmin=1 if user define lambda
      lambda.min=as.double(1)
      if(any(lambda<0)) stop("lambdas should be non-negative")
      userlam=as.double(rev(sort(lambda))) #lambda is declining
      nlambda=as.integer(length(lambda))
   }

   ##treatment of normalization
   istd=as.integer(x.normalize)
   if(!istd == 0 & !istd ==1)stop('x.normalize must be either 0, 1 or boolean')

   ##treatment of reg
   if(missing(reg)) reg <- 'L2'
   reg = match.arg(reg)
   if(reg == 'L2') norm <- 2
   else if(reg == 'Linf') norm <- 0

   ##treatment of alpha
   if(alpha < 0 || alpha >1) stop('alpha must be between 0 and 1')

   ##treatment of eps
   if(missing(eps)){
      if(norm == 2) eps <- 1e-3
      if(norm == 0) eps <- 1e-3
   }
   eps=as.double(eps)

   if(! eps > 0 | is.na(eps))stop('eps bust be a positive number')

   ##treatment of sr
   sr=as.integer(sr)
   if(!sr == 0 & !sr ==1)stop('sr must be either 0, 1 or boolean')

   ##treatment of kktstop
   kktstop=as.integer(kktstop)
   if(!kktstop == 0 & !kktstop ==1)stop('kktstop must be either 0, 1 or boolean')

   ##treatment of dfmax
   dfmax=as.integer(dfmax)
   if(! dfmax > 0 | is.na(dfmax))stop('dfmax bust be positive integer')

   ##treatment of pmax
   pmax=as.integer(pmax)
   if(! pmax > 0 | is.na(pmax))stop('pmax bust be positive integer')

   ##treatment of exclude
   iex <- rep(1, nvars)
   # if(!missing(exclude))
   # {
   #    if(any(match(exclude,seq(nvars),0) == 0)) stop("Some excluded variables out of range; must be variable number.")
   #    iex[exclude] <- 0
   # }

   ##treatment of pf
   if(length(pf)==1) pf <- rep(pf,nvars)
   if(length(pf)!= nvars) stop("The size of penalty factor must be same as the number of input variables")
   if(any(pf<0)) stop('pf must be non-negative')
   pf=as.double(pf)

   ##treatment of maxit
   maxit=as.integer(maxit)
   if(! maxit > 0 | is.na(maxit))stop('maxit bust be positive integer')

   #################################################################################
   ### FORTRAN CALL ################################################################
   #################################################################################
   fit = .Fortran("f_mstweedie", PACKAGE="MSTweedie",
                  ntasks = ntasks,
                  nobs = nobs,
                  nobsmax = nobsmax,
                  nvars = nvars,
                  w = wmat,
                  x = xmat,
                  y = ymat,
                  pf = pf,
                  iex = as.integer(iex),
                  sr = sr,
                  kktstop = kktstop,
                  reg = as.integer(norm),
                  alpha = as.double(alpha),
                  istd = istd,
                  dfmax = dfmax,
                  pmax = pmax,
                  nlam = nlambda,
                  minlam = lambda.min,
                  userlam = userlam,
                  maxit = maxit,
                  rho = as.double(rho),
                  eps = eps,
                  gam = double((nvars+1)*nlambda),
                  nalam = integer(1),
                  alam = double(nlambda),
                  beta0 = double(ntasks*nlambda),
                  beta = double(ntasks*nvars*nlambda),
                  nbeta = integer(nlambda),
                  idvars = integer(pmax),
                  npass = integer(1),
                  jerr = integer(2),
                  kkt = double(nvars*nlambda*ntasks),
                  bnorm = double(nvars*nlambda),
                  M = integer(nvars*ntasks*nlambda),
                  stdesc = double(maxit),
                  iter = integer(maxit)

   )
   #################################################################################
   ### OUTPUT ######################################################################
   #################################################################################
   nalam=fit$nalam
   nbeta=fit$nbeta[seq(nalam)]
   nbetamax=max(nbeta)
   lam=fit$alam[seq(nalam)]

   if(missing(lambda)) {##first lambda is infinity; changed to entry point
      llam=log(lam)
      lam[1]=exp(2*llam[2]-llam[3])
   }

   stepnames=paste("s",seq(nalam),sep="")

   errmsg=err(fit$jerr[1], fit$jerr[2], maxit,pmax, dfmax)### error messages from fortran

   switch(paste(errmsg$n),
          "1"=stop(errmsg$msg,call.=FALSE),
          "-1"=warning(errmsg$msg,call.=FALSE)
   )

   dd=c(nvars,ntasks)
   df=rep(0,nalam)

   beta <- lapply(seq(nalam), function(l){
      matrix(fit$beta[seq(from = (nvars*ntasks)*(l-1)+1, to = (nvars*ntasks)*l, by = 1)],
                   nvars, ntasks, byrow= TRUE, dimnames = list(xnames,ynames))
   })

   kkt <- array(fit$kkt, c(nvars, ntasks, nlambda))[,,seq(nalam)]
   norm <- matrix(fit$bnorm, nvars, nlambda)[,seq(nalam)]
   #norm <- coefnorm(beta, ifelse(fit$reg == 0, 'inf',2))[,seq(nalam)]
   M <- array(fit$M, dim=c(nvars, ntasks, nlambda))[,,seq(nalam)]
   df <- apply(norm, 2, function(col){
      sum(col>0)
   })
   #gam <- matrix(fit$gam, nvars+1, nlambda)[,seq(nalam)]
   b0=fit$beta0[seq(ntasks*nalam)]
   b0=matrix(b0,ntasks,nalam,dimnames=list(ynames,stepnames))
   outlist=list(beta0=b0,beta=beta,df=df,
                lambda=lam,npasses=fit$npass, idvars = fit$idvars,
                #jerr=fit$jerr,
                dim=dd,call=this.call, pf = fit$pf, eps=eps,
                #gam = gam,
                #iex = fit$iex,
                kkt = kkt, norm = norm, reg = reg, alpha =alpha,
                y=y, x=x, w=w, rho=rho, M=M, time = proc.time() - t0
                #stdesc = fit$stdesc[seq(fit$npass)],
                #iter = fit$iter[seq(fit$npass)]
                )
   class(outlist) <- c("MSTweedie","bmd")
   outlist
}

coef.MSTweedie <- function(fit, s = NULL, tasks) {
   if(missing(s)) s<-as.integer(seq(length(fit$lambda)))
   if(missing(tasks)) tasks = seq(fit$dim[2])
   if(max(tasks) > length(fit$beta)) stop('tasks out of range')
   if(is.integer(s))i<-s
   else i<-match(s,fit$lambda)
   outlist <- lapply(tasks, function(k){
      b0 <- t(as.matrix(fit$beta0[k,]))
      rownames(b0) <- "(Intercept)"
      beta <- sapply(fit$beta, function(bl) bl[,k])
      nbeta <- rbind2(b0, beta)
      nbeta[,s]
   })
   names(outlist) <- tasks
   outlist
}

predict.MSTweedie <- function(fit, newx, tasks,  s,
                              type = c("response", "link")) {
   if(is.integer(s))l<-s
   else l<-match(s,fit$lambda)
   if(missing(tasks))tasks <- seq(fit$dim[2])
   if(missing(newx))newx <- lapply(tasks,function(k) fit$x[[k]])
   if(max(tasks) > length(fit$beta)) stop('tasks out of range')
   type = match.arg(type)
   lapply(seq(length(tasks)), function(i){
      k <- tasks[i]
      b0 <- fit$beta0[k,l]
      names(b0) <- "(Intercept)"
      beta <- fit$beta[[l]][,k]
      nbeta <- t(c(b0, beta))
      xx0 <- cbind2(1, newx[[i]])
      xx0 <- as.matrix(noquote(xx0))
      nfit <- xx0 %*% t(nbeta)
      switch(type, link = nfit, response = exp(nfit))
   })
}


print.MSTweedie <- function(x, digits = max(3, getOption("digits") -
                                               3)) {
   cat("\nCall: ", deparse(x$call), "\n\n")
   n <- min(length(x$df), length(x$lambda))
   print(cbind(Df = x$df[seq(n)], Lambda = signif(x$lambda, digits)[seq(n)]))
}
