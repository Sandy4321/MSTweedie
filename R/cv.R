cv.MSTweedie <- function(x, y, w, source, rho, nlambda = 100, lambda, lambda.min = 1e-3,
                         nfolds = 10, kktstop=F, foldid, adaptive = 0, x.normalize = TRUE,
                         reg = c('L2','Linf'),eps = 5e-4, sr = TRUE,
                         maxit = 1e6, pf = rep(1,nvars), alpha = 0, ...) {
   t0 <- proc.time()
   if (missing(rho)) rho <- 1.5
   pred.loss <- "deviance"
   ##adaptive lasso
      if(adaptive > 0){
         cv.fit <- cv.MSTweedie(x=x, y=y, w=w, rho=rho, nlambda = 10, lambda.min = 1e-3,
                                x.normalize=x.normalize, eps = eps, sr = sr, kktstop = FALSE,
                                reg = reg, dfmax = nvars+1, pmax = nvars+1,
                                pf = rep(1,nvars), maxit = maxit, adaptive = 0, alpha = alpha)
         # min
         id.min <- match(cv.fit$lambda.min, cv.fit$lambda)
         norm <- cv.fit$MSTweedie.fit$norm[,id.min]
         #new penalty factors
         pf <- pmax(norm, .Machine$double.eps)^(-adaptive)
         adaptive <- 0
      }
   # not in target format, do transformations
      x <- as.data.frame(x)
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
         if(!missing(w)){
            if(length(w)>1)stop("only one w index is accepted")
            idw<-w
            w <- lapply(sources, function(src) subset(x, SRC == src , select = w ))
            x <- subset(x, select =-idw )
         }
         if(missing(w)) w <- lapply(seq(ntasks), function(k) rep(1,nobs[k] ))
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
         if(!missing(w)){
            if(length(w)>1)stop("only one w index is accepted")
            idw <- w
            w <- lapply(idy, function(src) subset(x, T, select = w ))
            x <- subset(x, select =-idw )
         }
         if(missing(w)) w <- lapply(seq(ntasks), function(k) rep(1,nobs[k]))
         # extract x, drop y nidex, w index and source index
         x <- lapply(idy, function(src) subset(x, T, select =-c(idy) ))
         nvars <- ncol(x[[1]])
      }

   ##treatment of x
   nobs=sapply(x,nrow)
   nobsmax=max(nobs)
   xnames=colnames(x[[1]])
   if(is.null(xnames))xnames=paste("X",seq(nvars),sep="")

   ##treatment of y
   ynames <- sapply(y,names)
   if(is.null(ynames))ynames=paste("Y",seq(ntasks),sep="")

   ##treatment of reg
      if(missing(reg)) reg <- 'L2'
      reg = match.arg(reg)

   ##treatment of alpha
      if(alpha < 0 || alpha >1) stop('alpha must be between 0 and 1')

   ##treatment of lambda
      nlambda=as.integer(nlambda)
      if(missing(lambda))
      {
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

   ## folds preparation
      if(!missing(foldid)){
         ##user supplied folds
         nfolds <- max(foldid[[1]])
         if(class(foldid) == 'vector') foldid <- list(foldid)
         for(k in seq(ntasks)){
            if(length(foldid[[k]] != nobs(k))) stop('foldid does not have the same number of observations as the data')
         }
      }else{
         if(nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")
         ##we create folds
         foldid <- lapply(seq(ntasks), function(k){
            ids <- sample(seq(nobs[k]))
            folds <- cut(ids,breaks=nfolds,labels=FALSE)
         })
      }

   ###Fit the model once
      MSTweedie.obj <- MSTweedie(x=x , y=y, w=w, rho=rho, lambda = lambda,
                                 nlambda = nlambda, lambda.min = lambda.min,
                                 kktstop = kktstop, reg=reg, pf=pf, alpha=alpha, maxit=maxit,...)
      lambda <- MSTweedie.obj$lambda

      outlist <- as.list(seq(nfolds))
   ###Now fit the nfold models and store them
      dev <- matrix(NA, nrow=nlambda, ncol = nfolds)
      for (i in seq(nfolds)) {
         ## subset
         x.train <- list(seq(ntasks))
         y.train <- list(seq(ntasks))
         w.train <- list(seq(ntasks))
         x.test <- list(seq(ntasks))
         y.test <- list(seq(ntasks))
         w.test <- list(seq(ntasks))
         for(k in seq(ntasks)){
            wh <- foldid[[k]] == i
            x.train[[k]] <- x[[k]][!wh,]
            y.train[[k]] <- y[[k]][!wh,1]
            w.train[[k]] <- w[[k]][!wh]
            x.test[[k]] <- x[[k]][wh,]
            y.test[[k]] <- y[[k]][wh,1]
            w.test[[k]] <- w[[k]][wh]
         }
         ## fit the model on train data
         fit <- MSTweedie(x = x.train, y = y.train, w=w.train,
                                   rho = rho, lambda = lambda, kktstop = FALSE,
                                   reg=reg, pf=pf,alpha=alpha, maxit=maxit,...)
         ## get prediction on test data
         pred <- lapply(seq(length(fit$lambda)), function(l){
            predict.MSTweedie(fit, x.test, seq(ntasks),l, type = 'response')
            })
         ## compute deviance
         dev[,i] <- sapply(pred, function(pl) deviance.MSTweedie(y.test, pl, w.test, rho))
      }
      ## compute mean an se
      cvm <- apply(dev, 1, mean)
      cvsd <- apply(dev, 1, sd)/sqrt(nfolds)
      out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, reg=reg, alpha=alpha,
                  cvlo = cvm - cvsd, name = "Tweedie Deviance", MSTweedie.fit = MSTweedie.obj,
                  time = proc.time()-t0)
      ## find lambda min and 1se
         cvmin <- min(cvm, na.rm=TRUE)
         idmin <- cvm <= cvmin
         lambda.min <- max(lambda[idmin], na.rm=TRUE)
         idmin <- match(lambda.min, lambda)
         semin <- (cvm + cvsd)[idmin]
         idmin <- cvm <= semin
         lambda.1se <- max(lambda[idmin], na.rm=TRUE)
         lamin <- list(lambda.min = lambda.min, lambda.1se = lambda.1se)
      ## output
      obj <- c(out, as.list(lamin))
      class(obj) <- "cv.MSTweedie"
      obj
}




##wrappers
coef.cv.MSTweedie <- function(cv.MSTweedie.obj, s = c("lambda.1se", "lambda.min")) {
   if (missing(s))s<-"lambda.1se"
   if (is.numeric(s)){
      l <- match(s,cv.MSTweedie.obj$lambda)
   }else if(is.integer(s)){
      l<-s
   }else if (is.character(s)) {
         s <- match.arg(s)
         l <- match(cv.MSTweedie.obj[[s]],cv.MSTweedie.obj$lambda)
   }else stop("Invalid form for s")
   coef.MSTweedie(cv.MSTweedie.obj$MSTweedie.fit, s = l )
}
predict.cv.MSTweedie <- function(object, newx, s = c("lambda.1se",
                                                     "lambda.min")) {

   if (missing(s))s<-"lambda.1se"
   if (is.numeric(s))
      lambda <- s else if (is.character(s)) {
         s <- match.arg(s)
         lambda <- object[[s]]
      } else stop("Invalid form for s")
   predict.MSTweedie(object$MSTweedie.fit, newx, s = lambda)
}
