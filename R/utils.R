#internal
devik<- function(y, mu, rho=1.5){
   y <- as.vector(y)
   ((y^(2-rho)-y*mu^(1-rho))/(1-rho) - (y^(2-rho)-mu^(2-rho))/(2-rho))*2
}

deviance.MSTweedie <- function(y, mu, w, rho = 1.5){
   if(!class(y)=='list')y<-list(y)
   if(!class(mu)=='list')mu<-list(mu)
   if(missing(w)) w <- lapply(y, function(yy) rep(1, length(unlist(yy))))
   w <- lapply(w, function(wk) wk/sum(wk))
   sum(unlist(w) * devik(unlist(y), unlist(mu), rho))
}

#internal
error.bars <- function(x, upper, lower, width = 0.02, ...) {
   xlim <- range(x)
   barw <- diff(xlim) * width
   segments(x, upper, x, lower, ...)
   segments(x - barw, upper, x + barw, upper, ...)
   segments(x - barw, lower, x + barw, lower, ...)
   range(upper, lower)
}

#internal
err <- function(n, info, maxit, pmax, dfmax){
   if(n==0) msg=""
   if(n>0)
   {#fatal error
      if(n==1) msg="All penalty factors are <= 0"
      if(n==2) msg="weights must be non-negative"
      if(n==3) msg=paste("standardization fails at predictor ", info, sep="")
      if(n==4) msg=paste("majoration is 0 at predictor ", info%%10000, " at lambda ", (info-info%%10000)/10000, sep="")
      n=1
      msg=paste("in cmd fortran code -",msg)
   }
   if(n<0)
   {#non fatal error
      if(n==-1)msg=paste("Convergence for ",info,"th lambda value not reached after maxit=",maxit," iterations; solutions for larger lambdas returned",sep="")
      if(n==-2)msg=paste("Number of nonzero coefficients along the path exceeds pmax=",pmax, " at ",info,"th lambda value; solutions for larger lambdas returned",sep="")
      if(n==-3)msg=paste("Number of nonzero coefficients along the path exceeds dfmax=",dfmax, " at ",info,"th lambda value; solutions for larger lambdas returned",sep="")
      if(n==-4)msg=paste("Strict descent violation at pass ",info,".",sep="")
      n=-1
      msg=paste("from cmd fortran code -",msg)
   }
   list(n=n,msg=msg)
}

#internal
nonzeroCoef <- function (beta, bystep = FALSE) {
   ### bystep = FALSE means which variables were ever nonzero
   ### bystep = TRUE means which variables are nonzero for each step
   nr=nrow(beta)
   if (nr == 1) {#degenerate case
      if (bystep)
         apply(beta, 2, function(x) if (abs(x) > 0)
            1
            else NULL)
      else {
         if (any(abs(beta) > 0))
            1
         else NULL
      }
   }
   else {
      beta=abs(beta)>0 # this is sparse
      which=seq(nr)
      ones=rep(1,ncol(beta))
      nz=as.vector((beta%*%ones)>0)
      which=which[nz]
      if (bystep) {
         if(length(which)>0){
            beta=as.matrix(beta[which,,drop=FALSE])
            nzel = function(x, which) if (any(x))
               which[x]
            else NULL
            which=apply(beta, 2, nzel, which)
            if(!is.list(which))which=data.frame(which)# apply can return a matrix!!
            which
         }
         else{
            dn=dimnames(beta)[[2]]
            which=vector("list",length(dn))
            names(which)=dn
            which
         }

      }
      else which
   }
}

#internal
coefnorm <- function(list,q=2){
   if(q == 'inf'){
      tmp <- sapply(list, function(listl){
         apply(listl, 1, function(row){
            max(abs(row))
         })
      })
      tmp
   }
   else{
      tmp <- sapply(list, function(listl){
         apply(listl, 1, function(row){
            sum(abs(row)^q)
         })
      })
      tmp^(1/q)
   }

}


