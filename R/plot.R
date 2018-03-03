plot.MSTweedie=function(x, log.lambda = TRUE, type.coef=c("coef","norm"),
                        lambda.min=NA, lambda.1se=NA, ...){
   type.coef=match.arg(type.coef)
   if(type.coef=="coef"){
      ntasks <- ncol(x$beta[[1]])
      beta <- matrix(0, nrow(x$beta[[1]]), length(x$lambda))
      for( k in seq(ntasks)){
         for(l in seq(length(x$lambda))){
            beta[,l] <- x$beta[[l]][,k]
         }
         plotCoef(beta,
                  x$lambda, main = paste('Coef. source',k),
                  log.lambda = log.lambda ,
                  ylab= "",
                  ...)
         if(!is.na(lambda.min))abline(v = log(lambda.min), lty=1, col = 'gray')
         if(!is.na(lambda.1se)) abline(v = log(lambda.1se), lty=1)
      }
   }
   else {
      plotCoef(coefnorm(x$beta,switch(x$reg, 'L2' = 2, 'Linf' = 'inf')),x$lambda,
               log.lambda = log.lambda ,
               ylab="", main = 'Coef. norm',...)
      if(!is.na(lambda.min))abline(v = log(lambda.min), lty=1, col = 'gray')
      if(!is.na(lambda.1se)) abline(v = log(lambda.1se), lty=1)
   }
}



#internal
plotCoef=function(mat, lambda,
                  log.lambda = TRUE ,xlab = xlab, ylab="Coefficients",...){
   which <- nonzeroCoef(mat)
   nwhich <- length(which)
   p <- nrow(mat)
   switch(nwhich+1,#we add one to make switch work
          "0"={
             warning("No plot produced since all coefficients zero")
             return()
          },
          "1"=warning("1 or less nonzero coefficients; plot is not meaningful")
   )
   mat  <- as.matrix(mat[which,,drop=FALSE])

   if (log.lambda) {
      l <- log(lambda)
      xlab <- expression(log(lambda))
   } else {
      l <- lambda
      xlab <- expression(lambda)
   }

   plot.args <- list(x = l, y = 1:length(l), ylim = range(mat), xlab = xlab,
                     ylab = ylab, type = "n", xlim = range(l))
   new.args <- list(...)
   if (length(new.args)) {
      new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
      plot.args[names(new.plot.args)] <- new.plot.args
   }
   do.call("plot", plot.args)
   #line.args <- list(col = gray.colors(nwhich + 1, start = 0.05, end = 0.7, gamma = 2.2)[1:nwhich],
   #                  lwd = 1 + 1.2^(-p/20), lty = 1)
   line.args <- list(col = 1:nwhich,
                     lwd = 1 + 1.2^(-p/20), lty = 1)

   if (length(new.args))
      line.args[names(new.args)] <- new.args
   line.args$x <- l
   line.args$y <- t(mat)
   line.args$col <- rep(line.args$col, table(which))
   do.call("matlines", line.args)

   abline(h = 0, lwd = line.args$lwd)

}

plot.cv.MSTweedie=function(x,sign.lambda=1,...){
   cvobj=x
   xlab=expression(log(lambda))
   if(sign.lambda<0)xlab=paste("-",xlab,sep="")
   plot.args=list(x=sign.lambda*log(cvobj$lambda),
                  y=cvobj$cvm,ylim=range(cvobj$cvup,cvobj$cvlo),
                  xlab=xlab,ylab=cvobj$name,type="n")
   new.args=list(...)
   if(length(new.args))plot.args[names(new.args)]=new.args
   do.call("plot",plot.args)
   error.bars(sign.lambda*log(cvobj$lambda),cvobj$cvup,cvobj$cvlo,width=0.01,col="darkgrey")
   points(sign.lambda*log(cvobj$lambda),cvobj$cvm,pch=20,col="red")
   abline(v=sign.lambda*log(cvobj$lambda.min),lty=1, col = 'grey')
   abline(v=sign.lambda*log(cvobj$lambda.1se),lty=1)
   invisible()
}

