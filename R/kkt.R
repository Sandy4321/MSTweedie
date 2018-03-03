
kkt.check <- function(fit, eps = fit$eps, cond = c(1,2),
                      from = 1, to = length(fit$lambda),...){
   ids <- seq(from, to)
   cond <- cond[order(cond)]
   reg <- fit$reg
   # 1st plot : nonzero conditions
   if(reg == 'Linf'){
      norms <- apply(fit$kkt,c(1,3) , function(x) sum(abs(x)))
      ylab <- expression(paste(group('|',group('|',U[j],'|'),'|')[1]/v[j]-lambda))
   }
   if(reg == 'L2'){
      norms <- apply(fit$kkt,c(1,3) , function(x) sqrt(sum(x^2)))
      ylab <- expression(paste(group('|',group('|',U[j],'|'),'|')[2]/v[j]-lambda))
   }
   kkt <- norms / fit$pf -
         matrix(rep(fit$lambda, nrow(norms)),
               nrow=nrow(norms), ncol=length(fit$lambda), byrow = T)

   if(match(1,cond,nomatch=0)>0){
      range <- range(kkt[abs(fit$M[,1,]) == 1])*2
      plot(NA, xlim = range(log(fit$lambda[ids])),
           ylim = range,
           ylab = ylab,
           xlab = expression(log(lambda)),
           main = "KKT condition on non-zero coefficients"
      )
      abline(h=0)
      for(j in seq(nrow(kkt))){
         tmp <- kkt[j,ids]
         tmp[abs(fit$M[j,1,ids]) != 1] <- NA
         lines(y = tmp, x =log(fit$lambda[ids]), col=j)
      }
   }
   # 2nd plot : zero condition
   if(reg == 'Linf'){
      ylab <- expression( paste(group('|',group('|',U[j],'|'),'|')[1]/v[j]-lambda))
   }
   if(reg == 'L2'){
      ylab <- expression( paste(group('|',group('|',U[j],'|'),'|')[2]/v[j]-lambda))
   }
   kkt <- kkt[,ids]
   if(match(2,cond,nomatch=0)>0){
      max <- max(kkt[abs(fit$M[,1,ids]) != 1])
      if(max<=0) max <- 1
      range <- c(-3*max, max)
      plot(NA, xlim = range(log(fit$lambda[ids])), ylim = range,
           ylab = ylab,
           xlab = expression(log(lambda)),
           main = "KKT condition on null coefficients"
      )
      abline(h=0)
      for(j in seq(nrow(kkt))){
         #tmp <- kkt[j,ids]
         tmp <- kkt[j,]
         tmp[abs(fit$M[j,1,ids]) == 1] <- NA
         lines(y = tmp, x =log(fit$lambda[ids]), col=j)
      }
   }
   # 3rd plot : detailed non-zero conditions
   pts <- t(c(log(fit$lambda[1]),0))
   #kkt3 <- rep(0,length(fit$lambda))
   if(reg == 'Linf'){
      for(l in seq(length(fit$lambda))){
         kktM <- fit$kkt[,,l]
         kktM <- as.vector(kktM[fit$M[,,l] == -1])
         #kkt3[l] <- ifelse(length(kktM)>0,max(abs(kktM)),-Inf)
         pts <- rbind(pts, cbind(rep(log(fit$lambda[l]), length(kktM)),kktM))
      }
      ylab <- expression(U[j]^(k))
   }
   if(reg == 'L2'){
      for(l in seq(length(fit$lambda))){
         kktN <- fit$kkt[,,l] + fit$pf * fit$lambda[l] * fit$beta[[l]] /
            apply(fit$beta[[l]], 1, norm) %*% t(rep(1, ncol(fit$beta[[l]])))
         kktN <- kktN[fit$M[,,l] == 1]
         #kkt3[l] <- ifelse(length(kktN)>0,max(abs(kktN)),-Inf)
         pts <- rbind(pts, cbind(rep(log(fit$lambda[l]), length(kktN)), kktN))
      }
      ylab <- expression(paste(lambda, v[j], beta^(k)/
                                  group('|',group('|',beta[j],'|'),'|')[2]+U^(k)[j],
                               ' (non-zero)'))
   }

   if(match(3,cond,nomatch=0)>0){
      plot(NA, xlim = range(log(fit$lambda)), ylim = range(pts[,2]),
           ylab = ylab,
           xlab = expression(log(lambda)),
           main = "KKT condition on non-zero coefficients (by task)")
      abline(h=0)
      abline(h=eps,lty=3)
      abline(h=-eps,lty=3)
      points(x = pts[,1], y= pts[,2])
   }

   # 4th plot : full condition
   ylab <- ''
   #full
   kkt4 <- sapply(seq(length(fit$lambda)), function(l){
      -sum(sapply(seq(fit$dim[1]), function(j){
         sum(fit$kkt[j,,l] * fit$beta[[l]][j,])
      }))/
         sum(fit$pf * fit$norm[,l])
   }) - fit$lambda


   if(match(4,cond,nomatch=0)>0){
      m <- max(abs(kkt4[!is.na(kkt4)]))
      range <- c(-m,m)
      plot(NA, xlim = range(log(fit$lambda)), ylim =range,
           ylab = ylab,
           xlab = expression(log(lambda)),
           main = "KKT condition : aggregate"
      )

      abline(h = c(-eps,0,eps), lty = c(3,1,3))
      #lines(x=log(fit$lambda), y=fit$lambda)
      lines(x = log(fit$lambda),
            y = kkt4, lty=1 , col=2, cex = 2)
      for(j in seq(fit$dim[1])){
         kkt5 <- sapply(seq(length(fit$lambda)), function(l){
            ifelse(fit$norm[j,l] == 0, NA,
                   - sum(fit$kkt[j,,l] * fit$beta[[l]][j,])/
                      (fit$pf[j] * fit$norm[j,l])
            )}) - fit$lambda
         lines(x = log(fit$lambda),
               y = kkt5, lty = 2, col = grey.colors(j))
      }

   }
}

