! --------------------------------------------------
SUBROUTINE f_mstweedie(ntasks, nobs, nobsmax, nvars, w, x, y, pf, iex, sr, kktstop, reg, alpha,&
            istd, dfmax, pmax, nlam, minlam, userlam, maxit, rho, eps, &
            gam, nalam, alam, beta0, beta, nbeta, idvars, npass, jerr, kkt, bnorm, M, stdesc, iter)
! --------------------------------------------------

! - VARIABLES DECLARATIONS -------------------------
    IMPLICIT NONE
    ! - PARAMETERS ---------------------------------
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: small=1.0D-12
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
    ! - INPUT VARIABLES ----------------------------
        INTEGER :: ntasks
        INTEGER :: nobs(ntasks)
        INTEGER :: nobsmax
        INTEGER :: nvars
        DOUBLE PRECISION :: w(ntasks, nobsmax)
        DOUBLE PRECISION :: x(ntasks, nobsmax, nvars)
        DOUBLE PRECISION :: y(ntasks, nobsmax)
        DOUBLE PRECISION :: pf(nvars)
        INTEGER :: iex(nvars)
        INTEGER :: sr
        INTEGER :: kktstop
        INTEGER :: reg
        DOUBLE PRECISION :: alpha
        INTEGER :: istd
        INTEGER :: dfmax
        INTEGER :: pmax
        INTEGER :: nlam
        DOUBLE PRECISION :: minlam
        DOUBLE PRECISION :: userlam(nlam)
        INTEGER :: maxit
        DOUBLE PRECISION :: rho
        DOUBLE PRECISION :: eps
    ! - OUTPUT VARIABLES ---------------------------
        DOUBLE PRECISION :: gam(0:nvars,nlam)
        INTEGER :: nalam
        DOUBLE PRECISION :: alam(nlam)
        DOUBLE PRECISION :: beta0(ntasks, nlam)
        DOUBLE PRECISION :: beta(ntasks, nvars, nlam)
        INTEGER :: nbeta(nlam)
        INTEGER :: idvars(pmax)
        INTEGER :: npass
        DOUBLE PRECISION :: kkt(nvars, ntasks, nalam)
        DOUBLE PRECISION :: bnorm(nvars, nalam)
        INTEGER :: M(nvars, ntasks, nalam)
        INTEGER :: jerr(2)
    ! - LOCAL VARIABLES ----------------------------
        ! - COUNTERS -------------------------------
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
        INTEGER :: l
        INTEGER :: g
        INTEGER :: me
        ! - FOR STANDARDIZATION --------------------
        DOUBLE PRECISION :: xmean(ntasks, nvars)
        DOUBLE PRECISION :: xnorm(ntasks, nvars)
        INTEGER :: inullnorm = 0
        DOUBLE PRECISION :: maj
        ! - BETAS, RESIDUALS, ETC. -----------------
        DOUBLE PRECISION :: b(ntasks,0:nvars)
        DOUBLE PRECISION :: oldbeta(ntasks, 0:nvars)
        DOUBLE PRECISION :: eta(ntasks, nobsmax)
        DOUBLE PRECISION :: oldeta(ntasks, nobsmax)
        DOUBLE PRECISION :: r1(nobsmax)
        DOUBLE PRECISION :: r2(nobsmax)
        ! - OTHERS ---------------------------------
        INTEGER :: mnl
        DOUBLE PRECISION :: al
        DOUBLE PRECISION :: alf
        DOUBLE PRECISION :: al0
        DOUBLE PRECISION :: tlam
        INTEGER :: isr(nvars)
        INTEGER :: jx
        DOUBLE PRECISION :: u(ntasks, nvars)
        DOUBLE PRECISION :: ga(nvars)
        DOUBLE PRECISION :: ht
        DOUBLE PRECISION :: norm
        DOUBLE PRECISION :: old_kkt_norm
        INTEGER :: stop_kkt
        DOUBLE PRECISION :: v(ntasks)
        ! - IRLS WORKING VARIABLES -----------------
        DOUBLE PRECISION :: wt(ntasks, nobsmax)
        DOUBLE PRECISION :: zt(ntasks, nobsmax)
        DOUBLE PRECISION :: obj
        DOUBLE PRECISION :: stdesc(maxit)
        INTEGER :: iter(maxit)
        ! - CMD ITERATION VARIABLES ----------------
        DOUBLE PRECISION :: dif
        DOUBLE PRECISION :: oldb(ntasks)
        DOUBLE PRECISION :: sj(ntasks)
        DOUBLE PRECISION :: proj(ntasks)
        DOUBLE PRECISION :: t
        DOUBLE PRECISION :: d(ntasks)
        INTEGER :: oidvars(nvars)
        INTEGER :: ni
        DOUBLE PRECISION :: maxB
! - END VARIABLE DECLARATIONS ----------------------

! - BEGIN PREPARATION ------------------------------
    jerr = 0
    ! - CHECK IF ONE POSITIVE PENALTY FACTOR -------
    IF(maxval(pf) <= 0.0D0) THEN
        jerr(1)= 1
        RETURN
    ENDIF
    ! - NEGATIVE PF ARE SET TO ZERO ----------------
    pf=max(0.0D0,pf)
    ! - CHECK IF POSITIVE WEIGHTS ------------------
    IF (minval(w) < 0) THEN
        jerr(1) = 2
        RETURN
    ENDIF
    ! - STANDARDIZE WEIGHTS ------------------------
    DO k=1,ntasks
       w(k, 1:nobs(k))=w(k, 1:nobs(k))/sum(w)
       !w(k, 1:nobs(k))=w(k, 1:nobs(k))/sum(w(k,1:nobs(k)))
    ENDDO
    ! - STANDARDIZE X WITH WEIGHTS -----------------
    DO k=1,ntasks
        DO j=1,nvars
            IF (iex(j) == 0) CYCLE
            IF(istd==1) THEN
                xmean(k,j) = sum(x(k, 1:nobs(k),j)*w(k, 1:nobs(k)))
                x(k, 1:nobs(k) ,j) = x(k, 1:nobs(k) ,j)-xmean(k,j)
                maj = sum(x(k, 1:nobs(k),j)*w(k, 1:nobs(k))*x(k, 1:nobs(k),j))
                IF(maj < 1.0D0-16) THEN
                    inullnorm = j
                    RETURN
                ENDIF
                xnorm(k,j) = sqrt(maj)
                IF(xnorm(k,j)>small) x(k, 1:nobs(k),j) = x(k, 1:nobs(k),j)/xnorm(k,j)
            ENDIF
        ENDDO
    ENDDO
    IF(inullnorm > 0) THEN
        jerr(1) = 3
        jerr(2) = inullnorm
        RETURN
    ENDIF
    ! - INITIALIZATIONS ----------------------------
        al = 0.0D0
        al0 = 0.0D0
        alf = 0.0D0
        mnl = Min (mnlam, nlam)
        eta = 0.0D0
        b = 0.0D0
        oldbeta = 0.0D0
        old_kkt_norm = 0.0D0
        stop_kkt = 0
        idvars = 0
        oidvars = 0
        npass = 0
        ni = 0
        gam = 0.0D0
        ga = 0.0D0
        kkt = 0.0D0
        bnorm = 0.0D0
        M = 0
        isr = 1
        obj = 0.0D0
        stdesc = 0.0D0
        iter = 0
        IF(minlam < 1.0D0) THEN
            minlam = Max(mfl, minlam)
        ENDIF
! - END PREPARATION --------------------------------
! - BEGIN LAMBDA LOOP ------------------------------
    DO l=1,nlam
    ! - FIND THE VALUE OF LAMBDA ---------------
        ! - store previous lambda -
        al0 = al
        ! - is user defined -
        IF(minlam>=1.0D0) THEN
            al=userlam(l)
        ! - not user defined -
        ELSE
            IF (l == 1) THEN
                al = big
                DO k=1,ntasks
                    ! - initial beta 0 -
                    b(k,0) = log(sum(w(k, 1:nobs(k)) * y(k, 1:nobs(k)))) - log(sum(w(k, 1:nobs(k))))
                    eta(k, 1:nobs(k)) = b(k,0)
                ENDDO
            ELSE IF(l == 2) THEN
            ! - compute lambda max -
                al0 = 0.0D0
                DO j = 1,nvars
                    IF (iex(j) == 0) CYCLE
                    ! -------------------------------------------------TBC - done
                    IF(reg == 0) THEN
                        ga(j) = sum(abs(u(:,j)))
                        IF(pf(j)>1.0D-16) THEN
                            al0 = max(al0, ga(j)/pf(j))
                        ENDIF
                    ELSE IF (reg == 2) THEN
                        ga(j) = sum(u(:,j)*u(:,j))
                        IF(pf(j)>1.0D-16) THEN
                            al0 = max(al0, sqrt(ga(j))/pf(j))
                        ENDIF
                    ENDIF
                ENDDO
                alf=minlam**(1.0D0/(nlam-1.0D0))
                al = al0
            ELSE IF(l > 2) THEN
                al = al * alf
            ENDIF
        ! - end not user defined -
        ENDIF
    ! - STRONG RULE SET IDENTIFICATION -
        IF (sr == 1) THEN
                ! - compute gradients -
            DO k=1,ntasks
                r1(1:nobs(k)) = w(k,1:nobs(k))*y(k,1:nobs(k))*exp((1.0D0-rho)*eta(k,1:nobs(k)))
                r2(1:nobs(k)) = w(k,1:nobs(k))*exp((2.0D0-rho)*eta(k,1:nobs(k)))
                DO j=1,nvars
                    IF (iex(j) == 0) CYCLE
                    u(k,j) = sum( (-r1(1:nobs(k))+r2(1:nobs(k))) * x(k,1:nobs(k),j) )
                ENDDO
            ENDDO
            ! - compute norms and check-
            tlam = 2*al-al0
            DO j=1,nvars
                IF (iex(j) == 0) CYCLE
                ! -------------------------------------------------TBC - done
                IF(reg == 0)THEN
                    norm = sum(abs(u(:,j)))
                ELSE IF (reg == 2) THEN
                    norm = sqrt(sum(u(:,j)*u(:,j)))
                ENDIF
                IF(norm < pf(j)*tlam) isr(j) = 0
            ENDDO
        ENDIF
        iter(npass + 1 ) = 1
    ! - STRONG RULE LOOP ----------------------------------
        DO
    ! - THEN, FULL IRLS-CMD UPDATE -------------------
        ! - BEGIN OUTER LOOP (IRLS) -----------------------
            DO
                ! - SETUP WLS PROBLEM -------------------------
                ! - save previous beta_ o -
                oldbeta(:,0)=b(:,0)
                ! - save previous beta -
                IF(ni>0) THEN
                    DO j=1,ni
                        g=idvars(j)
                        oldbeta(:,g)=b(:,g)
                    ENDDO
                ENDIF
                ! - compute working weights and responses -
                wt = 0.0D0
                zt = 0.0D0
                DO k=1,ntasks
                    r1(1:nobs(k)) = w(k,1:nobs(k))*y(k,1:nobs(k))*exp((1.0D0-rho)*eta(k,1:nobs(k)))
                    r2(1:nobs(k)) = w(k,1:nobs(k))*exp((2.0D0-rho)*eta(k,1:nobs(k)))
                    wt(k,1:nobs(k)) = -(1.0D0-rho)*r1(1:nobs(k))+(2.0D0-rho)*r2(1:nobs(k))
                    zt(k,1:nobs(k)) = eta(k,1:nobs(k)) + ( r1(1:nobs(k))-r2(1:nobs(k)) ) / wt(k,1:nobs(k))
                ENDDO

                ! - compute majorations -
                gam = 0.0D0
                DO j=1,nvars
                    IF (iex(j) == 0) CYCLE
                    IF (isr(j) == 0) CYCLE
                    DO k=1,ntasks
                        !ht = dot_product(x(k,1:nobs(k),j)*wt(k,1:nobs(k)),x(k,1:nobs(k),j))
                        ht = sum(x(k,1:nobs(k),j)*x(k,1:nobs(k),j)*wt(k,1:nobs(k)))
                        gam(j,l) = max(gam(j,l), ht) * (1.0D0 + 1.0D-6)
                    ENDDO
                    IF (gam(j,l) < small) THEN
                        jerr(1) = 4
                        jerr(2) = j+10000*l
                        RETURN
                    ENDIF
                ENDDO
                ! - obselete by exact solution -
                !gam(0,l)=maxval(sum(wt,2)) * (1.0D0 + 1.0D-6)
                !IF (gam(0,l) <  small) THEN
                !    jerr(1) = 4
                !    jerr(2) = 10000*l
                !    RETURN
                !ENDIF

                ! - BEGIN MIDDLE LOOP (CMD CYCLE) ---------
                IF (iter(npass + 1 ) == 0) iter(npass + 1 ) = 2
                DO
                    npass=npass+1
                    dif=0.0D0
                    oldeta = eta
                    ! - update beta once for all j -
                    DO j=1,nvars
                        IF (isr(j) == 0) CYCLE
                        IF ( l == 1 ) CYCLE
                        ! - store current solution -
                        oldb=b(:,j)
                        ! - compute update -
                        DO k=1,ntasks
                            sj(k) = dot_product(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-zt(k,1:nobs(k))),x(k,1:nobs(k),j))
                        ENDDO
                        sj = b(:,j) - sj/gam(j,l)
                        IF(alpha > 0.0D0) THEN
                            ! get sgn
                            proj = 0
                            DO k=1,ntasks
                                IF(sj(k)>0.0D0) proj(k)=1.0D0
                                IF(sj(k)<0.0D0) proj(k)=-1.0D0
                            ENDDO
                            sj = abs(sj) - alpha * al*pf(j)/gam(j,l)
                            DO k=1,ntasks
                               IF(sj(k) < 0.0D0) sj(k) = 0.0D0
                            ENDDO
                            sj = sj*proj
                        ENDIF
                        !Linf
                        IF(reg == 0)THEN
                            IF (sum(abs(sj)) <= al*pf(j)/gam(j,l)) THEN
                                b(:,j) = 0.0D0
                            ELSE
                                CALL ProjB1(sj ,ntasks ,al*pf(j)*(1.0D0-alpha)/gam(j,l) , proj)
                                b(:,j) = sj - proj
                            ENDIF
                        !L2
                        ELSE IF(reg == 2) THEN
                            norm=sqrt(dot_product(sj,sj))
                            t=norm-pf(j)*al*(1.0D0-alpha)/gam(j,l)
                            ! - check positivity and update -
                            IF(t>0.0D0) THEN
                                b(:,j)=sj*t/norm
                            ELSE
                                b(:,j)=0.0D0
                            ENDIF
                        ENDIF


                        ! - update difference -
                        d=b(:,j)-oldb
                        ! - if there was update -
                        IF(any(abs(d)>small)) THEN
                            ! - update overall update difference -
                            ! -------------------------------------------------TBC - done
                            IF(reg == 0 )THEN
                                dif=max(dif,maxval(abs(d)))
                            ELSE IF(reg == 2) THEN
                                dif=max(dif,sqrt(dot_product(d,d)))
                            ENDIF
                            ! - udate linear predictor -
                            DO k=1,ntasks
                                eta(k,1:nobs(k))=eta(k,1:nobs(k))+x(k,1:nobs(k),j)*d(k)
                            ENDDO
                            ! - update active set -
                            IF(oidvars(j)==0) THEN
                                ! - we include a new predictor -
                                ni=ni+1
                                IF(ni>pmax) EXIT
                                oidvars(j)=ni
                                idvars(ni)=j
                            ENDIF
                        ENDIF
                    ENDDO
                    ! - update beta0 once -
                    DO k=1,ntasks
                        d(k) = sum(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-zt(k,1:nobs(k)))) / sum(wt(k,1:nobs(k)))
                    ENDDO

                    IF(any(abs(d) > small)) THEN
                        b(:,0)=b(:,0)-d
                        DO k=1,ntasks
                            eta(k,1:nobs(k))=eta(k,1:nobs(k))-d(k)
                        ENDDO
                        ! -------------------------------------------------TBC - done
                        IF(reg == 0 )THEN
                            dif=max(dif,maxval(abs(d)))
                        ELSE IF(reg == 2) THEN
                            dif=max(dif,sqrt(dot_product(d,d)))
                        ENDIF

                        ! - descent difference update of penalty -
                        ! - disabled feature, i think the index j was the problem...
                        ! - let future simon deal with this
                        IF ( reg == 0 ) THEN
                            !norm =  maxval(abs(oldb)) - maxval(abs(b(:,j)))
                        ELSEIF ( reg == 2 ) THEN
                            !norm = sqrt(sum(oldb*oldb)) - sqrt(sum(b(:,j)*b(:,j)))
                        ENDIF
                        stdesc(npass) = stdesc(npass) + al * pf(j) * norm
                    ENDIF
                    ! - descent difference of objective fn -
                    DO k = 1,ntasks
                        norm = sum(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-oldeta(k,1:nobs(k)))&
                            *(2*zt(k,1:nobs(k))-eta(k,1:nobs(k))-oldeta(k,1:nobs(k))))
                        stdesc(npass) = stdesc(npass) + 0.5D0*norm
                    ENDDO
                    ! - first one is irrelevent
                    stdesc(npass) = 0.0D0
                    ! - check convergence -
                    IF (ni > pmax) EXIT
                    IF (dif < eps) EXIT
                    IF(npass > maxit) THEN
                        jerr(1) = -1
                        jerr(2) = l
                        RETURN
                    ENDIF
                    ! - check strict descent of objective function -
!                        obj = - sum(u(:,k) * (b(:,j) - oldb))
!                        obj = obj - 0.5D0 * gam(j,l) * sum((b(:,j) - oldb)*(b(:,j) - oldb))
!                        IF ( reg == 0 ) THEN
!                            norm = maxval(b(:,j)) - maxval(oldb)
!                        ELSEIF ( reg == 2 ) THEN
!                        norm = sqrt(sum(b(:,j)*b(:,j))) - sqrt(sum(oldb*oldb))
!                        ENDIF
!                        obj = obj + al*pf(j)*norm
!                        IF (obj < -small) stdesc(j) = stdesc(j)+1


                    ! - cycle on new active set until convergence -
                    ! - BEGIN INNER LOOP ------------------
                    DO
                        npass=npass+1
                        dif=0.0D0
                        oldeta = eta
                        ! - update beta once for all j -
                        DO g=1,ni
                            j=idvars(g)
                            IF (iex(j) == 0) CYCLE
                            IF (isr(j) == 0) CYCLE
                            ! - store current solution -
                            oldb=b(:,j)

                            ! - compute update -
                            DO k=1,ntasks
                                sj(k) = dot_product(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-zt(k,1:nobs(k))),x(k,1:nobs(k),j))
                            ENDDO
                            sj = b(:,j) - sj/gam(j,l)
                            IF(alpha > 0.0D0) THEN
                                ! get sgn
                                proj = 0
                                DO k=1,ntasks
                                    IF(sj(k)>0.0D0) proj(k)=1.0D0
                                    IF(sj(k)<0.0D0) proj(k)=-1.0D0
                                ENDDO
                                sj = abs(sj) - alpha * al*pf(j)/gam(j,l)
                                DO k=1,ntasks
                                   IF(sj(k) < 0.0D0) sj(k) = 0.0D0
                                ENDDO
                                sj = sj*proj
                            ENDIF
                            !Linf
                            IF(reg == 0)THEN
                                IF (sum(abs(sj)) <= al*pf(j)/gam(j,l)) THEN
                                    b(:,j) = 0.0D0
                                ELSE
                                    CALL ProjB1(sj ,ntasks ,al*pf(j)*(1.0D0-alpha)/gam(j,l) , proj)
                                    b(:,j) = sj - proj
                                ENDIF
                            !L2
                            ELSE IF(reg == 2) THEN
                                norm=sqrt(dot_product(sj,sj))
                                t=norm-pf(j)*al*(1.0D0-alpha)/gam(j,l)
                                ! - check positivity and update -
                                IF(t>0.0D0) THEN
                                    b(:,j)=sj*t/norm
                                ELSE
                                    b(:,j)=0.0D0
                                ENDIF
                            ENDIF

                            ! - update difference -
                            d=b(:,j)-oldb
                            ! - if there was update -
                            IF(any(abs(d)>small)) THEN
                                ! - update overall update difference -
                                ! -------------------------------------------------TBC - done
                                IF(reg == 0 )THEN
                                    dif=max(dif,maxval(abs(d)))
                                ELSE IF(reg == 2) THEN
                                    dif=max(dif,sqrt(dot_product(d,d)))
                                ENDIF
                                ! - udate linear predictor -
                                DO k=1,ntasks
                                    eta(k,1:nobs(k))=eta(k,1:nobs(k))+x(k,1:nobs(k),j)*d(k)
                                ENDDO
                                ! - dont update active set -
                                ! - check strict descent of objective function -
!                                    obj = - sum(u(:,k) * (b(:,j) - oldb))
!                                    obj = obj - 0.5D0 * gam(j,l) * sum((b(:,j) - oldb)*(b(:,j) - oldb))
!                                    IF ( reg == 0 ) THEN
!                                        norm = maxval(b(:,j)) - maxval(oldb)
!                                    ELSEIF ( reg == 2 ) THEN
!                                        norm = sqrt(sum(b(:,j)*b(:,j))) - sqrt(sum(oldb*oldb))
!                                    ENDIF
!                                    obj = obj + al*pf(j)*norm
!                                    IF (obj < -small) stdesc(j) = stdesc(j)+1
                            ENDIF
                        ENDDO
                        ! - update beta0 once -
                        DO k=1,ntasks
                            d(k) = sum(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-zt(k,1:nobs(k)))) / sum(wt(k,1:nobs(k)))
                        ENDDO

                        IF(any(abs(d) > small)) THEN
                            b(:,0)=b(:,0)-d
                            DO k=1,ntasks
                                eta(k,1:nobs(k))=eta(k,1:nobs(k))-d(k)
                            ENDDO
                            ! -------------------------------------------------TBC - done
                            IF(reg == 0 )THEN
                                dif=max(dif,maxval(abs(d)))
                            ELSE IF(reg == 2) THEN
                                dif=max(dif,sqrt(dot_product(d,d)))
                            ENDIF

                            ! - descent difference update of penalty -
                            IF ( reg == 0 ) THEN
                                norm = maxval(abs(oldb)) - maxval(abs(b(:,j)))
                            ELSEIF ( reg == 2 ) THEN
                                norm = sqrt(sum(oldb*oldb)) - sqrt(sum(b(:,j)*b(:,j)))
                            ENDIF
                            stdesc(npass) = stdesc(npass) + al * pf(j) * norm
                        ENDIF

                        ! - descent difference of objective fn -
                        DO k = 1,ntasks
                            norm = sum(wt(k,1:nobs(k))*(eta(k,1:nobs(k))-oldeta(k,1:nobs(k)))&
                                *(2.0D0*zt(k,1:nobs(k))-eta(k,1:nobs(k))-oldeta(k,1:nobs(k))))
                            stdesc(npass) = stdesc(npass) + 0.5D0*norm
                        ENDDO

                        ! - check convergence -
                        IF (dif < eps) EXIT
                        IF(npass > maxit) THEN
                            jerr(1)=-1
                            jerr(2)=l
                            RETURN
                        ENDIF
                    ! - END INNER LOOP
                    ENDDO
                ! - END MIDDLE LOOP (CMD CYCLE) -----------
                ENDDO
                ! - FINAL CHECK FOR WLS -----------------------
                IF(ni>pmax) EXIT
                jx = 0
                ! - check if b0 changed enough
                IF(reg == 0 )THEN
                    norm = maxval(abs(b(:,0)-oldbeta(:,0)))
                ELSE IF(reg == 2) THEN
                    norm = sqrt(sum((b(:,0)-oldbeta(:,0)) * (b(:,0)-oldbeta(:,0))))
                ENDIF

                IF(norm >= eps) jx =  1

                ! - check if any Bj changed enough -------
                DO j=1,nvars
                    IF (iex(j) == 0) CYCLE
                    IF (isr(j) == 0) CYCLE
                    ! -------------------------------------------------TBC - done
                    IF(reg == 0 )THEN
                        norm = maxval(abs(b(:,j)-oldbeta(:,j)))
                    ELSE IF(reg == 2) THEN
                        norm = sqrt(sum((b(:,j)-oldbeta(:,j)) * (b(:,j)-oldbeta(:,j))))
                    ENDIF

                    IF(norm >= eps) jx =  1
                ENDDO
                IF(jx == 1) CYCLE
                ! - otherwise we have convergence -
                EXIT
        ! - END OUTER LOOP (IRLS) -------------------------
            ENDDO
        ! - STRONG RULE CHECK -
            jx = 0
            ! - compute gradients -
            DO k=1,ntasks
                r1(1:nobs(k)) = w(k,1:nobs(k))*y(k,1:nobs(k))*exp((1.0D0-rho)*eta(k,1:nobs(k)))
                r2(1:nobs(k)) = w(k,1:nobs(k))*exp((2.0D0-rho)*eta(k,1:nobs(k)))
                DO j=1,nvars
                    IF (iex(j) == 0) CYCLE
                    u(k,j) = sum((-r1(1:nobs(k))+r2(1:nobs(k)))*x(k,1:nobs(k),j))
                ENDDO
            ENDDO
            ! - compute norms and check-

            DO j=1,nvars
                IF (iex(j) == 0) CYCLE
                IF(reg == 0 )THEN
                    norm = sum(abs(u(:,j)))
                ELSE IF(reg == 2) THEN
                    norm = sqrt(sum(u(:,j)*u(:,j)))
                ENDIF
                kkt(j ,:, l) =  u(:,j)
                IF (isr(j) == 1) CYCLE
                ! - check excluded only, M=0 -
                IF(norm <= al*pf(j)) THEN
                    ! - failed test -
                    jx = 1
                    isr(j) = 1
                ENDIF
            ENDDO
            !bug(:,l) = isr(:)
            ! - restart IRLS with new active set -
            IF(jx == 1) CYCLE
            ! - all strong rule ok for beta_j = 0 -
            EXIT
        ! - END STRONG RULE LOOP -------------------
        ENDDO
        ! - KKT CHECK ------------------------------
            norm = 0.0D0
            ! - zero
            DO j=1,nvars
                IF (iex(j) == 0) CYCLE
                IF (maxval(b(:,j)) >= small) CYCLE
                M(j,:,l) = 0

                IF(reg == 0 )THEN
                    !  the condition ||Uj||1 < lambda * v_j
                    norm = max(norm, sum(abs(u(:,j))) /pf(j) - al)
                ELSE IF(reg == 2) THEN
                    ! the condition ||Uj||2 < lambda * v_j
                    norm = max(norm, sqrt(sum(u(:,j)*u(:,j))) /pf(j) -al)
                ENDIF
            ENDDO

            ! - non-zero
            DO j=1,nvars
                IF (iex(j) == 0) CYCLE
                ! - check included only -
                IF (isr(j) == 0) CYCLE
                IF (maxval(b(:,j)) < small) CYCLE

                IF(reg == 0 )THEN
                    ! -------------------------------------------------TBC - done
                    ! - we identify the set of indices -
                    maxB = maxval(abs(b(:,j)))
                    DO k=1,ntasks
                        IF(abs(maxB - abs(b(k,j))) < small) THEN
                            M(j,k,l) = 1
                        ELSE
                            M(j,k,l) = -1
                            norm = max(norm, abs(u(k,j))/pf(j) )
                        ENDIF
                    ENDDO
                    !  the condition on ||Uj|| = lambda * v_j
                    norm = max(norm, abs( sum(abs(u(:,j))) /pf(j) - al ) )
                ELSE IF(reg == 2) THEN
                    M(j,:,l) = 1
                    ! the condition on ||Uj||=lambda * v_j
                    v =  al*pf(j) * b(:,j) / sqrt(sum(b(:,j)*b(:,j)) + u(:,j))
                    norm =  max(norm, sqrt(sum(v*v)))
                ENDIF
            ENDDO

            stop_kkt = 0
            IF( norm < eps *1.0D1 ) stop_kkt = 1
            IF( abs(old_kkt_norm-norm) >= eps/1.0D1 ) stop_kkt =0

            DO j=1,nvars
                kkt(j,:,l) = u(:,j)
            ENDDO
            old_kkt_norm = norm
        ! - END KKT CHECK --------------------------
        ! - after first update or irls converged or terminated -
        ! - if we have too many predictors in -
        IF(ni>pmax) THEN
            jerr(1) = -2
            jerr(2) = l
            EXIT
        ENDIF
        ! - we save estimates and some data for output -
        IF(ni>0) THEN
            DO g=1,ni
                j=idvars(g)
                IF (iex(j) == 0) CYCLE
                beta(:,j,l)=b(:,j)
            ENDDO
        ENDIF
        nbeta(l)=ni
        beta0(:,l)=b(:,0)
        alam(l)=al
        nalam=l
        ! - check if enough lambdas done -
        IF (l < mnl) CYCLE
        ! - count number of non-zero variables -
        me = 0
        DO g=1,ni
            j=idvars(g)
            IF (iex(j) == 0) CYCLE
            IF(any(abs(beta(:,j,l))>small)) me=me+1
        ENDDO
        IF(me>dfmax) THEN
            jerr(1) = -3
            jerr(2) = l
            EXIT
        ENDIF
        ! - KKT stop -
        IF(kktstop == 1 .AND. l> nlam /2.0D0)THEN
            IF (stop_kkt == 1) EXIT
        ENDIF
! - END LAMBDA LOOP
    ENDDO
! - DESTANDARDIZE -----------------------------------
    DO l=1,nalam
        DO k=1,ntasks
            DO j=1,nvars
                IF(iex(j)==0) CYCLE
                IF(reg == 2) THEN
                    bnorm(j,l) = sqrt(sum(beta(:,j,l)*beta(:,j,l)))
                ELSE IF(reg == 0) THEN
                    bnorm(j,l) = maxval(abs(beta(:,j,l)))
                ENDIF
                IF(istd==1 .AND. xnorm(k,j)>small) beta(k,j,l) = beta(k,j,l)/xnorm(k,j)
            ENDDO
            IF(istd==1) THEN
                beta0(k,l) = beta0(k,l) - dot_product(beta(k,:,l),xmean(k,:))
            ENDIF
        ENDDO
    ENDDO

    RETURN
END SUBROUTINE f_mstweedie


! --------------------------------------------------
SUBROUTINE ProjB1(v, p, z, w)
! --------------------------------------------------

! - VARIABLES DECLARATIONS -------------------------
    ! - INPUTS -
    INTEGER :: p
    DOUBLE PRECISION :: v(p)
    DOUBLE PRECISION :: z
    ! - OUTPUTS -
    DOUBLE PRECISION :: w(p)
    ! - LOCAL VARS -
    DOUBLE PRECISION :: vp(p)
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: s,ds
    DOUBLE PRECISION :: theta
    INTEGER :: U(p),G(p),L(p)
    INTEGER :: rho, drho
    INTEGER :: i,j,k,m,fl,cl
    INTEGER :: sgn(p)

! - END VARIABLE DECLARATIONS ----------------------

    ! - absolute components -
    vp = abs(v)
    ! - unit ball test -
    !z = 1.0D0
    ! - if already within the ball
    IF(sum(vp) <= z) THEN
        w = v
    ! - otherwise we have to project -
    ELSE
    ! - retrieve sgn -
        DO j=1,p
            IF(v(j) == 0.0D0) THEN
                sgn(j) = 0
            ELSEIF(v(j) > 0.0D0)THEN
                sgn(j) = 1
            ELSE
                sgn(j) = -1
            ENDIF
        ENDDO
    ! - initializations -
        U = 1
        s = 0.0D0
        rho = 0.0D0
    ! - while loop -
        DO
        ! - produce random index -
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(r)
            cl = CEILING(sum(U) * r)
            i = 0
            DO j=1,p
                IF (U(j) == 1) i = i + 1
                IF(i == cl) THEN
                    k = j
                    i = i+1
                    EXIT
                ENDIF
            ENDDO
        ! - partition -
            G = 0
            L = 0
            DO j=1,p
                IF(U(j) == 0) CYCLE
                IF(vp(j) >= vp(k)) THEN
                    G(j) = 1
                ELSE
                    L(j) = 1
                ENDIF
            ENDDO
        ! - update rule -
            drho = sum(G)
            ds = sum(vp * G)
            IF ((s+ds)-(rho+drho)*vp(k) < z ) THEN
                s = s + ds
                rho = rho + drho
                U = L
            ELSE
                U = G
                U(k) = 0
            ENDIF
        ! - while A non empty -
            IF (sum(U)==0) EXIT
        ENDDO
        ! - output -
        theta = (s-z) / rho
        DO j=1,p
            w(j) = max(vp(j) - theta, 0.0D0)*sgn(j)
        ENDDO
    ENDIF
    z=theta
END SUBROUTINE ProjB1


! --------------------------------------------------
SUBROUTINE ProjB1Mich(v, p, z, w)
! --------------------------------------------------

! - VARIABLES DECLARATIONS -------------------------
    ! - INPUTS -
    INTEGER :: p
    DOUBLE PRECISION :: v(p)
    DOUBLE PRECISION :: z
    ! - OUTPUTS -
    DOUBLE PRECISION :: w(p)
    ! - LOCAL VARS -
    DOUBLE PRECISION :: vp(p)
    DOUBLE PRECISION :: tau, rho
    INTEGER :: j,ch
    INTEGER :: sgn(p)
    INTEGER :: A(p)
! - END VARIABLE DECLARATIONS ----------------------

    ! - absolute components -
    vp = abs(v)
    ! - if already within the ball
    IF(sum(vp) <= z) THEN
        w = v
    ! - otherwise we have to project -
    ELSE
    ! - retrieve sgn -
        DO j=1,p
            IF(v(j) == 0.0D0) THEN
                sgn(j) = 0
            ELSEIF(v(j) > 0.0D0)THEN
                sgn(j) = 1
            ELSE
                sgn(j) = -1
            ENDIF
        ENDDO
    ! - cycle -
    A = 1
    rho = (sum(vp) - z) / sum(A)
    DO
        ch = 0
        DO j=1,p
            IF ( A(j) == 0 ) CYCLE
            IF (vp(j) <= rho) THEN
                A(j) = 0
                ch = 1
            ENDIF
        ENDDO
        rho = (sum(vp*A) - z) / sum(A)
        IF ( ch == 0 ) EXIT
    ENDDO
    tau = rho
    ! - produce projection on simplex and ajust sign -
        DO j=1,p
            w(j) = max(vp(j) - tau, 0.0D0)*sgn(j)
        ENDDO
    ENDIF
    z=tau
END SUBROUTINE ProjB1Mich


