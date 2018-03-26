# MSTweedie
Extending the **HDtweedie** package to multi-source data, this package implements **multi-source sparse Tweedie modelling**. Given non-negative data with excess zeroes (typically insurance claim data), the algorithm fits the Tweedie distribution while penalizing the coefficient. Flexible regularization is achieved by an adjustable penalty factor that can balance between-feature sparsity and between-sources sparsity within features.

## Installation
You can install the R package directly from GitHub using the **devtools** package function `install_github` using the following command :
    
    devtools::install_github("fontaine618/MSTweedie")
    
## Examples
Load package via 

    library(MSTweedie)

### AutoClaim dataset
The package contains a dataset of auto insurance claims, `AutoClaim`, which is a preprocessed version of the dataset of the same name from the **cplm** package. 

    data(AutoClaim)

The manual contains information on the dataset; you can display the first rows of it by

    head(AutoClaim)

The sources are defined by the `REVOLKED` variable (i.e. whether or not the policyholder had its license revoked in the past.) Here we classify the data by sources and by whether there was a claim or not.

    table(AutoClaim$REVOLKED, AutoClaim$CLM_AMT5 > 0)

### Fit the model and analyse
We fit the penalized Tweedie model with L1/Linf-regularization. The response is the variable `CLM_AMT5`, the claim amount.

    fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

We can extract the estimated coefficients at, say, the 34th regularization parameter value.

    coef.MSTweedie(fit, s=3)

We can predict the response using the coefficients from the 34th regularization parameter value.

    head(predict.MSTweedie(fit, s=34L)[[1]])

Print the regularization parameter value with the number of variables in the model along the regularization path.

    print.MSTweedie(fit)

We can plot the regularization pth of the norm of the coefficients

    plot.MSTweedie(fit, type.coef='norm')

or the coefficients themselves in each source.

    par(mfrow=c(2,1))
    plot.MSTweedie(fit, type.coef='coef')

We can check the convergence of the algorithm through the KKT conditions.

    par(mfrow=c(2,1))
    kkt.check(fit)

The deviance of the model can be computed (at the 34th regularization paramter for example) by

    deviance.MSTweedie(y= fit$y, mu = predict.MSTweedie(fit, s=34L))

### Cross-validation
We can perform 10-folds cross-validation.

    cv<-cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

Cross-validation is used to perform model selection. We can extract the estimated coefficients at le selected value of the regularization paramter according to the one-standard-error rule.

    coef.cv.MSTweedie(cv)

Also, we can use the selected coefficients to predict the responses.

    head(predict.cv.MSTweedie(cv)[[1]])

Finally, we can plot the CV model deviance with its standard error along the regularization path.

    plot.cv.MSTweedie(cv)
    
## Built with
This **R** package was built using the **devtools** package. The optimization routine is written in **FORTRAN90**.

## Author
Simon Fontaine

## License
GPL-2

## Acknowledgments
The code for the **R** end of the developpment was based in part on the **HDtweedie** package and the **glmnet** package. The **FORTRAN90** subroutine was inspired from that of the **HDtweedie** pakckage.
