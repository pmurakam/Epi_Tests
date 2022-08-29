# Epi_Tests
epi.tests.bayes(), epi.tests.robust(), and epi.tests.adjust() functions

Description:  
These 3 functions return point and interval estimates for sensitivity, specificity, and diagnostic 
accuracy.  Diagnostic accuracy is defined as the proportion of all tests that give a correct result.

- epi.tests.bayes estimates a posterior distribution for each value using a beta(a,b) prior combined 
  with the binomial likelihood. beta(1,1), the default, is the same as a uniform(0,1) prior.

- epi.tests.robust obtains estimates using Generalized Estimating Equations (GEE).  GEE uses the
  Huber-White Sandwich estimator of variance, which is robust to violation of the typical assumption 
  of uncorrelated data.  I.e., the confidence intervals reported from this function are valid in the 
  presence of correlated responses from cluster samples.  (In linear models the sandwich estimator is 
  also robust to heteroscedasticity.)

- epi.tests.adjust estimates sensitivity, specificity, and diagnostic accuracy adjusting for covariates.


Usage:
```
epi.tests.bayes(x, y, level=0.95, dig=3, method="hpd", a=1, b=1)

epi.tests.robust(x, y, id, level=0.95, dig=3, corstr="independence", tog=TRUE)

epi.tests.adjust(x, y, d, level=0.95, dig=3)
```

Arguments:
```
x      - logical (or 0/1) vector with predictions
y      - logical (or 0/1) vector with results
level  - confidence (or for epi.tests.bayes, "credible") interval level desired.
dig    - integer indicating the number of decimal places to be used in the output.
method - "hpd" for the obtaining the bayesian credible interval from the posterior distribution 
         using the highest posterior density interval.  Anything else will result in simply using 
         the (1-level)/2 and 1-(1-level)/2 quantiles.
a      - first shape parameter of the prior beta distribution to use.
b      - second shape parameter of the prior beta distribution to use.
id     - vector identifying the clusters
corstr - a character string specifying the correlation structure (see ?gee in the gee package)
tog    - if possible, estimate sensitivity and specificity in one model (TRUE) or don't (FALSE).
d      - data frame with columns for all the covariates to adjust for, and only those covariates.
         They will be adjusted for as is, so be sure that categorical variables are not of numeric 
         or integer class or they will be treated as continuous.
```

Notes:
All observations with any missing data are deleted (i.e., these functions perform complete-case analysis).
Also, the epi.tests.bayes is exact.  It is appropriate for either small or large samples.
epi.tests.robust and epi.tests.adjust, on the other hand, rely on asymptotic theory.  Thus they are not
appropriate for small samples.

See also:  
epi.tests in the epiR package.

Examples:
```
### Generate some fake data:
## This data is consistent with no predictive ability:
set.seed(834)
predictor = sample(c(FALSE,TRUE), 1000, replace=TRUE) # 0,1 data is acceptable also.
result = sample(c(FALSE,TRUE), 1000 , replace=TRUE)   # 0,1 data is acceptable also.

###epi.tests.bayes:
epi.tests.bayes(x=predictor, y = result, level=.9)

###epi.tests.robust:
id = rep(1:10,each=100) # variable identifying the clusters
epi.tests.robust(x=predictor, y = result, id=id, dig=3)

###epi.tests.adjust:
###Balanced design:
##Estimate is equal to the overall mean:
d2 = data.frame(iv = rep(letters[1:2], each=500))
epi.tests.adjust(x=predictor, y = rep(1,1000), d=d2) #see row for sensitivity
mean(predictor)

###Unbalanced design: 
##Each estimate is no longer the overall mean but the mean of the cell means:
d3 = data.frame(iv = rep(letters[1:2], times=c(300,700)))
epi.tests.adjust(x=predictor, y = rep(1,1000), d=d3) #see row for sensitivity
mean(c(mean(predictor[1:300]), mean(predictor[301:1000])))

##If continuous covariates are included, the estimates are the values estimated at the baseline 
##levels of those continuous covariates (here, when iv=0):
d1 = data.frame(iv = rnorm(1000))
epi.tests.adjust(x=predictor, y = result, d=d1)
##If you want estimates at different levels of the continuous covariates, you could modify the function 
##code to estimate linear combination of model parameters, or (easier) you could just re-center the 
##continuous covariate(s) that you give to the function.
```
