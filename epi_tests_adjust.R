################################################################################
## Copyright (C) 2010 Peter Murakami <peter.murakami@gmail.com>
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

epi.tests.adjust <- function(x,y,d,level=0.95,dig=3) {
    ## Checking and preprocessing:
    if(!inherits(level,c("numeric","integer"))) stop("level must be number between 0 and 1.")
    if(level<=0 | level>=1) stop("level must be number between 0 and 1.")  
    if(!is.vector(x) | !is.vector(y)) stop("x and y must be logical vectors or a numeric vectors of zeros for FALSE and 1 for TRUE.")
    if(!is.data.frame(d)) stop("d must be a data frame with all covariates.")    
    if(length(x)!=length(y))  stop("lengths of x and y differ.")
    if(length(x)!=nrow(d)) stop("lengths of x and y differ from number of rows in d.")
    if(length(x)<5) stop("Too few observations to do inference") #even 4 is too few.
    if(inherits(x,c("character","factor"))) stop("x must be 0,1 or FALSE,TRUE.")
    if(inherits(y,c("character","factor"))) stop("y must be 0,1 or FALSE,TRUE.")
    rm = is.na(x) | is.na(y) | apply(d,1,function(x) any(is.na(x)))
    x  = x[!rm]
    y  = y[!rm]
    d = d[!rm,,drop=FALSE]    
    if(is.numeric(x) | is.integer(x)){ if(!all(x%in%c(0,1))) stop("x must be 0,1 or FALSE,TRUE.")}
    if(is.integer(y) | is.integer(y)){ if(!all(y%in%c(0,1))) stop("y must be 0,1 or FALSE,TRUE.")}
    y = as.logical(y)
    x = as.numeric(x)

    ## Calculations:
    sens =  x[y]              #mean of this is the unadjusted sensitivity estimate
    sens_d = d[y,,drop=FALSE]
    spec = !x[!y]             #mean of this is the unadjusted specificity estimate
    spec_d = d[!y,,drop=FALSE]
    dacc = c(x[y],!x[!y])     #mean of this is the unadjusted diagnostic accuracy estimate
    dacc_d = rbind(d[y,,drop=FALSE],d[!y,,drop=FALSE]) 

    op <- options()
    options(contrasts=c("contr.sum","contr.poly"), warn=2)
    ##Recall that, say you have model response~a1+b1+a1:b1, then with balanced data and sum-to-zero contrasts, the intercept is the overall mean of the response; the coefficient of a1 is the mean of the response in category a1 minus the overall mean; the coefficient of a1:b1 is the mean of the response in cell a1, b1 minus the overall mean and the coefficients of a1 and b1; etc. For unbalanced data (and balanced data) the intercept is the mean of the cell means; the coefficient of a1 is the mean of cell means at level a1 minus the intercept; etc.

    tab <- function(mod, level=level) {
        ci = try(confint(mod, level=level, parm="(Intercept)"), silent=TRUE)
        if(inherits(mod,"try-error")|inherits(ci,"try-error")) mytab = rep(NA,3) else {
            est = c(ci[1], coef(mod)["(Intercept)"], ci[2])
            mytab = exp(est)/(1+exp(est))
        }
        names(mytab) = c("lower","estimate","upper")
        mytab
    }
    
    ##Sensitivity:
    mod = try(glm(sens~.,data=sens_d,family=binomial(link="logit")), silent=TRUE)
    sensit = tab(mod,level=level)
    sensit = c(sensit,length(sens))
    names(sensit)[4] = "n"

    ##Specificity:
    mod = try(glm(spec~.,data=spec_d,family=binomial(link="logit")), silent=TRUE)
    specif = tab(mod,level=level)
    specif = c(specif,length(spec))
    names(specif)[4] = "n"    

    ##Diagnostic accuracy:
    mod = try(glm(dacc~.,data=dacc_d,family=binomial(link="logit")), silent=TRUE)
    diagacc = tab(mod,level=level)
    diagacc = c(diagacc,length(dacc))
    names(diagacc)[4] = "n"    

    options(op)
    
    out = round(rbind(sensit,specif,diagacc), dig)
    rownames(out) = c("sensitivity","specificity","diagnostic_accuracy")
    out
}
