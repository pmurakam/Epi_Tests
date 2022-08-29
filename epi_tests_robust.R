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

epi.tests.robust <- function(x,y,id,level=0.95,dig=3,corstr="independence",tog=TRUE) {
    require(gee)
    ## Checking and preprocessing:
    if(!inherits(level,c("numeric","integer"))) stop("level must be number between 0 and 1.")
    if(level<=0 | level>=1) stop("level must be number between 0 and 1.") 
    if(!is.vector(x) | !is.vector(y)) stop("x and y must be logical vectors or a numeric vectors of zeros for FALSE and ones for TRUE.")
    if(!is.vector(id)) stop("id must be a vector identifying the clusters.")
    if(length(x)!=length(y))  stop("lengths of x and y differ.")
    if(length(x)!=length(id)) stop("lengths of x and y differ from length of id.")
    if(length(x)<5) stop("Too few observations to do inference") #even 4 is too few.
    if(inherits(x,c("character","factor"))) stop("x must be 0,1 or FALSE,TRUE.")
    if(inherits(y,c("character","factor"))) stop("y must be 0,1 or FALSE,TRUE.")
    rm = is.na(x) | is.na(y) | is.na(id)
    x  = x[!rm]
    y  = y[!rm]
    id = id[!rm]
    if(is.numeric(x) | is.integer(x)){ if(!all(x%in%c(0,1))) stop("x must be 0,1 or FALSE,TRUE.")}
    if(is.integer(y) | is.integer(y)){ if(!all(y%in%c(0,1))) stop("y must be 0,1 or FALSE,TRUE.")}
    y = as.logical(y)
    x = as.numeric(x)
    #mean of x[y] is the raw sensitivity estimate
    #mean of !x[!y] is the raw sensitivity estimate
    #mean of c(x[y],!x[!y]) is the raw sensitivity estimate

    ##Sensitivity and Specificity:
    mod = try(gee(x~ -1+y,id=id,family=binomial(link="logit"),corstr=corstr), silent=TRUE)
    if(!inherits(mod,"try-error") & tog) {
        mod = summary(mod)$coef
        ci = matrix(NA,2,2)
        ci[1,] = mod["yFALSE","Estimate"]+c(-1,1)*qnorm(1-(1-level)/2)*mod["yFALSE","Robust S.E."]
        ci[2,] = mod["yTRUE", "Estimate"]+c(-1,1)*qnorm(1-(1-level)/2)*mod["yTRUE", "Robust S.E."]
        est = cbind(ci[,1], mod[c("yFALSE","yTRUE"),"Estimate"], ci[,2])
        senspec = exp(est)/(1+exp(est))
        senspec["yFALSE",] = rev(1-senspec["yFALSE",])
        senspec = cbind(senspec,table(y)[c("FALSE","TRUE")])
        colnames(senspec) = c("lower","estimate","upper","n")
    } else {
        ##Specificity:
        if(sum(!y)<2) specif=c(rep(NA,3),sum(!y)) else {
            spec = !x[!y]
            spec_id = id[!y]
            mod = try(gee(spec~1,id=spec_id,family=binomial(link="logit"),
                          corstr=corstr), silent=TRUE)
            if(inherits(mod,"try-error")) specif=rep(NA,3) else {
                mod = summary(mod)$coef
                ci = mod["(Intercept)","Estimate"]+
                     c(-1,1)*qnorm(1-(1-level)/2)*mod["(Intercept)","Robust S.E."]
                est = c(ci[1], mod["(Intercept)","Estimate"], ci[2])
                specif = exp(est)/(1+exp(est))
            }
            specif = c(specif,length(spec))
        }
        names(specif) = c("lower","estimate","upper","n")
        
        ##Sensitivity:
        if(sum(y)<2) sensit=c(rep(NA,3),sum(y)) else {
            sens =  x[y]
            sens_id = id[y]
            mod = try(gee(sens~1,id=sens_id,family=binomial(link="logit"),
                          corstr=corstr), silent=TRUE)
            if(inherits(mod,"try-error")) sensit=rep(NA,3) else {  
                mod = summary(mod)$coef
                ci = mod["(Intercept)","Estimate"]+
                     c(-1,1)*qnorm(1-(1-level)/2)*mod["(Intercept)","Robust S.E."]
                est = c(ci[1], mod["(Intercept)","Estimate"], ci[2])
                sensit = exp(est)/(1+exp(est))
            }
            sensit = c(sensit,length(sens))
        }
        names(sensit) = c("lower","estimate","upper","n")

        senspec = rbind(specif,sensit)
    }

    ##Diagnostic accuracy:
    dacc = c(x[y],!x[!y])
    dacc_id = c(id[y],id[!y])
    mod = try(gee(dacc~1,id=dacc_id,family=binomial(link="logit"),corstr=corstr), silent=TRUE)
    if(inherits(mod,"try-error")) diagacc=rep(NA,3) else {   
        mod = summary(mod)$coef
        ci = mod["(Intercept)","Estimate"]+
             c(-1,1)*qnorm(1-(1-level)/2)*mod["(Intercept)","Robust S.E."]
        est = c(ci[1], mod["(Intercept)","Estimate"], ci[2])
        diagacc = exp(est)/(1+exp(est))
    }
    diagacc = c(diagacc,length(dacc))
    names(diagacc) = c("lower","estimate","upper","n")

    out = round(rbind(senspec,diagacc), dig)
    rownames(out) = c("specificity","sensitivity","diagnostic_accuracy")
    ##Put rows in same order as the other epi.tests functions:
    out[c("sensitivity","specificity","diagnostic_accuracy"),]
}
