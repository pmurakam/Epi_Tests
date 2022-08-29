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

epi.tests.bayes <- function(x,y,level=0.95,dig=3, method="hpd", a=1, b=1) { #,plotDensities=TRUE
    ## Checking and preprocessing:
    if(!inherits(level,c("numeric","integer"))) stop("level must be number between 0 and 1.")
    if(level<=0 | level>=1) stop("level must be number between 0 and 1.")
    if(!is.vector(x) | !is.vector(y)) stop("x and y must be logical vectors or a numeric vectors of zeros for FALSE and 1 for TRUE.")  
    if(length(x)!=length(y))  stop("lengths of x and y differ.")
    if(inherits(x,c("character","factor"))) stop("x must be 0,1 or FALSE,TRUE.")
    if(inherits(y,c("character","factor"))) stop("y must be 0,1 or FALSE,TRUE.")
    rm = is.na(x) | is.na(y)
    x  = x[!rm]
    y  = y[!rm]    
    if(is.numeric(x) | is.integer(x)){ if(!all(x%in%c(0,1))) stop("x must be 0,1 or FALSE,TRUE.")}
    if(is.integer(y) | is.integer(y)){ if(!all(y%in%c(0,1))) stop("y must be 0,1 or FALSE,TRUE.")}
    y = as.logical(y)
    x = as.numeric(x)

    ## Calculations:
    sens =  x[y]           #mean of this is the sensitivity estimate
    spec = !x[!y]          #mean of this is the specificity estimate
    dacc = c(x[y],!x[!y])  #mean of this is the diagnostic accuracy estimate

    #if(plotDensities) {
    #    par(mfrow=c(1,3), mar=c(5,4,4,.6))
    #    plotDensities(v=sens,a=a,b=b,s=s,f=f,main="sensitivity",ylab="Density")
    #    plotDensities(v=spec,a=a,b=b,s=s,f=f,main="specificity",ylab="")
    #    plotDensities(v=dacc,a=a,b=b,s=s,f=f,main="diagnostic accuracy",ylab="")        
    #}

    if(method=="hpd") {
        require(TeachingDemos)
        sensit = hpd(qbeta, shape1=sum(sens)+1, shape2=sum(!sens)+1, conf=level)
        specif = hpd(qbeta, shape1=sum(spec)+1, shape2=sum(!spec)+1, conf=level)
        diagac = hpd(qbeta, shape1=sum(dacc)+1, shape2=sum(!dacc)+1, conf=level)        
    } else {
        sensit = qbeta(p=c((1-level)/2,1-(1-level)/2), shape1=sum(sens)+1, shape2=sum(!sens)+1)
        specif = qbeta(p=c((1-level)/2,1-(1-level)/2), shape1=sum(spec)+1, shape2=sum(!spec)+1)
        diagac = qbeta(p=c((1-level)/2,1-(1-level)/2), shape1=sum(dacc)+1, shape2=sum(!dacc)+1)
    }

    sensit = c(sensit[1],mean(sens),sensit[2],length(sens))
    specif = c(specif[1],mean(spec),specif[2],length(spec))
    diagac = c(diagac[1],mean(dacc),diagac[2],length(dacc))
    out = round(rbind(sensit,specif,diagac), dig)
    colnames(out) = c("lower","estimate","upper","n")
    rownames(out) = c("sensitivity","specificity","diagnostic_accuracy")
    out
}

#plotDensities <- function(v,a,b,s,f,...) {
#    p = seq(0,1,length=600)
#    s = sum(v)
#    f = sum(!v)
#    prior = dbeta(p,a,b)
#    like  = dbeta(p,s+1,f+1)
#    post  = dbeta(p,a+s,b+f)
#    plot(p,like,type="l",lty=1,lwd=3,...)
#    lines(p,prior,lty=3,lwd=3)
#    lines(p,post,lty=2,lwd=3,col="red")
#    rug(mean(v))
#    legend("topright",c("Prior","Likelihood","Posterior"), lty=c(3,1,2), lwd=c(3,3,3), col=c("black","black","red"))
#}
