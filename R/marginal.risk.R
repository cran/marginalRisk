# only pass ph2 data to these functions
marginal.risk=function(fit.risk, fit.s, data, categorical.s, weights=rep(1, nrow(data)), t=NULL, ss=NULL, verbose=FALSE) {
    if(categorical.s) {
        marginal.risk.cat  (fit.risk, fit.s, data, weights=weights, t=t, verbose=verbose) 
    } else {
        marginal.risk.cont (fit.risk, fit.s, data, weights=weights, t=t, ss=ss, verbose=verbose) 
    }
}

# categorical markers
marginal.risk.cat=function(fit.risk, fit.s, data, weights=rep(1, nrow(data)), t=NULL, verbose=FALSE) {  
    
    if("coxph" %in% class(fit.risk)) {
        time.var=as.character(fit.risk$terms[[2]][[2]])
        y.var=as.character(fit.risk$terms[[2]][[3]])
    }
    
    marker.name= as.character(fit.s$terms[[2]])
    ss=unique(data[[marker.name]]); ss=sort(ss[!is.na(ss)])
    probs=predict(fit.s, newdata=data, type="probs")
    stopifnot(as.character(ss)==colnames(probs))
    
    if(any(is.na(probs))) stop("NA's found in fit.s")

    if (!"coxph" %in% class(fit.risk)) {
        # logistic regression
        dat.tmp.mrc=data
        risks=sapply(ss, function(s) {
            f.s.zi = probs[,s]
            dat.tmp.mrc[[marker.name]]=s    
            risks = predict(fit.risk, newdata=dat.tmp.mrc, type="response") # glm
            sum(weights * risks  * f.s.zi) / sum(weights*f.s.zi)    
        })
        names(risks)=levels(ss)
        risks        
            
    } else {
        # coxph or svycoxph
        if (is.null(t)) {
            # return risk versus time
            tt=sort(unique(data[[time.var]][data[[y.var]]==1]))        
            risks=sapply(tt, function (t) {
                dat.tmp.mrc=data
                dat.tmp.mrc[[time.var]]=t
                risks=sapply(ss, function(s) {        
                    f.s.zi = probs[,s]
                    dat.tmp.mrc[[marker.name]]=s    
                    risks = 1 - exp(-predict(fit.risk, newdata=dat.tmp.mrc, type="expected"))# coxph survival prob
                    sum(weights * risks  * f.s.zi) / sum(weights*f.s.zi)    
                })
            })
            risks=t(risks)
            colnames(risks)=as.character(ss)        
            list(time=tt, risk=risks)
            
        } else {
            if (verbose) print("return risk at time t")
            dat.tmp.mrc=data
            time.var=fit.risk$terms[[2]][[2]]
            dat.tmp.mrc[[time.var]]=t        
            risks=sapply(ss, function(s) {
                f.s.zi = probs[,s]
                dat.tmp.mrc[[marker.name]]=s    
                risks = 1 - exp(-predict(fit.risk, newdata=dat.tmp.mrc, type="expected")) # coxph survival prob
                sum(weights * risks  * f.s.zi) / sum(weights*f.s.zi)    
            })
            names(risks)=levels(ss)
            risks        
        }
    }
}

# continuous markers
marginal.risk.cont=function(fit.risk, fit.s, data, weights=rep(1, nrow(data)), t=NULL, ss=NULL, verbose=FALSE) {
    marker.name= as.character(fit.s$terms[[2]])
    ss.is.null=is.null(ss) 
    if (ss.is.null) ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01))
        
    dat.tmp.mri=data
    if (!is.null(t)) {
        time.var=fit.risk$terms[[2]][[2]]
        dat.tmp.mri[[time.var]]=t
    }
    
    risks=sapply(ss, function(s) {
        if (class(fit.s)[1]=="lm") {
            f.s.zi = dnorm(s, mean=predict(fit.s, newdata=dat.tmp.mri), sd=summary(fit.s)$sigma)
        } else if (class(fit.s)[1]=="svyglm") {
            f.s.zi = dnorm(s, mean=predict(fit.s, newdata=dat.tmp.mri), sd=sqrt(summary(fit.s)$dispersion))
        } else stop(paste0("this class of fit.s is not supported: ", class(fit.s)[1]))
        
        if(any(is.na(f.s.zi))) stop("NA's found in fit.s")
        
        dat.tmp.mri[[marker.name]]=s    
        risks = if(is.null(t)) {
            # glm
            predict(fit.risk, newdata=dat.tmp.mri, type="response")
        } else {
            # coxph survival prob
            1 - exp(-predict(fit.risk, newdata=dat.tmp.mri, type="expected"))
        }
        #if(any(is.na(risks))) stop("NA's found in fit.risk")
        
        sum(weights * risks  * f.s.zi) / sum(weights*f.s.zi)    
    })
    
    if (ss.is.null) cbind(marker=ss, prob=risks) else risks
}
