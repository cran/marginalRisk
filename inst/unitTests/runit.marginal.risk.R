test.marginal.risk <- function() {

library("RUnit")
library("marginalRisk")
tolerance=1e-5
RNGkind("Mersenne-Twister", "Inversion")


library(survival)
dat=subset(lung, !is.na(wt.loss) & !is.na(ph.ecog))

f1=Surv(time, status) ~ wt.loss + ph.ecog + age + sex
f2=wt.loss ~ ph.ecog + age + sex

fit.risk = coxph(f1, data=dat)
   fit.s = lm(f2, data=dat)
     
ss=quantile(dat$wt.loss, seq(.05,.95,by=0.1))
t0=1000
prob = marginal.risk(fit.risk, fit.s, dat, categorical.s=FALSE, t = t0, ss=ss)

checkEqualsNumeric (prob, c(0.9121199, 0.9113016, 0.9113016, 0.9110252, 0.9106587, 0.9103430, 0.9100705, 0.9097611, 0.9092597, 0.9086514), tolerance=tolerance)

}
