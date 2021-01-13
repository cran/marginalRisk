# equation 4
E.value=function(rr) { rr1=1/rr; rr1 + sqrt(rr1 * (rr1-1)) }

# equation above equation 5
bias.factor=function(RRud, RReu) RRud*RReu/(RRud+RReu-1)


controlled.risk.bias.factor=function(ss, s.cent, s1, s2, RRud) {
    # compute RR(s.cent, s) according to the display after equation (6) in the manuscript
    RR.cent.s = exp( (ss-s.cent) / (s2-s1) * log(RRud) )
    RR.s.cent = exp( (s.cent-ss) / (s2-s1) * log(RRud) )
    # compute B accoring to the display before equation (5) in the manuscript, assuming RRud and RReu are the same
    B.s.cent = RR.s.cent*RR.s.cent/(RR.s.cent+RR.s.cent-1)
    B.cent.s = RR.cent.s*RR.cent.s/(RR.cent.s+RR.cent.s-1)
    # compute bias factor according to equation (6)
    Bias = ifelse(ss>=s.cent, B.cent.s, 1/B.s.cent)
}
# ss is res[,"marker",trial]; s2 is res[s2,"marker",trial]; s1 is res[s1,"marker",trial]
#tmp=controlled.risk.bias.factor(ss=res[,"marker",trial], s.cent=s.ref, s1=res[s1,"marker",trial], s2=res[s2,"marker",trial], RRud) 
