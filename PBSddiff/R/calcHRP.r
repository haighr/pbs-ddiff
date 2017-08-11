#calcHRP--------------------------------2017-08-02
# Calculate Historical Ref Pts (standardise)
#-----------------------------------------------RH
calcHRP = function(A, aveYr, minYr=NULL, Burn=0, Thin=1)
{
	yr  = A$yr;  nyear = length(yr)
	yrs = A$yrs; nyrs  = length(yrs)
	HRP.list = list(); collectHRP = character()
	
	## Get MPD biomass
	mpd.bt = A$sbt ## loadScenario now replaces final year pseudo-projection, but just in case:
	names(mpd.bt) = yrs
	mpd.bt[nyrs]  = subset(A$mpdproj,tac==0)[,paste0("B",A$currYr)]
	mpd.abt       = mean(mpd.bt[is.element(yrs,aveYr)])
	mpd.dt        = mpd.bt/mpd.abt
	collectHRP = c(collectHRP, "mpd.bt", "mpd.abt", "mpd.dt")

	## Get the real current year data from projection file
	mcproj     = subset(A$mcproj, TAC==0)
	post.mcproj= as.data.frame(window(mcmc(mcproj),start=Burn+1,thin=Thin))  ## 1000 samples of Bt over N years for averaging
	post.cbt   = post.mcproj[,c(paste0("B",A$currYr),paste0(c("F"),A$ lastYr))]
	post.cbt[,paste0("U",A$lastYr)] = 1. - exp(-post.cbt[,paste0("F",A$lastYr)])
	collectHRP = c(collectHRP, "post.cbt")

	## Calculate average biomass as the median 1000 MCMC means
	mcavebt    = A$mc.sbt[,is.element(yrs,aveYr)]
	colnames(mcavebt) = yrs[is.element(yrs,aveYr)]
	post.avebt = as.data.frame(window(mcmc(mcavebt),start=Burn+1,thin=Thin))  ## 1000 samples of Bt over N years for averaging
	post.abt   = apply(post.avebt,1,function(x){mean(x)})                     ## 1000 averages over the averaging period
	avebt      = median(post.abt)                                             ## median of the 1000 averages
	collectHRP = c(collectHRP, "post.avebt", "post.abt", "avebt")

	## Calculate the depletion over the entire period
	mcbt    = A$mc.sbt[,1:nyear]
	colnames(mcbt) = yr
	post.bt = as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
	## Add the current year biomass
	post.bt = cbind(post.bt,post.cbt[,paste0("B",A$currYr)])
	colnames(post.bt) = yrs
	med.bt  = sapply(post.bt,median)
	## For each 1000 MCMCs, divide Bt (1967-2016) by the corresponding mean Bt (1967-2016 or shorter interval)
	post.dt = t(apply(post.bt,1,function(x){x/mean(x[is.element(yrs,aveYr)])}))
	post.dt2 = sweep(post.bt,1,post.abt,"/")  ## equivalent to previous line
	collectHRP = c(collectHRP, "post.bt", "med.bt", "post.dt")

	if (is.null(minYr))
		minYr  = yrs[is.element(med.bt,min(med.bt))]  ## overrides GUI value or user's value

	## The following takes minimum biomass (Bt) across year(s) designated as minimum;
	##     perhaps should be the median though if only one minimum year is specified, then it makes no difference.
	bLRPs    = apply(post.bt[,match(minYr,yrs),drop=FALSE],1,min); ## across years therefore 1000 mins
	bLRPci   = quantile(bLRPs,probs=c(0.025,0.5,0.975))
	bLRP     = median(bLRPs)
#browser();return()
	bUSRs    = apply(post.bt[,match(minYr,yrs),drop=FALSE],1,function(x){2*min(x)}); ## across years therefore 1000 mins
	bUSRci   = quantile(bUSRs,probs=c(0.025,0.5,0.975))
	bUSR     = median(bUSRs)
	collectHRP = c(collectHRP, "minYr", "bLRPs", "bLRP", "bLRPci", "bUSRs", "bUSR", "bUSRci")

	## The following takes minimum depletion (BT/Bavg) across year(s) designated as minimum;
	##     perhaps should be the median though if only one minimum year is specified, then it makes no difference.
	dLRPs    = apply(post.dt[,match(minYr,yrs),drop=FALSE],1,min); ## across years therefore 1000 mins
	dLRPci   = quantile(dLRPs,probs=c(0.025,0.5,0.975))
	dLRP     = median(dLRPs)
#browser();return()
	dUSRs    = apply(post.dt[,match(minYr,yrs),drop=FALSE],1,function(x){2*min(x)}); ## across years therefore 1000 mins
	dUSRci   = quantile(dUSRs,probs=c(0.025,0.5,0.975))
	dUSR     = median(dUSRs)
	collectHRP = c(collectHRP, "dLRPs", "dLRP", "dLRPci", "dUSRs", "dUSR", "dUSRci")

	## Calculate maximum fishing mortlaity rate
	## Years already indicate the mid-year point so there will be one less f_t than b_t.
	mcft       = A$mc.ft[,1:nyear]; colnames(mcft) = yr
	post.ft    = as.data.frame(window(mcmc(mcft),start=Burn+1,thin=Thin))     ## 1000 samples of Ft over N years
	post.ut    = 1. - exp(-post.ft)
	post.maxft = apply(post.ft,1,function(x){max(x)})                         ## 1000 maxima over N years
	maxft      = median(post.maxft)                                           ## median of the 1000 maxima
#browser();return()
	collectHRP = c(collectHRP, "post.ft", "post.ut", "post.maxft", "maxft")

	## Calculate average harvest rate as the median 1000 MCMC means
	mcaveft    = A$mc.ft[,is.element(yr,aveYr)]
	colnames(mcaveft) = yr[is.element(yr,aveYr)]
	post.aveft = as.data.frame(window(mcmc(mcaveft),start=Burn+1,thin=Thin))  ## 1000 samples of Ft over N years for averaging
	post.aveut = 1. - exp(-post.aveft)                                        ## 1000 samples of ut over N years for averaging
	post.aut   = apply(post.aveut,1,function(x){mean(x)})                     ## 1000 averages over the averaging period
	collectHRP = c(collectHRP, "post.aveft", "post.aveut", "post.aut")

	for (i in collectHRP) 
		eval(parse(text=paste0("HRP.list[[\"",i,"\"]] = ",i)))
	return(HRP.list)
#browser();return()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcHRP
#out = calcHRP(A=A, aveYr=1967:2016, Burn=200)