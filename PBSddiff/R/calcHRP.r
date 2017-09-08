#calcHRP--------------------------------2017-09-07
# Calculate Historical Ref Pts (standardise)
#-----------------------------------------------RH
calcHRP = function(A, aveYr, minYr=NULL, Burn=0, Thin=1)
{
	yr  = A$yr;  nyear = length(yr)    ## Reconstruction years
	yrs = A$yrs; nyrs  = length(yrs)   ## Reconstruction years + 1 pseudo-projection year
	HRP.list = list(); collectHRP = character()

	## Get MPD biomass
	mpd.bt = A$sbt ## loadScenario now replaces final year pseudo-projection, but just in case:
	names(mpd.bt) = yrs
	mpd.bt[nyrs]  = subset(A$mpdproj,tac==0)[,paste0("B",A$currYr)]   ## replace final year
	mpd.abt       = mean(mpd.bt[is.element(yrs,aveYr)])               ## mean biomass of the years identified for averaging
	mpd.dt        = mpd.bt/mpd.abt                                    ## MPD depletion Bt/Bavg
	collectHRP = c(collectHRP, "mpd.bt", "mpd.abt", "mpd.dt")

	## MCMC -- get the real current year data from projection file
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

	## For each row (# MCMCs) find the year of Bmin, Bmin, and Dmin (depletion)
	minList = t(apply(post.bt, 1, findBmin, aveYr=aveYr, yrs=yrs))
#browser();return()
	minList = as.data.frame(minList)
	#names(minList)=c("minYr", "Bmin", "Dmin", "Bavg")  ## no need to do this if 'findMin' purges original element names
	Bmin    = minList$Bmin
	Dmin    = minList$Dmin
	#Bavg    = minList$Bavg (no need to gather again as 'post.abt' = 'Bavg'
	if (is.null(minYr)) {  ## this makes no sense any more
		## Override GUI value or user's value
		minYrs = minList$minYr
		minYr  = findMode(minYrs)
	}

	## The following takes minimum biomass (Bt) across year(s) designated as minimum;
	#bLRPs    = apply(post.bt[,match(minYr,yrs),drop=FALSE],1,min); ## across years therefore 1000 mins
	bLRPs    = Bmin
	bLRPci   = quantile(bLRPs,probs=c(0.025,0.5,0.975))
	bLRP     = median(bLRPs)

	#bUSRs    = apply(post.bt[,match(minYr,yrs),drop=FALSE],1,function(x){2*min(x)}); ## across years therefore 1000 mins
	bUSRs    = 2 * bLRPs
	bUSRci   = quantile(bUSRs,probs=c(0.025,0.5,0.975))
	bUSR     = median(bUSRs)
	collectHRP = c(collectHRP, "minYr", "bLRPs", "bLRP", "bLRPci", "bUSRs", "bUSR", "bUSRci")

	## The following takes minimum depletion (BT/Bavg) across year(s) designated as minimum;
	#dLRPs    = apply(post.dt[,match(minYr,yrs),drop=FALSE],1,min); ## across years therefore 1000 mins (old method)
	dLRPs    = Dmin
	dLRPci   = quantile(dLRPs,probs=c(0.025,0.5,0.975))
	dLRP     = median(dLRPs)

	#dUSRs    = apply(post.dt[,match(minYr,yrs),drop=FALSE],1,function(x){2*min(x)}); ## across years therefore 1000 mins (old method)
	dUSRs    = 2 * dLRPs
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

#findBmin-------------------------------2017-09-07
# Find the minimum B (biomass) in a vector of B over time.
# Vector B must have names that denote years.
# PJS requirement: B_t trajectory after minYr has to reach B_avg
#-----------------------------------------------RH
findBmin = function(B, aveYr, yrs, Q=c(0.005,0.02,0.02))  ## Q[1]=first Qmin, Q[2]=Qmin series base, Q[3]=Qmin series increment
{
	#nS <<- nS + 1  ## for debugging to determine row being evaluated
	if (is.data.frame(B))
		B = unlist(B)
	if (is.null(names(B)) && (missing(yrs) || length(B)!=length(yrs)))
		stop("Input 'B' must either have names denoting years OR\n\tlength of input 'yrs' must equal length of 'B'.")
	if (is.null(names(B)))
		names(B) = yrs
	else if (!missing(yrs))
		B = B[is.element(names(B),as.character(yrs))]  ## in case  user wants to subset years
	else
		yrs = as.numeric(names(B))
	currYr = rev(yrs)[1]
	## Calculate a Bavg for each MCMC sample (exactly the same as 'post.abt' calculated elsewhere ~ line 71)
	Bavg =  mean(B[is.element(names(B),as.character(aveYr))])
	## Identify a set of candidate Bmins
	zout = FALSE; Qmin = 0  ## initialize
	## If no Bmins led to a recovery (all zout = FALSE) then increase Qmin and try again
	while (sum(zout)==0) {
		if (Qmin==0)         Qmin = Q[1]
		else if (Qmin==Q[1]) Qmin = Q[2]
		else                 Qmin = Qmin + Q[3]
		qmin = B <= quantile(B,Qmin)
		qmin[length(qmin)] = FALSE  ## Final year cannot be a minimum from which recovery is possible
		if (!any(qmin)) next
		ymin = yrs[qmin]
		## Process consecutive lows to find one real low
		consecs=split(B[qmin],cumsum(c(1,diff(ymin) != 1))) ## https://stackoverflow.com/questions/24837401/find-consecutive-values-in-vector-in-r
		names(consecs) = NULL
		ymin = as.numeric(names(sapply(consecs,function(x){z=min(x);x[x%in%z]})))
		## Create vectors from break points
		zbrk = c(ymin, currYr)
		## debugging: .flush.cat(nS, ": ", zbrk, "\n")
		zout = sapply(1:length(ymin), function(i,x,y,z){
			zvec = x[i]:x[i+1]
			if (length(zvec)==1) return(FALSE) ## cannot have a single year in z-vector
			zchr = as.character(zvec)
			any(y[zchr] > z)
		}, x=zbrk, y=B, z=Bavg, simplify=TRUE)
#if(length(zbrk)==3) {browser();return()}
	}
	## If no Bmins led to a recovery (all FALSE) then set to all TRUE and choose the minimum
	#if (sum(zout)==0){
#browser();return()
		#stop ("No minimum with recovery was detected in at least one of the MCMC samples.\n\tIncrease the quantile 'Qmin' to select more candidate Bmins.")
		#nomin <<- nomin + 1
		#zout = rep(TRUE,length(zout))
	#}
	if (sum(zout)==1)
		minYr = ymin[zout]
	else {
		## Find the true minimum from various Bmin options that satisfy recovery
		yymin = as.character(ymin[zout])
		minymin = B[yymin][is.element(B[yymin],min(B[yymin]))]
		minYr = as.numeric(names(minymin))
	}
#if (minYr==1989) {browser();return()}
	## As long as elements contain no names, then vector elements retain names assigned
	Omin = c(minYr=minYr, Bmin=as.vector(B[as.character(minYr)]), Dmin=as.vector(B[as.character(minYr)]/Bavg), Bavg=Bavg, Qmin=Qmin)
#{browser();return()}
	return(Omin)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findBmin

#findMode-------------------------------2017-09-06
## Function to find the mode
## https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
#-----------------------------------------------KW
findMode <- function(x)
{
	ux <- unique(x)
	ux[which.max(tabulate(match(x, ux)))]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findMode

##==========================================================
## Testing
#source("C:/Users/haighr/Files/GFish/PSARC17/WAP/DDiff/rscripts-pbs/calcHRP.r")
#out = calcHRP(A=A, aveYr=1967:2016, Burn=200)

## The following two methods yield the same median Dmin (LRP=0.199659)

## Method 1: Bmin recalculated from 8000 pooled Bt samples using 'findBmin' function
#nS = 0
#minListN2 = t(apply(A.aveN$Bt.mcmc[,1:51], 1, findBmin, aveYr=1967:2016, Q=c(0.005,0.02,0.02)))
#minListN2 = as.data.frame(minListN2)
#median(minList$Dmin)
#minListS2 = t(apply(A.aveS$Bt.mcmc[,1:51], 1, findBmin, aveYr=1967:2016, Q=c(0.005,0.02,0.02)))
#minListS2 = as.data.frame(minListS2)


## Method 2: 8000 pooled Bmin divided by 8000 pooled Bavg
#median(A.aveN$Bmin.mcmc/A.aveN$Bavg.mcmc)
