saveFig <- function(filename){
  # Save the currently plotted figure to disk
  # in the 'figDir' folder with the name
  # filename.plotType
  if(saveon){
    filename <- paste(figDir,filename,".",plotType,sep="")
    savePlot(filename,type=plotType)
    cat(paste("Saved figure ",filename,"...\n",sep=""))
  }
}

fig.fishingMortality <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	plot(A$yr,A$ft[1,],type="l",xlab="Year",ylab="Fishing Mortality (/yr) and Harvest Rate (%)", lwd=2, 
		ylim=c(0,1.2*max(A$ft[1,])), xlim=c(A$yr[1],A$yr[length(A$yr)]), las=1)
	Z <- A$ft[1,]+A$m
	HR <- (A$ft[1,]/Z)*(1-exp(-Z))
	lines(A$yr, HR, lwd=2, col="blue")
	legend("topright", c("F", "HR"), col=c(1,"blue"), lwd=2, lty=1,bty="n")
	addScenLab()
	saveFig("fig.fishing.mortality.mpd")
}

fig.RecAnomMPD <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	plot(A$yr,A$log_rec_devs,type="p",pch=20, main="Log recruitment anomalies", col=1, xlab="Year",
		ylab="MPD Recruitment anomalies",ylim=c(-5,5), xlim=c(A$yr[1],A$yr[length(A$yr)]), cex.lab=1.2, cex.axis=1.2,las=1)
	#points(A$yr,A$anom_obs,pch=4, col="darkblue")
	#points(A$yr,A$anom_obs*A$dt,pch=4, col=2)
	abline(h=0, lty=2, lwd=0.5)
	#legend("topleft", legend=c("Predicted log rec deviations", "Additional log anomalies (obs)", "Additional log anomalies (obs) * dt (pred)"), lty=0, pch=c(20, 4,4), col=c(1,"darkblue",2), bty="n", cex=1.2)
	addScenLab()
	saveFig("fig.recruit.anomalies.mpd")
}

fig.RecAnomMCMC <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.RecAnomMCMC): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	mc <- log(A$mc.RecDevs)
	mc.recdev <- as.data.frame(window(mcmc(mc),start=Burn+1,thin=Thin)) 
	recdev <- apply(mc.recdev,2,quantile,probs=c(0.025,0.5,0.975)) #gets quantiles 
	
	#Get the observed anomalies scaled by estimates of dt
	#mc_ap <- matrix(ncol=length(A$anom_obs), nrow=length(A$mc$dt))
	#for(i in 1:length(A$mc$dt)) mc_ap[i,] <- A$anom_obs*A$mc$dt[i] #multiply the observed anomalies by the estimated mcmc values of the scalar dt
	
	#mc.anompred <- as.data.frame(window(mcmc(mc_ap),start=Burn+1,thin=Thin))
	#anompred <- apply(mc.anompred,2,quantile,probs=c(0.025,0.5,0.975)) #gets quantiles 
	
	ylim = range(recdev)
	plot(yr, recdev[2,], type="p", pch=20, lwd=2, ylim=ylim, ylab="log recruitment deviations", xlab="Year",las=1)
	arrows(yr, recdev[1, ],yr,recdev[3,],code=3,angle=90,length=0.01)
	points(yr,A$anom_obs,pch=4, col="darkblue")
	#points(yr,anompred[2, ],pch=4, col=2)
	#arrows(yr, anompred[1, ],yr,anompred[3,],code=3,angle=90,length=0.01, col=2)
	abline(h=0, lty=2, lwd=0.5)
	# legend("topleft", legend=c("Predicted log rec deviations", "Additional log anomalies (obs)", "Additional log anomalies (obs) * dt (pred)"), lty=0, pch=c(20, 4,4), col=c(1,"darkblue",2), bty="n", cex=1.2)
	addScenLab()
	saveFig("fig.recruit.anomalies.mcmc")
}

fig.weightFit <- function(){ ## modified by RH (161222)
	simplify = usingSweave
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())
	wyrs  <- A$wyrs
	Amnwt <- A$obs_annual_mean_wt
	if (is.vector(Amnwt))
		obsMW = Amnwt
	else {
#browser();return()
		wyrs  <- Amnwt[,1]
		obsMW <- Amnwt[,2]
	}
	wsig  <- A$weight_sig
	ylim=c(0,max(c(obsMW,A$annual_mean_wt[1:length(A$yr)], obsMW + wsig*obsMW)))
	if (!simplify) par(mfrow=c(1,1), mar=c(4,4,2,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.5,.5,0))
	plot(wyrs,obsMW,type="p",pch=19,xlab=ifelse(simplify,"","Year"),ylab=ifelse(simplify,"","Annual mean weight (kg)"), ylim=ylim,
		#main=ifelse(simplify,"","Mean weight"), xlim=c(A$yr[1],A$yr[length(A$yr)]), ylim=c(0,1.5*max(obsMW)), 
		main=ifelse(simplify,"","Mean weight"), xlim=c(A$yr[1],A$yr[length(A$yr)]), , 
		las=1, cex.axis=1.2, cex.lab=1.5, cex.main=1.5)
	arrows(wyrs, obsMW + wsig*obsMW, wyrs, obsMW - wsig*obsMW, code=3, angle=90, length=0.04)
	lines(A$yr, A$annual_mean_wt[1:length(A$yr)], col="red", lwd=2)
	if (!simplify){
		legend("bottomleft", legend=c("Observed (+/- CV)","Predicted"), lty=c(0,1), lwd=c(0,2), pch=c(20, NA_integer_), col=c(1,"red"), bty="n", cex=1.2, inset=0.025)
		addScenLab(0.95, ifelse(simplify,0.9,0.95), adj=c(1,0))
	}
	if (simplify && par()$mfg[1]==1 && par()$mfg[2]==1){
		legend("bottomleft", legend=c("Obs. (+/- CV)","Pred."), lty=c(0,1), lwd=c(0,2), pch=c(20, NA_integer_), col=c(1,"red"), bty="n", cex=0.9, inset=0.025)
		#addScenLab(0.95, ifelse(simplify,0.9,0.95), adj=c(1,0))
	}
	saveFig("fig.annual.mean.weight")
}

fig.catchFit <- function() {
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	crows = (1:nrow(A$obs_ct))[!apply(A$obs_ct,1,function(x){all(x==0)})]
	par(mfrow=c(1,1), mai=c(0.75,0.75,0.3,0.2), omi=c(0.1,0.1,0.1,0.1), mgp=c(2.5,.5,0))
	plot(A$yr,A$ct[1,],type="n", xlab="Year", ylab="Catch (t)", ylim=c(0,1.5*max(A$obs_ct,A$ct)), xlim=c(A$yr[1],A$yr[length(A$yr)]),main="Catch", las=1)
	for (i in crows) {
		points(A$yr,A$obs_ct[i,],pch=16,cex=1,col="red")
		lines(A$yr,A$ct[i,], col="dodgerblue", lwd=2)
		points(A$yr,A$ct[i,],pch=1,cex=1.5,col="dodgerblue")
	}
	addLegend(0.05,0.95, pch=c(16,1), col=c("red","dodgerblue"), lty=c(0,1), legend=c("Observed","Fitted"), bty="n", cex=1.2)
	addScenLab(0.95, 0.95, adj=c(1,0))
	saveFig("fig.catch.fit")
}

fig.survey.age.residuals1 <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(AgeLikelihood!=3){
		#counters
		ii <- 1
		ij<-A$na_nobs[1]
		for(k in 1:A$na_gears)
		{
			Res <- Asurv_res[ii:ij,]
			iyr=Res[,1]
			if(length(iyr)<=4) par(mfcol=c(2,2), oma=c(2,3,1,1)) #4 graphs
			if(length(iyr)>4) par(mfcol=c(2,3), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) # 5 or 6
			if(length(iyr)>6) par(mfcol=c(3,3), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #7,8 or 9
			if(length(iyr)>9) par(mfcol=c(3,4), oma=c(2,3,1,1), mai=c(0.25,0.2,0.2,0.2)) #10, 11 or 12
			if(length(iyr)>12) par(mfcol=c(4,4), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #13-16
			if(length(iyr)>16) par(mfcol=c(5,5), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #17-25

			for(i in 1:length(iyr)){
				r <- Res[i,3:(nage-age[1]+3)]
				iclr <- r
				iclr[r>=0] <- 1
				iclr[iclr!=1] <- 2

				plot(age,r,type="h",ylim=c(-6,6),main=iyr[i],col=iclr, las=1, xlab="", ylab="")
				points(age,r,pch=19,col=iclr)
			}
			mtext("Residual",2,outer=T,las=3,line=1, cex=1.3)
			mtext("Age",1,outer=T,line=1, cex=1.3)
			addScenLab()
			saveFig(paste("fig.survey.age.residuals1.survey",k,sep=""))

			#update counters
			ii<-ij+1
			if(A$na_gears>1) {
				ij<-ij+A$na_nobs[k+1]
				if(k!=A$na_gears) 	windows()
			}
		}
	} else cat("WARNING (fig.survey.age.residuals1): No plot for ageless model\n")  
}

fig.survey.age.residuals2 <- function() {
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(AgeLikelihood!=3){
		#counters
		ii <- 1
		ij<-A$na_nobs[1]

		for(k in 1:A$na_gears)
		{	
			
			Res <- Asurv_res[ii:ij,]
			iyr=Res[,1]
			 r <- Res
			bubble.plot(r[,1],age,r[,3:(nage-age[1]+3)],scale=0.3,xlab="Year",ylab="Age",add=F,log.scale=T, las=1,cex.lab=1.3)
			addScenLab()
			saveFig(paste("fig.survey.age.residuals2.survey",k,sep=""))
				
			#update counters
			ii<-ij+1
			if(A$na_gears>1) {
				ij<-ij+A$na_nobs[k+1]
				if(k!=A$na_gears) 	windows()
			}	
		}	
		
	} else cat("WARNING (fig.survey.age.residuals2): No plot for ageless model\n")  
}

fig.survey.age.props <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(AgeLikelihood!=3){
		#counters
		ii <- 1
		ij<-A$na_nobs[1]

		for(k in 1:A$na_gears)
		{	
			
			par(mfcol=c(1,1))
			Prop <- Asurv_obs[ii:ij,]
			iyr=Prop[,1]
			r <- Prop

			bubble.plot(r[,1],age,r[,3:(nage-age[1]+3)],scale=0.3,xlab="Year",ylab="Age",add=F,log.scale=T, las=1)
			addScenLab()
			saveFig(paste("fig.survey.age.props.survey",k,sep=""))
			
			#update counters
			ii<-ij+1
			if(A$na_gears>1) {
				ij<-ij+A$na_nobs[k+1]
				if(k!=A$na_gears) 	windows()
			}
		}
	} else cat("WARNING (fig.survey.age.props): No plot for ageless model\n")  
}

fig.survey.age.props.fit <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(AgeLikelihood!=3){
		#Plots the collapsed conditional age-length keys to proportions at age in the survey
		#counters
		ii <- 1
		ij<-A$na_nobs[1]

		for(k in 1:A$na_gears)
		{	
			
			Obs<- Asurv_obs[ii:ij,]
			Est<- Asurv_est[ii:ij,]
			Gear <- Asurv_obs[ii,2]
			iyr=Obs[,1]
			if(length(iyr)<=4) par(mfcol=c(2,2), oma=c(2,3,1,1)) #4 graphs
			if(length(iyr)>4) par(mfcol=c(2,3), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) # 5 or 6
			if(length(iyr)>6) par(mfcol=c(3,3), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #7,8 or 9
			if(length(iyr)>9) par(mfcol=c(3,4), oma=c(2,3,1,1), mai=c(0.25,0.2,0.2,0.2)) #10, 11 or 12
			if(length(iyr)>12) par(mfcol=c(4,4), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #13-16
			if(length(iyr)>16) par(mfcol=c(5,5), oma=c(2,3,1,1), mai=c(0.2,0.2,0.2,0.2)) #17-25

		for(i in 1:length(iyr)){
				mp <- barplot(Obs[i,3:(nage-age[1]+3)],names.arg=paste(age),ylim=c(0.,0.8),main=iyr[i], las=1,legend.text="Observed", args.legend=c(bty="n"))
				lines(mp,Est[i,3:(nage-age[1]+3)],type="o",pch=20,cex=2,col=2)
				legend("topright", legend="Predicted",col=2, pch=19, bty="n")
		}
		mtext("Proportion",2,outer=T,las=3,line=2, cex=1.2)
		mtext("Age",1,outer=T,line=-1, cex=1.2)
		mtext(paste("Gear",Gear),3,outer=T, line=-1,cex=1.2)
		addScenLab()
		saveFig(paste("fig.survey.age.props.fit.survey",k,sep=""))
		
		#update counters
		ii<-ij+1
		if(A$na_gears>1) {
			ij<-ij+A$na_nobs[k+1]
			if(k!=A$na_gears) 	windows()
		}
	}
	} else cat("WARNING (fig.survey.age.props.fit): No plot for ageless model\n")  
}

### SURVEY BIOMASS FITS
#RF NOTE - ADD AN MCMC VERSION OF THIS
fig.surveybiomass.fit <- function(n=NULL)
{
	simplify = usingSweave
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())

	if (exists("cribtab",envir=.GlobalEnv)){
		if (simplify && exists("CRIBTAB",envir=.GlobalEnv)) cribtab = CRIBTAB
		Iseries  = attributes(cribtab)$indexSeries
		#c("NMFS Triennial","WCVI Synoptic","QCS Synoptic","Hecate St. Synoptic","WCHG Synoptic","Commercial Trawl CPUE") ## ONLY relevant to 2015 SST
		#http://stackoverflow.com/questions/17009628/extracting-numbers-from-string-in-r
		runNo = unlist(regmatches(A$ControlFile,gregexpr("[0-9]+", A$ControlFile)))
		runNm = paste0(getOptions(wapmod,"fdPrefix"),runNo)
		runI  = eval(parse(text=cribtab[runNm,"I"]))
		SurvName = Iseries[runI]
		if (all(is.na(runI))) runI=1:A$nits
	} else {
		runI = 1:A$nits
		SurvName = paste0("Index Series ",runI)
	}

	#get list of cvs from data file (horrible name format, will fix one day)
	Survdat = as.matrix(A$dat$"#Survey")  ## survey indices always start at 2 regardless of skipped surveys
	Survtab = table(Survdat[,3])
	Survgrs = unique(Survdat[,3])
	names(Survtab) = runI

	## Attempt to calculate the right survey from info in cribtab
	Survdat[,3] =  as.numeric(rep(names(Survtab),Survtab))
	Survgears <- unique(Survdat[,3])
#browser();return()
	if (!all(Survgears==runI)) stop("Survey indexing has gone wacky")

	if (simplify && !is.element(n,Survgears)) {
		frame(); addLabel(.5,.5,"No index",col="red",cex=2); return()
	}
	#print(Survdat)
	if(A$nits==1){
		A$iyr <- t(as.matrix(A$iyr))
		A$it  <- t(as.matrix(A$it))
		A$pit <- t(as.matrix(A$pit))
	}
	if (!simplify) {
		rc = .findSquare(A$nits)
		expandGraph(mfrow=rc, mar=c(3.5,4,2,1), oma=c(0,0,0,0), mgp=c(2.5,.5,0))
	}

	### RH: The results matrix did not have ncol = max # survey years (at least when there is commercial CPUE)
	### RH: OK, I changed the reptoRlist function and now the matrix has one row for each survey (150908)
	survindex = Survdat[,3]
	iyrvec = as.vector(t(A$iyr)); iyrvec = iyrvec[!is.na(iyrvec)]
	pitvec = as.vector(t(A$pit)); pitvec = pitvec[!is.na(pitvec)]
	itvec  = as.vector(t(A$it));  itvec  = itvec[!is.na(itvec)]

	N = if (simplify) n else 1:A$nits
	for (i in N){
		sgear = ifelse(simplify,i,Survgears[i])
		z   = is.element(survindex,sgear)  
		iyr = iyrvec[z]
		pit = pitvec[z]
		it  = itvec[z]
		#iyr <- A[["survBfits"]][[paste0("surv",i)]][["iyr"]] <<- A$iyr[i,which(A$iyr[i,]>0)] #get rid of NAs  ## doesn't save to master A list unless use "<<-"
		#pit <- A[["survBfits"]][[paste0("surv",i)]][["pit"]] <<- A$pit[i,which(A$pit[i,]>0)] #get rid of NAs
		#it  <- A[["survBfits"]][[paste0("surv",i)]][["it"]] <<- A$it[i,which(A$it[i,]>0)] #get rid of NAs
#if (i==7) {browser(); return()}

		#Get annual CVs from survey1
		srows <- which(Survdat[,3]==sgear)
		CV    <- 1/(Survdat[srows,4] )
		ylim = extendrange(c(it+CV*it,pit)); ylim[1]=0
		kt   = median(ylim) > 1000
		if (kt) {
			ylim = ylim / 1000.
			it   = it   / 1000.
			pit  = pit  / 1000.
		}

		plot(iyr, pit, type="n", col="slateblue", lwd=2, xlab=ifelse(simplify,"","Year"), ylab=ifelse(simplify,"",paste0("Index", ifelse(kt, " (,000)", ""))), main=ifelse(simplify,"",SurvName[i]), ylim=ylim, xlim=c(iyr[1],rev(iyr)[1]), las=ifelse(simplify,0,1), cex.lab=1.5, cex.axis=ifelse(simplify,1,1.2))
#browser();return()
		arrows(iyr, it+CV*it, iyr, it-CV*it, code=3,angle=90,length=0.03,lwd=2,col="grey70")
		points(iyr,it, pch=21, col="blue", bg="gold", cex=ifelse(simplify,1.5,2))
		points(iyr, pit, col="slateblue", pch=25, lwd=2, cex=ifelse(simplify,1.25,1.5))
		if (!simplify && par()$mfg[1]==1 && par()$mfg[2]==1){
			addScenLab(adj=c(0,0.5),cex=1.2)
		}
	}
	saveFig("fig.survey.biomass.fit")
}

fig.selectivity <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	
	if(A$delaydiff==0){
		sel<-A$selectivity
		vva<-sel[which(sel[,1]==1),] 
		va<-vva[length(vva[,1]),2:length(vva[1,])] #pick final year
		leg <- "Gear 1" #legend text
		cha <- 20 #plot character for legend
		
		plot(age, va,type="o",xlab="Age",ylab="Selectivity in final year", las=1, col=1, pch=cha) #ylim=c(0,1),
		for(i in 2:ngear) {	
			vva2<-sel[which(sel[,1]==i),] 
			va2<-vva2[length(vva2[,1]),2:length(vva2[1,])] #first column is gear id
			lines(age,va2,type="o",pch=i, lty=i, col=i)
			abline(h=0.5, lty=4, lwd=2, col="darkgray")
			leg <- c(leg,paste("Gear",i))
			cha <- c(cha,i)
		}
		leg <- c(leg,"50%")
		cha <- c(cha, -1)
		
		legend("bottomright",leg,lty=c(1:ngear,4),pch=cha,bty="n", cex=1.5, col=c(1:ngear,"darkgray"), lwd=c(rep(1,ngear),2))
		addScenLab()
		saveFig("fig.selectivities")
	} else cat("WARNING (fig.selectivity): No estimated selectivity for delay difference model\n")
}

fig.phase <- function(yUpperLimit=2){
	# The phase plot of Bt/Bmsy vs (1-SPR)/(1-SPRmsy)
	# yUpperLimit is the upper limit of the y-axis
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	spr <- A$sprmsy_status
	if (is.null(spr)) {
		cat("WARNING (fig.phase): No spawner-per-recruit data available\n"); return(invisible("No SPR data")) }
	sbstatus <- A$sbtmsy_status
	sbstatus <- sbstatus[1:length(A$yr)]
	maxX <- max(sbstatus)
	maxY <- yUpperLimit

	plot(sbstatus, spr, type="n", col=1, ylim=c(0, maxY), xlim=c(0, maxX), xlab="B/BTarget", ylab="(1-SPR)/(1-SPRTarget)")
	#dum <- fried.egg(sbstatus, spr)
	lines(sbstatus, spr, type="o",col="blue")
	colVector <- vector(length=length(A$yr))
	colVector[1] <- "green"
	colVector[length(colVector)] <- "red"
	for(i in 2:(length(colVector)-1)){
		colVector[i] <- 0
	}
	points(sbstatus, spr, pch=20, col=colVector)
	abline(h=1, v=1, lty=3,col="red")
	addScenLab()
	saveFig("fig.phase")
}

### Pairs with histograms
fig.estimated.params.pairs <- function(ghatGear=1){
  # ghatGear: 1 is survey 2 is commercial
	#Pair plots of estimated parameters
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.estimated.params.pairs): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc",paste0("q",1:A$nits)) ## RH not all estimated
	nPar <- length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

	j = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	doLog = grep("^q[1-9]",estPar)
	j[,doLog] = log(j[,doLog])
	newPar = estPar
	newPar[zLog]=paste0("ln(",sub("m","M",gsub("log\\.","",newPar[zLog])),")")
	names(j) = newPar

	mle <- vector(mode="numeric",length=nPar)
	names(mle) = newPar

	zPar = sapply(j,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
#browser();return()

	#jmess=j[,plotPar]; attr(jmess,"mpd")=mle[plotPar]
	pairs(j[,plotPar],gap=0,pch=20,cex=.2,upper.panel=panel.smoothie,diag.panel=panel.histie, lower.panel=panel.smoothie,cex.axis=1.2,cex.labels=1.5)
	addScenLab()
	saveFig("fig.estimated.params.pairs")
}

panel.smoothie = function (x, y, col = "gainsboro", bg = NA, pch = 0.2, 
   cex = 1, col.smooth = "blue", span = 0.5, iter = 3, ...) 
{
	points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	ok <- is.finite(x) & is.finite(y)
	if (any(ok)) 
		lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, lwd=2, ...)
}
panel.histie <- function(x, mpd, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="moccasin", ...)
	#if (!is.null(attributes(x)$mpd)) abline(v=attributes(x)$mpd,lty=2,col="red")
}

### Pairs with histograms (no logs)
fig.estimated.params.pairs2 <- function(ghatGear=1){
	# ghatGear: 1 is survey 2 is commercial
	#Pair plots of estimated parameters
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.estimated.params.pairs2): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc",paste0("q",1:A$nits)) ## RH not all estimated
	nPar <- length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

	j = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	noLog = grep("^[log]",estPar)
	j[,noLog] = exp(j[,noLog])
	newPar = estPar
	newPar[zLog]=sub("m","M",gsub("log\\.","",newPar[zLog]))
	names(j) = newPar

	mle <- vector(mode="numeric",length=nPar)
	names(mle) = newPar

	zPar = sapply(j,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
#browser();return()

	#jmess=j[,plotPar]; attr(jmess,"mpd")=mle[plotPar]
	pairs(j[,plotPar],gap=0,pch=20,cex=.2,upper.panel=panel.smoothie,diag.panel=panel.histie, lower.panel=panel.smoothie,cex.axis=1.2,cex.labels=1.5)
	addScenLab()
	saveFig("fig.estimated.params.pairs.no.log")
}

### KEY PAIRS ONLY
fig.estimated.params.pairs.key <- function(ghatGear=1)
{
	#Pair plots of estimated parameters
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.estimated.params.pairs.key): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc") #,paste0("q",1:A$nits)) ## RH not all estimated
	nPar <- length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

	j = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	doLog = grep("^q[1-9]",estPar)
	j[,doLog] = log(j[,doLog])
	newPar = estPar
	newPar[zLog]=paste0("ln(",sub("m","M",gsub("log\\.","",newPar[zLog])),")")
	names(j) = newPar

	mle <- vector(mode="numeric",length=nPar)
	names(mle) = newPar

	zPar = sapply(j,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)

	#jmess=j[,plotPar]; attr(jmess,"mpd")=mle[plotPar]
	pairs(j[,plotPar],gap=0,pch=20,cex=1,upper.panel=panel.smoothie,diag.panel=panel.histie, lower.panel=panel.smoothie,cex.axis=1.2,cex.labels=2)
	addScenLab()
	saveFig("fig.key.estimated.params.pairs")
}

### KEY PAIRS ONLY (NO LOG)
fig.estimated.params.pairs.no.log.key <- function(ghatGear=1)
{
	#Pair plots of estimated parameters
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.estimated.params.pairs.no.log.key): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc")#,paste0("q",1:A$nits)) ## RH not all estimated
	nPar <- length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

	j = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	noLog = grep("^[log]",estPar)
	j[,noLog] = exp(j[,noLog])
	newPar = estPar
	newPar[zLog]=sub("m","M",gsub("log\\.","",newPar[zLog]))
	names(j) = newPar

	mle <- vector(mode="numeric",length=nPar)
	names(mle) = newPar

	zPar = sapply(j,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)

	pairs(j[,plotPar],gap=0,pch=20,cex=1,upper.panel=panel.smoothie,diag.panel=panel.histie, lower.panel=panel.smoothie,cex.axis=1.2,cex.labels=2)
	addScenLab()
	saveFig("fig.key.estimated.params.pairs.no.log")
}

fig.mcmc.priors.vs.posts <- function(exFactor=1.0, qPriorFunction=4,ghatGear=1,ghatPriorFunction=4,showEntirePrior=T){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.priors.vs.posts): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
  # exFactor is a multiplier for the minimum and maximum xlims.
  # qPriorFunction: 4 for gamma
  # ghat=1 for survey selectivity, 2 for commercial selectivity
  # showEntirePrior, if T then plot the entire prior function to its limits
  #  and ignore posterior distribution limits
  #  plus one in code to compensate for the transfer from 0-based to 1-based
  qPriorFunction <- qPriorFunction + 1
 
  # npar is the number of parameters with priors + 1 for q which is seperate from the others
	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc",paste0("q",1:A$nits)) ## RH not all estimated
	nPar <- A$npar + A$nits - 1
	#zPar = sapply(A$mc,function(x){(length(sort(unique(x)))>1)})
	#nPar = sum(zPar)
	#estPar = names(A$mc)[zPar]
	zMCMC  = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC  = length(zMCMC)
   #if(A$nits==1) nPar <- nPar  # ( for q)
   #if(A$nits==2) nPar <- nPar + 1 # (+1 for q)
   #if(A$nits==3) nPar <- nPar + 2 # (+2 for q)
   #if(A$nits==4) nPar <- nPar + 3 # (+3 for q)
 
#browser();return()
	d = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	doLog = grep("^q[1-9]",estPar)
	d[,doLog] = log(d[,doLog])
	newPar = estPar
	newPar[zLog]=paste0("ln(",sub("m","M",gsub("log\\.","",newPar[zLog])),")")
	names(d) = newPar
  #d <- matrix(nrow=nMCMC,ncol=nPar, dimnames=list(mcmc=zMCMC,par=estPar))
  mle <- vector(mode="numeric",length=nPar)
  names(mle) = newPar
  # The values in the REPORT file for each of priorN are:
  # 1. ival  = initial value
  # 2. lb    = lower bound
  # 3. ub    = upper bound
  # 4. phz   = phase
  # 5. prior = prior distribution funnction
  #             0 = Uniform
  #             1 = normal    (p1=mu,p2=sig)
  #             2 = lognormal (p1=log(mu),p2=sig)
  #             3 = beta      (p1=alpha,p2=beta)
  #             4 = gamma     (p1=alpha,p2=beta)
  # 6. p1 (defined by 5 above)
  # 7. p2 (defined by 5 above)
  functionNames <- c(dunif,dnorm,dlnorm,dbeta,dgamma)
  functionNamesR <- c(runif,rnorm,rlnorm,rbeta,rgamma)
  # Set up column names depending on switches
  #  d[,1] <- A$mc$log.ro[Burn:nrow(A$mc)]
  #  d[,2] <- A$mc$h[Burn:nrow(A$mc)]
  #d[,3] <- A$mc$log.m[Burn:nrow(A$mc)]
  #d[,4] <- A$mc$log.rbar[Burn:nrow(A$mc)]
  #d[,5] <- A$mc$log.rinit[Burn:nrow(A$mc)]
  #d[,6] <- A$mc$rho[Burn:nrow(A$mc)]
  #d[,7] <- A$mc$varphi[Burn:nrow(A$mc)]
  #d[,8] <- A$mc$qc[Burn:nrow(A$mc)]
  # d[,9] <- log(A$mc$q1[Burn:nrow(A$mc)])
  #if(A$nits>1) d[,10] <- log(A$mc$q2[Burn:nrow(A$mc)])
  #if(A$nits>2) d[,11] <- log(A$mc$q3[Burn:nrow(A$mc)])
  #if(A$nits==4) d[,12] <- log(A$mc$q4[Burn:nrow(A$mc)])
  #if(A$nits == 1) colnames(d) <- c("ln(r0)","h","ln(M)","ln(rbar)","ln(rinit)","rho","varphi","qc","ln(q1)")
  #if(A$nits == 2) colnames(d) <- c("ln(r0)","h","ln(M)","ln(rbar)","ln(rinit)","rho","varphi","qc","ln(q1)","ln(q2)")
  #if(A$nits == 3) colnames(d) <- c("ln(r0)","h","ln(M)","ln(rbar)","ln(rinit)","rho","varphi","qc","ln(q1)","ln(q2)","ln(q3)")
  #if(A$nits == 4) colnames(d) <- c("ln(r0)","h","ln(M)","ln(rbar)","ln(rinit)","rho","varphi","qc","ln(q1)","ln(q2)","ln(q3)","ln(q4)")

   # Set up prior function names 
  #mle[1] <- A$theta[1]
  #mle[2] <- A$theta[2]
  #mle[3] <- A$theta[3]
  #mle[4] <- A$theta[4]
  #mle[5] <- A$theta[5]
  #mle[6] <- A$theta[6]
  #mle[7] <- A$theta[7]
  #mle[8] <- A$theta[8]
  #if(A$nits == 1) mle[9] <- log(A$q[1])
  #if(A$nits > 1) mle[10] <- log(A$q[2])
  #if(A$nits > 2) mle[11] <- log(A$q[3])
  #if(A$nits == 4) mle[12] <- log(A$q[4])
	mle[1:8] = A$theta[1:8]
	mle[9:nPar] = log(A$q)

  #post.samp <- window(mcmc(d),start=Burn+1,thin=Thin)  ## already adjusted above using zMCMC
  #colnames(post.samp) <- colnames(d)
  nm <- colnames(d)
  # the following are +1 because they are 0-based in the CTL file,
  # and R's vectors are 1-based
  fn <- c(functionNames[A$priorPars_theta1[5]+1],
          functionNames[A$priorPars_theta2[5]+1],
          functionNames[A$priorPars_theta3[5]+1],
          functionNames[A$priorPars_theta4[5]+1],
          functionNames[A$priorPars_theta5[5]+1],
          functionNames[A$priorPars_theta6[5]+1],
          functionNames[A$priorPars_theta7[5]+1],
          functionNames[A$priorPars_theta8[5]+1],
          functionNames[A$q_prior[1:A$nits]+1])
  mu <- c(A$priorPars_theta1[6],
          A$priorPars_theta2[6],
          A$priorPars_theta3[6],
          A$priorPars_theta4[6],
          A$priorPars_theta5[6],
          A$priorPars_theta6[6],
          A$priorPars_theta7[6],
          A$priorPars_theta8[6],
          rep(log(GT0(0)),A$nits) )
          #A$q_mu)
  sig <- c(A$priorPars_theta1[7],
           A$priorPars_theta2[7],
           A$priorPars_theta3[7],
           A$priorPars_theta4[7],
           A$priorPars_theta5[7],
           A$priorPars_theta6[7],
           A$priorPars_theta7[7],
           A$priorPars_theta8[7],
           #A$q_sd)
           rep(log(1),A$nits) )
  fnr <- c(functionNamesR[A$priorPars_theta1[5]+1],
          functionNamesR[A$priorPars_theta2[5]+1],
          functionNamesR[A$priorPars_theta3[5]+1],
          functionNamesR[A$priorPars_theta4[5]+1],
          functionNamesR[A$priorPars_theta5[5]+1],
          functionNamesR[A$priorPars_theta6[5]+1],
          functionNamesR[A$priorPars_theta7[5]+1],
          functionNamesR[A$priorPars_theta8[5]+1],
          functionNamesR[A$q_prior[1:A$nits]+1])         

	zPar = sapply(d,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
	rc = .findSquare(length(plotPar))
	par(mfrow=rc, mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))
#browser();return()

  #plotPar = c(1:5,9:nPar)
  #if(length(plotPar<=9)) par(mfrow=c(3,3),mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))  # show all parameters both fixed and estimated
  #else	 par(mfrow=c(4,4),mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))

	for(param in plotPar){
		x <- list(p=d[,param],mu=mu[param],sig=sig[param],fn=fn[[param]],nm=nm[param],mle=mle[param])
#if(param==9) {browser();return()}
		plot.marg(x,breaks="sturges",col="moccasin",exFactor=exFactor)
		if (param == plotPar[1]) addScenLab(cex=1.2)
	}
	saveFig("fig.mcmc.priors.vs.posts")
return(invisible(plotPar))  ## RH: I don't think we need to see prior plots alone
  windows()
  expandGraph()

 #Plot priors on their own
#  if(length(plotPar<=9)) par(mfrow=c(3,3),mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))  # show all parameters both fixed and estimated
#	else	 par(mfrow=c(4,4),mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))

	n<-0
	exfactor=c(1,1.1,2.75,1,1,10,10,1,rep(1,A$nits))
	for(param in plotPar){
		n=n+1
		px <- list(p=d[,param],mu=mu[param],sig=sig[param],fn=fn[[param]],nm=nm[param],mle=mle[param],fnr=fnr[[param]])
#if(n==2) {browser();return()}
		plot.prior(px,breaks="sturges",col="blue",exFactor=exfactor[n]) 
	}
	addScenLab()
	saveFig("fig.priors")
}


plot.marg <- function(xx,breaks="sturges",exFactor=1.0,...)
{
	# xx is a list(p=samples, mu=prior mean, s=prior varian, fn=prior distribution)
	# exFactor is a multiplier for the minimum and maximum xlims.
	# showEntirePrior, if T then plot the entire prior function to its limits
	#  and ignore posterior distribution limits
	ssNoPlot <- hist(xx$p,breaks=breaks,plot=FALSE)
	xl <- seq(min(ssNoPlot$breaks)/exFactor,max(ssNoPlot$breaks)*exFactor,length=250)
	pd <- xx$fn(xl,xx$mu,xx$sig)
	z <- cbind(xl,pd)
	Xlim <- range(c(xl,xx$mle)) #c(min(xl),max(xl))
	ss <- hist(xx$p,prob=T,breaks=breaks,main=xx$nm,xlab="", cex.axis=1.2, cex.main=1.5, xlim=Xlim,...)#,
	lines(xl,pd,col="blue",lwd=2)     
	abline(v=xx$mle, lwd=2, lty=2, col=2)
}

plot.prior <- function(xx,breaks="sturges",exFactor=1.0,...){
	#xx is a list(p=samples, mu=prior mean, s=prior varian, fn=prior distribution)
	randsamp<-xx$fnr(10000,  xx$mu, xx$sig) #get random sample
	xl <- seq(min(randsamp),max(randsamp),length=250)
	pd <- xx$fn(xl,xx$mu,xx$sig)
	plot(xl,pd,col=1,lwd=2, type="l", ,main=xx$nm, las=1, cex.axis=1.2)     
}

#NO LOGS OR PRIORS 
# exFactor is a multiplier for the minimum and maximum xlims.
fig.mcmc.priors.vs.posts2 <- function(exFactor=1.0,qPriorFunction=4,ghatGear=1,ghatPriorFunction=4){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.priors.vs.posts2): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

  # npar is the number of parameters with priors + 1 for q which is seperate from the others
  #nPar <- A$npar +4 # plotting 4 ref points
	basePar = estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc",paste0("q",1:A$nits)) ## RH not all estimated
	msyPar  = c("bo","bmsy","msy","fmsy")
	if(A$nits>=1)
		estPar = c(basePar,msyPar)
	bPar  = length(basePar)
	nPar  = length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

#  d <- matrix(ncol=nPar,nrow=length(Burn:nrow(A$mc)))
#  mle <- vector(mode="numeric",length=nPar)

	d = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	noLog = grep("^[log]",estPar)
	d[,noLog] = exp(d[,noLog])
	newPar = estPar
	newPar[zLog]=sub("m","M",gsub("log\\.","",newPar[zLog]))
	names(d) = newPar
	mle = vector(mode="numeric",length=bPar)
	names(mle) = newPar[1:bPar]

	# Set up mpd
	mle[1:8] = A$theta[1:8]
	mle[9:bPar] = A$q
	mle[noLog] = exp(mle[noLog])

	#post.samp <- window(mcmc(d),start=Burn+1,thin=Thin)   ## already adjusted above using zMCMC
	#colnames(post.samp) <- colnames(d)
	nm <- colnames(d)

	zPar = sapply(d,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
	rc = .findSquare(length(plotPar))
	par(mfrow=rc, mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))

	for(param in plotPar){
		x <- list(p=d[,param],nm=nm[param],mle=mle[param])
		plot.marg2(x,breaks=17,col="moccasin",exFactor=exFactor,showEntirePrior,xlim=extendrange(r=range(c(x$p,x$mle),na.rm=T)))
		if (param == plotPar[1]) addScenLab(x=0.95, adj=c(1,0.5), cex=1.2)
	}
	saveFig("fig.mcmc.priors.vs.posts2")
#browser();return()
}

#NO LOGS OR PRIORS 
fig.mcmc.priors.vs.postskey <- function(exFactor=1.0,qPriorFunction=4,ghatGear=1,ghatPriorFunction=4)
{
# exFactor is a multiplier for the minimum and maximum xlims.
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.priors.vs.postskey): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

  # npar is the number of parameters with priors + 1 for q which is seperate from the others
	basePar = estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc") #,paste0("q",1:A$nits)) ## RH not all estimated
	msyPar  = c("bo","bmsy","msy","fmsy")
	if(A$nits>=1)
		estPar = c(basePar,msyPar)
	bPar  = length(basePar)
	nPar  = length(estPar)
	zMCMC = (Burn+1):nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC = length(zMCMC)

#  d <- matrix(ncol=nPar,nrow=length(Burn:nrow(A$mc)))
#  mle <- vector(mode="numeric",length=nPar)

	d = A$mc[zMCMC,estPar]
	zLog = grep("^[log]|^q[1-9]",estPar)
	noLog = grep("^[log]",estPar)
	d[,noLog] = exp(d[,noLog])
	newPar = estPar
	newPar[zLog]=sub("m","M",gsub("log\\.","",newPar[zLog]))
	names(d) = newPar
	mle = vector(mode="numeric",length=bPar)
	names(mle) = newPar[1:bPar]

	# Set up mpd
	mle[1:8] = A$theta[1:8]
	#mle[9:bPar] = A$q
	mle[noLog] = exp(mle[noLog])


	#post.samp <- window(mcmc(d),start=Burn+1,thin=Thin)  ## already adjusted above using zMCMC
	#colnames(post.samp) <- colnames(d)
	nm <- colnames(d)

	zPar = sapply(d,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
	rc = .findSquare(length(plotPar))
	par(mfrow=rc, mai=c(0.3,0.2,0.3,0.2), omi=c(0.2,0.2,0.2,0.2))

	for(param in plotPar){
		x <- list(p=d[,param],nm=nm[param],mle=mle[param])
		plot.marg2(x,breaks=17,col="moccasin",exFactor=exFactor,showEntirePrior,
			xlim=extendrange(r=range(c(x$p,x$mle),na.rm=T)))
		if (param == plotPar[1]) addScenLab(x=0.95, adj=c(1,0.5), cex=1.2)
	}
	saveFig("fig.mcmc.priors.vs.posts2")
}


plot.marg2 <- function(xx,breaks="sturges",exFactor=1.0,showEntirePrior=F,...){
  #xx is a list(p=samples, mu=prior mean, s=prior varian, fn=prior distribution)
  # exFactor is a multiplier for the minimum and maximum xlims.
  # showEntirePrior, if T then plot the entire prior function to its limits
  #  and ignore posterior distribution limits
	ssNoPlot <- hist(xx$p,breaks=breaks,plot=F)
	xl <- seq(min(ssNoPlot$breaks),max(ssNoPlot$breaks),length=250)
	ss <- hist(xx$p,prob=T,breaks=breaks,main=xx$nm,xlab="", cex.axis=1.2, cex.main=1.5,...)
	abline(v=xx$mle, lwd=2, lty=2, col=2)
}

fig.mcmc.trace <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.trace): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	estPar = parEstimated()
	nPar   = length(estPar)
	Iseries = attributes(cribtab)$indexSeries[eval(parse(text=(cribtab[getScenLab(),"I"])))]
	#mcmcData <- window(mcmc(A$mc[,estPar]), start=Burn+1, thin=Thin)
	mcmcData <- mcmc2(A$mc[,estPar], start=Burn+1, thin=Thin)
	names(mcmcData)[grep("^q[1-9]",names(mcmcData))] = paste0("q (",Iseries,")")
#browser();return()

	rc = .findSquare(nPar)
	expandGraph(mfrow=rc, mar=c(2,3,2,2))
	for(param in 1:length(estPar)){
		mcmcTrace <- as.matrix(mcmcData[,param])
		plot(mcmcTrace,main=colnames(mcmcData)[param],type="l",ylab="",xlab="",axes=F,col="gainsboro")
		lines(cquantile.vec(mcmcTrace,0.05),col="blue",lty=2)
		lines(cquantile.vec(mcmcTrace,0.95),col="blue",lty=2)
		lines(cquantile.vec(mcmcTrace,0.5),col="red")
		points(mcmcTrace[1],col="red",bg="pink",pch=21,cex=1.5)
		box()
		#at <- seq(0,end(mcmcData)-start(mcmcData),200)
		#labels <- seq(start(mcmcData),end(mcmcData),200)
		at <- seq(0,nrow(mcmcData),200)
		labels <- seq(0,nrow(mcmcData),200)
		axis(1,at=at,labels=labels)
		axis(2)
		if (param==1) addScenLab()
	}
	saveFig("fig.mcmc.trace")
}

### MCMC DENSITY
fig.mcmc.density <- function(color="blue",opacity="20")
{
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.density): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	shade <- getShade(color,opacity)
	d = A$mc
	zPar = sapply(d,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
	rc = .findSquare(length(plotPar))
	par(mfrow=rc, mai=c(0.4,0.4,0.3,0.2), omi=c(0.1,0.1,0.1,0.1))

	mcmcData <- window(mcmc(d[,plotPar]), start=Burn+1, thin=Thin)
	for(param in 1:ncol(mcmcData)){
		ii = colnames(mcmcData)[param]
		err = try(dens <- density(mcmcData[,param], bw="SJ", kernel="gaussian"), silent=T)
		if (!inherits(err,"try-error")) {
			plot(dens, main=colnames(mcmcData)[param])
			xx <- c(dens$x,rev(dens$x))
			yy <- c(rep(min(dens$y),length(dens$y)),rev(dens$y))
			polygon(xx, yy, density=NA, col=shade)
		} else {
			hist(mcmcData[,param], nclass=10,main=ii,col=shade,xlab="Sparse sample")
		}
		if (param==1) addScenLab(0.95, 0.95, adj=c(1,1))
#browser();return()
	}
	saveFig("fig.mcmc.density")
}

fig.mcmc.autocor <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.autocor): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	estPar = parEstimated()
	nPar   = length(estPar)
	Iseries = attributes(cribtab)$indexSeries[eval(parse(text=(cribtab[getScenLab(),"I"])))]
	mcmcData <- window(mcmc(A$mc[,estPar]), start=Burn+1, thin=Thin)
	colnames(mcmcData)[grep("^q[1-9]",colnames(mcmcData))] = paste0("q (",Iseries,")")
	estPar[grep("^q[1-9]",estPar)] = paste0("q (",Iseries,")")

	Lags=  c(0, 1, 5, 10, 15,20,30,40,50)
	lags <- matrix(nrow=length(Lags), ncol=length(estPar))
#browser();return()
	colnames(lags) = colnames(mcmcData[,estPar])

	minmax = apply(mcmcData[,estPar],2,function(x){range(acf(x,lag=100,plot=F)[[1]][-1],na.rm=TRUE)})
	ylim   = c(min(-0.15,minmax[1,]), max(minmax[2,]))

	rc = .findSquare(nPar)
	expandGraph(mfrow=rc, mar=c(3,3,2,0.5))

	n <- 0
	for(i in estPar){
		n <- n+1
		#x <- window(mcmc(A$mc[,i]),start=Burn+1,thin=Thin)
		x <- mcmcData[,i]
		ac <- autocorr(x, lags = Lags, relative=TRUE)
		lags[,n] <- ac
		# autocorr.plot(x, lag.max=100, auto.layout=F, main=i) #paste(colnames(lags)[n]))  ## does not add CI
		acf(x, lag=100, ylim=ylim, main="",cex.axis=1.2,cex.lab=1.5)
		mtext(i,side=3,cex=1,line=0.25)
#browser();return()
		if (n==1) addScenLab(0.95, 0.95, adj=c(1,1))
	}
	saveFig("fig.mcmc.autocor")
}

fig.mcmc.gelman <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.gelman): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	x <- window(mcmc(A$mc[,1:12]), start=Burn+1, thin=Thin)
	gelman.plot(x,bin.width = 10,max.bins = 50,confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)
	addScenLab()
	saveFig("fig.mcmc.gelman")
}

fig.mcmc.geweke <- function(frac1=0.1,frac2=0.5,nbins=20,pvalue=0.05,silent=F, useShades=F, ...){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.geweke): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	zPar = sapply(A$mc,function(x){(length(sort(unique(x)))>1)})
	plotPar = grep(TRUE,zPar)
	rc = .findSquare(sum(zPar))

	x <- window(mcmc(A$mc[,plotPar]),start=Burn+1,thin=Thin)
	#x <- mcmc(read.table("gewekeExample.csv",sep=",",header=T)[,1:10])

  ystart <- seq(from = start(x), to = (start(x) + end(x))/2, length = nbins)
  gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)), 
               dimnames = c(ystart, varnames(x), chanames(x)))
  for (i in 1:length(ystart)) {
    geweke.out <- try(geweke.diag(window(x, start = ystart[i]),frac1 = frac1, frac2 = frac2), silent=silent)
    for (k in 1:nchain(x)){
      gcd[i, , k] <- geweke.out[[k]]
    }
  }
#browser();return()

  expandGraph(mfrow=rc, mar=c(3,2.75,2,0), las=0)
  climit <- qnorm(1 - pvalue/2)
  for (k in 1:nchain(x)){
    for (j in 1:nvar(x)) {
      ylimit <- max(c(climit, abs(gcd[, j, k])))
      plot(ystart, gcd[, j, k], type = "p", xlab = "Iteration", 
           ylab = "Z-score", pch = 1, ylim = c(-ylimit, ylimit), col="blue", 
           ...)
      lines(ystart,gcd[,,1][,j],lwd=2, col="blue")
      if(useShades){
        xx <- c(ystart,rev(ystart))
        yy <- c(rep(100*max(gcd[,,1][,j]),length(gcd[,,1][,j])),rev(gcd[,,1][,j]))
        shade <- getShade("green","15")
        polygon(xx,yy,density=NA,col=shade)

        xx <- c(ystart,rev(ystart))
        yy <- c(rep(100*min(gcd[,,1][,j]),length(gcd[,,1][,j])),rev(gcd[,,1][,j]))
        shade <- getShade("blue","15")
        polygon(xx,yy,density=NA,col=shade)

      }
      abline(h = c(climit, -climit), lty = 2)
      if (nchain(x) > 1) {
        title(main = paste(varnames(x, allow.null = FALSE)[j], 
                " (", chanames(x, allow.null = FALSE)[k], ")", 
                sep = ""))
      }
      else {
        title(main = paste(varnames(x, allow.null = FALSE)[j], 
                sep = ""))
      }
    }
  }
	addScenLab()
  saveFig("fig.mcmc.geweke")
}

### PARTITION VARIANCE RHO
fig.variance.partitions <- function()
{
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.variance.partitions): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	rho <- A$mc$rho[(Burn+1):nrow(A$mc)]
	if (length(sort(unique(rho)))==1) {
		cat("WARNING (fig.variance.partitions): Rho is singular (probably not estimated)\n"); return(invisible("No MCMC data")) }
	varphi <- 1/A$mc$varphi[(Burn+1):nrow(A$mc)]
	sig <- rho*varphi
	tau <- (1-rho)*varphi
	d <- cbind(vp=varphi,sig,tau)
	pairs(d, pch=".", upper.panel=NULL, gap=0)
	addScenLab()
	saveFig("fig.variance.partitions")
}

fig.time.varying.selectivity <- function(gear=1){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(A$delaydiff==0){
			with(A, {
				plot.sel<-function(x, y, z, ...){
					z <- z/max(z)
					z0 <- 0 #min(z) - 20
					z <- rbind(z0, cbind(z0, z, z0), z0)
					x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
					y <- c(min(y) - 1e-10, y, max(y) + 1e-10)
					clr=colorRampPalette(c("honeydew","lawngreen"))
					nbcol=50
					iclr=clr(nbcol)
					nrz <- nrow(z)
					ncz <- ncol(z)
					zfacet <- z[-1, -1]+z[-1, -ncz]+z[-nrz, -1]+z[-nrz, -ncz]
					facetcol <- cut(zfacet, nbcol)
					fill <- matrix(iclr[facetcol],nr=nrow(z)-1,nc=ncol(z)-1)
					fill[ , i2 <- c(1,ncol(fill))] <- "white"
					fill[i1 <- c(1,nrow(fill)) , ] <- "white"

					par(bg = "transparent")
					persp(x, y, z, theta = 35, phi = 25, col = fill, expand=5, 
						shade=0.75,ltheta=45 , scale = FALSE, axes = TRUE, d=1,  
						xlab="Year",ylab="Age",zlab="Selectivity", 
						ticktype="detailed", ...)
				}
				ix <- 1:length(yr)
			plot.sel(yr, age, exp(log_sel[log_sel[,1]==gear,-1]),main=paste("Gear", gear))
		})
		addScenLab()
		if(gear==1){ # WARNING (fig.time.varying.selectivity) - assumes gears laid out a certain way !!!
			saveFig("fig.time.varying.comm.sel")
		} else if(gear==2){
			saveFig("fig.time.varying.surv.sel")
		} 
	}else cat("WARNING (fig.time.varying.selectivity): No estimated selectivity for delay difference model\n")
}

### RUN TIME STATS
plotRuntimeStats <- function(type=0,ylab=""){
	# plots runtime stats for all scenarios.
	# assumes all scenarios have the same number of projected years
	# assumes the opList structure is used.
	# types:
	# 1 = objFun, 2 = max gradient, 3 = number of function evals, 4 = hangcode, any other value = exit code
	op <- par(no.readonly=T)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))

	if(type==1 | type==2 | type==3){ # use PLOTBUBBLES from PBSModelling for these ones
		if(type==1){
			# Objective function values
			# GREEN means value is positive GOOD
			# RED means a bad objective function value, i.e. returned -1,#IND (which is represented as 0.0 from GrMPE)
			if (is.null(opList[[1]][[4]]$ObjectiveFunction)){
				frame(); addLabel(0.5,0.5,"No ObjectiveFunction value available",cex=2,col="red"); return() }
			dat = numeric()
			for(scenario in 1:length(opList)){
					dat <- rbind(dat,abs(opList[[scenario]][[4]]$ObjectiveFunction))
			}
			dat <- t(dat)
			colnames(dat) <- 1:length(opList)
			rownames(dat) <- ""
			plotBubbles(dat,dnam=F,cpro=F,ylab=ylab,clrs=c("green4","red","black"),xaxt='n',yaxt='n',lwd=2)
			text(1:length(opList),1.04,dat,srt=-45,adj=1)
			text(1:length(opList),1,1:length(opList))
			title("Objective function values")
		} else if(type==2){
			# Maximum Gradient values
			# GREEN represents a good gradient, i.e. one that is smaller than .maxGrad
			# RED represents anything greater than .MAXGRAD
			if (is.null(opList[[1]][[4]]$MaxGrad)){
				frame(); addLabel(0.5,0.5,"No MaxGrad value available",cex=2,col="red"); return() }
			dat = numeric()
#browser();return()
			for(scenario in 1:length(opList)){
				dat <- rbind(dat,opList[[scenario]][[4]]$MaxGrad)
			}
			dat <- t(dat)
			colnames(dat) <- 1:length(opList)
			rownames(dat) <- ""
			dat <- ifelse(dat>.MAXGRAD,0,dat)
			dat <- ifelse(dat<.MAXGRAD,dat,-dat)
			plotBubbles(dat,dnam=F,cpro=F,ylab=ylab,clrs=c("green4","red","red"),xaxt='n',yaxt='n',lwd=2) 
			text(1:length(opList),1.04,dat,srt=-45,adj=1)
			text(1:length(opList),1,1:length(opList))
			title(paste("Maximum gradient values (<",.MAXGRAD,")"))
		} else if(type==3){
			# Number of function evaluations
			# GREEN means the number of function evaluations was greater than .FUNEVALS
			# RED means the number of function evaluations was less than .FUNEVALS
			if (is.null(opList[[1]][[4]]$nf)){
				frame(); addLabel(0.5,0.5,"No nf value available",cex=2,col="red"); return() }
			dat = numeric()
			for(rep in 1:length(opList)){
				dat <- rbind(dat,opList[[rep]][[4]]$nf)
			}
			dat <- t(dat)
			colnames(dat) <- 1:length(opList)
			rownames(dat) <- ""
			dat <- ifelse(dat<.FUNEVALS,-dat,dat)
			plotBubbles(dat,dnam=F,cpro=F,ylab=ylab,clrs=c("green","red","black"),xaxt='n',yaxt='n')
			text(1:length(opList),1.04,dat,srt=-45,adj=1)
			text(1:length(opList),1,1:length(opList))
			title(paste("Number of function evaluations (>",.FUNEVALS,")"))
		}
	} else {
	# NOT USING PLOTBUBBLES!!
		if(type==4){
			# Hang codes
			if (is.null(opList[[1]][[4]]$HangCode)){
				frame(); addLabel(0.5,0.5,"No Hang Code value available",cex=2,col="red"); return() }
			plotcharCol <- ifelse(opList[[1]][[4]]$HangCode==1,"red","green")
			# GREEN means no error condition
			# RED means no improvement in function value when 10th to last value compared with
			#		 current value.
		}else{
			# Exit codes
			if (is.null(opList[[1]][[4]]$ExitCode)){
				frame(); addLabel(0.5,0.5,"No Exit Code value available",cex=2,col="red"); return() }
			plotcharCol <- ifelse(opList[[1]][[4]]$ExitCode==1,"green","red")
			# GREEN for normal exit - i.e. all derivatives satisfy conditions
			# RED for problem with the initial estimate for the Hessian matrix.
			#		 - The hessian matrix must be positive definite
			plotcharCol <- ifelse(opList[[1]][[4]]$ExitCode==2,"blue",plotcharCol)
			# BLUE for problem with the derivatives, either:
			# a) There is an error in the derivatives or
			# b) function does not decrease in direction of search, perhaps due to numerical
			#		round off error, or too stringent a convergence criterion
			plotcharCol <- ifelse(opList[[1]][[4]]$ExitCode==3,"purple",plotcharCol)
			# PURPLE for Maximum number of function calls exceeded
		}
#	par( oma=c(2,2,4,1), mar=c(3,3,3,1), mfrow=c(1,1) )
		plot(1,1,
			 pch=.PCHCODE,
			 xlab="Scenario",
			 ylab=ylab,
			 col=plotcharCol,
			 xlim=c(1,length(opList)),
			 ylim=c(1,1))
		for(rep in 2:length(opList)){
			if(type==4){
				plotcharCol <- ifelse(opList[[rep]][[4]]$HangCode==1,"red","green")
				title("Hang code values")
			}else{
				plotcharCol <- ifelse(opList[[rep]][[4]]$ExitCode==1,"green","red")
				plotcharCol <- ifelse(opList[[rep]][[4]]$ExitCode==2,"blue",plotcharCol)
				plotcharCol <- ifelse(opList[[rep]][[4]]$ExitCode==3,"purple",plotcharCol)
				title("Exit code values")
			}
			points(rep,1,pch=.PCHCODE,col=plotcharCol)
		}
	}
	
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#New figures for Pacific cod performance measures
### CONTROL POINTS
fig.Allcontrol.pts.Box <- function(tac.use=0, useHRP=FALSE, minYr, aveYr) ## WAP Outside TAC=3110
{
	#This function takes the output from ADMB projections and plots the distribution of alternative control points
	#Control points are benchmark variables that do not change with tac
	#The mcmc data files has nrow =ntac x nctlpts	 only take the first tac level (tac=0) as these control points are insensitive to future tac (all based on historical catch/biomass/F)
	simplify = usingSweave
	sideways = TRUE
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())
	if (!isScenMcmc()) {
		if (simplify){
			frame(); addLabel(.5,.5,"No MCMC",col="red",cex=2); return()
		} else {
			cat("WARNING (fig.Allcontrol.pts.Box): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data"))
		}
	}
	ddmcmc <- A$mcproj
	tacs = sort(unique(ddmcmc$TAC))
	ntac = which(abs(tacs-tac.use)==min(abs(tacs-tac.use))) ## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	tac.use = tacs[ntac]
#browser();return()
	dmcmc <- subset(ddmcmc, TAC==tac.use)   #take only the first tac from each posterior sample - the control points do not change with tac
	nsamp <- length(dmcmc[,1])
	dmcmc <- as.data.frame(window(mcmc(dmcmc),start=Burn+1,thin=Thin))
	mcnames = names(ddmcmc)

	if (useHRP) {
		##********** change to using calcHRP() eventually
		btmcmc  = A$mc.sbt[,1:nyear]
		ftmcmc  = A$mc.ft[,1:nyear]
		post.bt = as.data.frame(window(mcmc(btmcmc),start=Burn+1,thin=Thin))
		post.ft = as.data.frame(window(mcmc(ftmcmc),start=Burn+1,thin=Thin))
		post.dt = t(apply(post.bt,1,function(x){x/mean(x[match(aveYr,yr)])}))
		ddmpd   = A$mpdproj

		unpackList(renderVals(minYr=minYr, aveYr=aveYr))
#browser();return()
		ctls.hrp = list()
		dmcmc[,"Bmin"]  = apply(post.bt[,match(minYr,yr),drop=FALSE],1,min)
		dmcmc[,"2Bmin"] = 2 * dmcmc[,"Bmin"]
		dmcmc[,"Bavg"] = apply(post.bt[,match(aveYr,yr),drop=FALSE],1,mean)
		dmcmc[,"Favg"] = apply(post.ft[,match(aveYr,yr),drop=FALSE],1,mean)
		dmcmc[,"Fmin"] = apply(post.ft[,match(minYr,yr),drop=FALSE],1,min)
		dmcmc[,"Fmax"] = apply(post.ft,1,function(x){max(x)})                     ## 1000 maxima over N years
		#ctls.hrp[["Davg"]] = apply(post.dt[,match(aveYr,yr),drop=FALSE],1,mean)  ## all values = 1
		mcnames = names(dmcmc)
		nmB = c(paste0("B",A$currYr+c(0,3)), "Bavg", "Bmin", "2Bmin")
		Bctl = match(nmB,mcnames)
		nmF = c(paste0("F",A$currYr+c(-1)), "Favg", "Fmin", "Fmax")
		Fctl = match(nmF,mcnames)
	} else {
		dmcmc[,"BMSY80"] = dmcmc[,match("BMSY",mcnames)]*0.80
		dmcmc[,"BMSY40"] = dmcmc[,match("BMSY",mcnames)]*0.40
		dmcmc[,"FMSY80"] = dmcmc[,match("FMSY",mcnames)]*0.80
		dmcmc[,"FMSY40"] = dmcmc[,match("FMSY",mcnames)]*0.40
		mcnames = names(dmcmc)
		nmB = c(paste0("B",A$currYr+c(0,3)),paste0("BMSY",c("","80","40")),paste0("B0",c("","40","20")))
		Bctl = match(nmB,mcnames)
		nmF = c(paste0("F",A$currYr+c(-1)),paste0("FMSY",c("","80","40")))
		Fctl = match(nmF,mcnames)
	}
	
	mcmcDataB <- dmcmc[,Bctl]
	mcmcDataF <- dmcmc[,Fctl]
#browser();return()

	if (!sideways) par(mfrow=c(2,1), mar=c(2,5,2,2), mgp=c(2.5,0.5,0))
	for(i in 1:2){
		if(i==1) {
			#if (simplify) par(plt=c(0.15,0.99,0.55,0.99), new=FALSE)
			if (simplify) par(plt=c(0.10,0.55,0.20,0.90), new=FALSE)
			else if (sideways) par(plt=c(0.10,0.55,0.15,0.95), mgp=c(2.25,0.5,0), new=FALSE)
			if (useHRP) colnames(mcmcDataB) = gsub("min",minYr,colnames(mcmcDataB))
			quantBox(mcmcDataB/1000, names=colnames(mcmcDataB), range=0.95, pch=".", col=ifelse(sideways,"green","honeydew"), las=1, xlab="", outline=FALSE, ylab=ifelse(simplify,"","Biomass (1000 t)"), main=ifelse(sideways,"","Biomass-based control points"),
				ylim=c(0,0.001*max(sapply(mcmcDataB,quantile,0.95))), boxwex=0.5, cex.lab=1.2, xaxt=ifelse(sideways,"n","s"))#,cex.axis=0.8) #
			if(simplify || sideways) {
				if(par()$mfg[1]==par()$mfg[3]){
					if (useHRP)
						Blab = paste0("expression(", paste0(c(rep("",4), "2*"), paste0("italic(B)", c(paste0("[",A$currYr,"]"), paste0("[",A$currYr+3,"]"), rep("[avg]",1), rep(paste0("[",minYr,"]"),2)) ) ), ")")
					else
						Blab = paste0("expression(", paste0(c(rep("",3), "0.8~", "0.4~", "", "0.4~", "0.2~"), paste0("italic(B)", c(paste0("[",A$currYr,"]"), paste0("[",A$currYr+3,"]"), rep("[MSY]",3), rep("[0]",3)))), ")")
					for (l in 1:length(Blab))
						eval(parse(text=paste0("axis(1,at=",l,",",Blab[l],",las=2,cex.axis=1.5,padj=0.3)")))
				}
				addLabel(0.05,0.05,"B",adj=c(0,0),cex=1.2)
			}
		} else {
			#if (simplify) par(plt=c(0.15,0.99,0.07,0.47), new=TRUE)
			if (simplify) par(plt=c(0.65,0.99,0.20,0.90), new=TRUE)
			else if (sideways) par(plt=c(0.65,0.99,0.15,0.95), new=TRUE)
			if (useHRP) colnames(mcmcDataF) = gsub("min",minYr,colnames(mcmcDataF))
			quantBox(mcmcDataF, names=colnames(mcmcDataF), pch=".", range=0.95, col=ifelse(sideways,"cadetblue1","aliceblue"), las=1, xlab="", outline=FALSE, ylab=ifelse(simplify,"","Fishing mortality (/y)"), main=ifelse(sideways,"","F-based control points"),
				ylim=c(0,max(sapply(mcmcDataF,quantile,0.95))), boxwex=ifelse(sideways,0.5,0.3), cex.lab=1.2, xaxt=ifelse(sideways,"n","s"))
#browser();return()
			if(simplify || sideways) {
				if(par()$mfg[1]==par()$mfg[3]){
					if (useHRP)
						Flab = paste0("expression(", paste0(c(rep("",3)), paste0("italic(F)", c(paste0("[",A$currYr-1,"]"), "[avg]", paste0("[",minYr,"]"), "[max]") ) ), ")")
					else
						Flab = paste0("expression(", paste0(c(rep("",2), "0.8~", "0.4~"), paste0("italic(F)", c(paste0("[",A$currYr-1,"]"), rep("[MSY]",3)))), ")")
					for (l in 1:length(Flab))
						eval(parse(text=paste0("axis(1,at=",l,",",Flab[l],",las=2,cex.axis=1.5,padj=0.3)")))
				}
				#addLabel(0.05,0.05,"F",adj=c(0,0),cex=1.2)
				addLabel(0.95,0.05,"F",adj=c(1,0),cex=1.2)
			}
		}
	}
	if(!simplify) saveFig("fig.ctlpt.box")
}

### MSY-BASED CONTROL POINTS
fig.MSYcontrol.pts <- function()
{
	#This function takes the output from ADMB projections and plots the distribution of alternative control points
	#Control points are benchmark variables that do not change with tac
	#The mcmc data files has nrow =ntac x nctlpts	 only take the first tac level (tac=0) as these control points are insensitive to future tac (all based on historical catch/biomass/F)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.MSYcontrol.pts): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	ddmcmc <- A$mcproj
	ddmpd  <- A$mpdproj
	#ctlpts <-c(8,12) ## columns of projection file that contain control points
	dmcmc  <- subset(ddmcmc, TAC==0)	   #take only the first tac from each posterior sample - the control points do not change with tac
	ctlpts = match(c("BMSY","FMSY"), names(dmcmc))
	ctlpts.mpd = match(c("BMSY","FMSY"), names(ddmpd)); names(ctlpts.mpd)=ctlpts
	#colour = 0
	colours= c("green4","blue")
	rc = .findSquare(length(ctlpts))
	par(mfrow=rc, mar=c(3,5,2,1), mgp=c(3.2,0.5,0))

	for(i in ctlpts) {
		ii = as.character(i)
		iii = grep(i,ctlpts)
		#colour = colour+1
		colour = colours[iii]
		shade  = getShade(colour,20)

		mcmcData <- window(mcmc(dmcmc[,i]), start=Burn+1, thin=Thin)
#browser();return()
		#dens<-density(mcmcData)
		dens = hist(mcmcData, nclass=100, col=shade, main=colnames(ddmcmc)[i])
		Med <- quantile(mcmcData, na.rm=TRUE, 0.5) ## Use quantiles to get the median. The median calculated from density is a bit off
		mpdCtlPts <- ddmpd[1,ctlpts.mpd[ii]] ## take only the first tac for the MPD
		abline(v=mpdCtlPts, col=2,lwd=2,lty=2)
		text(mpdCtlPts,0.95*par()$usr[4],"MPD",cex=0.7,col="red",pos=4)
	}
	for(i in length(ctlpts)) saveFig(paste("fig.ctlpt.msy.freq",i, sep=""))
}

### HISTORICAL CONTROL POINTS
fig.Histcontrol.pts <- function(minYr, aveYr)
{
	#This function takes the output from ADMB projections and plots the distribution of alternative control points
	#Control points are benchmark variables that do not change with tac
	#The mcmc data files has nrow =ntac x nctlpts	 only take the first tac level (tac=0) as these control points are insensitive to future tac (all based on historical catch/biomass/F)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.Histcontrol.pts): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	btmcmc  = A$mc.sbt[,1:nyear]
	ftmcmc  = A$mc.ft[,1:nyear]
	post.bt = as.data.frame(window(mcmc(btmcmc),start=Burn+1,thin=Thin))
	post.ft = as.data.frame(window(mcmc(ftmcmc),start=Burn+1,thin=Thin))
	post.dt = t(apply(post.bt,1,function(x){x/mean(x[match(aveYr,yr)])}))
	ddmpd   = A$mpdproj

	unpackList(renderVals(minYr=minYr, aveYr=aveYr))
	ctls.hrp = list()
#browser(); return()
	ctls.hrp[["Bmin"]] = apply(post.bt[,match(minYr,yr),drop=FALSE],1,min)
	ctls.hrp[["Bavg"]] = apply(post.bt[,match(aveYr,yr),drop=FALSE],1,mean)
	ctls.hrp[["Fmin"]] = apply(post.ft[,match(minYr,yr),drop=FALSE],1,min)
	ctls.hrp[["Favg"]] = apply(post.ft[,match(aveYr,yr),drop=FALSE],1,mean)
	#ctls.hrp[["Davg"]] = apply(post.dt[,match(aveYr,yr),drop=FALSE],1,mean)  ## all values = 1
	ctlpts = match(names(ctls.hrp), names(ctls.hrp))
	mmY =  median(minYr) ## median year of minimum biomass


	colours= c("green4","blue","red","gold")
	rc = .findSquare(length(ctlpts))
	par(mfrow=rc, mar=c(3,3,2,1), mgp=c(3.2,0.5,0))

	for(i in ctlpts) {
		ii = names(ctls.hrp)[i]
		iexp  = sub("D","italic(B[.(mmY)]/B)", sub("F","italic(F)", sub("B","italic(B)", sub("avg","[avg]", sub("min","[min]",ii)))))
		iii = grep(i,ctlpts)
		#colour = colour+1
		colour = colours[iii]
		shade  = getShade(colour,20)

		mcmcData <- ctls.hrp[[i]]
		hist(mcmcData, nclass=100, col=shade, main=eval(parse(text=paste0("expression(",iexp,")"))))
		medCtlPts <- median(mcmcData)
		abline(v=medCtlPts, col=2,lwd=2,lty=2)
		text(medCtlPts,0.95*par()$usr[4],as.expression(bquote(Q[50]==.(round(medCtlPts,2)))),cex=0.8,col="red",pos=4)
	}
#browser();return()
	for(i in length(ctlpts)) saveFig(paste("fig.ctlpt.hrp.freq",i, sep=""))
}
fig.Histcontrol.pts.old <- function()
{
	#This function takes the output from ADMB projections and plots the distribution of alternative control points
	#Control points are benchmark variables that do not change with tac
	#The mcmc data files has nrow =ntac x nctlpts	 only take the first tac level (tac=0) as these control points are insensitive to future tac (all based on historical catch/biomass/F)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.Histcontrol.pts): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	ddmcmc <- A$mcproj
	if (!is.element("Bmin",names(ddmcmc))) {
		cat("WARNING (fig.Histcontrol.pts): MCMCs do not contain historical reference points\n"); return(invisible("No historical ref pts")) }
	ddmpd  <- A$mpdproj
	#ctlpts <-c(14,16,18,20,22) #columns of projection file that contain control points
	dmcmc  <- subset(ddmcmc, TAC==0)	   #take only the first tac from each posterior sample - the control points do not change with tac
	ctlpts = match(c("Bmin","BAvg_S","FAvg_S","BAvg_L","FAvg_L"), names(dmcmc))
	ctlpts.mpd = match(c("BMSY","FMSY"), names(ddmpd)); names(ctlpts.mpd)=ctlpts
	colour = 0
	rc = .findSquare(length(ctlpts))
	par(mfrow=rc, mar=c(3,5,2,1), mgp=c(4,0.5,0))
	
	for(i in ctlpts) {
		ii = as.character(i)
		colour=colour+1
		mcmcData <- window(mcmc(dmcmc[,i]), start=Burn+1, thin=Thin)
		dens<-density(mcmcData)
		Med <- quantile(mcmcData, na.rm=T,0.5)  #Use quantiles to get the median. The median calculated from density is a bit off
		#Median <- median(mcmcData, na.rm=T)    #test
		mpdCtlPts <- ddmpd[1,ctlpts.mpd[ii]] #take only the first tac for the MPD

		plot(dens,main=colnames(ddmcmc)[i], xlab=colnames(ddmcmc)[i], ylab="Density",las=1)
		 xx <- c(dens$x,rev(dens$x))
		 yy <- c(rep(min(dens$y),length(dens$y)),rev(dens$y))
		 shade <- getShade(1,20)
		polygon(xx,yy,density=NA,col=shade)
		abline(v=mpdCtlPts, col=2,lwd=2,lty=2)
	}
	for(i in length(ctlpts)) saveFig(paste("fig.ctlpt.hist.density",i, sep=""))
}

### BENCHMARKS
fig.Benchmarks <- function()
{
	#This function takes the output from ADMB projections and plots the distribution of alternative control points
	#Control points are benchmark variables that do not change with tac
	#The mcmc data files has nrow =ntac x nctlpts	 only take the first tac level (tac=0) as these control points are insensitive to future tac (all based on historical catch/biomass/F)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.Benchmarks): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	ddmcmc <- A$mcproj 
	ddmpd  <- A$mpdproj
	ctlpts <-c(2,5) #columns of projection file that contain control points
	dmcmc  <- subset(ddmcmc, TAC==0)	   #take only the first tac from each posterior sample - the control points do not change with tac
	ctlpts = match(c(paste0("B",A$currYr),paste0("F",A$currYr-1)), names(dmcmc))
#browser();return()
	ctlpts.mpd = match(c(paste0("B",A$currYr),paste0("F",A$currYr-1)), names(ddmpd)); names(ctlpts.mpd)=ctlpts
	#colour = 0
	colours= c("green4","blue")
	rc = .findSquare(length(ctlpts))
	par(mfrow=rc, mar=c(3,5,2,1), mgp=c(4,0.5,0))

	for(i in ctlpts) {
		ii = as.character(i)
		iii = grep(i,ctlpts)
		#colour = colour+1
		colour = colours[iii]
		shade  = getShade(colour,20)

		mcmcData <- window(mcmc(dmcmc[,i]), start=Burn+1, thin=Thin)
		#dens<-density(mcmcData)
		dens = hist(mcmcData, nclass=100, col=shade, main=colnames(ddmcmc)[i])
		Med <- quantile(mcmcData, na.rm=T,0.5)	  	 #Use quantiles to get the median. The median calculated from density is a bit off
		#Median <- median(mcmcData, na.rm=T)	 #test
		mpdCtlPts <- ddmpd[1,ctlpts.mpd[ii]] #take only the first tac for the MPD
		#plot(dens,main=colnames(ddmcmc)[i], xlab=colnames(ddmcmc)[i], ylab="Density",las=1)
		#xx <- c(dens$x,rev(dens$x))
		#yy <- c(rep(min(dens$y),length(dens$y)),rev(dens$y))
		#shade <- getShade(colour,20)
		#polygon(xx,yy,density=NA,col=shade)
		abline(v=mpdCtlPts, col=2,lwd=2,lty=2)
		text(mpdCtlPts,0.95*par()$usr[4],"MPD",cex=0.7,col="red",pos=4)
	}
	for(i in length(ctlpts)) saveFig(paste("fig.ctlpt.benchmarks.density",i, sep=""))
}


######~~~~~~~~Unused in 2012 so far~~~~~~~~~~~###################################

fig.mcmc.diagnostics <- function(){
  #This function runs diagnostic plots for the posterior samples
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.mcmc.diagnostics): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
  if(length(unique(A$mc$abar))==1){
    # abar was not estimated, so don't include it
    if(length(unique(A$mc$gbar))==1){
      # gbar was not estimated, so don't include it
      np <- c(1:5,9)
    }else{
      np <- c(1:6,9)
    }
  }else{
    np <- c(1:7,9)
  }

	xyplot(post.samp[,np])
	acfplot(post.samp[,1:np])
	addScenLab()
	saveFig("fig.mcmc.diagnostics")
}

fig.fmsy.steepness <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.fmsy.steepness): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	par(mfrow=c(2,2), mai=c(0.75,0.751,0.15,0.15), omi=c(0.35,0.45,0.15,0.15), las=1)
	fmc=A$mc[,2]
	hmc=A$mc[,11]
	CRmc=4*hmc/(1-hmc)
	Mmc=A$mc[,3]

	hist(fmc,prob=T,breaks=30,xlab="Fmsy",ylab="Posterior density Fmsy",main="",col="moccasin", cex.axis=1.1, cex.lab=1.1)
	lines(density(fmc,na.rm=T),col="red")
	mtext(paste("Mean =",round(mean(fmc),2), " Median =", round(median(fmc),2), "SD =", round(sd(fmc),2)), side=3, line=-0.5, cex=0.8)
	addScenLab()

	hist(hmc,prob=T,breaks=30,xlab="h",ylab="Posterior density Steepness",main="", xlim=c(0.2,1),col="moccasin", cex.axis=1.1, cex.lab=1.1)
	lines(density(hmc,na.rm=T),col="red")
	mtext(paste("Mean =",round(mean(hmc),2), " Median =", round(median(hmc),2), "SD =", round(sd(hmc),2)), side=3, line=-0.5, cex=0.8)

	hist(CRmc,prob=T,breaks=30,xlab="CR",ylab="Posterior density Compensation Ratio",main="",col="moccasin", cex.axis=1.1, cex.lab=1.1)
	lines(density(CRmc,na.rm=T),col="red")
	mtext(paste("Mean =",round(mean(CRmc),2), " Median =", round(median(CRmc),2), "SD =", round(sd(CRmc),2)), side=3, line=-0.5, cex=0.8)
	
	hist(Mmc,prob=T,breaks=30,xlab="M",ylab="Posterior density M",main="",col="moccasin", cex.axis=1.1, cex.lab=1.1)
	lines(density(Mmc,na.rm=T),col="red")
	mtext(paste("Mean =",round(mean(Mmc),2), " Median =", round(median(Mmc),2), "SD =", round(sd(Mmc),2)), side=3, line=-0.5, cex=0.8)
	saveFig("fig.fmsy.steepness")
}

fig.spr.vs.management.target <- function(){
	#The relative spawning potential ratio (1-spr)/(1-spr.at.msy)  
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.spr.vs.management.target): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	spr <- A$mc.f40spr #read.table("tinss.f40spr",h=F)
	if (is.null(spr)) {
		cat("WARNING (fig.spr.vs.management.target): No spawner-per-recruit data available\n"); return(invisible("No SPR data")) }
	post.spr <- as.data.frame(window(mcmc(spr),start=Burn+1,thin=thin))
	sprci <- apply(post.spr,2,quantile,probs=c(0.025,0.5,0.975))
	matplot(A$yr,t(sprci),type="l",col=c(2,1,2),lty=c(3,1,3), lwd=2, pch=c(-1, 0, 1),ylim=c(0,max(sprci))
		,xlab="Year",ylab="SPR")
	abline(h=1)
	text(1980, 1, "Management target", pos=3)
	addScenLab()
}

fig.yields.4panel <- function(A,type){
	#plot the equilibrium yield curves  
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (is.null(A$equil)) {
		cat("WARNING fig.yields.4panel): Equilibrium data not available\n"); return(invisible("No equlibrium data")) }
	par(mfcol=c(2,2))
	#A$equil comes from the TINSS.rep file
	fe <- A$equil[, 1]
	ye <- A$equil[, 2]
	de <- A$equil[, 3]
	spr <- A$equil[, 4]
	
	plot(fe, ye, type="l", xlab="Fishing mortality (Fe)", ylab="Equilibrium yield"); gletter(1)
	addScenLab()
	plot(de, ye, type="l", xlab="Spawning depletion", ylab="Equilibrium yield", lty=2, col=2);gletter(2)
	plot(spr,ye, type="l", xlab="Spawning potential ratio", ylab="Equilibrium yield", lty=3, col=3);gletter(3)
	matplot(cbind(fe, de, spr), ye/max(ye)*100, type="l",xlab="Fe, depletion,  SPR",  ylab="Relative equilibrium yield")
	gletter(4)
	saveFig("fig.yields.4panel")
}

fig.yield.depletion.relrecuitment.spr <- function(){
	#Relationship between fishing mortlaity ~ yield,  recruitment,  SBe,  SPR
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (is.null(A$equil)) {
		cat("WARNING (fig.yield.depletion.relrecuitment.spr): Equilibrium data not available\n"); return(invisible("No equlibrium data")) }
	par(mfcol=c(2,2))
	fe <- A$equil[, 1]
	ye <- A$equil[, 2]
	de <- A$equil[, 3]
	spr <- A$equil[, 4]
	re <- A$equil[, 5]  
	ix <- c(min(which(ye==max(ye))),min(which(de<=0.4)) , min(which(spr<=0.4)))

	plot(fe, ye, type="l",xlab="", ylab="Equilibrium yield (million mt)", lwd=2)
	segments(fe[ix],0,fe[ix],ye[ix],lty=c(1, 2, 3))
	segments(0,ye[ix],fe[ix],ye[ix],lty=c(1, 2, 3))
	addScenLab()
	
	re <- re/re[1]
	plot(fe, re, type="l",xlab="", ylab="Relative recruitment", lwd=2) 
	segments(fe[ix],0,fe[ix],re[ix],lty=c(1, 2, 3)) 
	segments(0,re[ix],fe[ix],re[ix],lty=c(1, 2, 3)) 
	
	plot(fe, de, type="l",xlab="", ylab="Spawning depletion", lwd=2)
	segments(fe[ix],0,fe[ix],de[ix],lty=c(1, 2, 3))
	segments(0,de[ix],fe[ix],de[ix],lty=c(1, 2, 3))

	plot(fe, spr, type="l",xlab="", ylab="Spawning Potential Ratio", ylim=c(0, 1), lwd=2)
	segments(fe[ix],0,fe[ix],spr[ix],lty=c(1, 2, 3))
	segments(0,spr[ix],fe[ix],spr[ix],lty=c(1, 2, 3))   
	
	legend("topright", c("MSY", "SB40", "SPR40"), lty=1:3, bty="n")
	
	mtext("Equilibrium fishing mortality rate", 1, outer=T, line=-1)                               
	saveFig("fig.yield.depletion.rel.recruitment.spr")
}

gweke.chain <- function (x, frac1 = 0.1, frac2 = 0.5, nbins = Nbin, pvalue = 0.05, auto.layout = TRUE, ...){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (gweke.chain): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
    x <- as.mcmc.list(x)
    if (auto.layout) 
        oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x), 
            Nparms = nvar(x)))
    ystart <- seq(from = start(x), to = (start(x) + end(x))/2, 
        length = nbins)
    if (is.R()) 
        gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)), 
            dimnames = c(ystart, varnames(x), chanames(x)))
    else gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)), 
        dimnames = list(ystart, varnames(x), chanames(x)))
    for (n in 1:length(ystart)) {
        geweke.out <- geweke.diag(window(x, start = ystart[n]), 
            frac1 = frac1, frac2 = frac2)
        for (k in 1:nchain(x)) gcd[n, , k] <- geweke.out[[k]]$z
    }
    climit <- qnorm(1 - pvalue/2)
    for (k in 1:nchain(x)) for (j in 1:nvar(x)) {
        ylimit <- max(c(climit, abs(gcd[, j, k])))
        plot(ystart, gcd[, j, k], type = "p", xlab = "First iteration in segment", 
            ylab = "Z-score", pch = 4, ylim = c(-ylimit, ylimit), 
            ...)
        abline(h = c(climit, -climit), lty = 2)
        if (nchain(x) > 1) {
            title(main = paste(varnames(x, allow.null = FALSE)[j], 
                " (", chanames(x, allow.null = FALSE)[k], ")", 
                sep = ""))
        }
        else {
            title(main = paste(varnames(x, allow.null = FALSE)[j], 
                sep = ""))
        }
        if (k == 1 && j == 1) 
            oldpar <- c(oldpar, par(ask = ask))
    }
 }

getScenLab = function(scenario) {
	if (missing(scenario))
		scenlab = getWinVal()$scenarioHeader[getWinVal()$entryScenario,1]
	else
		scenlab = opList[[scenario]][[1]]
#	if (is.null(scenlab)) 
#		scenlab = paste0(getOptions(wapmod)$fdPrefix,pad0(baseRep,2))
#	else
	scenlab = basename(scenlab)
	#scenlab = gsub("/","",gsub("\\.\\./Scenarios/","",scenlab))
	return(scenlab)
}
isScenMcmc = function (){
	runstr = getWinVal()$scenarioHeader[getWinVal()$entryScenario,1]
	if (is.null(runstr)) 
		runstr = pad0(baseRep,2)
	else {
		runstr = basename(runstr)
		#runstr = gsub("/","",gsub("\\.\\./Scenarios/","",runstr))
		prefix = getOptions(wapmod,"fdPrefix")
		runstr = sub(prefix,"",runstr)
	}
#browser();return()
	is.mcmc = getOptions(wapmod,"is.mcmc")
	return(is.mcmc[runstr])
}
parEstimated = function(){
	gwval = getWinVal()
	if (length(gwval)!=0) {
		burn=gwval$burn
		zMCMC  = burn:nrow(A$mc)
	} else zMCMC = 1:nrow(A$mc)
	zMCMC = seq(zMCMC[1],rev(zMCMC)[1],Thin)
	nMCMC  = length(zMCMC)
#browser();return()
	estPar = c("log.ro","h","log.m","log.rbar","log.rinit","rho","varphi","qc",paste0("q",1:A$nits)) ## RH not all estimated
	j = A$mc[zMCMC,estPar]
	zPar = sapply(j,function(x){(length(sort(unique(x)))>1)})
	return(names(zPar)[zPar])
}
### AME's running quantile function
cquantile.vec = function (z, prob)
{
	cquant <- rep(NA, length(z))
	if (length(prob) != 1) 
		stop("length prob should be 1")
	for (i in 1:length(z)) {
		cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
	}
	return(cquant)
}
addScenLab = function(x=0.05,y=0.95, adj=c(0,0), scenario, ...) {
	addLabel(x, y, getScenLab(scenario), adj=adj, col="grey40", ...)
}

#plotSnail------------------------------2017-08-03
# Plot snail-trail plots for MCMC analysis.
# Use B_t and U_t-1
#-------------------------------------------AME/RH
plotSnail=function (BoverBmsy, UoverUmsy, p=c(0.1,0.9), xLim=NULL, yLim=NULL, Lwd=c(1,2,2), ## c(trace,refpts,final)
   ngear=1, Cnames="Commercial", png=FALSE, plotname="fig.snail", PIN=c(7,7), path="./",
   useHRP=FALSE, RPs=c(0.4,0.8), minYr=NULL)  ## only for plot labelling
{
	simplify = usingSweave && !png
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())
	BUlist = as.list(0:ngear); names(BUlist)=c("Spawning Biomass",Cnames[1:ngear])
	## In POP, we matched the current Bt with previous Ut (e.g., B2017 and U2016)
	#BUlist[[1]] = BoverBmsy[,-length(BoverBmsy)]
	BUlist[[1]] = BoverBmsy[,-1]  ## get rid of first year and keep last
#browser();return()
	if (ngear==1)
		BUlist[[2]] = UoverUmsy[,-length(UoverUmsy)]  ## get rid of last year and keep first (prob. should use 'ncol' but keep 'length' for now)
	else {
		for (g in 1:ngear) {
			gfile = UoverUmsy[,grep(paste0("_",g),names(UoverUmsy))]
			names(gfile) = substring(names(gfile),1,4)
			BUlist[[g+1]] = gfile[,-length(gfile)]
		}
	}
	# Calculate medians to be plotted
	BUmed  = sapply(BUlist,function(x){apply(x,2,median)},simplify=FALSE)  # median each year
	colPal = colorRampPalette(c("grey90", "grey30"))
	colSlime = rep(c("ivory3","slategray3"),ngear)[1:ngear]
	colStart = rep(c("dodgerblue","steelblue"),ngear)[1:ngear]
	colStop = rep(c("purple","darkred"),ngear)[1:ngear]
	nB = length(BUmed[[1]])
	if (is.null(xLim))
		xLim=c(0, max(c(BUmed[[1]], apply(BUlist[[1]],2,quantile,p[2]), 1)))
		#xLim=c(0, max(c(BUmed[[1]], quantile(apply(BUlist[[1]],2,quantile,p[2]),1), 1)))
	if (is.null(yLim))
		yLim=c(0, max(c(sapply(BUmed[(1:ngear)+1],max), sapply(BUlist[(1:ngear)+1],function(x,p){allyrs=apply(x,2,quantile,p); rev(allyrs)[1]},p=p[2]), 1)))
		#yLim=c(0, max(c(sapply(BUmed[(1:ngear)+1],max), quantile(sapply(BUlist[(1:ngear)+1],function(x,p){apply(x,2,quantile,p)},p=p[2]),1), 1)))
	if (!simplify) {
		if (png) png(filename=paste(path,plotname,".png",sep=""),width=PIN,height=PIN[2],units="in",res=300)
		expandGraph(mfrow=c(1,1),mar=c(3.2,3.5,0.5,0.5),mgp=c(2,0.5,0))
	}
	Lwd = rep(Lwd,3)[1:3]
	plot(0,0, xlim=xLim, ylim=yLim, type="n", 
		xlab = ifelse(simplify,"",ifelse(useHRP,expression(italic(B[t])/italic(B)[avg]),expression(italic(B[t])/italic(B)[MSY]))), 
		ylab = ifelse(simplify,"",ifelse(useHRP,expression(italic(u[t-1])/italic(u)[avg]),expression(italic(u[t-1])/italic(u)[MSY]))),
		cex.lab=1.5,cex.axis=1.2,las=1)
	abline(h=1, col=c("grey20"), lwd=Lwd[2], lty=3)
	abline(v=RPs, col=c("red","green4"), lwd=Lwd[2], lty=2)
	for (i in ngear:1) {
		lines(BUmed[[1]], BUmed[[i+1]], col=colSlime[i], lwd=Lwd[1])
		points(BUmed[[1]], BUmed[[i+1]], type="p", pch=19, col=colPal(nB),cex=ifelse(saveon,0.8,1.2))
		points(BUmed[[1]][1], BUmed[[i+1]][1], pch=19, col=colStart[i],cex=1.5)
		zmin = match(minYr,names(BUmed[[1]]))
		points(BUmed[[1]][zmin], BUmed[[i+1]][zmin], type="p", pch=19, col="red",cex=1.5)
#browser();return()
		segments(quantile(BUlist[[1]][,as.character(currYear)],p[1]),
			BUmed[[i+1]][as.character(currYear-1)],
			quantile(BUlist[[1]][,as.character(currYear)], p[2]), 
			BUmed[[i+1]][as.character(currYear-1)], col=colStop[i], lwd=Lwd[3])
		segments(BUmed[[1]][as.character(currYear)], 
			quantile(BUlist[[i+1]][, as.character(currYear-1)], p[1]),
			BUmed[[1]][as.character(currYear)], 
			quantile(BUlist[[i+1]][, as.character(currYear-1)], p[2]), col=colStop[i], lwd=Lwd[3])
		points(rev(BUmed[[1]])[1], rev(BUmed[[i+1]])[1], pch=19, col=colStop[i],cex=1.5)
	}
	if (ngear>1)  addLegend(0.88,0.97,legend=Cnames,lty=1,lwd=Lwd[2],col=colSlime,seg.len=4,xjust=1,bty="n",cex=0.8)
	addLegend(0.95,0.90,
		legend=c(paste0(c("LRP","USR")," = ",round(RPs,2)), paste0(c("start","end","min B")," year = ",c(names(BUmed[[1]])[c(1,length(BUmed[[1]]))],minYr))),
		lty=c(2,2,NA,NA,NA), lwd=c(Lwd[2],Lwd[2],1,1,1), pch=c(NA,NA,19,19,19), pt.cex=c(NA,NA,1.5,1.5,1.5),
		col=c("red","green4",colStart[1],colStop[1],"red"), seg.len=4, xjust=1, bty="n", cex=1)
	addScenLab(x=0.95,adj=1)
	box()
	if (!simplify && png) dev.off()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSnail

#fig.snail------------------------------2017-08-03
# Wrapper for function `plotSnail' to choose
# between MSY-based or HRP-based references.
#-------------------------------------------AME/RH
fig.snail = function(useHRP=FALSE, minYr, aveYr, ...)
{
	if (!isScenMcmc()) {
		cat("WARNING (fig.snail): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	val = getWinVal()
	tac.use = val$currTAC
	ddmcmc <- A$mcproj
	tacs = sort(unique(ddmcmc$TAC))
	ntac = which(abs(tacs-tac.use)==min(abs(tacs-tac.use))) ## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	tac.use = tacs[ntac]
	assign("currYear", A$currYr, envir=.GlobalEnv)
	assign("lastYear", A$currYr - 1, envir=.GlobalEnv)
	mcmc.prj  = mcmc2(subset(A$mcproj, TAC==tac.use),start=(Burn+1),thin=Thin)
	mcmc.sbt  = mcmc2(A$mc.sbt[,1:length(A$yr)],start=(Burn+1),thin=Thin)  ## RH: final year in A$mc.sbt is a pseudo-projection
	Bt.mcmc   = cbind(mcmc.sbt,Bcurr=mcmc.prj[,paste0("B",currYear)])
	names(Bt.mcmc) = A$yrs
	## RH: likely not appropriate for non-equilibrium start:
	Ut.mcmc   = mcmc2(cbind(U0=rep(0,nrow(A$mc.ft)),1.-exp(-A$mc.ft[,1:(length(A$yr)-1)])),start=(Burn+1),thin=Thin)
	Ut.mcmc   = cbind(Ut.mcmc,Ucurr=1.-exp(-mcmc.prj[,paste0("F",lastYear)]))
	## RH: Name with year of Bt even though it is behind 0.5 years
	names(Ut.mcmc) = A$yrs
	if (useHRP) {
		if (missing(minYr) || missing(aveYr))
			unpackList(renderVals(aveYr=getWinVal()$aveYr,minYr=getWinVal()$minYr))
		else
			unpackList(renderVals(minYr=minYr, aveYr=aveYr))

		## Calculate the depletion over the entire period
		HRPs = calcHRP(A, aveYr, minYr, Burn)
		Bmsy.mcmc  = HRPs$post.abt
		Umsy.mcmc  = HRPs$post.aut
		RPs    = c(HRPs$dLRP, HRPs$dUSR)
	} else {
		Bmsy.mcmc = mcmc.prj[,"BMSY"]
		Umsy.mcmc = mcmc.prj[,"UMSY"]
		RPs = c(0.4, 0.8)
	}
	BoverBmsy = Bt.mcmc/Bmsy.mcmc
	UoverUmsy = Ut.mcmc/Umsy.mcmc
#browser();return()
	plotSnail(BoverBmsy, UoverUmsy, useHRP=useHRP, RPs=RPs, minYr=minYr, ...)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fig.snail

#plotChains-----------------------------2015-10-30
plotChains=function (mcmc, nchains=3, pdisc=0, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
   tck=0.5, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
   lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
   col.quant="black", plot=TRUE, probs=c(0.025, 0.5, 0.975), ...)  # AME probs
{
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (plotChains): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	panel.trace <- function(x, y, ...)
	{
		panel.xyplot(x, y, type="n")
		chainlink=rep(1:nchains,f)
		for (i in 1:nchains) {
			z=is.element(chainlink,i)
			panel.xyplot(x[z], y[z], type="l", lty=lty.trace, lwd=2, col=rep(col.trace,nchains)[i])
			panel.text(min(x)+0.04*diff(range(x)),1,labels=ifelse(panel.number()==1,getScenLab(),""),adj=c(0,1),cex=0.8,col="grey40")
		}
		#panel.xyplot(x, y, type="l", lty=lty.trace, lwd=lwd.trace, col=col.trace)
	}

	mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	Iseries = attributes(cribtab)$indexSeries[eval(parse(text=(cribtab[getScenLab(),"I"])))]
	names(mcmc)[grep("^q[1-9]",names(mcmc))] = paste0("q (",Iseries,")")

#browser();return()
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	f=rep(round(n/nchains),nchains-1); f=c(f,n-sum(f))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,f),p), Value=as.vector(as.matrix(mcmc)))
	mess = c(
	"require(grid, quietly=TRUE, warn.conflicts=FALSE)",
	"require(lattice, quietly=TRUE, warn.conflicts=FALSE)"
	)
	eval(parse(text=mess))
	if (trellis.par.get()$background$col == "#909090") {
		for (d in dev.list()) dev.off()
		trellis.device(color=FALSE)
	}
	mymain <- list(label=main, cex=cex.main)
	myxlab <- list(label=xlab, cex=cex.lab)
	myylab <- list(label=ylab, cex=cex.lab)
	myrot <- switch(as.character(las), `0`=0, `1`=0, `2`=0, `3`=90)     # AME changed '0'=90 to 0
	myscales <- list(x=list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot), 
		y= list(draw=axes, relation=relation, cex=cex.axis, tck=tck, tick.number=tick.number, rot=myrot))
	mystrip <- list(cex=cex.strip)

	dat$Index=paste(dat$Factor,dat$Chain,sep="-")
	vList=split(dat$Value,dat$Index)
	qList=sapply(vList,function(x){
		xsort=sort(x)
		xscal=xsort - min(xsort)
		ycumu=cumsum(xscal)/sum(xscal)
		out=cbind(x=xsort,y=ycumu)
		return(out) }, simplify=FALSE )
	dat$CumFreq=dat$ValueSort=NA
	for (i in names(qList)) {
		z=is.element(dat$Index,i)
		dat$ValueSort[z]=qList[[i]][,"x"]
		dat$CumFreq[z]  =qList[[i]][,"y"]
	}
	graph <- xyplot(CumFreq ~ ValueSort  | Factor, panel=panel.trace, 
		data=dat, as.table=TRUE, between=between, main=mymain, 
		xlab=myxlab, ylab=myylab, par.strip.text=mystrip, 
		scales=myscales, ylim=c(0,1), ...)
#browser();return()
	if (plot) {
		print(graph)
		#invisible(dat)
	}
	else {
		invisible(graph)
	}
	saveFig("fig.three.chains")
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotChains

