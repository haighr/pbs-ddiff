#**********************************************************************************
# ccamSens.r
# This file contains the code necessary to plot sensitivities overlaid on top
# of one another by sensitivity group which is set in the file 'SensitivityGroup'
# in each Scenario folder.  Scenarios with the same sensitivity group will be plotted
# together.  The 'SensitivityGroup' file should have a single number in it, nothing more.
# The Ref model will be SensitivityGroup=0 and will be plotted against all other
# Sensitivity groups.
# 
# Assumes the opList has been built and contains the SensitivityGroup element:
# opList[[scenario]][[4]]$SensitivityGroup
#
# Author            : Chris Grandin
# Development Date  : December 2011 - January 2012
# Modified          : Rowan Haigh (2016-12-23)
#
#**********************************************************************************

.writeSensitivityGroups <- function(){
  val <- getWinVal()
  for(scenario in 1:length(opList)){
    filename <- paste(opList[[scenario]][[1]],fnSensitivityGroup,sep="")
    opList[[scenario]][[4]]$SensitivityGroup <<- val$scenarioHeader[scenario,2]
    write(val$scenarioHeader[scenario,2],filename)
  }
  cat("Saved the SensitivityGroup information.\n")
}

fig.base.vs.sens <- function(sensitivityGroup=NULL, whichPlot="biomass", useHRP=FALSE, minYr, aveYr, ylimit=6, useMaxYlim=T, lty=2, lwd=2, pch=20, offset=0.3, opacity="20")
{
	# plots Spawning stock biomiass, depletion, and recruitment for a given sensitivity group
	#  - plot the MCMC posterior data for each, with confidence limits
	#  - offset is the number of years to offset each vertical bar from each other.
	#  - opacity is a two digit string from 00-99 
	# whichPlot can be:
	# 1. "biomass"
	# 2. "depletion"
	# 3. "recruits"
	if (!isScenMcmc()) {
		.flush.cat("WARNING (fig.base.vs.sens): No MCMCs generated for this scenario\n\n"); return(invisible("No MCMC data")) }
	on.exit(expandGraph(mfrow=c(1,1), new=FALSE))
	base   <- 0
	color  <- 1
	colors <- color
	nocrib = !exists("cribtab",envir=.GlobalEnv)

	if (exists("usingSweave",envir=.GlobalEnv) && usingSweave && !nocrib){
		SGs = CRIBTAB$SG; names(SGs) = rownames(CRIBTAB)  ## need the original for proper indexing
	} else {
		val = getWinVal()
		SGs = val$scenarioHeader[,2]
		if (is.null(SGs) && nocrib){
			.flush.cat("WARNING (fig.base.vs.sens): No Sensitivtiy Groups found\n\n"); return(invisible("No Sensitivity Groups"))
		}
		names(SGs) = basename(val$scenarioHeader[,1])
		if (is.null(sensitivityGroup)) sensitivityGroup=val$entrySensitivityGroup
	}
	if (!any(is.element(SGs,sensitivityGroup))) {
		.flush.cat("WARNING (fig.base.vs.sens):\n\tChoose a sensitivity 'Group number' that appears in the Sensitivity Group column (upper right in GUI)\n\n")
		return(invisible("No matching Group number"))
	}
	#SGs = SGs[is.element(SGs,c(base,sensitivityGroup))]
	if (!any(SGs==0)) {
		.flush.cat(paste0("WARNING (fig.base.vs.sens):\n\tNo Reference Case (0) found in Sensitivtiy Groups\n\tComparisons are among scenarios in Sensitivity Group ",sensitivityGroup),"\n\n")
		SGs = SGs[is.element(SGs,sensitivityGroup)]
		ibase = grep(sensitivityGroup,SGs)[1]
		isens = grep(sensitivityGroup,SGs)[-1]
	} else {
		## Get index numbers for the base and sensitivity group
		ibase = grep(0,SGs)
		if (length(ibase)>1) {
			.flush.cat("WARNING (fig.base.vs.sens):\n\tMore than one Reference Case (0) specified.\n\tOnly the first will be used and the others ignored.\n\n")
			ibase = ibase[1]
		}
		isens = grep(sensitivityGroup,SGs)  ## numerical index of sensitivityGroup in SGs
	}
#browser();return()
	isens = isens[getOptions(wapmod,"is.mcmc")[isens]] ## use only MCMC sensitivities
	Burn  = opList[[ibase]][[4]]$mc.burn
	plotR <- TRUE
	baseRep <- opList[[ibase]][[4]]
	baseName = basename(opList[[ibase]][[1]]) #strsplit(opList[[ibase]][[1]],"/")[[1]][3]
	names(ibase) = baseName
	Nyears = length(baseRep$yr)
	if (nocrib) { ## not revised
		runNames <- paste0("Ref: ",baseName) ## Gets the run's folder name out of the reletive path
		sensNum = 0:31
		oruns   = c(5,12,1,6,13:18,20,19,21,2:4,11,7,10,9,8,22:31)        ## order of the runs starting with base/ref (make sure this matches `oruns' in `modelResults.Rnw')
		names(sensNum) = paste0("assess",pad0(oruns,2))
		senscols = rep("white",length(oruns))
		ncols = grep("^0|99",sensNum,invert=TRUE)
		scols = rainbow(length(ncols),end=0.75)
		senscols[ncols] = scols
		senscols[1]= "gainsboro"
		names(senscols) = c(runNames,names(sensNum))
	} else {
		#cribtab = cribtab[!is.na(cribtab$sens),]
		runNames = paste0(cribtab[baseName,c("Slab","label")],collapse=": ")
		sensOrd  = order(cribtab$sens)
		sensNum  = cribtab$sens; names(sensNum) = row.names(cribtab); sensNum = sensNum[sensOrd]
		oruns    = cribtab$run; names(oruns) = row.names(cribtab); oruns = oruns[sensOrd]
		sensNum[is.na(sensNum)] = 99
		senscols = cribtab$col[sensOrd]; names(senscols) = row.names(cribtab)[sensOrd]
	}
	if (any(grepl(0,SGs)))
		senscols[names(SGs)[grep(0,SGs)]] = "gainsboro"
	#names(senscols) = oruns     ## label sensitivity colours with actual run number
	osens = NULL
	for (i in isens)
		osens = c(osens, basename(opList[[i]][[1]])) #strsplit(opList[[i]][[1]],"/")[[1]][3] )
	names(isens) = osens
#whichPlot="refpts"
#browser();return()
###########################
### BIOMASS  
	if (whichPlot == "biomass") {
		#Get the maximum year from the longest time series for xlim - i.e. might be comparing results from different lengths of data
		#Need to do this here so the xlim is long enough
		sensyr <- baseRep$yr
		mcbt <- baseRep$mc.sbt[,1:Nyears]/1000
		post.bt <- as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
		MaxBt <- apply(post.bt,2,quantile,probs=c(0.025,0.5,0.975))
		MaxBt <- max(MaxBt)
#browser();return()
		for(scenario in isens){
			#Get max xlim
			maxnyr <- length(baseRep$yr)
			yrtest <- opList[[scenario]][[4]]$yr
			if(length(yrtest) > maxnyr) {
				maxnyr <- length(opList[[scenario]][[4]]$yr)
				sensyr <- opList[[scenario]][[4]]$yr
			}
			#Now get max ylim
			iBurn  = opList[[i]][[4]]$mc.burn
			mcbt <- opList[[scenario]][[4]]$mc.sbt[,1:Nyears]/1000
			post.bt <- as.data.frame(window(mcmc(mcbt),start=iBurn+1,thin=Thin))
			ScBt <- apply(post.bt,2,quantile,probs=c(0.025,0.5,0.975))
			if(max(ScBt) > MaxBt) MaxBt <- max(ScBt)
		}
		mcbo <- baseRep$mc$bo/1000
		post.bo <- as.data.frame(window(mcmc(mcbo),start=Burn+1,thin=Thin))
		boci <- apply(post.bo,2,quantile,probs=c(0.025,0.5,0.975))

		mcbt <- baseRep$mc.sbt[,1:Nyears]/1000
		post.bt <- as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
		btci <- apply(post.bt,2,quantile,probs=c(0.025,0.5,0.975))
#browser();return()
		if(usingSweave || val$maxBiomassSensYlim){
			yUpperLimit <- 1.2*MaxBt
		} else {
			yUpperLimit <- val$biomassSensYlim
		}
		matplot(baseRep$yr,
			t(btci),
			type="l",
			col=color,
			lty=c(2,1,2),
			lwd=2,
			ylim=c(0,yUpperLimit),
			xlim=c(min(sensyr),max(sensyr)),
			xlab="Year",
			ylab="Biomass",
			main="Biomass",
			las=1)

		# Shade the confidence interval
		xx <- c(baseRep$yr,rev(baseRep$yr))
		yy <- c(btci[1,],rev(btci[3,]))
		shade <- getShade(color,opacity)
		polygon(xx,yy,density=NA,col=shade)
		# End shade the confidence interval

		points(baseRep$yr[1]-0.8,boci[2],col=color,pch=1)
		arrows(baseRep$yr[1]-0.8,boci[1],baseRep$yr[1]-0.8,boci[3],col=color, code=0, lwd=1.5)
		currOffset <- offset

		for(i in isens){
			iBurn  = opList[[i]][[4]]$mc.burn
			mcbo <- opList[[i]][[4]]$mc$bo/1000
			post.bo <- as.data.frame(window(mcmc(mcbo),start=iBurn+1,thin=Thin))
			boci <- apply(post.bo,2,quantile,probs=c(0.025,0.5,0.975))
			mcbt <- opList[[i]][[4]]$mc.sbt[,1:Nyears]/1000
			post.bt <- as.data.frame(window(mcmc(mcbt),start=iBurn+1,thin=Thin))
			btci <- apply(post.bt,2,quantile,probs=c(0.025,0.5,0.975))
			color <- color + 1
			colors <- c(colors,color)

			par(new=T)
			matplot(opList[[i]][[4]]$yr,
				t(btci),
				type="l",
				col=color,
				lty=c(2,1,2),
				lwd=2,
				xlim=c(min(sensyr),max(sensyr)),
				ylim=c(0,yUpperLimit),
				xlab="",
				ylab="",
				las=1)

			# Shade the confidence interval
			xx <- c(opList[[i]][[4]]$yr,rev(opList[[i]][[4]]$yr))
			yy <- c(btci[1,],rev(btci[3,]))
			shade <- getShade(color,opacity)
			polygon(xx,yy,density=NA,col=shade)
			# End shade the confidence interval

			par(new=F)
			points(opList[[i]][[4]]$yr[1]-0.8+currOffset,boci[2],col=color,pch=1)
			arrows(opList[[i]][[4]]$yr[1]-0.8+currOffset,boci[1],opList[[i]][[4]]$yr[1]-0.8+currOffset,boci[3],col=color, code=0, lwd=1.5)
			currOffset <- currOffset + offset
			#runNames <- c(runNames,paste0("Sens ",i,": ",basename(opList[[i]][[1]]))) #strsplit(opList[[i]][[1]],"/")[[1]][3]))
			runNames <- c(runNames,paste0(cribtab[i,c("Slab","label")],collapse=": "))
		}
#browser();return()
		lty <- c(rep(1,(length(runNames))))
		legend("topright", runNames, lty=lty, col=colors, bty="n", lwd=2)
		
		filename <- paste("fig.sensi.group",sensitivityGroup,".biomass",sep="")
		saveFig(filename)
###########################
### DEPLETION
	} else if (whichPlot == "depletion") {
#browser();return()
		## Get the maximum year from the longest time series for xlim --
		## i.e. might be comparing results from different lengths of data
		## Need to do this here so the xlim is long enough
		for(i in isens){
			maxnyr <- length(baseRep$yr)
			yrtest <- opList[[i]][[4]]$yr
			if(length(yrtest) > maxnyr) {
				maxnyr <- length(opList[[i]][[4]]$yr)
				sensyr <- opList[[i]][[4]]$yr
			} else sensyr <- baseRep$yr
		}
		if (useHRP){
			unpackList(renderVals(minYr=minYr, aveYr=aveYr))
			mcbt = baseRep$mc.sbt[,1:Nyears]
			post.bt = as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
			post.dt = t(apply(post.bt,1,function(x){x/mean(x[match(aveYr,yr)])}))

			med.bt = sapply(post.bt,median)
			minYr  = yr[match(min(med.bt),med.bt)]  ## overrides GUI value or user's value

			LRPs = apply(post.dt[,match(minYr,yr),drop=FALSE],1,min); ## across years therefore 1000 mins
			LRPci = quantile(LRPs,probs=c(0.025,0.5,0.975))
			LRP = median(LRPs)
#browser();return()
			USR = 2 * LRP
			USRci = 2 * LRPci
			bS = names(getOptions(wapmod)$is.mcmc)[1]
			mmY =  median(minYr) ## median year of minimum biomass
			legtxt =  as.expression(c(bquote(USR[.(bS)]:2%*%LRP==.(round(USR,2))), bquote(LRP[.(bS)]:~italic(B[.(mmY)])/italic(B)[avg]==.(round(LRP,2)))))
		} else {
			mcdt    <- baseRep$mc.sbdepletion[,1:Nyears]
			post.dt <- as.data.frame(window(mcmc(mcdt),start=Burn+1,thin=Thin))
			legtxt =  as.expression(c(bquote(Ref:~0.4*italic(B)[0])))
		}
		dtci    <- apply(post.dt,2,quantile,probs=c(0.025,0.5,0.975))
		if(usingSweave || val$maxDepletionSensYlim){
			yUpperLimit <- 1.2*max(dtci)
		} else {
			yUpperLimit <- val$depletionSensYlim
		}
		matplot(baseRep$yr,
			t(dtci),
			type="l",
			col=color,
			lty=c(2,1,2), 
			lwd=2,
			ylim=c(0,yUpperLimit),
			xlim=c(min(sensyr),max(sensyr)),
			xlab="Year",
			ylab="Depletion",
			main="Depletion")

		if (useHRP) {
			mtLineColor = c("green4","orange")
			mtLineType =c (2,2); mtLineWidth =c (2,2)
			abline(h=c(USR,LRP), lwd=mtLineWidth, col=mtLineColor, lty=mtLineType)
		} else {
			abline(h=1,col="gainsboro")
			abline(h=0.40, lwd=mtLineWidth, col=mtLineColor, lty=mtLineType)
		}
		
		# Shade the confidence interval
		xx <- c(baseRep$yr,rev(baseRep$yr))
		yy <- c(dtci[1,],rev(dtci[3,]))
		shade <- getShade(color,opacity)
		polygon(xx,yy,density=NA,col=shade)
		# End shade the confidence interval

		for(i in isens){
			iBurn  = opList[[i]][[4]]$mc.burn
			if (useHRP){
				unpackList(renderVals(minYr=minYr, aveYr=aveYr))
				mcbt = opList[[i]][[4]]$mc.sbt[,1:Nyears]
				post.bt = as.data.frame(window(mcmc(mcbt),start=iBurn+1,thin=Thin))
				post.dt = t(apply(post.bt,1,function(x){x/mean(x[match(aveYr,yr)])}))
				dtci    <- apply(post.dt,2,quantile,probs=c(0.025,0.5,0.975))
			} else {
				mcdt <- opList[[i]][[4]]$mc.sbdepletion[,1:Nyears]
				post.dt  <- as.data.frame(window(mcmc(mcdt),start=iBurn+1,thin=Thin))
			}
			dtci <- apply(post.dt,2,quantile,probs=c(0.025,0.5,0.975))
			color <- color + 1				  
			colors <- c(colors,color)

			par(new=T)
			matplot(opList[[i]][[4]]$yr,
				t(dtci),
				type="l",
				col=color,
				lty=c(2,1,2), 
				lwd=2,
				xlim=c(min(sensyr),max(sensyr)),
				ylim=c(0,yUpperLimit),
				xlab="",
				ylab="")

			# Shade the confidence interval
			xx <- c(opList[[i]][[4]]$yr,rev(opList[[i]][[4]]$yr))
			yy <- c(dtci[1,],rev(dtci[3,]))
			shade <- getShade(color,opacity)
			polygon(xx,yy,density=NA,col=shade)
			par(new=F)
			#runNames <- c(runNames,paste0("Sens ",i,": ", basename(opList[[i]][[1]]))) #strsplit(opList[[i]][[1]],"/")[[1]][3]))
			runNames <- c(runNames,paste0(cribtab[i,c("Slab","label")],collapse=": "))
		}
		lty <- c(rep(1,(length(runNames))),mtLineType)
		lwd <- c(rep(2,(length(runNames))),mtLineWidth)
		colors <- c(colors,mtLineColor)
		runNames <- c(runNames,legtxt)
		legend("topright",runNames,lty=lty,col=colors,bty="n", lwd=lwd)
		filename <- paste("fig.sensi.group",sensitivityGroup,".depletion",sep="")
		saveFig(filename)

###########################
### RECRUITS  
	} else if (whichPlot == "recruits") {
		#Get the maximum year from the longest time series for xlim - i.e. might be comparing results from different lengths of data
		#Need to do this here so the xlim is long enough
		sage    = baseRep$sage
		rnyr    = length(baseRep$yr)
		baseryr = baseRep$yr[(1+sage):rnyr]
		ymax    = max(baseRep$mc.rt/1000)*1.2
#print(ymax)
		for(i in c(ibase,isens)){
			ryrtest <- opList[[i]][[4]]$yr
			if(length(ryrtest) >length(baseRep$yr)) {
				rnyr <- length(opList[[i]][[4]]$yr)
				sensryr <- opList[[i]][[4]]$yr[(1+sage):rnyr]
			} else sensryr <- baseryr
			ymax2 = max(sapply(opList[[1]][[4]]$mc.rt/1000,quantile,0.975))
			ymax  = max(ymax,ymax2)
#print(ymax2)
		}

		mc <- baseRep$mc.rt/1000
		mc.rt <- as.data.frame(window(mcmc(mc),start=Burn+1,thin=Thin)) 
		rt <- apply(mc.rt,2,quantile,probs=c(0.025,0.5,0.975))

		if(usingSweave || val$maxRecruitmentSensYlim){
			yUpperLimit <- ymax
		} else {
			yUpperLimit <- val$recruitmentSensYlim
		}

		xp <- plot(baseryr,
			rt[2,],
			type="p",
			pch=20,
			col=color,
			xlim=c(min(sensryr),max(sensryr)),
			ylim=c(0,yUpperLimit),
			xlab="Year", 
			ylab="Recruits (millions)",
			main="Recruits",
			las=1)
		arrows(baseryr, rt[1, ],baseryr,rt[3,],code=3,angle=90,length=0.01,col=color)
		abline(h=median(as.matrix(mc.rt)),col=2,lty=2)
		abline(h=mean(as.matrix(mc.rt)),col=3,lty=2)
		currOffset <- offset

		for(i in isens){
			iBurn  = opList[[i]][[4]]$mc.burn
			sage  <-opList[[i]][[4]]$sage
			rnyr  <- length(opList[[i]][[4]]$yr)
			ryr   <- opList[[i]][[4]]$yr[(1+sage):rnyr]
			mc    <- opList[[i]][[4]]$mc.rt/1000
			mc.rt <- as.data.frame(window(mcmc(mc),start=iBurn+1,thin=Thin)) 
			rt    <- apply(mc.rt,2,quantile,probs=c(0.025,0.5,0.975))
			color <- color + 1
			colors<- c(colors,color)

			par(new=T)
			plot(ryr+currOffset,
				rt[2,],
				pch=20,
				#xlim=c(min(ryr),max(ryr)),
				xlim=c(min(sensryr),max(sensryr)),
				ylim=c(0,yUpperLimit),
				col=color,
				xlab="",
				ylab="",
				las=1,
				axes=F)
			arrows(ryr+currOffset, rt[1, ],ryr+currOffset,rt[3,],code=3,angle=90,length=0.01,col=color)
			par(new=F)
			currOffset <- currOffset + offset
			#runNames <- c(runNames,paste0("Sens ",i,": ", basename(opList[[i]][[1]]))) #strsplit(opList[[i]][[1]],"/")[[1]][3]))
			runNames <- c(runNames,paste0(cribtab[i,c("Slab","label")],collapse=": "))
		}
		lty <- c(rep(1,(length(runNames))),2,2)
		runNames <- c(runNames,"base long-term median","base long-term mean")
		legend("topright",
				runNames,
				lty=lty,
				col=c(colors,2,3),
				bty="n",
				lwd=2) 
		filename <- paste("fig.sensi.group",sensitivityGroup,".recruits",sep="")
		saveFig(filename)

###########################
### REFERENCE POINTS BOXPLOTS  
	} else if (whichPlot == "refpts") {
		if (usingSweave){
			outline = useLog = FALSE
		} else {
			outline  = val$outline  ## for boxplot outliers (or not)
			if (is.null(outline)) outline = FALSE
			useLog   = val$useLog
			if (is.null(useLog)) useLog = FALSE
		}
		if (useLog) {
				scaleFn  = log10; scaleNm = "log[10]~bgroup(\"[\","
		} else {
			scaleFn = function(x){x}; scaleNm = ""
		}
		## Cannot use cbind when MCMC samples are different sizes
		post.rfps = list()

		mcBt    <- subset(baseRep$mcproj,TAC==0)[,paste0("B",baseRep$currYr)]/1000.
		post.Bt <- as.vector(window(mcmc(mcBt),start=Burn+1,thin=Thin))
		mcut    <- 1-exp(-subset(baseRep$mcproj,TAC==0)[,paste0("F",baseRep$lastYr)])
		post.ut <- as.vector(window(mcmc(mcut),start=Burn+1,thin=Thin))

		if (useHRP){
			if(!usingSweave)
				unpackList(renderVals(minYr=minYr, aveYr=aveYr))
			mcbt    = baseRep$mc.sbt[,1:Nyears]/1000.
			post.bt = as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
			post.bo = post.rfps[["bo"]][[baseName]] = apply(post.bt,1,function(x){mean(x[match(aveYr,yr)])})  ## B0 = Bavg
			post.dt = post.rfps[["dt"]][[baseName]] = post.Bt/post.bo

			post.rfps[["ut"]][[baseName]] = post.ut
			mcft    = baseRep$mc.ft[,1:Nyears]
			post.ft = as.data.frame(window(mcmc(mcft),start=Burn+1,thin=Thin))
			post.ht = 1 - exp(-post.ft)  ## harvest rate
			post.ha = post.ht[,match(aveYr,yr)]
			post.uo = post.rfps[["uo"]][[baseName]] = apply(post.ha,1,function(x){mean(x)}) ##  uo = Uavg
			post.ut.uo = post.rfps[["ut.uo"]][[baseName]] = post.ut/post.uo

			#mmY =  median(minYr) ## median year of minimum biomass ## minYr not need for Ref Pt boxplots
#browser();return()

		} else {
			mcbo   <- baseRep$mc$bo/1000
			mcbmsy <- baseRep$mc$bmsy/1000
			mcmsy  <- baseRep$mc$msy #/1000
			mcfmsy <- baseRep$mc$fmsy
			#mcdt   <- baseRep$mc.sbdepletion[,Nyears]
			mcumsy <- 1-exp(-baseRep$mc$fmsy)

			post.bo   = post.rfps[["bo"]][[baseName]] = as.vector(window(mcmc(mcbo),start=Burn+1,thin=Thin))
			post.bmsy = post.rfps[["bmsy"]][[baseName]] = as.vector(window(mcmc(mcbmsy),start=Burn+1,thin=Thin))
			post.msy  = post.rfps[["msy"]][[baseName]] = as.vector(window(mcmc(mcmsy),start=Burn+1,thin=Thin))
			post.fmsy = post.rfps[["fmsy"]][[baseName]] = as.vector(window(mcmc(mcfmsy),start=Burn+1,thin=Thin))
			post.umsy = post.rfps[["umsy"]][[baseName]] = as.vector(window(mcmc(mcumsy),start=Burn+1,thin=Thin))
			
			post.dt   = post.rfps[["dt"]][[baseName]] = post.Bt/post.bo
			post.rfps[["ut"]][[baseName]] = post.ut
			post.Bt.bmsy = post.rfps[["Bt.bmsy"]][[baseName]] = post.Bt/post.bmsy
			post.ut.umsy = post.rfps[["ut.umsy"]][[baseName]] = post.ut/post.umsy
		}

		boxCols   <- senscols[names(ibase)]
#browser();return()

		for(scenario in ibase){
			#runNames <- paste0("Ref: ",strsplit(opList[[scenario]][[1]],"/")[[1]][3])
			AxisName = if (SGs[1]==0) "S00"else paste0("S",pad0(sensNum[baseName],2))
		}
		#for(scenario in isens){
		for (scenarioName in names(sort(sensNum[osens]))){
			i       = isens[scenarioName]
			iBurn   = opList[[i]][[4]]$mc.burn
			
			sensRep = opList[[i]][[4]]
			boxCols = c(boxCols,senscols[names(i)])

			## Current Biomass and Harvest rate
			mcBt     <- subset(sensRep$mcproj,TAC==0)[,paste0("B",sensRep$currYr)]/1000.
			post.Bt2 <- as.vector(window(mcmc(mcBt),start=iBurn+1,thin=Thin))
			mcut     <- 1-exp(-subset(sensRep$mcproj,TAC==0)[,paste0("F",sensRep$lastYr)])
			post.ut2 <- as.vector(window(mcmc(mcut),start=iBurn+1,thin=Thin))

			if (useHRP){
				mcbt     = sensRep$mc.sbt[,1:Nyears]/1000.
				post.bt2 = as.data.frame(window(mcmc(mcbt),start=iBurn+1,thin=Thin))
				post.bo2 = post.rfps[["bo"]][[scenarioName]] = apply(post.bt2,1,function(x){mean(x[match(aveYr,yr)])}) ## B0 = Bavg
				post.dt2 = post.rfps[["dt"]][[scenarioName]] = post.Bt2/post.bo2

				post.rfps[["ut"]][[scenarioName]] = post.ut2
				mcft     = sensRep$mc.ft[,1:Nyears]
				post.ft2 = as.data.frame(window(mcmc(mcft),start=iBurn+1,thin=Thin))
				post.ht2 = 1 - exp(-post.ft2)  ## harvest rate
				post.ha2 = post.ht2[,match(aveYr,yr)]
				post.uo2 = post.rfps[["uo"]][[scenarioName]] = apply(post.ha2,1,function(x){mean(x)}) ##  uo = Uavg
				post.rfps[["ut.uo"]][[scenarioName]] = post.ut.uo2 = post.ut2/post.uo2
#browser();return()

				## Cannot use cbind when MCMC samples are different sizes
				#post.bo   <- cbind(post.bo,post.bo2)
				#post.dt   <- cbind(post.dt,post.dt2)
				#post.ut   <- cbind(post.ut,post.ut2)
				#post.uo   <- cbind(post.uo,post.uo2)
				#post.ut.uo <- cbind(post.ut.uo, post.ut.uo2)

			} else {
				mcbo   <- sensRep$mc$bo/1000
				mcbmsy <- sensRep$mc$bmsy/1000
				mcmsy  <- sensRep$mc$msy #/1000
				mcfmsy <- sensRep$mc$fmsy
				#mcdt   <- opList[[i]][[4]]$mc.sbdepletion[,Nyears]  ## subset the correct final year
				mcumsy <- 1-exp(-sensRep$mc$fmsy)

				post.bo2   = post.rfps[["bo"]][[scenarioName]] = as.vector(window(mcmc(mcbo),start=iBurn+1,thin=Thin))
				post.bmsy2 = post.rfps[["bmsy"]][[scenarioName]] = as.vector(window(mcmc(mcbmsy),start=iBurn+1,thin=Thin))
				post.msy2  = post.rfps[["msy"]][[scenarioName]] = as.vector(window(mcmc(mcmsy),start=iBurn+1,thin=Thin))
				post.fmsy2 = post.rfps[["fmsy"]][[scenarioName]] = as.vector(window(mcmc(mcfmsy),start=iBurn+1,thin=Thin))
				post.umsy2 = post.rfps[["umsy"]][[scenarioName]] = as.vector(window(mcmc(mcumsy),start=iBurn+1,thin=Thin))

				post.dt2   = post.rfps[["dt"]][[scenarioName]] = post.Bt2/post.bo2
				post.rfps[["ut"]][[scenarioName]] = post.ut2
				post.Bt.bmsy2 = post.rfps[["Bt.bmsy"]][[scenarioName]] = post.Bt2/post.bmsy2
				post.ut.umsy2 = post.rfps[["ut.umsy"]][[scenarioName]] = post.ut2/post.umsy2

				## Cannot use cbind when MCMC samples are different sizes
				#post.bo   <- cbind(post.bo,post.bo2)
				#post.bmsy <- cbind(post.bmsy,post.bmsy2)
				#post.msy  <- cbind(post.msy,post.msy2)
				#post.fmsy <- cbind(post.fmsy,post.fmsy2)
				#post.dt   <- cbind(post.dt,post.dt2)
				#post.ut   <- cbind(post.ut,post.ut2)
				#post.umsy <- cbind(post.umsy,post.umsy2)
				#post.Bt.bmsy <- cbind(post.Bt.bmsy,post.Bt.bmsy2)
				#post.ut.umsy <- cbind(post.ut.umsy,post.ut.umsy2)
			}
#browser();return()
			sensName  = basename(opList[[i]][[1]]) #strsplit(opList[[i]][[1]],"/")[[1]][3]
			if (nocrib)
				runNamesi = paste0("Sens ",sensNum[sensName],": ",sensName)
			else
				runNamesi = paste0(cribtab[names(i),c("Slab","label")],collapse=": ")
			#runNamesi <- paste0("Sens ",i,": ",strsplit(opList[[i]][[1]],"/")[[1]][3])
			runNames  <- c(runNames, runNamesi)
			AxisName  <- c(AxisName,paste0(substr(runNamesi,1,1),pad0(sensNum[sensName],2)))
		}
#browser();return()

		par(mfrow=c(2,2), mar=c(2,3,2,0.5), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		#par(mfrow=c(2,2), mai=c(0.3,0.5,0.4,0.2), oma=c(1.,1.2,0.2,0.1))

		if (useHRP) {
			Bt = baseRep$currYr; B0 = "avg"; ut = baseRep$lastYr; MSY = "avg"; yield = "italic(u)[avg]*italic(B)[avg]~(t)" }
		else {
			Bt = baseRep$currYr; B0 = 0; ut = baseRep$lastYr; MSY = "MSY"; yield = "MSY~(t)" }

		## -------------------------------
		## Virgin biomass plot (B0 or Bavg)
		## -------------------------------
		#Set limits for Bo plot
		#yy = scaleFn(post.bo)
		yy = sapply(post.rfps[["bo"]],scaleFn,simplify=FALSE)
		maxBSY = ifelse(usingSweave, TRUE, val$maxBiomassSensYlim)
#browser();return()
		if (is.null(maxBSY)) maxBSY = TRUE
		if(maxBSY){
			#yUpperLimit = max(apply(yy,2,quantile,ifelse(outline,1,0.95)))
			yUpperLimit = max(sapply(yy,quantile,ifelse(outline,1,0.95)))
		} else {
			yUpperLimit <- scaleFn(val$biomassSensYlim)
		}
		#cex.axis = ifelse(dim(yy)[[2]]>5, 0.8, 1.2)
		cex.axis = ifelse(length(yy)>5, 0.8, 1.2)
		main = paste0("as.expression(bquote(",scaleNm,"italic(B)[.(B0)]~(1000~t)",ifelse(useLog,",\"]\")",""),"))")
		quantBox(yy, pch=".", col=boxCols, names=AxisName, main=eval(parse(text=main)), las=1,
			cex.main=1.5, cex.axis=1.2, cex=1.2, ylim=c(0,yUpperLimit), outline=outline, boxwex=0.5, xaxt="n")
		axis(1,at=1:length(AxisName),labels=AxisName,cex.axis=cex.axis)
		if (length(sort(unique(boxCols)))==1)
			legend("bottomleft",runNames, bty="n")
		else
			legend("topright",runNames, lty=1, col=boxCols, bty="n", lwd=2)

		## -------------------------------
		## Depletion plot Bt/B0 or Bt/Bavg
		## -------------------------------
		post.dt = post.rfps[["dt"]]
		maxDSY = ifelse(usingSweave, TRUE, val$maxDepletionSensYlim)
		if (is.null(maxDSY)) maxDSY=TRUE
		if(maxDSY){
			#yUpperLimit = max(apply(post.dt,2,quantile,ifelse(outline,1,0.95)))
			yUpperLimit = max(sapply(post.dt,quantile,ifelse(outline,1,0.95)))
		} else {
			yUpperLimit <- val$depletionSensYlim
		}
		quantBox(post.dt, pch=".", col=boxCols, names=AxisName, main=as.expression(bquote(italic(B[.(Bt)])/italic(B)[.(B0)])), las=1, cex.main=1.5, cex.axis=1.2, cex=1.2, ylim=c(0,yUpperLimit), outline=outline, boxwex=0.5, xaxt="n")
		axis(1,at=1:length(AxisName),labels=AxisName,cex.axis=cex.axis)

		## -------------------------------
		## Harvest rate ut vs. umsy or uavg
		## -------------------------------
		#post.ut.ubase = if (useHRP) post.ut.uo else post.ut.umsy
		post.ut.ubase = if (useHRP) post.rfps[["ut.uo"]] else post.rfps[["ut.umsy"]]
		#yUpperLimit = max(apply(post.ut.ubase,2,quantile,ifelse(outline,1,0.95)))
		yUpperLimit = max(sapply(post.ut.ubase,quantile,ifelse(outline,1,0.95)))
		#yUpperLimit = max(apply(post.fmsy,2,quantile,ifelse(outline,1,0.95)))
		#boxplot(post.fmsy, pch=".", range=0.95, col=boxCols, names=AxisName, main="FMSY (/y)", las=1, cex.axis=1.2, cex=1.2, ylim=c(0,1)) #1.1*max(post.fmsy)
		#quantBox(post.fmsy, pch=".", col=boxCols, names=AxisName, main="FMSY (/y)", las=1, 
		quantBox(post.ut.ubase, pch=".", col=boxCols, names=AxisName, main=as.expression(bquote(italic(u[.(ut)])/italic(u)[.(MSY)])), las=1, 
			cex.main=1.5,cex.axis=1.2, cex=1.2, ylim=c(0,yUpperLimit), outline=outline, boxwex=0.5, xaxt="n")
		axis(1,at=1:length(AxisName), labels=AxisName, cex.axis=cex.axis)

#browser();return()
		#legend("topleft",runNames,lty=1, col=boxCols,bty="n", lwd=2)

		## -------------------------------
		## MSY or Uavg  Yield plots
		## -------------------------------
		#post.yield = if (useHRP) post.uo * post.bo else post.umsy
		post.yield = 
			if (useHRP) sapply(names(post.rfps[["uo"]]),function(i,x,y){x[[i]]*y[[i]]},x=post.rfps[["uo"]],y=post.rfps[["bo"]], simplify=FALSE)
			else post.rfps[["umsy"]]
		#yy = scaleFn(post.yield)
		yy = sapply(post.yield,scaleFn,simplify=FALSE)
		maxRSY = ifelse(usingSweave, TRUE, val$maxRefptSensYlim)
		if (is.null(maxRSY)) maxRSY=TRUE
		if(maxRSY){
			#yUpperLimit = max(apply(yy,2,quantile,ifelse(outline,1,0.95)))
			yUpperLimit = max(sapply(yy,quantile,ifelse(outline,1,0.95)))
		} else {
			yUpperLimit <- scaleFn(val$RefptSensYlim)   #set ylimit high on GUI and assume that MSY will be smaller than Bo or Bmsy
		}
		#yLowerLimit = min(apply(yy,2,quantile,ifelse(outline,0,0.05)))
		yLowerLimit = min(sapply(yy,quantile,ifelse(outline,0,0.05)))
		#boxplot(scaleFn(post.msy), pch=".", range=0.95 , col=boxCols, names=AxisName, main=paste(scaleNm,"MSY (1000 t)"), las=1, cex.axis=1.2, cex=1.2, ylim=c(yLowerLimit,yUpperLimit))
		main = paste0("as.expression(bquote(",scaleNm,yield,ifelse(useLog,",\"]\")",""),"))")
		#quantBox(yy, pch=".", col=boxCols, names=AxisName, main=paste(scaleNm,"MSY (1000 t)"), las=1,
#browser();return()
		quantBox(yy, pch=".", col=boxCols, names=AxisName, main=eval(parse(text=main)), las=1,
			cex.main=1.5, cex.axis=1.2, cex=1.2, ylim=c(yLowerLimit,yUpperLimit), outline=outline, boxwex=0.5, xaxt="n")
		axis(1,at=1:length(AxisName),labels=AxisName,cex.axis=cex.axis)

		save("post.rfps",file=paste0(fdScenarios,baseName[1],"/",fdTables,"post.rfps.rda"))
		assign("post.rfps",post.rfps, envir=.GlobalEnv)

		figDir.curr = figDir
		assign("figDir",paste0(fdScenarios,baseName[1],"/",fdFigures),envir=.GlobalEnv) ## RH temporarily override the global in case user is not in the right scenario
		filename <- paste("fig.sensi.group",sensitivityGroup,".refpts",sep="")
		saveFig(filename)
		assign("figDir",figDir.curr,envir=.GlobalEnv)
	}
}

### Redefine boxplot to show quantiles (RH 150910)
### http://r.789695.n4.nabble.com/Box-plot-with-5th-and-95th-percentiles-instead-of-1-5-IQR-problems-implementing-an-existing-solution-td3456123.html
myboxplot.stats <- function (x, coef=NULL, do.conf=TRUE, do.out=TRUE)
{
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- quantile(x, c(.05,.25,.5,.75,.95), na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  out <- x < stats[1] | x > stats[5]
  conf <- if (do.conf)
    stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)
  list(stats = stats, n = n, conf = conf, out = x[out & nna])
} 

boxcode = deparse(boxplot.default)
boxcode = gsub("boxplot\\.stats","myboxplot.stats",boxcode)
eval(parse(text=c("qboxplot=",boxcode)))

quantBox = function (x, use.cols = TRUE, ...) ## taken from boxplot.matrix
{
	if (rev(class(x))[1]=="matrix") {
		groups <- if (use.cols) 
			split(x, rep.int(1L:ncol(x), rep.int(nrow(x), ncol(x))))
		else split(x, seq(nrow(x)))
		if (length(nam <- dimnames(x)[[1 + use.cols]])) 
		names(groups) <- nam
		qboxplot(groups, ...)
	}
	else qboxplot(x, ...)
}

updateSG = function(){
	if (!exists("cribtab",envir=.GlobalEnv)) return()
	scenarioHeader = getWinVal()$scenarioHeader
	SGlab  = basename(scenarioHeader[,1])
	scenarioHeader[,2] = cribtab[SGlab,"SG"]
	updateGUI()
}
