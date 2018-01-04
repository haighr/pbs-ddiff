fig.a <- function(){
	# Total landings
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	barplot(A$obs_ct, names.arg=A$yr, xlab="Year", ylab="Landings", las=1, space=0, col="moccasin")  #must plot observed catch
	addScenLab()
	saveFig("fig.catch")
}

### FISHING EFFORT
fig.effort <- function() {
	# Total effort
	if (is.null(A$effort) || is.na(A$effort) || all(A$effort==0)) {
		cat("WARNING (fig.effort): Effort data all zeroes or not available\n"); return(invisible("No effort data")) }
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	plot(A$yr[1:length(A$effort)], A$effort/1000,  xlab="Year", ylab="Fishing effort (1000 hr)", type="l", lty=1, lwd=2, col="darkblue",main="Effort",las=1)  
	addScenLab()
	saveFig("fig.effort")
}

### SPAWNING BIOMASS
fig.b <- function(includeMPD=FALSE, ylimit=6, useMaxYlim=T, useHRP=F, minYr, aveYr, opacity="20", includeCatch=T, tac.use=0, ...)
{
	simplify = usingSweave
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())
	if (!isScenMcmc()) {
		if (simplify){
			frame(); addLabel(.5,.5,"No MCMC",col="red",cex=2); return()
		} else {
			cat("WARNING (fig.b): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data"))
		}
	}
	if(A$delaydiff==1){
		ddmcmc  <- A$mcproj 
		tacs = sort(unique(ddmcmc$TAC))
		ntac = which(abs(tacs-tac.use)==min(abs(tacs-tac.use))) ## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
		tac.use = tacs[ntac]
		dmcmc   <- subset(ddmcmc, TAC==tac.use)   #take only the first tac from each posterior sample - the control points do not change with tac
		post.dmcmc <- as.data.frame(window(mcmc(dmcmc),start=(Burn+1),thin=Thin))
		if ("Bmin" %in% names(post.dmcmc)) {
			Bmin    <-  post.dmcmc$Bmin
			Bmean   <-  post.dmcmc$BAvg_S
			bmin    <- median(Bmin)
			bmean   <- median(Bmean)
			Bcurr   <- post.dmcmc[,paste0("B",A$currYr),drop=FALSE]/1000.
		} else { 
			Bmin=Bmean=bmin=bmean=NULL
			Bcurr = post.dmcmc[,paste0("B",c(A$currYr,A$projYr[1:Nproj]))]/1000. ## current year is now in post.bt
		}
		Bcurr.ci <-  sapply(Bcurr,quantile,quants3,na.rm=T)  ## current and projected biomass
	}
	
	# Spawning stock virgin biomass
	mcbo <- (A$mc$bo)/1000.
	post.bo <- as.data.frame(window(mcmc(mcbo),start=Burn+1,thin=Thin))
	boci <- apply(post.bo,2,quantile,probs=quants3)

	if (useHRP){
		if (missing(minYr) || missing(aveYr))
			unpackList(renderVals(aveYr=getWinVal()$aveYr,minYr=getWinVal()$minYr))
		else
			unpackList(renderVals(minYr=minYr, aveYr=aveYr))

		HRPs = calcHRP(A=A, aveYr=aveYr, Burn=Burn)  ## one function to collect all of the stuff below
		unpackList(HRPs)
#browser();return()
	}
	else {
		avebt=NULL
		unpackList(renderVals(aveYr=getWinVal()$aveYr,minYr=getWinVal()$minYr))
	}

	if(A$delaydiff==1 && !useHRP) mcbt <- (A$mc.sbt[,1:nyear])/1000
	if(A$delaydiff==0) mcbt <- (A$mc.sbt[,1:nyrs])/1000 #projection code not written yet for ASM so just use projection year from main model
	if (useHRP){
		post.bt = post.bt/1000.
		avebt   = avebt/1000.
	}
	else
		post.bt <- as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
	btci <- apply(post.bt,2,quantile,probs=quants3)
	if(useMaxYlim){
		yUpperLimit <- max(btci,(A$sbt)/1000)
	} else {
		yUpperLimit <- ylimit
	}
	if(A$delaydiff==1) 
		btci <- cbind(btci,Bcurr.ci[,-1]) # add the projection year from the projection model (exclude current year)
	projYrs = c(A$currYr,A$projYr[1:Nproj])
	allYrs  = A$yr[1]:rev(A$projYr[1:Nproj])[1]; nyrs = length(allYrs)
	ryear   = nyear + 1 ## number of reconstruction (not projection) years
	xlim = range(c(A$yr,projYrs))

	if (!simplify)  par(mar=c(3,4,2,1), mgp=c(2.5,0.5,0))
	matplot(A$yrs,t(btci[,1:ryear]),type="l",col=1,lty=c(2,1,2), lwd=2,ylim=c(0,yUpperLimit), xlim=xlim,
		ylab=ifelse(simplify,"","Biomass (thousand tonnes)"), xlab=ifelse(simplify,"",""), main=ifelse(simplify,"","Biomass"), las=1, cex.axis=1.2, cex.lab=1.5)
#browser();return()
	if (!simplify) mtext("Year",side=1,line=2,cex=1.5)
	abline(h=avebt, col=1, lty=3, lwd=2)
	xx <- c(A$yrs,rev(A$yrs))
	yy <- c(btci[1,1:ryear],rev(btci[3,1:ryear]))
	shade <- getShade(1,opacity)
	polygon(xx,yy,density=NA,col=shade)
	arrows(A$yrs[1]-0.8,boci[1],A$yrs[1]-0.8,boci[3],col="green4", code=0, lwd=2)
	points(A$yrs[1]-0.8,boci[2],cex=1.5,col=1,pch=21,bg="green")
#browser();return()

	#projection year
	lines(projYrs,  btci[2,ryear:nyrs], col=2, lwd=2)
	lines(projYrs,  btci[1,ryear:nyrs], col=2, lwd=2, lty=2)
	lines(projYrs,  btci[3,ryear:nyrs], col=2, lwd=2, lty=2)
	xxp <- c(projYrs,rev(projYrs))
	yyp <- c(btci[1,ryear:nyrs],rev(btci[3,ryear:nyrs]))
	shade <- getShade(2,opacity)
	polygon(xxp,yyp,density=NA,col=shade)
	#if(A$delaydiff==1) legend("topright",c("Median Bmin","Median BAvg 1956-2004"),lty=c(2,2),pch=c(-1,-1),lwd=c(1,1),col=c(2,3),bty="n")

	legtxt =c(
		paste0("Bt (", A$yrs[1],"-",rev(A$yrs)[1], ")"),
		paste0("Bt proj. (", projYrs[1],"-",rev(projYrs)[1], ")"),
		paste0("Bavg (", aveYr[1],"-",rev(aveYr)[1], ")"),
		paste0(ifelse(useHRP,"Binit","B0")," (", A$yr[1], ")"),
		paste0("Bcurr (", A$currYr, ")")
		)
	lty = c(1,1,3,NA,NA); pch=c(-1,-1,-1,21,21); lwd=c(2,2,2,1,1); col=c(1,2,1,1,2); pt.bg=c(NA,NA,NA,"green","yellow")

	if(includeMPD){
		Bt<-A$sbt 
		lines(A$yrs,Bt/1000,type="l",col=mpdLineColor,lty=1, lwd=2,ylim=c(0,yUpperLimit), las=1, xlab="",ylab="")
		legtxt = c(legtxt, paste0("Bt MPD (", A$yrs[1],"-",rev(A$yrs)[1], ")"))
		lty = c(lty,1); pch=c(pch,-1); lwd=c(lwd,2); col=c(col,mpdLineColor); pt.bg=c(pt.bg,NA)
	}
	points(A$currYr,Bcurr.ci[2,paste0("B",A$currYr)],pch=21,cex=1.5,col="red",bg="yellow")

	if(includeCatch){
		removals = apply(A$obs_ct,2,sum)/1000.
		drawBars(A$yr,removals,width=1,col="red3",lwd=2)
		tacos = rep(tac.use/1000.,length(projYrs))  ## add on TAC removals
#browser();return()
		drawBars(projYrs,tacos,width=1,col="pink2",lwd=2)
	}
	legend("topright",inset=0.01, legend=legtxt, title="MCMC (posterior distribution)",
		lty=lty, pch=pch, lwd=lwd, col=col, pt.bg=pt.bg, bty="n")

#browser();return()
	if (!simplify) {
		addScenLab(x=0.25)
		saveFig("fig.spawning.biomass", fun.call=match.call())
	}
}

fig.biomass.mpd <- function(ylimit=45, useMaxYlim=TRUE,...){
	# Spawning stock biomass
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(A$delaydiff==1){
		ddmpd <- A$mpdproj
		## Bmin and Bavg are hard-wired in tpl code for P.cod.
		bmin <- ddmpd$Bmin[1]
		bmean <-  ddmpd$BAvg_S[1] }

	if(useMaxYlim){
		yUpperLimit <- 1.05*max(A$sbt)/1000 #max((A$sbt)/1000)
	} else {
		yUpperLimit <- ylimit
	}
	Bt    <- A$sbt[1:nyear]
	bmin  <- min(Bt)
	aveYr = renderVals(aveYr=getWinVal()$aveYr,simplify=T)
	bmean <- mean(Bt[is.element(A$yr,aveYr)],na.rm=TRUE)  #mean(Bt)
	yrmin <- max(A$yr[match(bmin,Bt)])  ## to avoid ties with early-period lows
	junk = try(setWinVal(list(minYr=yrmin)),silent=TRUE)

	plot(A$yr,Bt/1000,type="l",col=mpdLineColor,lty=1, lwd=2,ylim=c(0,yUpperLimit),ylab="Biomass (thousand t)", xlab="Year", main="Biomass", las=1) #for pcod testing - sbo is male + female
	abline(v=yrmin, col="red", lty=3)
	#points(A$yr[1]-0.8,(A$sbo)/1000,col=mpdLineColor,pch=1)
	points(A$yr[1]-0.8,(A$sbo)/1000,cex=1.5,col=1,pch=21,bg="green")
	if(A$delaydiff==1){
		abline(h=bmin/1000,col="red",lty=2)
		abline(h=bmean/1000,col="green4",lty=2)
		legend("topright",inset=0.025, legend=c(
			paste0("Bt (", A$yr[1],"-",rev(A$yr)[1], ")"),
			paste0("Bavg (", aveYr[1],"-",rev(aveYr)[1], ")"),
			paste0("Bmin (",paste0(yrmin,collapse=","),")"),
			paste0("B0 (", A$yr[1], ")")
			),
			title="MPD (mode of posterior distribution)",
			lty=c(1,2,2,NA), pch=c(-1,-1,-1,21), lwd=c(2,1,1,1), col=c(mpdLineColor,"green4","red","black"), pt.bg="green", bty="n")
	}
	#mtext("Effective SB0_2005_model = SB0 * brat = 21,078 x 0.675 = 14228 t",side=3, line=0, outer=F)
	addScenLab()
	saveFig("fig.spawning.biomass.mpd")
}

#total biomass (MCMC)
fig.bt <- function(includeMPD=F,ylimit=150,useMaxYlim=T,opacity="20",includeCatch=T,...){
	# Total stock biomass
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.bt): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	mcbt <- (A$mc.tbt)/1000
	post.bt <- as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
	btci <- apply(post.bt,2,quantile,probs=quants3)
	if(useMaxYlim){
		yUpperLimit <- max(btci,(A$tbt)/1000)
	}else{
		yUpperLimit <- ylimit
	}
	matplot(A$yr,t(btci[,1:nyear]),type="l",col=1,lty=c(2,1,2), lwd=2,ylim=c(0,yUpperLimit),ylab="Total biomass (thousand t)", las=1)
	xx <- c(A$yr,rev(A$yr))
	yy <- c(btci[1,1:nyear],rev(btci[3,1:nyear]))
	shade <- getShade(1,opacity)
	polygon(xx,yy,density=NA,col=shade)

	if(includeMPD){ ## A$tbt only available for nyear (not nyrs)
		lines(A$yr,(A$tbt[1:nyear])/1000,type="l",col=mpdLineColor,lty=1, lwd=2,ylim=c(0,yUpperLimit), las=1, xlab="",ylab="")
		legend("topright",c("MPD estimate"),lty=1,lwd=2,col=mpdLineColor,bty="n")
	}
	if(includeCatch){
		removals = apply(A$obs_ct,2,sum)/1000.
		drawBars(A$yr,removals,width=1,col="orange",lwd=2)
	}
	addScenLab()
	saveFig("fig.total.biomass")
}

fig.bt.mpd <- function(ylimit=150,useMaxYlim=T,...){
	# Total stock biomass
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if(useMaxYlim){
		yUpperLimit <- max(A$tbt/1000)
	} else {
		yUpperLimit <- ylimit
	}

	matplot(A$yr,(A$tbt[,1:nyear])/1000,type="l",col=mpdLineColor,lty=1, lwd=2,ylim=c(0,yUpperLimit),ylab="Total biomass (thousand t)", las=1) 
	abline(h=median(A$tbt[1:nyear])/1000,col=2,lty=2)
	legend("topright","MPD long-term median",lty=2,pch=-1,lwd=1,col=2,bty="n")
	addScenLab()
	saveFig("fig.total.biomass.mpd")
}


fig.biomass.recruits <- function(yBiomassYlim=45,
                                 useMaxBiomassYlim=T,
                                 yRecruitmentYlim=250,
                                 useMaxRecruitmentYlim=T){
  # make two-panel plot of biomass and recruitment
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	par(oma=c(1,1,1,1),mar=c(3,3,1,1))
	par(mfrow=c(2,1))
	fig.b(ylimit=yBiomassYlim,useMaxYlim=useMaxBiomassYlim, xlab="")
	fig.d(ylimit=yRecruitmentYlim,useMaxYlim=useMaxRecruitmentYlim,xlab="Year")
	saveFig("fig.biomass.recruitment")
}

fig.depletion.mpd <- function(ylimit=3.5, useMaxYlim=TRUE, useHRP=FALSE, minYr, aveYr)
{
	# Spawning stock biomass depletion
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (useHRP){
		unpackList(renderVals(minYr=minYr, aveYr=aveYr))
		dBt = A$sbt[1:nyear]/mean(A$sbt[1:nyear])
		LRP = dBt[match(minYr,yr)]
		USR = 2 * LRP
		#legtxt = c(expression(B[avg]),expression(2~B[avg]))
		mmY =  median(minYr) ## median year of minimum biomass
		legtxt =  as.expression(c(bquote(USR:2%*%LRP==.(round(USR,2))), bquote(LRP:~italic(B[.(mmY)])/italic(B)[avg]==.(round(LRP,2)))))
	} else {
		dBt = A$sbt[1:nyear]/A$sbo
		LRP = 0.2
		USR = 0.4
		legtxt = as.expression(lapply(c(USR,LRP), function(x) bquote(.(x)*italic(B)[0])))
	}
	if(useMaxYlim){
		yUpperLimit <- 1.05*max(dBt)
	} else {
		yUpperLimit <- ylimit
	}
	matplot(A$yr, dBt, type="l", col="blue", lty=1, lwd=2, ylim=c(0,yUpperLimit), xlab="Year", ylab="Spawning Depletion", main="Spawning Depletion", las=1, mgp=c(2,0.5,0))
	abline(h=c(LRP,USR), lwd=mtLineWidth, col=c(lrpLineColor,mtLineColor), lty=mtLineType)
	legend("topright", legend=legtxt, lty=c(mtLineType), lwd=c(mtLineWidth), col=c(mtLineColor,lrpLineColor), bty="n")
	addScenLab(x=0.4)
	saveFig("fig.depletion.mpd")
}

## Spawning stock depletion
fig.c <- function(includeMPD=F, ylimit=3.5, useMaxYlim=T, useHRP=F, minYr, aveYr, opacity="20", ...)
{
	simplify = usingSweave
	adieu = if (simplify) invisible else function(){expandGraph(mfrow=c(1,1), new=FALSE)}
	on.exit(adieu())
	if (!isScenMcmc()) {
		if (simplify){
			frame(); addLabel(.5,.5,"No MCMC",col="red",cex=2); return()
		} else {
			cat("WARNING (fig.c): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data"))
		}
	}
	legtxt =  as.expression(bquote(italic(B[t])/italic(B)[0]))
	if (useHRP){
		if (missing(minYr) || missing(aveYr))
			unpackList(renderVals(aveYr=getWinVal()$aveYr,minYr=getWinVal()$minYr))
		else
			unpackList(renderVals(minYr=minYr, aveYr=aveYr))

		## MPD estimate of depletion
		dBt = A$sbt[1:nyear]/mean(A$sbt[match(aveYr,yr)]) ## MPD

		HRPs = calcHRP(A=A, aveYr=aveYr, Burn=Burn)  ## one function to collect all of the stuff below
		unpackList(HRPs)
#browser();return()
		legtxt =  as.expression(bquote(italic(B[t])/italic(B)[avg]))
		legtxt  =  as.expression(c(legtxt, bquote(USR:2%*%LRP==.(show0(round(dUSR,2),2))), bquote(LRP:~italic(B)[min]/italic(B)[avg]==.(show0(round(dLRP,2),2)))))
	} else {
		dBt = A$sbt[1:nyear]/A$sbo ## MPD
		mcdt <- A$mc.sbdepletion[,1:nyear]
		post.dt  <- as.data.frame(window(mcmc(mcdt),start=Burn+1,thin=Thin))
		LRP = 0.4
		USR = 0.8
		legtxt = as.expression(c(legtxt,lapply(c(USR,LRP), function(x) bquote(.(x)*italic(B)[MSY]))))
	}
	legcol = c("black","green3","red")
	leglty = c(1,3,3)

	ddmcmc  <- A$mcproj 
	dmcmc   <- subset(ddmcmc, TAC==0)   #take only the first tac from each posterior sample - the control points do not change with tac
	post.dmcmc <- as.data.frame(window(mcmc(dmcmc),start=(Burn+1),thin=Thin))
	Bcurr = post.dmcmc[,paste0("B",c(A$currYr,A$projYr[1:Nproj]))]
	if (useHRP)
		Dcurr = apply(Bcurr,2,function(x,y){x/y}, y=post.abt)  ## divide Bcurr by Bavg
	else 
		Dcurr = apply(Bcurr,2,function(x,y){x/y}, y=post.dmcmc[,"B0"])  ## divide Bcurr by B0
	Dcurr.ci <- apply(Dcurr,2,quantile,quants3,na.rm=T)  ## current and projected depletion

#browser();return()

	dtci <- apply(post.dt,2,quantile,probs=quants3)
	dtci =  cbind(dtci,Dcurr.ci[,1]) ## add current year
	if(useMaxYlim){
		yUpperLimit <- max(dtci,dBt)
	}else{
		yUpperLimit <- ylimit
	}
	if (!simplify) par(mar=c(3.5,4,2,0.5),mgp=c(2.2,0.5,0))
	txt.main = if (useHRP) "Biomass relative to average biomass" else "Spawning biomass depletion" 
	txt.ylab = if (useHRP) as.expression(bquote(italic(B[t])/italic(B)[avg])) else as.expression(bquote(italic(B[t])/italic(B)[0]))
	matplot(A$yrs, t(dtci[,1:nyrs]), type="n", col=1, lty=c(2,1,2), lwd=2, ylim=c(0,yUpperLimit), main=ifelse(simplify,"",txt.main), xlab=ifelse(simplify,"","Year"), las=1, ylab=ifelse(simplify,"",txt.ylab), cex.lab=1.5, cex.axis=1.2)
	#abline(h=1); abline(v=2017)
	#mcproj <- mcmc2(subset(A$mcproj, TAC==0), start=Burn+1, thin=Thin)
	if (useHRP) {
		for (i in 1:2) {
			ihrp = switch(i, dLRPci, dUSRci)
			polygon(par()$usr[c(1,1,2,2)],ihrp[c(1,3,3,1)],col=ifelse(i==1,lucent("pink",0.3), lucent("greenyellow",0.5)),border="gainsboro")
			lines(par()$usr[c(1,2)],ihrp[c(2,2)],col=ifelse(i==1,"red","green3"),lty=3,lwd=2)
		}
	} else {
		for (i in c(USR,LRP)){
			imsy = quantile(post.dmcmc$BMSY*i/post.dmcmc$B0,c(.05,0.5,0.95))
			polygon(par()$usr[c(1,1,2,2)],imsy[c(1,3,3,1)],col=ifelse(i==0.4,lucent("pink",0.3), lucent("greenyellow",0.5)),border="gainsboro")
			lines(par()$usr[c(1,2)],imsy[c(2,2)],col=ifelse(i==0.4,"red","green3"),lty=3,lwd=2)
		}
	}
	xx <- c(A$yrs,rev(A$yrs))
	yy <- c(dtci[1,1:nyrs],rev(dtci[3,1:nyrs]))
	shade <- getShade(1,opacity)
#browser();return()
	polygon(xx,yy,density=NA,col=shade)
	matlines(A$yrs,t(dtci[,1:nyrs]),type="l",col=1,lty=c(2,1,2), lwd=2)
	#abline(h=c(0.2,0.4), lwd=mtLineWidth, col=c(lrpLineColor,mtLineColor), lty=mtLineType)

	if(includeMPD){
		lines(A$yr,dBt,type="l",col=mpdLineColor,lty=1, lwd=2,ylim=c(0,yUpperLimit), las=1, xlab="",ylab="")
		legtxt = as.expression(c(legtxt, ifelse(simplify,"MPD","MPD estimate")))
		legcol = c(legcol,mpdLineColor)
		leglty = c(leglty,1)
	}
	legend("topright", title="MCMC (posterior distribution)", legend=legtxt, lty=leglty, lwd=rep(2,length(legcol)), col=legcol, seg.len=3, bty="n")
#browser();return()
	box()
	if (!simplify) {
		addScenLab(x=0.4)
		saveFig("fig.depletion", fun.call=match.call())
	}
}

#RECRUITS
#this was figure e in 2010 assessment
fig.d <- function(ylimit=250,useMaxYlim=T, ...){   
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	sage<-A$sage
	ryr <- yr[(1+sage):nyear]
	#ryr<-yr
	
	plot(ryr, A$rt/1000, lty=1, col=1,type="o", pch=19,ylim=c(0,1.2*max(A$rt/1000)), xlab="Year",ylab="Recruits (million)", las=1, main="Recruits")
	abline(h=median(A$rt)/1000,col=2,lty=2)
	abline(h=mean(A$rt)/1000,col=3,lty=2)
	legend("topright",legend=c("MPD long-term median","MPD long-term mean"),lty=c(2,2),pch=c(-1,-1),lwd=c(1,1),col=c(2,3),bty="n")
	addScenLab()
	saveFig("fig.recruits.mpd")
}

#RECRUITS
#this was figure e in 2010 assessment
fig.dmcmc <- function(ylimit=250, useMaxYlim=TRUE, ...){
	#recruits
	simplify = usingSweave
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.dmcmc): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }

	## Divide by 1000 to get millions (assume input catch in tonnes and body weight in kg)
	dfac = 1000.
	lcol = c("#D55E00", "#009E73")  ## colour-blind vermillion, bluegreen
	sage <- A$sage
	ryr  <- yr[(1+sage):nyear]
	#ryr <- yr
	
	mc    <- A$mc.rt
	mc.rt <- as.data.frame(window(mcmc(mc), start=Burn+1, thin=Thin)) 
	rt    <- apply(mc.rt, 2, quantile, probs=quants3) #gets quantiles for number of age 1 recruits

	if(useMaxYlim){
		yUpperLimit <- max(rt)/dfac
	}else{
		yUpperLimit <- ylimit
	}
	if (!simplify)  par(mar=c(3,4,2,1), mgp=c(2.5,0.5,0))

	plot(ryr, rt[2,]/dfac, type="n", ylim=c(0,yUpperLimit), xlab="", ylab="Recruits (million)", main="Recruits", las=1, cex.lab=1.5, cex.axis=1.2, cex.main=1.75)
	mtext("Year", side=1, line=2, cex=1.5)
	axis(1, at=seq(1975,2015,10), tcl=-0.25, labels=FALSE)

	abline(h=median(as.matrix(mc.rt))/dfac, col=lcol[1], lty=2, lwd=2)  ## colour-blind vermillion
	abline(h=mean(as.matrix(mc.rt))/dfac,   col=lcol[2], lty=2, lwd=2)  ## colour-blind bluegreen
	arrows(ryr, rt[1, ]/dfac, ryr, rt[3,]/dfac, code=3, angle=90, length=0.01, lwd=1.5)
	points(ryr, rt[2,]/dfac, pch=20, cex=1.5, col="grey20")

	##points(xp,A$nt[,1],pch=19,cex=1)
	legend("topright", c("MCMC long-term median","MCMC long-term mean"), lty=c(2,2), pch=c(-1,-1), lwd=c(2,2), col=lcol, bty="n", seg.len=4)
	addScenLab()
	saveFig("fig.recruits.mcmc", fun.call=match.call())
}

#SPR status - fmsy
fig.e1 <- function(){
	#The relative spawning potential ratio (1-spr)/(1-spr.at.msy)  
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.e1): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	spr <- A$mc.sprstatus_fmsy #read.table("ccam.spr",h=F)
	if (is.null(spr)) {
		cat("WARNING (fig.e1): No spawner-per-recruit data available\n"); return(invisible("No SPR data")) }
	post.spr <- as.data.frame(window(mcmc(spr),start=Burn+1,thin=Thin))
	sprci <- apply(post.spr,2,quantile,probs=quants3)

	matplot(A$yr,t(sprci),type="l",col=c(2,1,2),lty=c(3,1,3), lwd=2, pch=c(-1, 0, 1),ylim=c(0,2)
		,xlab="Year",ylab="(1-SPR)/(1-SPR at fmsy)")
	abline(h=1, lwd=mtLineWidth, col=mtLineColor, lty=mtLineType)
	text(1980, 1, "Management target", pos=3)
	addScenLab()
	saveFig("fig.spr.fmsy")
}

#SPR status - f40
fig.e2 <- function(){
	#The relative spawning potential ratio (1-spr)/(1-spr.at.msy)  
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.e2): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	spr <- A$mc.sprstatus_f40 #read.table("ccam.spr",h=F)
	if (is.null(spr)) {
		cat("WARNING (fig.e2): No spawner-per-recruit data available\n"); return(invisible("No SPR data")) }
	post.spr <- as.data.frame(window(mcmc(spr),start=Burn+1,thin=Thin))
	sprci <- apply(post.spr,2,quantile,probs=quants3)
	
	matplot(A$yr,t(sprci),type="l",col=c(2,1,2),lty=c(3,1,3), lwd=2, pch=c(-1, 0, 1),ylim=c(0,2)
		,xlab="Year",ylab="(1-SPR)/(1-SPR at f40)")
	abline(h=1, lwd=mtLineWidth, col=mtLineColor, lty=mtLineType)
	text(1980, 1, "Management target", pos=3)
	addScenLab()
	saveFig("fig.spr.f40")
}

#FISHING MORTALITY   MCMC
fig.Fmcmc <- function(opacity="20", tac.use=600, ...){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.Fmcmc): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	ddmcmc <- A$mcproj
	tacs = sort(unique(ddmcmc$TAC))
	ntac = which(abs(tacs-tac.use)==min(abs(tacs-tac.use))) ## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	tac.use = tacs[ntac]
	dmcmc <- subset(ddmcmc, TAC==tac.use)   #take only the first tac from each posterior sample - the control points do not change with tac
	post.dmcmc <- as.data.frame(window(mcmc(dmcmc),start=Burn+1,thin=Thin))
	if ("Bmin" %in% names(post.dmcmc)) {
		Fmean <-  post.dmcmc$FAvg_S
		Fmean <-  median(Fmean)
	} else {
		Fmean=NULL; use.hist=F
	}
	#Fishing mortality - GEAR 1 ONLY
	ft <- A$mc.ft #read.table("ccam.bt2",h=F)
	post.ft <- as.data.frame(window(mcmc(ft),start=Burn+1,thin=Thin))
	ftci <- apply(post.ft,2,quantile,probs=quants3)

	expandGraph(mfrow=c(1,1))
	matplot(A$yr,t(ftci[, 1:length(A$yr)]),type="l",col=1,lty=c(2,1,2), lwd=2, pch=c(-1, 0, 1)
		,xlab="Year",ylab="Fishing mortality rate (/yr)",main="Fishing mortality", ylim=c(0,1.1*max(ftci)))#ylim=c(0,   1.1*max(ftci)
	xx <- c(A$yr,rev(A$yr))
	yy <- c(ftci[1,],rev(ftci[3,]))
	shade <- getShade(1,opacity)
	polygon(xx,yy,density=NA,col=shade)
	
	#abline(h=median(as.matrix(ft)),col=3,lty=2)
	if (use.hist) {
		abline(h=Fmean,col=2,lty=2)
		legend("topleft",c("MCMC Fishing mortality", "Median FAvg 1956-2004"),lty=c(1,2),lwd=c(2,1),col=c(1, 2),bty="n")
	} else
		legend("topleft",c("MCMC Fishing mortality"),lty=c(1),lwd=c(2),col=c(1),bty="n",inset=0.05)
	addScenLab(0.95, 0.95, adj=c(1,0))
	saveFig("fig.fishing.mortality.mcmc")
}

#FISHING MORTALITY   MPD
fig.Fmpd <- function(opacity="20",...){
	#Fishing mortality - GEAR 1 ONLY
	on.exit(expandGraph(mfrow=c(1,1), new=FALSE, mgp=c(3,0.5,0)))
	plot(A$yr, A$ft[1,1:length(A$yr)],type="l",col=mpdLineColor,lty=1, lwd=2, las=1,
		xlab="Year",ylab="Fishing mortality rate (/yr)",main="Fishing mortality", ylim=c(0,1.1*max(A$ft)))#ylim=c(0,
	lines(A$yr, A$ut[1:length(A$yr)], col="red",lty=1, lwd=2)
	abline(h=median(A$ft[1,]),col="green4",lty=2)
	abline(h=mean(A$ft[1,]),col=mpdLineColor,lty=2)
	abline(h=max(A$ft),col=mpdLineColor,lty=3)
	abline(h=max(A$ut),col="red",lty=3)
	legend("topleft",c(
		"Fishing mortality",
		"Median fishing mortality",
		"Mean fishing mortality",
		"Max fishing mortality",
		"Exploitation (harvest) rate",
		"Max exploitation rate"
		), title="MPD (mode of posterior distribution)", inset=0.02, 
		lty=c(1,2,2,3,1,3), lwd=c(2,1,1,1,2,1), col=c(mpdLineColor,"green4",rep(mpdLineColor,2),rep("red",2)), bty="n")
	addScenLab(0.95, 0.95, adj=c(1,0))
	saveFig("fig.fishing.mortality.mpd")
}

### PA Phase plots
fig.h <- function(){
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (!isScenMcmc()) {
		cat("WARNING (fig.h): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
	expandGraph(mfrow=c(2,1), mar=c(3,4,0.5,1))
	n <- length(A$yr)

	S   = Burn+1; N = dim(A$mc.ft)[1]
	ft  = mcmc2(A$mc.ft,S,N,Thin)   #read.table("ccam.ft",h=F)
	sbt = mcmc2(A$mc.sbt,S,N,Thin)  #read.table("ccam.sbt",h=F)
	mc  = mcmc2(A$mc,S,N,Thin)      #read.table("ccam.mcmc", h=T)
	sbstatus <- sbt/mc$bmsy
	fstatus <- ft/mc$fmsy
	
	sbci <- apply(sbstatus, 2, quantile, probs=c(0.5))
	ftci <- apply(fstatus, 2, quantile, probs=c(0.5))
	if (any(ftci>2)) {
		yfun = function(x){xx=log2(GT0(x)); names(xx)=x; return(xx)}
		ynuf = function(x){2^x}
		ylog = TRUE
	} else {
		yfun = ynuf = function(x){x; names(x)=x; return(x)}
		ylog = FALSE
	}
	ftci = yfun(ftci)
#browser();return()
	maxX = 1.1*max(sbci[1:n])
	maxY = 1.1*max(ftci,1)

	plot(sbci[1:n], ftci[1:n], type="n",xlim=c(0,maxX), ylim=c(min(min(ftci),0),maxY), xlab="Median Bt / Bmsy", ylab="", yaxt="n", mgp=c(1.8,0.5,0), cex.lab=1.2)
	#ytcklab = if (ylog) 2^(0:6) else pretty(ftci)
	#ytcklab = pretty(as.numeric(names(ftci)))
	#if (ylog) { ytcklab = setdiff(ytcklab,0); ytcklab=sort(unique(c(0.1,0.2,0.5, 1,2, 5, ytcklab))) }
	ytcklab = sort(unique(rep(c(1,2,4),each=4)*10^c(-2:1)))
	ytck    = yfun(ytcklab)
	axis(2, at=ytck, labels=ytcklab, las=1, cex=0.8)
	mtext("Median Ft / Fmsy",side=2, line=2.5, cex=1.2)
	
#browser();return()
#	plot(sbci[1:n], log(ftci[1:n]), type="n",xlim=c(0,maxX), ylim=c(0.01,log(maxY)), xlab="Median Bt / Bmsy", ylab="Median Ft/Fmsy", log="y")
	abline(h=yfun(1),v=1,lty=2,col=2)
	lines(sbci[1:n], ftci[1:n],type="o")
	gletter(1)
	addScenLab(0.95, 0.95, adj=c(1,0.5))
	
	#2nd plot
	sbstatus <- sbt/mc$bo
	sbci <- apply(sbstatus, 2, quantile, probs=c(0.5))
	ftci <- apply(ft, 2, quantile, probs=c(0.5))
	ftci = yfun(ftci)
	ssb <- seq(0,max(sbci)*1.1,length=100)
	fp <- median(mc$fmsy)*(ssb-0.1)/0.3
	fp[ssb<=0.1] <- 0
	fp[ssb>0.4] <- median(mc$fmsy)
	fp=yfun(fp)
	maxX = 1.1*max(sbci[1:n])
	maxY = 1.1*max(ftci,fp)
	#maxY <- max(fp,median(A$ft))*1.1
#browser();return()

	plot(ssb, fp, type="l", xlim=c(0,maxX), ylim=c(min(min(ftci),0),maxY), lwd=2, xlab="Median Bt / B0",ylab="", yaxt="n", mgp=c(1.8,0.5,0), cex.lab=1.2)
	axis(2, at=ytck, labels=ytcklab, las=1, cex=0.8)
	mtext("Median Ft",side=2, line=2.5, cex=1.2)
	lines(ssb, fp, lwd=2,col="green4")
	#lines(sbci[1:n],ftci, type="o")
	lines(sbci[1:n],ftci,col="grey")
	cr = colorRampPalette(colors=c("cyan","blue","darkblue"))
	points(sbci[1:n],ftci,col=cr(n),pch=20,cex=1.2)
	points(sbci[n],ftci[n],col="darkblue",bg="pink",pch=21,cex=1.5)
	gletter(2)
	saveFig("fig.pa.phase")
}

### EQUILIBRIUM YIELD
fig.i <- function()
{
	#plot the equilibrium yield curves  
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))
	if (is.null(A$equil)) {
		cat("WARNING (fig.i): Equilibrium data not available\n"); return(invisible("No equlibrium data")) }
	#A$equil comes from the ccam.rep file
	fe <- A$equil[, 2]
	ye <- A$equil[, 3]
	sde <- A$equil[, 5]
	spr <- A$equil[, 7]

	par(mfcol=c(2,2))
	plot(fe, ye, type="l", xlab="Fishing mortality (Fe)", ylab="Equilibrium yield")
	addScenLab()
	gletter(1)
	plot(sde, ye, type="l", xlab="Spawning depletion", ylab="Equilibrium yield", lty=2, col=2)
	gletter(2)
	plot(spr,ye, type="l", xlab="Spawning potential ratio", ylab="Equilibrium yield", lty=3, col=3)
	gletter(3)
	matplot(cbind(fe, sde, spr), ye/max(ye)*100, type="l",xlab="Fe, depletion,  SPR",  ylab="Relative equilibrium yield")
	gletter(4)
	saveFig("fig.equil.yield")
}

### EQUILIBRIUM FISHING MORTALITY
fig.j <- function()
{
	#Relationship between fishing mortlaity ~ yield,  recruitment,  SBe,  SPR
	op <- par(no.readonly=T)
	on.exit(expandGraph(mfrow=c(1,1),new=FALSE))

	if (is.null(A$equil)) {
		cat("WARNING (fig.j): Equilibrium data not available\n"); return(invisible("No equlibrium data")) }
	fe <- A$equil[, 2]
	ye <- A$equil[, 3]
	sde <- A$equil[, 5]
	re <- A$equil[, 6] 
	spr <- A$equil[, 7]

	ix <- c(min(which(ye==max(ye))),min(which(sde<=0.4)) , min(which(spr<=0.4)))
	par(mfcol=c(2,2))

	plot(fe, ye, type="l",xlab="", ylab="Equilibrium yield (million mt)", lwd=2)
	segments(fe[ix],0,fe[ix],ye[ix],lty=c(1, 2, 3))
	segments(0,ye[ix],fe[ix],ye[ix],lty=c(1, 2, 3))
	addScenLab()

	re <- re/re[1]
	plot(fe, re, type="l",xlab="", ylab="Relative recruitment", lwd=2)
	segments(fe[ix],0,fe[ix],re[ix],lty=c(1, 2, 3)) 
	segments(0,re[ix],fe[ix],re[ix],lty=c(1, 2, 3)) 

	plot(fe, sde, type="l",xlab="", ylab="Spawning depletion", lwd=2)
	segments(fe[ix],0,fe[ix],sde[ix],lty=c(1, 2, 3))
	segments(0,sde[ix],fe[ix],sde[ix],lty=c(1, 2, 3))

	plot(fe, spr, type="l",xlab="", ylab="Spawning Potential Ratio", ylim=c(0, 1), lwd=2)
	segments(fe[ix],0,fe[ix],spr[ix],lty=c(1, 2, 3))
	segments(0,spr[ix],fe[ix],spr[ix],lty=c(1, 2, 3))   

	legend("topright", c("MSY", "SB40", "SPR40"), lty=1:3, bty="n")
	mtext("Equilibrium fishing mortality rate", 1, outer=T, line=-1)                               
	saveFig("fig.equilibrium")
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

table.b <- function() {
	mcbt <- A$mc.sbt
	post.bt <- as.data.frame(window(mcmc(mcbt),start=Burn+1,thin=Thin))
	btci <- t(apply(post.bt,2,quantile,probs=quants3))
	
	filename <- paste(tabDir,"tableb_sbt.tex",sep="")
	filenamecsv <- paste(tabDir,"tableb_sbt.csv",sep="")
	t.b <- cbind(A$yrs,btci[1:nyrs,])
	
	cap="Recent trends in estimated female spawning stock biomass (million mt) 
		based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","Female biomass")
	ncgrp	<- c(1,3)
	colnames(t.b)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.b, 10),file=filename,caption=cap,label="tableb",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.b,file=filenamecsv, row.names=F)
  cat(paste("Saved table ",filename,"...\n",sep=""))
  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
  
}	

table.c <- function() {
	
	mcdt <- A$mc.sbdepletion
	post.dt  <- as.data.frame(window(mcmc(mcdt),start=Burn+1,thin=Thin))
	dtci <- t(apply(post.dt,2,quantile,probs=quants3))
	
	filename<-paste(tabDir,"tablec_depletion.tex",sep="")
	filenamecsv<-paste(tabDir,"tablec_depletion.csv",sep="")
	t.c	<- cbind(A$yrs,dtci[1:nyrs,])
	
	cap="Recent trends in estimated spawning depletion level
		based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","Depletion")
	ncgrp	<- c(1,3)
	colnames(t.c)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.c, 10),file=filename,caption=cap,label="tablec",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.c,file=filenamecsv, row.names=F)
#  cat(paste("Saved table ",filename,"...\n",sep=""))
#  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
}	

table.d <- function() {	
	mcrt <- A$mc.rt 
	post.rt  <- as.data.frame(window(mcmc(mcrt),start=Burn+1,thin=Thin))
	rtci <- t(apply(post.rt,2,quantile,probs=quants3))
	
	filename<-paste(tabDir,"tabled_recruits.tex",sep="")
	filenamecsv<-paste(tabDir,"tabled_recruits.csv",sep="")
	t.d	<- cbind(ryr,rtci[1:length(ryr),])
	
	cap="Recent trends in estimated recruitment (billions of age 1 fish)
		based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","Recruits")
	ncgrp	<- c(1,3)
	colnames(t.d)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.d, 10),file=filename,caption=cap,label="tabled",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.d,file=filenamecsv, row.names=F)
  cat(paste("Saved table ",filename,"...\n",sep=""))
  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
}

table.e1 <- function() {
	spr <- A$mc.sprstatus_fmsy
	post.spr <- as.data.frame(window(mcmc(spr),start=Burn+1,thin=Thin))
	sprci <- t(apply(post.spr,2,quantile,probs=quants3))
	
	filename<-paste(tabDir,"tablee_sprfmsy_status.tex",sep="")
	filenamecsv<-paste(tabDir,"tablee_sprfmsy_status.csv",sep="")
	t.e	<- cbind(A$yr,sprci)
	
	cap="Recent trends in (1-spr)/(1-spr at fmsy) based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","(1-spr)/(1-spr at fmsy)")
	ncgrp	<- c(1,3)
	colnames(t.e)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.e, 10),file=filename,caption=cap,label="tablee",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.e,file=filenamecsv, row.names=F)
  cat(paste("Saved table ",filename,"...\n",sep=""))
  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
}

table.e2 <- function() {	
	spr <- A$mc.sprstatus_f40
	post.spr <- as.data.frame(window(mcmc(spr),start=Burn+1,thin=Thin))
	sprci <- t(apply(post.spr,2,quantile,probs=quants3))
	
	filename<-paste(tabDir,"tablee_sprf40_status.tex",sep="")
	filenamecsv<-paste(tabDir,"tablee_sprf40_status.csv",sep="")
	t.e	<- cbind(A$yr,sprci)
	
	cap="Recent trends in (1-spr)/(1-spr at f40) based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","(1-spr)/(1-spr at f40)")
	ncgrp	<- c(1,3)
	colnames(t.e)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.e, 10),file=filename,caption=cap,label="tablee",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.e,file=filenamecsv, row.names=F)
  cat(paste("Saved table ",filename,"...\n",sep=""))
  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
}	

table.f <- function() {
	
	bt3 <- A$mc.bt3
	ct <- A$ct[1,]
	cbt3 <- ct/bt3[, 1:length(A$yr)]
	post.cbt3 <- as.data.frame(window(mcmc(cbt3),start=Burn+1,thin=Thin))
	cbt3ci <- t(apply(post.cbt3,2,quantile,probs=quants3))
	
	filename<-paste(tabDir,"tablef_exploitationFraction.tex",sep="")
	filenamecsv<-paste(tabDir,"tablef_exploitationFraction.csv",sep="")
	t.f	<- cbind(A$yr, cbt3ci)
	
	cap="Recent trends in Ct/Bt3+ based on 5000 systematic samples from the joint posterior distribution."
	cgrp <- c(" ","Ct/Bt3+")
	ncgrp	<- c(1,3)
	colnames(t.f)	<- c("Year",c("2.5%","median","97.5%")) 
		
  dum <- latex(tail(t.f, 10),file=filename,caption=cap,label="tablef",rowname=NULL,cgroup=cgrp,n.cgroup=ncgrp)
  write.csv(t.f,file=filenamecsv, row.names=F)
  cat(paste("Saved table ",filename,"...\n",sep=""))
  cat(paste("Saved table ",filenamecsv,"...\n",sep=""))
}	

table.h <- function(mle=T,tableType="ssb",perc=c(0.25,0.75),stock="Female",catchFactor=1e6,writeCSV=T){
  # tableDir:  absolute directory where the files will be stored
  # A:         an objects returned from reptolist, i.e. list of the rep file entries
  # Assumes A$mc.for exists and has the case sensetive column names: Year, CtStream, Rt, Sbt, f40spr, depletion, OY
  # streams: a vector of catch weights to include in the table. If a stream is less than zero, the OY will be used as the catch stream
  #          value in the table.  These f-based catch streams will be appended to the table in the order they were
  #          entered in the catchStreams vector.
  #          -They are found in the mcmc forecast file (tinss.for) in the CtStream column
  # mle: if True, use mle outputs.  If false, use MCMC output.
  # tableType = "ssb" means make table using spawning stock biomass (sbt column in A$mc.for) - makes es.table.h.1.csv
  # tableType = "depletion" means make table using reletive deletion (depletion column in A$mc.for) - makes es.table.h.2.csv
  # tableType = "f40spr" means make table using reletive spawning potential ratio (1-spr)/(1-0.4) (f40spr column in A$mc.for) - makes es.table.h.3.csv
  # stock is either "Female" or "All".  If Female the biomass is divided by 2.
  # catchFactor:  number to multiply catch by so it appears in the table in a nicer format
  # perc is a vector of the locations of the biomass vector to split the catch stream data on.

  #splits <- c(0,sort(perc),1)
  splits <- c(0,1)
  sections <- 1:(length(splits)-1)
  tmpMeds <- vector("numeric",length=length(sections))
  if(mle){
    forc <- A$mlefor
  }else{
    forc <- A$mcfor
  }
  years <- unique(forc$Year)
  streams <- unique(forc$CtStream)
  
  #medians <- as.data.frame(matrix(nrow=0,ncol=length(splits)+1))  # medians for the fixed catch streams
  #fmsyMedians <- as.data.frame(matrix(nrow=0,ncol=length(splits)+1)) # medians for the f-based catch streams
  #f40Medians <- as.data.frame(matrix(nrow=0,ncol=length(splits)+1)) # medians for the f-based catch streams
  #ssMedians <- as.data.frame(matrix(nrow=0,ncol=length(splits)+1)) # medians for the SS-OY-based catch streams
  
  medians <- as.data.frame(matrix(nrow=0,ncol=4))  # medians for the fixed catch streams
  fmsyMedians <- as.data.frame(matrix(nrow=0,ncol=4)) # medians for the f-based catch streams
  f40Medians <- as.data.frame(matrix(nrow=0,ncol=4)) # medians for the f-based catch streams
  ssMedians <- as.data.frame(matrix(nrow=0,ncol=4))# medians for the SS-OY-based catch streams

  for(year in years){
    for(stream in streams){
      yStream <- subset(forc,forc$Year==year & forc$CtStream==stream)
      specialStream <- median(as.numeric(yStream$OY))
      ctstream<-stream

      #colnames(medians) <- c("Year","CtStream",paste(splits[1:(length(splits)-1)]," - ",splits[2:length(splits)]))
      #colnames(fmsyMedians) <- c("Year","CtStream",paste(splits[1:(length(splits)-1)]," - ",splits[2:length(splits)]))
      #colnames(f40Medians) <- c("Year","CtStream",paste(splits[1:(length(splits)-1)]," - ",splits[2:length(splits)]))
      #colnames(ssMedians) <- c("Year","CtStream",paste(splits[1:(length(splits)-1)]," - ",splits[2:length(splits)]))
      
      colnames(medians) <- c("Year","CtStream","ABC","Median")
      colnames(fmsyMedians) <- c("Year","CtStream","ABC","Median")
      colnames(f40Medians) <- c("Year","CtStream","ABC","Median")
      colnames(ssMedians) <- c("Year","CtStream","ABC","Median")
            
      if(mle){
        sorted <- yStream
      }else{
        nsamp <- length(yStream[,1])
        yStream <- yStream[(Burn+1):nsamp,]
        sorted <- yStream #[order(yStream$Rt_2009),]
      }
        
      lenSorted <- nrow(sorted)
      for(section in sections){
        lims <- splits[section:(section+1)]
        lowerInd <- floor(lenSorted*as.numeric(lims[1]))
        upperInd <- floor(lenSorted*as.numeric(lims[2]))
        if(tableType=="ssb"){
          if(mle){
            tmpMeds[section] <- sorted$SBt            
          }else{
            tmpMeds[section] <- median(as.numeric(sorted[lowerInd:upperInd,]$Sbt))
          }
        }
        if(tableType=="depletion") {
          if(mle){
            tmpMeds[section] <- sorted$depletion
          }else{
            tmpMeds[section] <- median(as.numeric(sorted[lowerInd:upperInd,]$depletion))
          }
        }
        if(tableType== "f40spr") {
	          if(mle){
	            tmpMeds[section] <- sorted$SPR40status
	          }else{
	            tmpMeds[section] <- median(as.numeric(sorted[lowerInd:upperInd,]$SPR40status))
	          }
        }
        if(section==length(sections)){
          if(stream>=0){
            medians <- rbind(medians,as.numeric(c(year,ctstream,as.numeric(stream)*catchFactor,tmpMeds)))
          }
          if(stream==.FMSYFORCFLAG){
            fmsyMedians <- rbind(fmsyMedians,as.numeric(c(year,ctstream,specialStream*catchFactor,tmpMeds)))
          }
          if(stream==.F40FORCFLAG){
            f40Medians <- rbind(f40Medians,as.numeric(c(year,ctstream,specialStream*catchFactor,tmpMeds)))
          }
         if(stream==.SSFORCFLAG){
           ssMedians <- rbind(ssMedians,c(year,ctstream,specialStream*catchFactor,tmpMeds))
          }
        }
      }
    }
  }
  
  sort.medians.for.output <- medians[order(medians$CtStream,medians$Year),]
  sort.medians.for.output <- rbind(sort.medians.for.output,fmsyMedians)  # add fmsy-based catch streams
  sort.medians.for.output <- rbind(sort.medians.for.output,f40Medians)  # add f40-based catch streams
  sort.medians.for.output <- rbind(sort.medians.for.output,ssMedians)  # add SS-OY-based catch streams

  if(writeCSV){
    if(tableType=="ssb")
      fn <- paste(tabDir,"table.h1.ssb",sep="")
    if(tableType=="depletion")
      fn <- paste(tabDir,"table.h2.depletion.",sep="")
    if(tableType=="f40spr")
      fn <- paste(tabDir,"table.h3.f40spr",sep="")
    if(mle){
      fn <- paste(fn,"_mle.csv",sep="")
    }else{
      fn <- paste(fn,"_mcmc.csv",sep="")
    }
    # remove NAs from output
    sort.medians.for.output <- sort.medians.for.output[!is.na(sort.medians.for.output$CtStream),]
    sort.medians.for.output <- sort.medians.for.output[!is.na(sort.medians.for.output$Year),]
    write.csv(sort.medians.for.output,file=fn,row.names=F)
    cat(paste("Saved table ",fn,"...\n",sep=""))
  }
}

table.i <- function(){
	#This is the final summary table for the executive summary
	lbl <- c("Unfished SBo (million mt)", 
           "Unfished age-1 recruits (billions)", 
           "\\underline{\\emph{\\textbf{REFERENCE POINTS based on SB$_{40\\%}$}}}", 
           "MSY proxy spawning biomass SB$_{40\\%}$", 
           "SPR resulting in SB$_{40\\%}$", 
           "Exploitation fraction (ct/Bt3) resulting in SB$_{40\\%}$", 
           "Yield with SB$_{40\\%}$", 
           "\\underline{\\emph{\\textbf{REFERENCE POINTS based on SPR$_{40\\%}$}}}", 
           "Spawning biomass at SPR$_{40\\%}$", 
           "SPR", 
           "Exploitation fraction (ct/Bt3) resulting in SPR$_{40\\%}$", 
           "Yield with SPR$_{40\\%}$", 
           "\\underline{\\emph{\\textbf{REFERENCE POINTS based on MSY}}}", 
           "Spawning biomass at MSY", 
           "SPR at MSY", 
           "Exploitation fraction (ct/Bt3) at MSY", 
           "MSY")
  
  mcs <- A$mcRefPoints
  post.mcs <- as.data.frame(window(mcmc(mcs),start=Burn+1,thin=Thin))
  mcsci <- t(apply(post.mcs,2,quantile,probs=quants3))
  ti <- mcsci
  ti <- rbind(ti[1:2, ], rep("", 3), ti[3:6, ], rep("", 3), ti[8:11, ], rep("", 3), ti[12:15, ])
  ti <- cbind(lbl, (ti))
  colnames(ti) <- c("Quantity", "2.5\\% percentile", "Median", "97.5\\% percentile")
  fn <- paste(tabDir,"table.i.csv",sep="")
  write.csv(ti,file=fn, row.names=F)
}

table.g1 <- function(){
  s <- read.csv(file="HakeManagement.csv")
  cap <- "Recent trend in Hake Management performance."
	filename <- paste(tabDir,"tablec.tex",sep="")
  dimnames(s)[[2]] <- c("Year","Landings","OY(mt)","ABC(mt)","Landings/OY(\\%)")
	latex(s,file=filename,caption=cap,label="tablec",rowname=NULL)  
}

table.g2<-function(){
	#This is a sideways table
  s <- read.csv(file="Summarytable.CSV")
	cap <- "Summary of recent trends in Pacific hake exploitation and stock levels."
	filename <- paste(tabDir,"tableg.tex",sep="")
  latex(s[, c(1, 2:11)],file=filename,caption=cap,label="tableg",rowname=NULL)
}

