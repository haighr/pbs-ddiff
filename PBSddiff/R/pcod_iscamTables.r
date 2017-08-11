#**********************************************************************************
# iscamTables.r
# This file contains the code for writing main body tables to disk.
# This file assumes that an object called 'opList' exists and is a valid opList as
#  described in loadScenarios.r
#
# Author            : Chris Grandin
# Development Date  : January 2012
# RF: Modified Nov 2013
# CG: Added table.mcmc.mpd November 2013
#
#**********************************************************************************

table.mcmc.mpd <- function(mcmcData,
                           burnin          = 0,
                           probs           = c(0.025,0.5,0.975),
                           mpdData         = NULL,
                           tableName       = NULL,
                           colLabels       = NULL,
                           formatOut       = "%1.6f",  # 4 decimal places shown regarless of what roundDec is
                           roundDec        = 6,        # Round to 4 decimal places
                           formatThousands = FALSE){
  # mcmcData is a years-column, mcmc-sample from posterior row matrix with the burnin not yet removed.
  # If the number of rows in mcmcData exceeds the number of years, a warning will be given and the
  #   extra rows (which are projection years) will be removed.
  # burnIn is the number of rows to remove from the beginning of the mcmcData matrix. If greater than number of
  #  rows of mcmcData, an error will be given.
  # mpdData is a vector, or a matrix (in which case only the first row will be used by this function.
  # probs are the probablities used in the quantile function through getMCMCQuantiles
  # tableName is the name of the table to be stored. If null, an error will be given. .csv will be appended to this.
  # years will be used as the names for the columns in the table. If NULL, an error will be given.
  # If th length of the years vector does not match with the number of columns of mcmcData or
  # number of elements in mpdDat, an error will be given.
  # roundDec is the number of decimal points to round to
  # formatThousands puts a comma seperator in thousands values if TRUE.
	print(tableName)

  if(is.null(tableName)){
    cat("table.mcmc.mpd: Error - You didn't supply a filename.\n")
    return()
  }
  if(is.null(colLabels)){
    cat("table.mcmc.mpd: Error - For ",tableName,", you didn't supply the column labels.\n",sep="")
    return()
  }
  if(burnin < 0 || burnin >= nrow(mcmcData)){
    cat("table.mcmc.mpd: Error - For ",tableName,", the burnin value is not in within the number of rows in mcmcData.\n",sep="")
    return()
  }
  if(ncol(mcmcData) < length(colLabels)){
    cat("table.mcmc.mpd: Error - For ",tableName,", the length of vector 'colLabels' is greater than the number of columns in mcmcData.\n",sep="")
    return()
  }
#if (tableName=="RecruitmentQuants") {browser();return()}
  if(ncol(mcmcData) > length(colLabels)){
    cat("table.mcmc.mpd: Warning - For ",tableName,", the length of vector 'colLabels' is less than the number of columns in mcmcData.\n",
        "                The extra rows are assumed to be early average recruitment samples so they have been removed.\n",sep="")
    #    "                The extra rows are assumed to be projections unwanted in the table so they have been removed.\n",sep="")
    #mcmcData <- mcmcData[,-((length(colLabels)+1):length(mcmcData))]
    mcmcData <- mcmcData[,rev(ncol(mcmcData)-(1:length(colLabels))+1)] ### discard early average recruitment data before kage
  }
  # Take first row of mpdData if it is a matrix
  if(is.matrix(mpdData)){
    mpdData <- mpdData[1,]
  }
  #if(length(mpdData) != length(colLabels)){
  #  cat("table.mcmc.mpd: Error - For ",tableName,", the length of vector 'colLabels' doesn't match the length of vector 'mpdData'.\n",sep="")
  #  return()
  #}
  if(length(mpdData) > length(colLabels)){
    cat("table.mcmc.mpd: Warning - For ",tableName,", the length of vector 'colLabels' is less than the number of records in mpdData.\n",
        "                The extra elements are assumed to be early average recruitment values so they have been removed.\n",sep="")
    mpdData <- mpdData[rev(length(mpdData)-(1:length(colLabels))+1)] ### discard early average recruitment data before kage
  }

  tableName <- paste0(tableName,".csv")
  fileName  <- paste0(tabDir,tableName)

  # Apply burnin, remove first 'burnin' rows from mcmcData
  if (burnin<1) mcmcDataBurnedIn = mcmcData
  else          mcmcDataBurnedIn = mcmcData[-(1:burnin),]
  # Get the mcmc quantiles
  quants <- sapply(as.data.frame(mcmcDataBurnedIn), quantile, probs = probs)

  # Add on the MPD output
  tableData <- rbind(quants, mpdData)
  # Make it look like a table
  tableData <- t(tableData)

  # Format the values
  tableOut <- NULL
  for(col in 1:ncol(tableData)){
    tableOut <- cbind(tableOut, rsprintf(tableData[,col], roundDec = roundDec, formatOut = formatOut, formatThousands = formatThousands))
  }

  colnames(tableOut) <- colnames(tableData)
  rownames(tableOut) <- colLabels
  write.csv(tableOut, fileName)
}

#This table calculates medians and quantiles for the control points (forecast and hindcast variables that are not affected by future (2014)TAC)
#	and performance measures (forecast variables that are affected by future (2014) TAC)
table.projections <- function(use.historic=FALSE){
	dd   <- A$mcproj #read.table(file=fn, head=T)
	tac  <- unique(dd$TAC)
	quan <- c(0.025, 0.25, 0.5, 0.75, 0.975)
	med  <- 0.5
	
	#Values that do not change with TAC
	zVal = sapply(dd[1:length(tac),],function(x){(length(sort(unique(x)))>1)})
	ctlpts = grep(FALSE,zVal)
	#ctlpts<-c(2,5,8,12,14,16,18,20,22)
	CtlPts_quan <- NULL
	
	#Values that change with TAC
	Perf_med    <- NULL
	Bproj_quan  <- NULL
	Fproj_quan  <- NULL
	
	#Values that do not change with TAC
	ctlnames=paste(colnames(dd[ctlpts]))
	d <- subset(dd, tac==0)  ## Take only values from one TAC - these control points are insensitive to TAC
	nsamp <- length(d[,1])
	d <- d[(Burn+1):nsamp,]
	jj=0
	for(i in ctlpts) {
		jj=jj+1
		CtlPts_quan=rbind(CtlPts_quan,c(ctlnames[jj],quantile(dd[,i], na.rm=T,quan)))
	}
	### RH Field names for projections
	fnBp = paste0("B",A$projYr)
	fnFp = paste0("F",A$projYr-1)
	fnBpBc = paste0("B",A$projYr[1],"B",A$currYr)
	fnFpFc = paste0("F",(A$projYr-1)[1],"F",A$currYr-1)
	fnBpBm = paste0("B",A$projYr[1],"BMSY")
	fnBpBm8 = paste0("B",A$projYr[1],"BMSY80")
	fnBpBm4 = paste0("B",A$projYr[1],"BMSY40")
	fnFpFm  = paste0("F",(A$projYr-1)[1],"FMSY")
	
	#Values that change with TAC
	for(i in tac){
		d <- subset(dd, tac==i)
		d <- d[(Burn+1):nsamp,]	  #Remove burn-in samples once subset has been done
		#print(nsamp)
		#print(length(d[,1]))
		
		#Medians
		Bproj       <- apply(d[,fnBp,drop=F],2,quantile,med,na.rm=T) # Biomass in projected year
		Fproj       <- apply(d[,fnFp,drop=F],2,quantile,med,na.rm=T) # F in projected year
		BprojBcurr  <- quantile(d[,fnBpBc],na.rm=T,med)  # Bproj/Bcurr
		FprojFcurr  <- quantile(d[,fnFpFc],na.rm=T,med)  # Fproj/Fcurr
		BprojBMSY   <- quantile(d[,fnBpBm],na.rm=T,med)  # Bproj/BMSY
		BprojBMSY80 <- quantile(d[,fnBpBm8],na.rm=T,med) # Bproj/08BMSY
		BprojBMSY40 <- quantile(d[,fnBpBm4],na.rm=T,med) # Bproj/04BMSY
		FprojFMSY   <- quantile(d[,fnFpFm],na.rm=T,med)  # Fproj/FMSY
		if (use.historic) {
			BprojBmin   <- quantile(d[,paste0("B",A$projYr[1],"Bmin")],na.rm=T,med)          # Bproj/B1971
			BprojBAvg_S <- quantile(d[,paste0("B",A$projYr[1],"BAvg_S")],na.rm=T,med)        # Bproj/BAvg 1956-2004
			BprojBAvg_L <- quantile(d[,paste0("B",A$projYr[1],"BAvg_L")],na.rm=T,med)        # Bproj/BAvg 1956-2012
			FprojFAvg_S <- quantile(d[,paste0("F",(A$projYr-1)[1],"FAvg_S")], na.rm=T,med)   # Fproj/FAvg 1956-2004
			FprojFAvg_L <- quantile(d[,paste0("F",(A$projYr-1)[1],"FAvg_L")],na.rm=T,med)    # Fproj/FAvg 1956-2012
		}
		#Quantiles
		Bprojq <- round(apply(d[,fnBp,drop=F],2,quantile,quan,na.rm=T),1)  # Biomass in projected year
 		Fprojq <- round(apply(d[,fnFp,drop=F],2,quantile,quan,na.rm=T),3)  # F in porjected year  #print(round(c(i, fspr), 2))
		
		Perf_med <- rbind(Perf_med, c(i,Bproj,Fproj,BprojBcurr,FprojFcurr,BprojBMSY,BprojBMSY80,BprojBMSY40,FprojFMSY))
		if (use.historic) 
			Perf_med = rbind(Perf_med,c(BprojBmin,BprojBAvg_S,BprojBAvg_L,FprojFAvg_S,FprojFAvg_L))
		Bproj_quan <-rbind(Bproj_quan,c(i,Bprojq))
		Fproj_quan <-rbind(Fproj_quan,c(i,Fprojq))
	}
#browser();return()
	colnames(Perf_med ) = c("TAC",fnBp,fnFp,fnBpBc,fnFpFc,fnBpBm,fnBpBm8,fnBpBm4,fnFpFm)
	if (use.historic)
		colnames(Perf_med) = c(colnames(Perf_med),"BprojBmin","BprojBAvg_S","BprojBAvg_L","FprojFAvg_S","FprojFAvg_L")
	colnames(Bproj_quan) = c("TAC",paste0(rep(fnBp,each=length(quan)),"Q",rep(quan,length(fnBp))))
	colnames(Fproj_quan) = c("TAC",paste0(rep(fnFp,each=length(quan)),"Q",rep(quan,length(fnFp))))

	write.table(Perf_med, file=paste(tabDir, "Median_Peformance_Measures.csv", sep=""), sep=",", row.names=F)
	write.table(Bproj_quan, file=paste(tabDir, "Bproj_Quantiles.csv", sep=""), sep=",", row.names=F)
	write.table(Fproj_quan, file=paste(tabDir, "Fproj_Quantiles.csv", sep=""), sep=",", row.names=F)
	write.table(CtlPts_quan, file=paste(tabDir, "Control_Points_Quantiles.csv", sep=""), sep=",", row.names=F)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#This table calculates probabilities of performance measures based on the posterior samples of performance measure quantities in pbs_pcod2013ddmcmc.proj
#	 For example: A performance measure is Bproj/0.8Bmsy. A posterior sample of length 1000 will have 1000 estimates of Bproj/0.8Bmsy
#       We are interested in P(Bproj < 0.8Bmsy, i.e. Bproj/0.8Bmsy <1
#	Therefore calculate the proportion of posterior samples for Bproj/0.8Bmsy that are less than 1
table.decision.old <- function(use.historic=FALSE){
	dd <- A$mcproj 
	tac <- unique(dd$TAC)

	#Values that change with TAC
	dtable <- NULL

	#Values that change with TAC
	for(i in tac){
		d <- subset(dd, tac==i)
		nsamp <- length(d[,1])
		d <- d[(Burn+1):nsamp,]  #Remove burn-in samples once subset has been done
		nd<-length(d[,1])
		#print(nsamp)
		#print(nd)
		#print(length(d[,1]))

		#Biomass-based performance measures
		#Want probability Bproj >  the control point  (Bcurr, Bmsy or one of the Bavgs)
		#When  Bproj > ctl, the performance measures below are > 1
		P_BprojBcurr  <- sum(d$BprojBcurr>1)/nd   #length(which(d$BprojBcurr<1))/nd   # Bproj/Bcurr
		P_BprojBMSY   <- sum(d$BprojBMSY>1)/nd    #length(which(d$BprojBMSY<1))/nd    # Bproj/BMSY
		P_BprojBMSY80 <- sum(d$BprojBMSY80>1)/nd  #length(which(d$BprojBMSY80<1))/nd  # Bproj/08BMSY
		P_BprojBMSY40 <- sum(d$BprojBMSY40>1)/nd  #length(which(d$BprojBMSY40<1))/nd  # Bproj/04BMSY
		if (use.historic) {
			P_BprojBmin   <- sum(d$BprojBmin>1)/nd    #length(which(d$BprojBmin<1))/nd    # Bproj/B1971
			P_BprojBAvg_S <- sum(d$BprojBAvg_S>1)/nd  #length(which(d$BprojBAvg_S<1))/nd  # Bproj/BAvg 1956-2004
			P_BprojBAvg_L <- sum(d$BprojBAvg_L>1)/nd  #length(which(d$BprojBAvg_L<1))/nd  # Bproj/BAvg 1956-2012
		}

		#Fishing mortality-based performance measures
		#Want probability Fproj >  the control point  (Fcurr, Fmsy or one of the Favgs)
		#When  F2015 > ctl, the performance measures below are >1
		P_FprojFcurr <- sum(d$FprojFcurr>1)/nd  #length(which(d$FprojFcurr>1))/nd  # Fproj/Fcurr
		P_FprojFMSY  <- sum(d$FprojFMSY>1)/nd   #length(which(d$FprojFMSY>1))/nd   # Fproj/FMSY
		if (use.historic) {
			P_FprojFAvg_S <- sum(d$FprojFAvg_S>1)/nd  #length(which(d$FprojFAvg_S>1))/nd  # Fproj/FAvg 1956-2004
			P_FprojFAvg_L <- sum(d$FprojFAvg_L>1)/nd  #length(which(d$FprojFAvg_L>1))/nd  # Fproj/FAvg 1956-2012
		}
		dtable <- rbind(dtable, c(i,P_BprojBcurr,P_FprojFcurr,P_BprojBMSY,P_BprojBMSY80,P_BprojBMSY40,P_FprojFMSY))
		if (use.historic)
			dtable <- rbind(dtable, c(P_BprojBmin,P_BprojBAvg_S,P_BprojBAvg_L,P_FprojFAvg_S,P_FprojFAvg_L))
	}
	colnames(dtable) <- c("TAC","P_BprojBcurr)","P_FprojFcurr","P_BprojBmsy","P_Bproj08Bmsy","P_Bproj04Bmsy","P_FprojFmsy")
	if (use.historic)
		colnames(dtable ) <- c(colnames(dtable),"P_BprojBmin","P_BprojBAvg_S","P_BprojBAvg_L","P_FprojFAvg_S","P_FprojFAvg_L")
	write.table(dtable, file=paste(tabDir, "Decision_Table_Probs.csv", sep=""), sep=",", row.names=F)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This table calculates probabilities of performance measures based on the
# posterior samples of performance measure quantities in pbs_pcod2013ddmcmc.proj
# For example: A performance measure is Bproj/0.8Bmsy. 
# A posterior sample of length 1000 will have 1000 estimates of Bproj/0.8Bmsy
# We are interested in P(Bproj < 0.8Bmsy, i.e. Bproj/0.8Bmsy <1
# Therefore calculate the proportion of posterior samples for Bproj/0.8Bmsy that are greater than 1
# Last modified -- RH 170412
#---------------------------------------------------------------------
table.decision <- function(use.historic=FALSE, useHRP=FALSE, minYr, aveYr){
	dd <- A$mcproj
	tac <- unique(dd$TAC)

	### RH Field names for projections
	fnBc = paste0("B",A$currYr)      ## Biomass -- current year
	fnBp = paste0("B",A$projYr)      ## Biomass -- projected years
	fnFc = paste0("F",A$currYr-1)    ## Fishing mortality -- current year
	fnFp = paste0("F",A$projYr-1)    ## Fishing mortality -- projected years
	fnUc = paste0("U",A$currYr-1)    ## Harvest rate -- current year
	fnUp = paste0("U",A$projYr-1)    ## Harvest rate -- projected years
	fnB0 = "B0"                      ## Equilibrium unfished biomass
	fnBm = "BMSY"                    ## Biomass at MSY
	fnUm = "UMSY"                    ## Harvest rate at MSY

	if (useHRP) Pout = c("P_Bt>LRP","P_Bt>USR","P_Bt>Bavg",paste0("P_Bt>B",A$currYr),"P_ut>uavg")
	else Pout = c("P_Bt>0.4Bmsy","P_Bt>0.8Bmsy","P_Bt>Bmsy","P_Bt>0.2B0","P_Bt>0.4B0",paste0("P_Bt>B",A$currYr),"P_ut>umsy")

	dtable = array(NA,dim=c(length(tac),length(c(A$currYr,A$projYr)),length(Pout)), 
		dimnames=list(tac,c(A$currYr,A$projYr),Pout))

	if (useHRP){
		if (missing(minYr) || missing(aveYr))
			unpackList(renderVals(aveYr=getWinVal()$aveYr,minYr=getWinVal()$minYr))
		else
			unpackList(renderVals(minYr=minYr, aveYr=aveYr))
		HRPs = calcHRP(A=A, aveYr=aveYr, Burn=Burn)  ## one function to collect all the HRP stuff (see "calcHRP.r")
		unpackList(HRPs)
	} else {
		#if (!is.element("BMSY40",names(dd)))
		dd$BMSY40 = 0.4 * dd$BMSY
		#if (!is.element("BMSY80",names(dd)))
		dd$BMSY80 = 0.8 * dd$BMSY
		#if (!is.element("B020",names(dd)))
		dd$B020 = 0.2 * dd$B0
		#if (!is.element("B040",names(dd)))
		dd$B040 = 0.4 * dd$B0
	}
	## Convert the fishing mortalities to harvest rates in the MCMC projection table
	for (i in c(fnFc,fnFp)) {
		ii = sub("F","U",i)
		dd[,ii] = (1. - exp(-dd[,i]))
	}
	## Values that change with TAC
	for(i in tac){
		ii  = as.character(i)
		j   = c(A$currYr,A$projYr)
		jj  = as.character(j)
		jjj = c(fnBc,fnBp)

		d <- subset(dd, tac==i)
		nsamp <- length(d[,1])
		d <- d[(Burn+1):nsamp,]  #Remove burn-in samples once subset has been done
		nd <- length(d[,1])
#browser();return()

		## Biomass-based performance measures
		## Want probability Bproj >  the control point  (Bcurr, Bmsy or one of the Bavgs)
		## When  Bproj > ctl, the performance measures below are > 1
		if (useHRP) {
			dtable[ii,jj,"P_Bt>LRP"]  <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=bLRPs)                   # Bproj > LRP (Bmin/Bavg)
			dtable[ii,jj,"P_Bt>USR"]  <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=bUSRs)                   # Bproj > USR (2*LRP)
			dtable[ii,jj,"P_Bt>Bavg"] <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=post.abt)               # Bproj > Bavg
			dtable[ii,jj,paste0("P_Bt>B",A$currYr)] <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,fnBc]) # Bproj > Bcurr
			dtable[ii,jj,"P_ut>uavg"]    <- apply(d[,c(fnUc,fnUp)],2,function(x,y){ sum(x>y)/length(x)},y=post.aut)   # Uproj > Uavg
#if (i==0) {browser();return()}
		} else {
			dtable[ii,jj,"P_Bt>0.4Bmsy"] <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,"BMSY40"]) # Bproj > 0.4BMSY
			dtable[ii,jj,"P_Bt>0.8Bmsy"] <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,"BMSY80"]) # Bproj > 0.8BMSY
			dtable[ii,jj,"P_Bt>Bmsy"]    <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,"BMSY"])   # Bproj > BMSY
			dtable[ii,jj,"P_Bt>0.2B0"]   <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,"B020"])   # Bproj > 0.2B0
			dtable[ii,jj,"P_Bt>0.4B0"]   <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,"B040"])   # Bproj > 0.4B0
			dtable[ii,jj,paste0("P_Bt>B",A$currYr)] <- apply(d[,jjj],2,function(x,y){ sum(x>y)/length(x)},y=d[,fnBc])   # Bproj > Bcurr
			dtable[ii,jj,"P_ut>umsy"]    <- apply(d[,c(fnUc,fnUp)],2,function(x,y){ sum(x>y)/length(x)},y=d[,"UMSY"])   # Uproj > UMSY
		}
#if (i==2000) {browser();return()}
	}
#browser();return()
	for (i in dimnames(dtable)[[3]]){
		itab = dtable[,,i]
		### RH: Note that years for ut>umsy will be incorrectly labelled one year forward, so subtract a year
		if (grepl("umsy", i) || grepl("uavg", i)) dimnames(itab)[[2]] = as.character(as.numeric(dimnames(itab)[[2]])-1)
		write.csv(itab, file=paste0(tabDir,"DTable_",gsub(">","_GT_",i),".csv"))
	}
}
#table.decision(useHRP=T)
