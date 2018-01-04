## PJS crib table (N scenarios)
## Last modified by RH (170324)
make.cribtab = function(fdPrefix)
{
	N=30
	cribtab = data.frame(
	run   = 1:N,
	sens  = 1:N,
	label =
		c("M.30+k3", "M.20+k3", "M.30+k2", "M.20+k2", "M.30+k4", 
		"M.30+k3+NvB", "M.30+k3+sMWi", "M.35+k3+sMWi", "M.30+k4+sMWi", "M.35+k4+sMWi",
		"M.35+k3+sMWi-CPUE", "M.30+k4+sMWi-CPUE", "M.35+k3+sMWi+DFOvB", "M.30+k4+sMWi+DFOvB", "M.30+k3+sMWu",
		"M.35+k3+sMWu", "M.30+k4+sMWu", "M.35+k4+sMWu", "M.30+k3+sMWs", "M.35+k3+sMWs",
		"M.30+k4+sMWs", "M.35+k4+sMWs", "M.35+k3+sMWu-CPUE", "M.30+k4+sMWu-CPUE", "M.35+k3+sMWu+DFOvB",
		"M.30+k4+sMWu+DFOvB", "M.35+k3+sMWu+th4", "M.35+k3+sMWu+th14", "M.30+k4+sMWs+th4", "M.30+k4+sMWs+th14"
		),
	## sMWi = standardised mean weight incomplete, sMWu = std mean wt unsorted, sMWs = std mean wt unsorted
	SG    = c(0, rep(1,5), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4)),  ## Sensitivity Group
	col   = c("gainsboro", paste0("red",1:4), "orange", paste0("gold",1:4), paste0("seagreen",1:4), paste0("skyblue",1:4), paste0("pink",1:4), paste0("orangered",1:4), paste0("plum",1:4))
	)
	row.names(cribtab)=paste0("assess",pad0(cribtab$run,2))
	cribtab$Slab[!is.na(cribtab$sens)] = paste0("S",pad0(cribtab$sens[!is.na(cribtab$sens)],2))
	kpos=regexpr("k",cribtab$label)+1
	cribtab$k = as.numeric(substring(cribtab$label,kpos,kpos))
	Mpos=regexpr("M",cribtab$label)+1
	cribtab$M = as.numeric(substring(cribtab$label,Mpos,Mpos+2))
	cribtab$sigO = c(rep(0.2,N))
	cribtab$sigR = c(rep(0.6,N))
	cribtab$sigW = c(rep(0.15,N))
	cribtab$t1   = c(rep(1967,N))
	cribtab$I    = c(rep("1:7",10),rep("1:6",2),rep("1:7",10),rep("1:6",2),rep("1:7",6))  ## Index for surveys/CPUE used
	cribtab$Rlab = cribtab$RG = rep(NA,nrow(cribtab))  ## Runs used for composite reference scenario
	RGI = paste0("assess",pad0(c(16),2))
	cribtab[RGI,"RG"] = 1:length(RGI)
	cribtab[RGI,"Rlab"] = paste0("R",pad0(1:length(RGI),2))
	cribtab$burn = rep(30,N)
	attr(cribtab, "indexSeries") = c("GB Reed", "HS Assemblage", "WCVI Synoptic", "QCS Synoptic", "HS Synoptic", "WCHG Synoptic", "WAP CPUE")

	## Northern WAP -------------------------------------------------------------
	clrs = c("blue","gold","seagreen","red","purple")
	if (fdPrefix=="Nassess") {
		row.names(cribtab)=paste0("Nassess",pad0(cribtab$run,2))
		cribtab$label =
			c("M.30+k3" , "M.35+k3", "M.30+k3-GIG", "M.30+k4", "M.30+k4-GIG",
			"M.30+k4-GIG-CPUE", "M.25+k4-GIG", "M.30+k3-GIG-CPUE", "M.25+k3-GIG", "M.35+k4-GIG",
			"M.35+k3-GIG", "M.35+k4", "M.25+k5-GIG", "M.30+k5-GIG", "M.35+k5-GIG",
			"M.25+k3", rep("TBA", 14)
			)
		cribtab$run  = c( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16, rep(NA,14))
		cribtab$sens = c( 0, 3, 8, 2, 9,15, 6,14, 5,12,11, 4, 7,10,13, 1, rep(NA,14)) ## Run 1 = base run
		cribtab$ssrd = c( 0, 3, 7, 2, 8,11, 6,10, 5, 9,NA, 4,NA,NA,NA, 1, rep(NA,14)) ## Sensitivity runs used in the RD
		cribtab$SG   = c( 0, 1, 2, 1, 2, 3, 2, 3, 2, 2, 2, 1, 2, 2, 2, 1, rep(9,14))  ## Sensitivity Group
		cribtab$RG   = c( 1, rep(0,15), rep(NA,14))  ## Reference Group -- in case we need to use a group for model averaging
		cribtab$M    = c(0.30, 0.35, rep(0.30,4), 0.25, 0.30, 0.25, rep(0.35,3), 0.25, 0.30, 0.35, 0.25, rep(NA,14))
		cribtab$k    = c(rep(3,3), rep(4,4), rep(3,2), 4, 3, 4, rep(5,3), 3, rep(NA,14))
		cribtab$I    = c(rep("1:5",2), "2:5", "1:5", "2:5", "2:4", "2:5", "2:4", rep("2:5",3), "1:5", rep("2:5",3), "1:5", rep("",14))  ## Index for surveys/CPUE used
		cribtab$col  = c("gainsboro", "lightblue3", "seagreen1","lightblue2", "seagreen2", "purple2", "gold2", "purple1", "gold1", "red2", "red1", "lightblue4", "gold3", "seagreen3", "red3", "lightblue1", rep("white",14))
		cribtab$burn = c(200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 0, 200, 0, 0, 0, 200, rep(0,14))
		## Rank quality of MCMCs (1=good, 2=acceptable, 3=poor)
		NRrank       = c(   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  12,  16) ## runs North
		Nrank        = c(1.25, 1.5,   3, 2.5, 1.5,   2,   2,   2,   3,1.25, 1.5,   3) ## MCMC ranking (subjective average RH + PJS)
		#Nrank        = c(1.25, 1.5,   3,   3,   3,   3,   3,   2,   3,   3,   3,   3) ## MCMC ranking (after RPR meeting)
		Nfmax        = c(0.71,0.51,8.01,19.4,18.9,19.9,19.0,18.4,9.97,17.5,16.1,1.69) ## median Fmax across MCMCs
		Nfbig        = c(   0,   0,   0,   7,   6,   8,   6,   1,   1,   5,   6,   0) ## no. years when median (annual) Ft > 2

		names(Nfbig) = names(Nfmax) = names(Nrank) = NRrank
		cribtab$Fbig = cribtab$Fmax = cribtab$rank = rep(NA,30)
		runchar      = as.character(cribtab$run)
		cribtab$rank = Nrank[runchar]  ## it works but how?
		cribtab$Fmax = Nfmax[runchar]
		cribtab$Fbig = Nfbig[runchar]
#browser();return()
		attr(cribtab, "indexSeries") = c("GB Reed", "HS Assemblage", "HS Synoptic", "WCHG Synoptic", "North CPUE")
	}
	## Southern WAP -------------------------------------------------------------
	if (fdPrefix=="Sassess") {  
		row.names(cribtab)=paste0("Sassess",pad0(cribtab$run,2))
		cribtab$label =
			c("M.30+k3+GoA", "M.30+k3+3CDvB", "M.30+k4+3CDvB", "M.30+k3+OSvB", "M.30+k4+OSvB",
			"M.30+k3+SGvB", "M.30+k4+SGvB", "M.30+k4+OSvB-CPUE", "M.25+k4+OSvB", "M.30+k3+OSvB-CPUE",
			"M.25+k3+OSvB", "M.35+k3+OSvB", "M.35+k4+OSvB", "M.25+k5+OSvB", "M.30+k5+OSvB", 
			"M.35+k5+OSvB", "M.30+k4+GoA",
			rep("TBA", 13)
			)
		cribtab$run  = c( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17, rep(NA,13))
		cribtab$sens = c(11,13,14, 0, 1,15,16,10, 4, 9, 3, 6, 7, 5, 2, 8,12, rep(NA,13)) ## Run 4 = base run
		cribtab$ssrd = c(NA,NA,NA, 0, 1,NA,NA,10, 4, 9, 3, 6, 7, 5, 2, 8,NA, rep(NA,13)) ## Sensitivity runs used in the RD
		cribtab$SG   = c( 4, 4, 4, 0, 1, 4, 4, 1, 2, 1, 2, 3, 3, 2, 1, 3, 4, rep(9,13))  ## Sensitivity Group
		cribtab$RG   = c( rep(0,3), 1, rep(0,13), rep(NA,13))  ## Reference Group -- in case we need to use a group for model averaging
		cribtab$M    = c(rep(0.3,8), 0.25, 0.30, 0.25, rep(0.35,2), 0.25, 0.30, 0.35, 0.30, rep(NA,13))
		cribtab$k    = c(rep(3,2), 4, 3, 4, 3, rep(4,3), rep(3,3), 4, rep(5,3), 4, rep(NA,13))
		cribtab$I    = c(rep("1:4",7), "1:3", "1:4", "1:3", rep("1:4",7), rep("1:4",13))  ## Index for surveys/CPUE used
		cribtab$col  = c(paste0("red",1:3),"gainsboro", "blue1", "red4", "orange1", "blue2", "gold1","blue3", "gold2", paste0("seagreen",1:2), "gold3", "blue4", "seagreen3","orange2", rep("white",13))
		cribtab$burn = c(333, 250, 0, 200, 200, 0, 0, 200, 200, 200, 200, 200, 200, 200, 200, 200, 0, rep(0,13))
		## Rank quality of MCMCs (1=good, 2=acceptable, 3=poor)
		#SRrank       = c( 4, 5, 15, 11,   9, 14,  12,  13, 16, 10, 8) ## runs South
		#Srank        = c( 2, 2,  2,  1, 1.5,  2, 2.5, 2.5,  2,  1, 3) ## MCMC ranking (subjective ave RH + PJS)
		SRrank       = c(   4,   5,   8,   9,  10,  11,  12,  13,  14,  15,  16) ## runs South
		Srank        = c(   2,   2,   3, 1.5,   1,   1,   2,   2,   2,   2,   2) ## MCMC ranking (subjective average RH + PJS)
		#Srank        = c(   2,   2,   3, 1.5,   3,   1,   2,   2,   3,   3,   3) ## MCMC temp. ranking for RPR meeting
		Sfmax        = c(0.28,18.3,18.4,18.3,19.2,0.49,0.12,14.2,19.4,17.9,19.7) ## median Fmax across MCMCs
		Sfbig        = c(   0,   1,   2,   1,   2,   0,   0,   1,   5,   4,   3) ## no. years when median (annual) Ft > 2

		names(Sfbig)  = names(Sfmax)  = names(Srank) = SRrank
		cribtab$Fbig = cribtab$Fmax = cribtab$rank = rep(NA,30)
		runchar      = as.character(cribtab$run)
		cribtab$rank = Srank[runchar]  ## it works but how?
		cribtab$Fmax = Sfmax[runchar]
		cribtab$Fbig = Sfbig[runchar]
#browser();return()
		attr(cribtab, "indexSeries") = c("GB Reed", "WCVI Synoptic", "QCS Synoptic", "South CPUE")
	}
	cribtab$Slab[!is.na(cribtab$sens)] = paste0("S",pad0(cribtab$sens[!is.na(cribtab$sens)],2))
	cribtab$Slab[is.na(cribtab$sens)] = NA
	cribtab = cribtab[!is.na(cribtab$run),]
	write.csv(cribtab,"cribtab.csv")
	return(cribtab)
}
#cribtab=make.cribtab("Sassess")


