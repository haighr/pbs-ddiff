#**********************************************************************************
# iscamLoadScenarios.r
# This file contains the code to load multiple ADMB scenarios into an 'opList'
# All filenames and foldernames for iscam are set here.
#
# Author            : Chris Grandin/Robyn Forrest
# Development Date  : December 2011 - January 2012
#
#**********************************************************************************

### Set up an PBSmodelling options manager (RH)
if (!exists("wapmod",envir=.GlobalEnv) || (!exists("usingSweave",envir=.GlobalEnv) || !usingSweave))  ## This is set to TRUE in `modelResults.Rnw'
wapmod = new("PBSoptions",filename="wapmod.txt",
	#initial.options=list(fdPrefix="assess", fnInput="wapassessZZZdd", model="iscamdelaydiff"))
	initial.options=list(fdPrefix="Nassess", fnInput="wapNassessZZZdd", model="iscamdelaydiff"))

## To change the wapmod options only if the Results Sweave is being run
if (exists("wapmod",envir=.GlobalEnv) && exists("fdPrefix",envir=.GlobalEnv) && exists("usingSweave",envir=.GlobalEnv) && usingSweave)  ## This is set to TRUE in `modelResults.Rnw'
	setOptions(wapmod, fdPrefix=fdPrefix, fnInput=paste0("wap",fdPrefix,"ZZZdd"))

cribtab = make.cribtab( getOptions(wapmod,"fdPrefix"))

# Set up folder structure - prefix fd means foldername - all folders here must have trailing '/'
fdADMB          <- "../src/"  # location of the TPL and EXE for iscam.
fdControlFiles  <- "../src/"  # location to which the admb .CTL and .DAT files are copied
errWin          <- try(focusWin("iscamGui"),silent=TRUE)
if (!exists("usingSweave",envir=.GlobalEnv) || !usingSweave)
	fdScenarios     <- if (!is.null(errWin)) "../Scenarios_North/" else getWinVal()$fdScenarios
fdFigures       <- "figs/"
fdTables        <- "tables/"
dirList         <- dir(fdScenarios)

# Model executable name
exModel                 <- "iscamdelaydiff" # don't use iscam.exe - it will not work on some machines!

# prefix fn means filename
# ADMB iscam input filenames
fnISCAM                 <- "iscamdelaydiff.dat"
fnCTL                   <- paste0(getOptions(wapmod,"fnInput"),".ctl")
fnDat                   <- paste0(getOptions(wapmod,"fnInput"),".dat")
fnPFC                   <- paste0(getOptions(wapmod,"fnInput"),".pfc")

# ADMB iscam output filenames written to the rscripts directory
#fnForecastMCMC         <- "iscam_forecast_mcmc.rep"
#fnForecastMLE          <- "iscam_forecast_mle.rep"
fnFit                   <- "iscamdelaydiff"
fnRecruitment           <- "rt.mcmc"
fnSBiomass              <- "sbt.mcmc"
fnFishingMort           <- "ft.mcmc"
#fnSpawnersPR           <- "1-spr.mcmc"
#fnSPRStatus            <- "sprmsy_status.mcmc"
#fnSPR40                <- "spr40_status.mcmc"
#fnBiomass2Plus         <- "bt2.mcmc"
fnBiomass               <- "tbt.mcmc"
fnSurveyResid           <- "surveyres.mcmc"

#GETRID OF THESE ONCE CONDITIONAL CODE BELOW IS WORKING
fnSurveyResid2 <- "surveyres2.mcmc" 
fnSurveyResid3 <- "surveyres3.mcmc"
fnSurveyResid4 <- "surveyres4.mcmc"
#only read subsequent survey residual files if they exist
#THIS CURRENTLY NOT WORKING - NEED TO WORK ON THE PATH IN FILE.EXISTS
#FOR NOW, ISCAM IS OUTPUTTING EMPTY FILES FOR THE 2,3 AND 4 SURVEY CASES
#mult_survey				<- file.exists("surveyres2.mcmc","surveyres3.mcmc","surveyres4.mcmc")
#print(mult_survey)
#break
#if(mult_survey[1]==T) fnSurveyResid2 <- "surveyres2.mcmc"
#if(mult_survey[2]==T) fnSurveyResid3 <- "surveyres3.mcmc"
#if(mult_survey[3]==T) fnSurveyResid4 <- "surveyres4.mcmc"

fnMCMC                  <- "iscamdelaydiff.mcmc"
#fnEquilibrium          <- "iscamdelaydiff.eql"
fnDepletion             <- "sbtdepletion.mcmc"
fnRecDevs               <- "recDevs.mcmc"
fnCov                   <- "admodel.cov"  # not used in the GUI
fnDep                   <- "admodel.dep"  # not used in the GUI
fnHes                   <- "admodel.hes"  # not used in the GUI
fnSims                  <- "sims"         # not used in the GUI
fnEigRpt                <- "eigv.rpt"     # not used in the GUI
fnFminlog               <- "fmin.log"     # not used in the GUI
fnVariance              <- "variance"     # not used in the GUI

# ADMB iscam output filenames written to the admb directory
fnRep                   <- "iscamdelaydiff.rep"
fnRetro                 <- "iscamdelaydiff.ret"    # REP file template for a retrospective run
                                          # these will have numbers appended for each
                                          # retorospective year, i.e. iscam.ret5 means
                                          # this is a REP file for a run using data from
                                          # beginning of time series to 5 years ago.
fnPar                   <- "iscamdelaydiff"
fnPsv                   <- "iscamdelaydiff.psv"    # not used in the GUI
fnStd                   <- "iscamdelaydiff.std"    # not used in the GUI
fnMcm                   <- "iscamdelaydiff.mcm"    # not used in the GUI
fnMc2                   <- "iscamdelaydiff.mc2"    # not used in the GUI
fnHst                   <- "iscamdelaydiff.hst"    # not used in the GUI
fnLog                   <- "iscamdelaydiff.log"    # not used in the GUI
fnEva                   <- "iscamdelaydiff.eva"    # not used in the GUI
fnEcm                   <- "iscamdelaydiff.ecm"    # not used in the GUI
fnCor                   <- "iscamdelaydiff.cor"    # not used in the GUI
fnBar                   <- "iscamdelaydiff.bar"    # not used in the GUI
fnB01                   <- "iscamdelaydiff.b01"    # not used in the GUI
fnB02                   <- "iscamdelaydiff.b02"    # not used in the GUI
fnB03                   <- "iscamdelaydiff.b03"    # not used in the GUI
fnR01                   <- "iscamdelaydiff.r01"    # not used in the GUI
fnR02                   <- "iscamdelaydiff.r02"    # not used in the GUI
fnR03                   <- "iscamdelaydiff.r03"    # not used in the GUI
fnP01                   <- "iscamdelaydiff.p01"    # not used in the GUI
fnP02                   <- "iscamdelaydiff.p02"    # not used in the GUI
fnP03                   <- "iscamdelaydiff.p03"    # not used in the GUI
fnSensitivityGroup      <- "SensitivityGroup.txt" 
#new objects for pcod
fnProjmcmc              <- "iscammcmc.proj"      # Use Chris' new function `write_proj_output' which creates a csv (see below)
fnProjmpd               <- "iscammpd.proj"

opFilesCurrDir <- c(
                    fnFit,
                    fnRecruitment,
                    fnSBiomass,
                    fnFishingMort,
                    fnBiomass,
                    fnSurveyResid,
                    fnMCMC,
                    #fnEquilibrium,
                    fnDepletion,
                    fnCov,
                    fnDep,
                    fnHes,
                    fnSims,
                    fnEigRpt,
                    fnFminlog,
                    fnVariance,
                    fnProjmcmc,
                    fnProjmpd
                     )
# Booleans for whether or not the opFilesCurrDir files are mcmc outputs
opFilesCurrDirIsMCMC <- c(T,F,F,T,T,T,T,T,T,T,F,T,T,F,F,F,F,F,F,F,T,F,T,F,F) #T,T,T, (7,8,16)

opFilesADMBDir <- c(fnRep,   
			fnDat,
                    fnPar,
                    fnPsv,
                    fnStd,
                    fnMcm,
                    fnMc2,
                    fnHst,
                    fnLog,
                    fnEva,
                    fnEcm,
                    fnCor,
                    fnBar,
                    fnB01,
                    fnB02,
                    fnB03,
                    fnR01,
                    fnR02,
                    fnR03,
                    fnP01,
                    fnP02,
                    fnP03  )

loadScenarios <- function(silent=FALSE){
  # Create empty opList object of correct dimensions
  # prefix op is short for output.
  # opList - is a list of length equal to the number of scenario runs
  #        - each scenario list contains the following (indexed by number)
  #          1. Folder name - taken directly from their folder name on the OS.
  #          2. Figure folder name
  #          3. Tables folder name
  #          4. Data Frame representing the ouput of the model (REP file from ADMB)
  #          5. list of Retrospective data (with filenames as $name) - REP file contents
  opList <<- rep(list(list("","","",NULL,NULL)),length(dirList))

  setOptions(wapmod,is.mcmc = rep(FALSE,length(dirList)))
  allPrefix = names(table(dirList))
  popPrefix = names(rev(sort(table(substring(allPrefix,1,nchar(allPrefix)-2))))[1])  # determine most used prefix
  setOptions(wapmod, fdPrefix=popPrefix, fnInput=paste0("wap",popPrefix,"ZZZdd"))
  assign("fnCTL", paste0(getOptions(wapmod,"fnInput"),".ctl"), envir=.GlobalEnv)
  assign("fnDat", paste0(getOptions(wapmod,"fnInput"),".dat"), envir=.GlobalEnv)
  assign("fnPFC", paste0(getOptions(wapmod,"fnInput"),".pfc"), envir=.GlobalEnv)
  
  assign("cribtab", make.cribtab(getOptions(wapmod,"fdPrefix")), envir=.GlobalEnv)

  for(scenario in 1:(length(dirList))){
    loadScenario(scenario=scenario, silent=silent)
  }
  if(!silent){
    .flush.cat("loadScenarios: opList now contains all model data.\n")
  }
}

loadScenario <- function(scenario, silent=FALSE){
  # Setup the directories for this scenario
  ## patch for SST which uses the fdPrefix "assess"
  prefix = getOptions(wapmod,"fdPrefix")
  runstr = substring(dirList[scenario],nchar(prefix)+1,)
  fnCTL  = sub("ZZZ",runstr,fnCTL)
  fnDat  = sub("ZZZ",runstr,fnDat)
  fnPFC  = sub("ZZZ",runstr,fnPFC)

  fdScenario <- paste0(fdScenarios,dirList[[scenario]],"/")
  tmp.mcmc = getOptions(wapmod,"is.mcmc")
  is.mcmc = tmp.mcmc[scenario] = file.exists(paste0(fdScenario,fnPsv))
  names(tmp.mcmc)[scenario] = runstr
  setOptions(wapmod,is.mcmc=tmp.mcmc)
#PBStools::.flush.cat(runstr)

  opList[[scenario]][[1]] <<- dirList[[scenario]] ## RH: why is this over-written by the following line?
  opList[[scenario]][[1]] <<- fdScenario
  opList[[scenario]][[2]] <<- paste(fdScenario,fdFigures,sep="")
  opList[[scenario]][[3]] <<- paste(fdScenario,fdTables,sep="")

  # Try creating the opList for this Scenario.  If the data are not there, it will silently
  # continue to the next scenario instead of throwing an error.
  opList[[scenario]][[4]]                  <<- try(reptoRlist(paste0(fdScenario,fnRep)),silent=silent)
  opList[[scenario]][[4]]$dat              <<- try(dattoRlist(paste0(fdScenario,fnDat)),silent=silent)
  opList[[scenario]][[4]]$pfc              <<- try(pfctoRlist(paste0(fdScenario,fnPFC)),silent=silent)
  # opList[[scenario]][[4]]$cntrl          <<- opList[[scenario]][[4]]$rep$cntrl
  opList[[scenario]][[4]]$fit              <<- try(read.fit(paste0(fdScenario,fnPar)),silent=silent)
  #opList[[scenario]][[4]]$yrs             <<- try(c(opList[[scenario]][[4]]$yr,max(opList[[scenario]][[4]]$yr)+1),silent=silent)  ## RH already loaded from rep file
  opList[[scenario]][[4]]$lastYr           <<- rev(opList[[scenario]][[4]]$yr)[1]
  currYr = opList[[scenario]][[4]]$currYr  <<- rev(opList[[scenario]][[4]]$yrs)[1]  ## RH `yrs' implicitly contains reconstructed years + 1
  ikage  = grep("#kage",names(opList[[scenario]][[4]]$dat))
  if (length(ikage)>0)  opList[[scenario]][[4]]$kage <<- as.numeric(opList[[scenario]][[4]]$dat[ikage][[1]])
  opList[[scenario]][[4]]$mc.burn              <<- 0
  if (is.mcmc) {
    #projYr = opList[[scenario]][[4]]$projYr    <<- rev(opList[[scenario]][[4]]$yrs)[1]  ## RH temporary, will be over-written
    opList[[scenario]][[4]]$mc.rt              <<- try(read.table(paste(fdScenario,fnRecruitment,sep="")),silent=silent)
    opList[[scenario]][[4]]$mc.sbt             <<- try(read.table(paste(fdScenario,fnSBiomass,sep="")),silent=silent)
    opList[[scenario]][[4]]$mc.ft              <<- try(read.table(paste(fdScenario,fnFishingMort,sep="")),silent=silent)
    opList[[scenario]][[4]]$mc.tbt             <<- try(read.table(paste(fdScenario,fnBiomass,sep="")),silent=silent)
    opList[[scenario]][[4]]$mc.resid           <<- try(read.table(paste(fdScenario,fnSurveyResid,sep="")),silent=silent)#
    opList[[scenario]][[4]]$mc.sbdepletion     <<- try(read.table(paste(fdScenario,fnDepletion,sep="")),silent=silent)
    opList[[scenario]][[4]]$mc                 <<- try(read.table(paste(fdScenario,fnMCMC,sep=""),header=TRUE),silent=silent)
    opList[[scenario]][[4]]$mc.burn            <<- cribtab[dirList[scenario],"burn"]  ## use cribtab burn entry (PJS prefers)
    #cribtab[dirList[scenario],"burn"]          <<- max(0,nrow(opList[[scenario]][[4]]$mc)-1000) ## overwrite burn entry in cribtab
#browser();return()
    if (file.exists(paste0(fdScenario,fnProjmcmc,".csv"))) {
      mcproj = try(read.csv(paste0(fdScenario,fnProjmcmc,".csv")),silent=silent)
      projnames = names(mcproj)
      ## RH executable now projects 5 years but we want only 3 years. 
      ## Start collecting proj yrs 2 fields right of `TAC' and end collecting 3 fields left of `B0'
      ## For 5 proj yrs, end collecting 1 field left of `B0'
      ## ***** maybe better to use all five and subset based on user's input of Nproj
      #projflds  = eval(parse(text=paste0("projnames[",paste0(grep("^TAC$|^B0$",projnames)+c(2,-(6-Nproj)),collapse=":"),"]")))
      projflds  = eval(parse(text=paste0("projnames[",paste0(grep("^TAC$|^B0$",projnames)+c(2,-1),collapse=":"),"]")))
      projYr = opList[[scenario]][[4]]$projYr  <<- as.numeric(gsub("B","",projflds))
      #opList[[scenario]][[4]]$yrs              <<- sort(unique(c(opList[[scenario]][[4]]$yrs,projYr)))
      #currYr = opList[[scenario]][[4]]$currYr  <<- rev(opList[[scenario]][[4]]$yr)[1]
      opList[[scenario]][[4]]$mcproj           <<- mcproj
#browser();return()
    } else {
      ### only for backwards compatibility with fixed 1-yr projection
      mcproj = try(read.table(paste(fdScenario,fnProjmcmc,sep=""),h=T),silent=silent)
      projnames = names(mcproj)
      projnames[1] = "TAC"
      #names(mcproj) = gsub("F2014","Fproj",gsub("F2013","Fcurr",gsub("B2015","Bproj",gsub("B2014","Bcurr",projnames))))
      names(mcproj) = gsub("F2014",paste0("F",currYr),gsub("F2013",paste0("F",currYr-1),gsub("B2015",paste0("B",currYr+1),gsub("B2014",paste0("B",currYr),projnames))))
      opList[[scenario]][[4]]$mcproj             <<- mcproj
    }
    if (file.exists(paste0(fdScenario,fnProjmpd,".csv"))) {
      opList[[scenario]][[4]]$mpdproj           <<- try(read.csv(paste0(fdScenario,fnProjmpd,".csv")),silent=silent)
      ## Replace the sbt current year pseudo-projection with real projection
      opList[[scenario]][[4]]$sbt.original <<- opList[[scenario]][[4]]$sbt
      opList[[scenario]][[4]]$sbt[length(opList[[scenario]][[4]]$yrs)] <<- subset(opList[[scenario]][[4]]$mpdproj,tac==0)[,paste0("B",currYr)]
    } else {
      ### only for backwards compatibility with fixed 1-yr projection
      mpdproj = try(read.table(paste(fdScenario,fnProjmpd,sep=""),h=T),silent=silent)
      projnames = names(mpdproj)
      projnames[1] = "TAC"
      #names(mpdproj) = gsub("F2014","Fproj",gsub("F2013","Fcurr",gsub("B2015","Bproj",gsub("B2014","Bcurr",projnames))))
      pnames = gsub("F2014",paste0("F",currYr),gsub("F2013",paste0("F",currYr-1),gsub("B2015",paste0("B",currYr+1),gsub("B2014",paste0("B",currYr),projnames))))
      names(mpdproj) = gsub("04BMSY","BMSY40",gsub("08BMSY","BMSY80",pnames))
      opList[[scenario]][[4]]$mpdproj             <<- mpdproj
#browser();return()
    }
		opList[[scenario]][[4]]$SensitivityGroup   <<- 
		if (file.exists(paste0(fdScenario,fnSensitivityGroup))) read.table(paste0(fdScenario,fnSensitivityGroup))[,1] else 99
    opList[[scenario]][[4]]$mc.RecDevs         <<- try(read.table(paste(fdScenario,fnRecDevs,sep=""),h=F),silent=silent)
  }
  # Try loading all the retrospective files for this scenario
  # get ret* filenames..
  dirScenario <- dir(fdScenario)
  retLoc <- grep(paste(fnRetro,"[0-9]+",sep=""),dirScenario,T)
  # retLoc gives the indexes for all ret* files
  # create a list structure for the ret* outputs
  opList[[scenario]][[5]] <<- rep(list(NULL),length(retLoc))
  if(length(retLoc>0)){
    for(retro in 1:length(retLoc)){
      opList[[scenario]][[5]][[retro]]      <<- try(reptoRlist(paste(fdScenario,dirScenario[retLoc[retro]],sep="")),silent=silent)
      opList[[scenario]][[5]][[retro]]$name <<- dirScenario[retLoc[retro]]
    }
  }
  # End of loading of retrospectives
  if(!silent){
    .flush.cat(opList[[scenario]][[1]]," - Successfully loaded.\n")
  }

  # Make sure the figures and tables directories exist, if not, create them
  tmpFdFigures <- sub("/$","",paste0(fdScenario,fdFigures))
  tmpFdTables  <- sub("/$","",paste0(fdScenario,fdTables))
  if(!file.exists(tmpFdFigures)){
    dir.create(tmpFdFigures,showWarnings=!silent)
    .flush.cat("Created the directory '",tmpFdFigures,"'\n")
  }
  if(!file.exists(tmpFdTables)){
    dir.create(tmpFdTables,showWarnings=!silent)
    .flush.cat("Created the directory '",tmpFdTables,"'\n")
  }
}
