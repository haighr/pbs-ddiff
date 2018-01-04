#**********************************************************************************
# pcod_iscam.r
# This file contains the code for a front end GUI controller for pcod_iscam using
# Tcl/Tk windows implemented using the R package 'PBSModelling'.  The data
# structure used is an opList, which is a list of lists, see loadScenarios.r for
# details on this structure. This file assumes that an object called 'opList'
# exists and is a valid opList.
#
# Author            : Chris Grandin/Robyn Forrest
# Development Date  : December 2011 - February 2012
#
# Source this file, then call pcod_iscam() with whatever arguments you need.
#
# pcod_iscam(reload=F,copyADMBExecutables=F,silent=F)
#
# To change/create global variables, find the function assignGlobals() in this file
# 
#**********************************************************************************

removeAllExcept <- function(vars="opList", envir = .GlobalEnv) {
  # removeAllExcept()
  # removes everything in the workspace except for what's in the vars list
  # Upon finishing, the workspace will contain whatever is in the vars list
  # plus the objects 'removeAllExcept' (this function) and 'modelLoaded'
  # which tells the software that the model has already been loaded.
  # - vars - A list of objects to keep, typically just 'opList'
  # - envir - environment to clear, typically .GlobalEnv
  
  vars <- c(vars, "removeAllExcept", "so","Rcode")
  keep <- match(x = vars, table = ls(all=T,envir = envir))
#browser();return()
  if(any(is.na(keep))){
    assign("modelLoaded",FALSE,envir=.GlobalEnv)
  }else{
    rm(list=ls(all=T,envir=envir)[-keep],envir=envir)
    assign("modelLoaded",TRUE,envir=.GlobalEnv)
  }
}
if (!exists("usingSweave",envir=.GlobalEnv) || !usingSweave){  ## This is set to TRUE in `modelResults.Rnw'
	removeAllExcept(vars=NULL)
	usingSweave = FALSE
}

#require(Hmisc)
require(MASS)
require(KernSmooth)
require(MCMCpack)
require(coda)
require(PBSmodelling)  ## changes from PBSmodelling (PBStools loads PBSmodelling)
require(grDevices)

modelDir = "C:/Users/haighr/Files/Archive/Bat/"
## IMPORTANT -- each user should define a directory that contains `iscamdelaydiff.exe';
##              this directory will be added to the path seen by R's Sys.getenv()["PATH"]

options(stringsAsFactors=FALSE)
## make warnings print out when they occur instead of at the end
options(warn=1) 
if (!exists("usingSweave",envir=.GlobalEnv) || !usingSweave) {
	## Standardise the quantiles to avoid confusion
	quants3 = c(0.05, 0.50, 0.95)
	quants5 = c(0.05, 0.25, 0.50, 0.75, 0.95)
	## Choose the number of projection years (max=5)
	Nproj = 2
}

# Runtime stats constants for plotting
.MAXGRAD <- 0.0001
.FUNEVALS <- 150
.PCHCODE <- 16

source("cribtab.r")
source("calcHRP.r")
source("reptolist.r")
#source("pcod_iscamFriedEgg.r")
source("pcod_iscamExecutiveSummary.r")
source("pcod_iscamFigs.r")
source("pcod_iscamTables.r")
source("pcod_iscamSens.r")
source("pcod_iscamRetro.r")
source("pcod_iscamADMBFileControl.r")
source("pcod_iscamUtils.r")
source("pcod_iscamLoadScenarios.r")

cat("Type pcod_iscam() to start GUI\n")
cat("Optional arguments: pcod_iscam(reload=T,silent=F,copyADMBExecutables=F)\n\n")

pcod_iscam <- function(reloadScenarios=TRUE, copyADMBExecutables=FALSE, silent=FALSE){
  # pcod_iscam()
  # launches the CCAM GUI.
  # - reloadScenarios T/F - reload the data from all model output files in all scenarios.
  # - copyADMBExecutables T/F copy the admb executable from admb folder to each scenario folder.
  # - silent T/F - show messages on command line
  
  #graphics.off()  # Destroy graphics window if it exists
  closeWin() ### GUI doesn't seem to refresh properly
  loadData(reloadScenarios=reloadScenarios,
           copyADMBExecutables=copyADMBExecutables,
           silent=silent)
  .pcod_iscamGUIsetup(win="pcod_iscamGui", silent=silent)
  invisible()
}

loadData <- function(reloadScenarios=TRUE, copyADMBExecutables=FALSE, silent=FALSE){
  # loads all model output data from all scenarios.
  # - reloadScenarios T/F - reload the data from all model output files in all scenarios.
  # - copyADMBExecutables T/T copy the admb executable from admb folder to each scenario folder.
  #   this is used for the case where you want to run the model from inside its scenario folder
  #   which makes it possible to run multiple models at one time on multi-processor machines.
  # - silent T/F - show messages on command line

  if(!exists("modelLoaded",envir=.GlobalEnv)){
    modelLoaded <<- FALSE
  }
  if(reloadScenarios || !modelLoaded){
    loadScenarios(silent=silent)
    modelLoaded <<- TRUE
    if(!silent){
      cat("loadData: Loading data from model output files.\n")
    }
  }else{
    if(!silent){
      cat("loadData: Warning! Using previously loaded data for GUI.\n\n")
    }
  }
  if(copyADMBExecutables){
    copyExecutableToScenariosFolders(silent=silent)
  }
}

.pcod_iscamGUIsetup <- function(win, silent=FALSE){
  if(win=="pcod_iscamGui"){
    viewHeader <- data.frame()
    viewSensitivityGroups <- data.frame()
    for(scenario in 1:length(opList)){
      viewHeader            <- rbind(viewHeader,opList[[scenario]][[1]])
      sengroup = opList[[scenario]][[4]]$SensitivityGroup
      if (is.null(sengroup)) sengroup = 0
      viewSensitivityGroups <- rbind(viewSensitivityGroups,sengroup)
    }
    colnames(viewHeader) <- "Scenario List"
    colnames(viewSensitivityGroups) <- "Sensitivity Group"

    scenarioHeader <- cbind(viewHeader,viewSensitivityGroups)
#browser();return()
    assign("viewHeader", viewHeader,envir=.GlobalEnv)
    assign("viewSensitivityGroups", viewSensitivityGroups,envir=.GlobalEnv)
    assign("scenarioHeader", scenarioHeader,envir=.GlobalEnv)

    assign("scenarioList", as.numeric(rownames(viewHeader)), envir=.GlobalEnv)
    assign("A", opList[[1]][[4]], envir=.GlobalEnv)  ## the 4th member in each model list is the data
    
    createWin(paste(getwd(),"/",win,"Win.txt",sep=""),env=.GlobalEnv)
    
    #winList <- c(entryScenario=1)
    winList <- list(entryScenario=1, fdScenarios=fdScenarios, burn=opList[[1]][[4]]$mc.burn)
    try(setWinVal(winList), silent=silent)

    .loadPlottingLimits()
    
    # Set Base as default start model, assumes Base 
    assignGlobals(scenario=1, gui=TRUE)
  }
}

assignGlobals <- function(scenario=1, silent=FALSE, gui=TRUE){
  # assignGlobals()
  # assigns global variables used in the model, directory names,  and the data object 'A'
  # used for plotting.
  # - scenario - which scenario to set up
  # - silent T/F - show messages on command line
  # A is the pcod_iscam model output object
#browser();return()
  assign("scenarioCurr",scenario,envir=.GlobalEnv)  ## primarily for writing ALL plot and tables
  assign("A",opList[[scenario]][[4]],envir=.GlobalEnv)
  assign("figDir",opList[[scenario]][[2]],envir=.GlobalEnv)
  assign("tabDir",opList[[scenario]][[3]],envir=.GlobalEnv)
  # saveon is a toggle for writing figures to disk
  assign("saveon",FALSE,envir=.GlobalEnv)
  # model variables - ideally these are read in not hardwired like this
  assign("age",opList[[scenario]][[4]]$age,envir=.GlobalEnv)
  assign("nage",age[length(age)],envir=.GlobalEnv)
  assign("yr",opList[[scenario]][[4]]$yr,envir=.GlobalEnv)
  assign("yrs",opList[[scenario]][[4]]$yrs,envir=.GlobalEnv)
  assign("nyrs",length(yrs),envir=.GlobalEnv)
  assign("nyear",length(yr),envir=.GlobalEnv)
  assign("nyr",yr[length(yr)],envir=.GlobalEnv)
  assign("ngear",opList[[scenario]][[4]]$ngear,envir=.GlobalEnv)
  assign("delaydiff",opList[[scenario]][[4]]$delaydiff,envir=.GlobalEnv)
  assign("Burn", opList[[scenario]][[4]]$mc.burn, envir=.GlobalEnv)
  # catch streams based on estimated parameters
  # assign(".FMSYFORCFLAG",-999,envir=.GlobalEnv)  # Needed by admb for catch streams (es.table.h)
  # assign(".F40FORCFLAG",-888,envir=.GlobalEnv)  # based on estimated values
  # assign(".SSFORCFLAG",-777,envir=.GlobalEnv)  # SS ABC catch stream
  # Set age comp data
  
  #@@@@@@@@@@@@@@ RF JULY 2012 additions @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #RF made this change so that the code doesn't try and read in age observations that aren't there in the ageless model (cntrl 14 == 3)
  assign("AgeLikelihood", opList[[scenario]][[4]]$cntrl[14], envir=.GlobalEnv)
#browser();return()

  if(AgeLikelihood!=3){  #do not do this if model is 'ageless'  3 is the only survey with age obs right now
    # assign("Acom_obs",opList[[scenario]][[4]]$A[which(opList[[scenario]][[4]]$A[,2]==1),],envir=.GlobalEnv)
    assign("Asurv_obs",opList[[scenario]][[4]]$A,envir=.GlobalEnv)  #[which(opList[[scenario]][[4]]$A[,2]==3),]
    # assign("Acom_est",opList[[scenario]][[4]]$Ahat[which(opList[[scenario]][[4]]$Ahat[,2]==1),],envir=.GlobalEnv)
    assign("Asurv_est",opList[[scenario]][[4]]$Ahat,envir=.GlobalEnv)
    #assign("Acom_res",opList[[scenario]][[4]]$A_nu[which(opList[[scenario]][[4]]$A_nu[,2]==1),],envir=.GlobalEnv)
    assign("Asurv_res",opList[[scenario]][[4]]$A_nu,envir=.GlobalEnv)
  }
  assign("nits",opList[[scenario]][[4]]$nits,envir=.GlobalEnv) # num survey indices

  # Set plot output type
  assign("plotType","png",envir=.GlobalEnv)
  # Set plot globals
  # mt variables control the management target line
  assign("mtLineColor","#009E73",envir=.GlobalEnv)   ## Upper Stock Reference (colour-blind bluegreen)
  assign("lrpLineColor","#CC79A7",envir=.GlobalEnv)  ## Lower/Limit Reference Point (colour-blind redpurple)
  assign("mtLineType",3,envir=.GlobalEnv)  
  assign("mtLineWidth",3,envir=.GlobalEnv)  
  #assign("mpdLineColor","#0072B2",envir=.GlobalEnv)  ## MPD lines (colour-blind blue) a bit too dark?
  assign("mpdLineColor","#56B4E9",envir=.GlobalEnv)  ## MPD lines (colour-blind skyblue) a bit too light?
  .setBurnThin(silent=silent, gui=gui)
}

.setBurnThin <- function(silent=FALSE, gui=TRUE){
  val <- getWinVal()
  if (gui) assign("Burn",val$burn,envir=.GlobalEnv)
  else     setWinVal(list(burn=Burn))
  assign("Thin",val$thin,envir=.GlobalEnv)
  assign("Nbin",val$nbin,envir=.GlobalEnv)
}

.subView <- function(silent=F){
  act <- getWinAct()
  val <- getWinVal()
  # scenarioList is a list of the scenario names (i.e. folder names)
  scenarioList <- as.numeric(rownames(viewHeader))
  if(length(act)>1)
    act <- act[1]
  if(act=="prevScenario"){
    prevScenario <- val$entryScenario-1
    if(prevScenario<as.numeric(min(scenarioList))){
      prevScenario <- as.numeric(min(scenarioList))
    }
    setWinVal(c(entryScenario=prevScenario))
    assignGlobals(prevScenario, gui=FALSE)
    .loadPlottingLimits()
  }else if(act=="nextScenario"){
    nextScenario <- val$entryScenario+1
    if(nextScenario>as.numeric(max(scenarioList))){
      nextScenario <- as.numeric(max(scenarioList))
    }
    setWinVal(c(entryScenario=nextScenario))
    assignGlobals(nextScenario, gui=FALSE)
    .loadPlottingLimits()
  }else if(act=="changeEntryScenario"){
    assignGlobals(val$entryScenario, gui=FALSE)
    .loadPlottingLimits()    
  }else if(act=="prevSens"){
    prevSens <- val$entrySensitivityGroup - 1
    setWinVal(c(entrySensitivityGroup=prevSens))  
  }else if(act=="nextSens"){
    nextSens <- val$entrySensitivityGroup + 1
    setWinVal(c(entrySensitivityGroup=nextSens))  
  }else if(act=="writePlots"){
    assignGlobals(getWinVal()$entryScenario, gui=TRUE) ## make sure you have selected Scenario
    .writePlots()
  }else if(act=="writeTables"){
    assignGlobals(getWinVal()$entryScenario, gui=TRUE) ## make sure you have selected Scenario
    .writeTables()
  }else if(act=="writeAllPlots"){
    .writeAllPlots(gui=FALSE)
  }else if(act=="writeAllTables"){
    .writeAllTables(gui=FALSE)
  }else if(act=="writeSensPlots"){
    .writeSensPlots()
  }else if(act=="writeRetroPlots"){
    .writeRetroPlots()
  }else if(act=="runCurrScenario"){
    .runCurrScenario()
  }else if(act=="changeBurnThin"){
    .setBurnThin()
  }else if(act=="changeSelectivityGroup"){

  }else if(act=="changeSensStatus"){
    .writeSensitivityGroups()
  }else if(act=="runRetros"){
    .runRetros()
  }else if(act=="runAllRetros"){
    .runAllRetros()
  }else if(act=="changeBiomassYlim"){
    # Set the sensitivity and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$biomassYlim <<- val$biomassYlim
    opList[[val$entryScenario]][[4]]$maxBiomassYlim <<- val$maxBiomassYlim
    winList <- c(biomassSensYlim=val$biomassYlim,
                 maxBiomassSensYlim=val$maxBiomassYlim,
                 biomassRetroYlim=val$biomassYlim,
                 maxBiomassRetroYlim=val$maxBiomassYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeBiomassSensYlim"){
    # Set the base and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$biomassYlim <<- val$biomassSensYlim
    opList[[val$entryScenario]][[4]]$maxBiomassYlim <<- val$maxBiomassSensYlim
    winList <- c(biomassYlim=val$biomassSensYlim,
                 maxBiomassYlim=val$maxBiomassSensYlim,
                 biomassRetroYlim=val$biomassSensYlim,
                 maxBiomassRetroYlim=val$maxBiomassSensYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeBiomassRetroYlim"){
    # Set the base and sensitivity ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$biomassYlim <<- val$biomassRetroYlim
    opList[[val$entryScenario]][[4]]$maxBiomassYlim <<- val$maxBiomassRetroYlim
    winList <- c(biomassYlim=val$biomassRetroYlim,
                 maxBiomassYlim=val$maxBiomassRetroYlim,
                 biomassSensYlim=val$biomassRetroYlim,
                 maxBiomassSensYlim=val$maxBiomassRetroYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeDepletionYlim"){
    # Set the sensitivity and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$depletionYlim <<- val$depletionYlim
    opList[[val$entryScenario]][[4]]$maxDepletionYlim <<- val$maxDepletionYlim
    winList <- c(depletionSensYlim=val$depletionYlim,
                 maxDepletionSensYlim=val$maxDepletionYlim,
                 depletionRetroYlim=val$depletionYlim,
                 maxDepletionRetroYlim=val$maxDepletionYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeDepletionSensYlim"){
    # Set the base and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$depletionYlim <<- val$depletionSensYlim
    opList[[val$entryScenario]][[4]]$maxDepletionYlim <<- val$maxDepletionSensYlim
    winList <- c(depletionYlim=val$depletionSensYlim,
                 maxDepletionYlim=val$maxDepletionSensYlim,
                 depletionRetroYlim=val$depletionSensYlim,
                 maxDepletionRetroYlim=val$maxDepletionSensYlim)
    try(setWinVal(winList), silent=silent)    
  }else if(act=="changeDepletionRetroYlim"){
    # Set the base and sensitivity ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$depletionYlim <<- val$depletionRetroYlim
    opList[[val$entryScenario]][[4]]$maxDepletionYlim <<- val$maxDepletionRetroYlim
    winList <- c(depletionYlim=val$depletionRetroYlim,
                 maxDepletionYlim=val$maxDepletionRetroYlim,
                 depletionSensYlim=val$depletionRetroYlim,
                 maxDepletionSensYlim=val$maxDepletionRetroYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeRecruitmentYlim"){
    # Set the sensitivity and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$recruitmentYlim <<- val$recruitmentYlim
    opList[[val$entryScenario]][[4]]$maxRecruitmentYlim <<- val$maxRecruitmentYlim
    winList <- c(recruitmentSensYlim=val$recruitmentYlim,
                 maxRecruitmentSensYlim=val$maxRecruitmentYlim,
                 recruitmentRetroYlim=val$recruitmentYlim,
                 maxRecruitmentRetroYlim=val$maxRecruitmentYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeRecruitmentSensYlim"){
    # Set the base and retro ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$recruitmentYlim <<- val$recruitmentYlim
    opList[[val$entryScenario]][[4]]$maxRecruitmentYlim <<- val$maxRecruitmentYlim
    winList <- c(recruitmentYlim=val$recruitmentSensYlim,
                 maxRecruitmentYlim=val$maxRecruitmentSensYlim,
                 recruitmentRetroYlim=val$recruitmentSensYlim,
                 maxRecruitmentRetroYlim=val$maxRecruitmentSensYlim)
    try(setWinVal(winList), silent=silent)
  }else if(act=="changeRecruitmentRetroYlim"){
    # Set the base and sensitivity ylimit entry boxes and check boxes
    opList[[val$entryScenario]][[4]]$recruitmentYlim <<- val$recruitmentRetroYlim
    opList[[val$entryScenario]][[4]]$maxRecruitmentYlim <<- val$maxRecruitmentRetroYlim
    winList <- c(recruitmentYlim=val$recruitmentRetroYlim,
                 maxRecruitmentYlim=val$maxRecruitmentRetroYlim,
                 recruitmentSensYlim=val$recruitmentRetroYlim,
                 maxRecruitmentSensYlim=val$maxRecruitmentRetroYlim)
    try(setWinVal(winList), silent=silent)
 }else if(act=="changeRefptSensYlim"){
      # Set the base and retro ylimit entry boxes and check boxes
      winList <- c(RefptSensYlim=val$RefptSensYlim,
                   maxRefptSensYlim=val$maxRefptSensYlim)
                   try(setWinVal(winList), silent=silent)
    }
  # Whichever radio button is selected will now be plotted for the scenario
  .doPlots()
}

.doPlots <- function(){
  # val is the value object from GetWinVal()
  val <- getWinVal()
  if (saveon) tput(val) ## store the current val in .PBSmodEnv for use by saveFig
  pType <- val$viewPlotType
  .setBurnThin()
  oldpar <- par( no.readonly=TRUE )
  if(.checkEntries()){
    #################################################    
    # Call figure code from pcod_iscamExecutiveSummary.r #
    #################################################
    if(pType=="sLandings"){
      fig.a()
    }else if(pType=="sEffort"){
       fig.effort()
   }else if(pType=="sBiomass"){
      fig.b(includeMPD=F,
            ylimit=val$biomassYlim,
            useMaxYlim=val$maxBiomassYlim,
            useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
            tac.use=val$currTAC,
            opacity=20,
            main="Spawning biomass",
            xlab="Year")
    }else if(pType=="sBiomassMPDOver"){
      fig.b(includeMPD=T,
            ylimit=val$biomassYlim,
            useMaxYlim=val$maxBiomassYlim,
            useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
            tac.use=val$currTAC,
            opacity=20,
            main="Spawning biomass",
            xlab="Year")
    }else if(pType=="sBiomassMPD"){
      fig.biomass.mpd(ylimit=val$biomassYlim,
                      useMaxYlim=val$maxBiomassYlim,
                      main="Spawning biomass",
                      xlab="Year")
    }else if(pType=="tBiomass"){
          fig.bt(includeMPD=F,
                ylimit=val$tbiomassYlim,
                useMaxYlim=val$maxtBiomassYlim,
                opacity=20,
                main="Total biomass",
                xlab="Year")
        }else if(pType=="tBiomassMPDOver"){
          fig.bt(includeMPD=T,
                ylimit=val$tbiomassYlim,
                useMaxYlim=val$maxtBiomassYlim,
                opacity=20,
                main="Total biomass",
                xlab="Year")
        }else if(pType=="tBiomassMPD"){
          fig.bt.mpd(ylimit=val$tbiomassYlim,
                          useMaxYlim=val$maxtBiomassYlim,
                          main="Total biomass",
                          xlab="Year")
    } else if(pType=="sBiomassRecruits"){
      fig.biomass.recruits(yBiomassYlim=val$biomassYlim,
                           useMaxBiomassYlim=val$maxBiomassYlim,
                           yRecruitmentYlim=val$recruitmentYlim,
                           useMaxRecruitmentYlim=val$maxRecruitmentYlim)

    }else if(pType=="sDepletion"){
      fig.c(includeMPD=F,
            ylimit=val$depletionYlim,
            useMaxYlim=val$maxDepletionYlim,
            useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
            xlab="Year",
            main="Spawning depletion")
    }else if(pType=="sDepletionMPDOver"){
      fig.c(includeMPD=T,
            ylimit=val$depletionYlim,
            useMaxYlim=val$maxDepletionYlim,
            useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
            xlab="Year",
            main="Spawning depletion")
    }else if(pType=="sDepletionMPD"){
      fig.depletion.mpd(ylimit=val$depletionYlim, useMaxYlim=val$maxDepletionYlim, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr)
    }else if(pType=="sRecruits"){
      fig.dmcmc(ylimit=val$recruitmentYlim,
            useMaxYlim=val$maxRecruitmentYlim,
            xlab="Year",
            main="Recruitment")
     }else if(pType=="sRecruitsMPD"){
       fig.d(ylimit=val$recruitmentYlim,
       useMaxYlim=val$maxRecruitmentYlim,
       xlab="Year",
    main="Recruitment")
    }else if(pType=="sSPRMSY"){
      fig.e1()
    }else if(pType=="sSPRf40"){
      fig.e2()
    }else if(pType=="sFMCMC"){
      fig.Fmcmc()
     }else if(pType=="sFMPD"){
      fig.Fmpd()
    }else if(pType=="sBtBmsy"){
      fig.h()
    }else if(pType=="sEquilYield"){
      fig.i()
    }else if(pType=="sEquilF"){
      fig.j()
    #####################################    
    # Call figure code from pcod_iscamFigs.r #
    #####################################
    }else if(pType=="sCommAgeResids1"){
      fig.comm.age.residuals1()
    }else if(pType=="sCommAgeResids2"){
      fig.comm.age.residuals2()
    }else if(pType=="sCommAgeProps"){
      fig.comm.age.props()
    }else if(pType=="sSurvAgeResids1"){
      fig.survey.age.residuals1()
    }else if(pType=="sSurvAgeResids2"){
      fig.survey.age.residuals2()
    }else if(pType=="sSurvAgeProps"){
      fig.survey.age.props()
    }else if(pType=="sCommAgePropsFit"){
      fig.comm.age.props.fit()
    }else if(pType=="sSurvAgePropsFit"){
      fig.survey.age.props.fit()
    }else if(pType=="sSurvBiomassFit"){
      fig.surveybiomass.fit()
    }else if(pType=="sSelectivities"){
      fig.selectivity()
    }else if(pType=="sCatchFit"){
      fig.catchFit()
     }else if(pType=="sWeightFit"){
      fig.weightFit()
    }else if(pType=="sFishingMortality"){
      fig.fishingMortality()
    }else if(pType=="sRecAnom"){ #RF July 2012
      fig.RecAnomMPD()
    }else if(pType=="sRecAnomMCMC"){ #RF July 2012
      fig.RecAnomMCMC()
    }else if(pType=="sannualMeanWt"){ #RF July 2012
      fig.annualMeanWt()
    }else if(pType=="sPhase"){
      fig.phase()
    }else if(pType=="sTimeVaryingSurvSel"){
      fig.time.varying.selectivity(2)
    }else if(pType=="sTimeVaryingCommSel"){
      fig.time.varying.selectivity(1)
    }else if(pType=="sCtlPtsBox" && delaydiff==1){
        if(nyr>=2012) {fig.Allcontrol.pts.Box(tac.use=val$currTAC, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr) 
      	 }else cat("No control point plot for 2005 bridging analyses\n")
    }else if(pType=="sMSYCtlPts" && delaydiff==1){
       if(nyr>=2012) {fig.MSYcontrol.pts()
       }else cat("No control point plot for 2005 bridging analyses\n")
    }else if(pType=="sHistCtlPts" && delaydiff==1){
        if(nyr>=2012) {fig.Histcontrol.pts(minYr=val$minYr, aveYr=val$aveYr)
    }else cat("No control point plot for 2005 bridging analyses\n")
    }else if(pType=="sBench" && delaydiff==1){
      if(nyr>=2012) {fig.Benchmarks()
      }else cat("No control point plot for 2005 bridging analyses\n")
    }else if(pType=="sFmax" && delaydiff==1){
      fig.Fdensity(scenario = scenarioCurr,
      opList  = opList,
      cribtab = cribtab,
      sensfld = "sens",
      type    = "box")
    }else if(pType=="sCtlPtsBox" && delaydiff==0){
         cat("No control point plot for age-structured model\n")
    }else if(pType=="sMSYCtlPts" && delaydiff==0){
      cat("No control point plot for age-structured model\n")
    }else if(pType=="sHistCtlPts" && delaydiff==0){
        cat("No control point plot for age-structured model\n")
    }else if(pType=="sBench" && delaydiff==0){
      cat("No control point plot for age-structured model\n")
    }else if(pType=="sSnailTrail"){
      fig.snail(useHRP=val$useHRP)  ##RH: taken from PBSawatea
    #################################################    
    # Call Parameter estimate code from pcod_iscamFigs.r #
    #################################################    
    }else if(pType=="sPosteriorParams"){
          fig.mcmc.priors.vs.posts(exFactor=1.0)
    }else if(pType=="sPosteriorParams2"){
          fig.mcmc.priors.vs.posts2(exFactor=1.0)
    }else if(pType=="sPosteriorParamskey"){
      fig.mcmc.priors.vs.postskey(exFactor=1.0)
    }else if(pType=="sParameterPairs"){
      fig.estimated.params.pairs()
    }else if(pType=="sParameterPairs2"){
      fig.estimated.params.pairs2()
    }else if(pType=="sParameterPairskey"){
      fig.estimated.params.pairs.key()
    }else if(pType=="sParameterPairsnologkey"){
      fig.estimated.params.pairs.no.log.key()
    }else if(pType=="sVariancePartitions"){
      fig.variance.partitions()
    }else if(pType=="sPriorsVsPosts"){
      fig.mcmc.priors.vs.posts(exFactor=1.0,showEntirePrior=T)
    }else if(pType=="sMCMCTrace"){
      fig.mcmc.trace()
    }else if(pType=="sMCMCChains"){
      dmcmc = mcmc2(A$mc[,parEstimated()],start=Burn+1,thin=Thin)
      plotChains(dmcmc,pdisc=0,axes=TRUE,between=list(x=0.15,y=0.2),col.trace=c("green","red","blue"),xlab="Sample",ylab="Cumulative Frequency")
    }else if(pType=="sMCMCAutocor"){
      fig.mcmc.autocor()
    }else if(pType=="sMCMCDensity"){
      fig.mcmc.density()
    }else if(pType=="sMCMCGeweke"){
      fig.mcmc.geweke()
    }else if(pType=="sMCMCGelman"){
      fig.mcmc.gelman()
    ###################################################    
    # Call Sensitivity plotting code from pcod_iscamSens.r #
    ###################################################
    }else if(pType=="sSensSB"){
      fig.base.vs.sens(sensitivityGroup=val$entrySensitivityGroup,
                       whichPlot="biomass",
                       useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                       ylimit=val$biomassYlim,
                       useMaxYlim=val$maxBiomassYlim,
                       offset=0.3)
    }else if(pType=="sSensD"){
      fig.base.vs.sens(sensitivityGroup=val$entrySensitivityGroup,
                       whichPlot="depletion",
                       useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                       ylimit=val$depletionYlim,
                       useMaxYlim=val$maxDepletionYlim)
    }else if(pType=="sSensRec"){
      fig.base.vs.sens(sensitivityGroup=val$entrySensitivityGroup,
                       whichPlot="recruits",
                       useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                       ylimit=val$recruitmentYlim,
                       useMaxYlim=val$maxRecruitmentYlim,
                       offset=0.3)
     }else if(pType=="sSensRefPts"){
      fig.base.vs.sens(sensitivityGroup=val$entrySensitivityGroup,
                       whichPlot="refpts",
                       useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                       ylimit=val$RefptSensYlim,
                       useMaxYlim=val$maxRefptSensYlim)                  
    ######################################################    
    # Call Retrospective plotting code from pcod_iscamRetro.r #
    ######################################################
    }else if(pType=="sRetroSB"){
      fig.retro(whichPlot="biomass",
                ylimit=val$biomassYlim,
                useMaxYlim=val$maxBiomassYlim)
    }else if(pType=="sRetroD"){
      fig.retro(whichPlot="depletion",
                ylimit=val$depletionYlim,
                useMaxYlim=val$maxDepletionYlim)
    }else if(pType=="sRetroRec"){
      fig.retro(whichPlot="recruits",
                ylimit=val$recruitmentYlim,
                useMaxYlim=val$maxRecruitmentYlim)      
    #############################################    
    # Call runtime Values code from pcod_iscamFigs.r #
    #############################################    
    }else if(pType=="sObjFuncVal"){
      plotRuntimeStats(1)
    }else if(pType=="sMaxGrad"){
      plotRuntimeStats(2)
    }else if(pType=="sFuncEvals"){
      plotRuntimeStats(3)
    }else if(pType=="sHangCodes"){
      plotRuntimeStats(4)
    }else if(pType=="sExitCodes"){
      plotRuntimeStats(5)
    }
  }
  par(oldpar)
  return(invisible())
}

.writeAllTables <- function(silent=FALSE, gui=FALSE){
  # write all tables for all scenarios to disk
  # scenarioList <- as.numeric(rownames(viewHeader))  ## already defined globally
  if(exists("scenarioList",envir=.GlobalEnv, inherits=FALSE))
    scenarioList = get("scenarioList",envir=.GlobalEnv)
  else stop("PROBLEM: scenarioList not found in Global environment")
  for(scenario in scenarioList){
    assignGlobals(scenario, gui=gui)
    .writeTables(silent=silent, gui=gui)
  }  
}

## --------------------------
## Last modified by RH 170906
## --------------------------
.writeTables <- function(silent=FALSE, gui=TRUE){
  ## write all tables for the given scenario to disk
  if(exists("scenarioCurr",envir=.GlobalEnv, inherits=FALSE))
    scenarioCurr = get("scenarioCurr",envir=.GlobalEnv)
  else stop("PROBLEM: scenarioCurr not found in Global environment")
  setWinVal(list(entryScenario=scenarioCurr))
  val <- getWinVal()
  if (!isScenMcmc()) {
    cat("WARNING (.writeTables): No MCMCs generated for this scenario\n"); return(invisible("No MCMC data")) }
  .setBurnThin(silent=silent, gui=gui)
  #Nyears = length(A$yr)

  cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
  #cat("No tables currently implemented\n")
  cat("Writing tables\n")
  cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
  if(delaydiff && nyr>=2012){
    quantProbs <- quants3 ## SST & WAP
    weightScale <- 1000.0
debug = F; if (!debug) {
    try(table.projections(), silent=silent)
    try(table.decision(useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr), silent=silent)
    ## Make Fishing mortality Quantiles and MPD table
    table.mcmc.mpd(mcmcData = A$mc.ft,
                  burnin    = Burn,
                  probs     = quantProbs,
                  mpdData   = A$ft,
                  colLabels = A$yr,
                  roundDec  = 6, # Number of decimal places
                  tableName = "FishingMortalityQuants")
    ## Make Biomass Quantiles and MPD table
    table.mcmc.mpd(mcmcData = A$mc.sbt / weightScale,
                  burnin    = Burn,
                  probs     = quantProbs,
                  mpdData   = A$sbt / weightScale,
                  colLabels = A$yrs,
                  formatOut = "%1.3f",
                  roundDec  = 6, # Number of decimal places
                  tableName = "BiomassQuants")
    ## Make Recruitment Quantiles and MPD table
    table.mcmc.mpd(mcmcData = A$mc.rt,
                  burnin    = Burn,
                  probs     = quantProbs,
                  mpdData   = A$rt,
                  #colLabels = A$yr[-(1:2)], # The 2 is because there is age-2 recruitment (P.cod)
                  colLabels = A$yr[-(1:ifelse(A$delaydiff==1,max(A$sage,A$kage-A$sage),A$sage))], # Use kage if delaydiff model for age-k recruitment
                  roundDec  = 6, # Number of decimal places
                  tableName = "RecruitmentQuants")

    ## Make Paramter Quantiles and MPD values table
    mcmcParamTable   <- cbind(exp(A$mc$log.ro), A$mc$h, exp(A$mc$log.m), exp(A$mc$log.rbar), exp(A$mc$log.rinit), A$mc$bo)
    paramNames       <- c("r0","steepness","m","rbar","rbar_init","b0")
    mpdParamVector   <- c(A$ro, A$steepness, A$m, A$rbar, A$rinit, A$sbo)
    # Add variable number of q's
    for(qInd in 1:length(A$q)){
      mcmcParamTable <- cbind(mcmcParamTable, eval(parse(text=paste0("A$mc$q",qInd))))
      mpdParamVector <- c(mpdParamVector, A$q[qInd])
      paramNames     <- c(paramNames, paste0("q",qInd))
    }
    colnames(mcmcParamTable) <- paramNames
    table.mcmc.mpd(mcmcData = mcmcParamTable,
                  burnin    = Burn,
                  probs     = quantProbs,
                  mpdData   = mpdParamVector,
                  colLabels = paramNames,
                  roundDec  = 6, # Number of decimal places
                  tableName = "ParameterQuants")
}

		### ========================================
		### Make one table with Parameters (P), Biomass-based quantities (B), and MST-based quantities (M)
		### The table is formatted to be latex-ready and sent to 'PBMQuants.csv' in the relevant Scenario's table folder.
		### When building the Model Results latex file, import the CSV and use the PBStools function `texArray'.
		### ========================================

		##  For compatibility with calcHRPs, adjust for burnin before sending to 'table.mcmc.mpd()'

		# 1. Make Paramter Quantiles and MPD values table
		mcmcParamTable <- cbind(rep(0,nrow(A$mc)),exp(A$mc$log.ro), A$mc$h, exp(A$mc$log.m), exp(A$mc$log.rbar), exp(A$mc$log.rinit))
		mpdParamVector <- c(0,A$ro, A$steepness, A$m, A$rbar, A$rinit)
		paramNames     <- c("Parameters","$R_0$","$h$","$M$","$\\bar{R}$","$\\bar{R}_\\mathrm{init}$")

		# Add variable number of q's
		for(qInd in 1:length(A$q)){
			mcmcParamTable <- cbind(mcmcParamTable, eval(parse(text=paste0("A$mc$q",qInd))))
			mpdParamVector <- c(mpdParamVector, A$q[qInd])
			paramNames     <- c(paramNames, paste0("$q_",qInd,"$"))
		}
		mcmcParamTable = mcmc2(mcmcParamTable, start=Burn+1) ## get rid of burnin now

		## Collect MPD values
		Bcurr.mpd = (rev(A$sbt))[1]  ## final year fixed in loadScenarios
		Ucurr.mpd = A$ut[length(A$yr)]

		## MSY-based MPD
		B0.mpd    = A$sbo
		msy.mpd   = A$msy
		Bmsy.mpd  = A$bmsy
		Umsy.mpd  = 1-exp(-A$fmsy)

		## HRP MPD
		aveYr     = renderVals(aveYr=val$aveYr,simplify=TRUE)
		sbt.mpd   = A$sbt; names(sbt.mpd) = A$yrs
		sbt.min   = findBmin(sbt.mpd,aveYr)
		#Bavg.mpd  = mean(A$sbt[is.element(yrs,aveYr)])
		#LRP.mpd   = min(A$sbt)
		Bavg.mpd  = sbt.min["Bavg"]
		LRP.mpd   = sbt.min["Bmin"]
		USR.mpd   = 2. * LRP.mpd
		Bcurr.Bavg.mpd = Bcurr.mpd/Bavg.mpd
		Bcurr.LRP.mpd  = Bcurr.mpd/LRP.mpd
		Bcurr.USR.mpd  = Bcurr.mpd/USR.mpd
		Uavg.mpd       = mean(A$ut[is.element(yr,aveYr)])
		Ucurr.Uavg.mpd = Ucurr.mpd/Uavg.mpd

		## Collect MCMC values
		#dmcmc = subset(A$mcproj, TAC==0)  # take the first tac (=0) from each posterior sample - the control points do not change with tac
		dmcmc = mcmc2(subset(A$mcproj, TAC==0), start=Burn+1)
		zeros = rep(0,nrow(dmcmc))
		Bcurr.mcmc = dmcmc[,paste0("B",A$currYr)]
		Ucurr.mcmc = 1-exp(-dmcmc[,paste0("F",A$lastYr)])
		## MSY-based
		B0.mcmc    = dmcmc[,"B0"]  ## same as A$mc$bo
		msy.mcmc   = dmcmc[,"MSY"]
		Bmsy.mcmc  = dmcmc[,"BMSY"]
		Umsy.mcmc  = dmcmc[,"UMSY"]

		## HRP MCMC -- use function `calcHRP' with Burn=0 in 'table.mcmc.mpd()'
		HRPs = calcHRP(A=A, aveYr=aveYr, Burn=Burn)  ## one function to collect all of the stuff below
		unpackList(HRPs)
#browser();return()
		#post.bt   = A$mc.sbt #[,1:Nyears]
		#post.avebt   = post.bt[,is.element(yrs,aveYr)]
		#Bavg.mcmc = apply(post.avebt,1,function(x){mean(x)}) ##  Bo = Bavg
		Bavg.mcmc = post.abt

		#med.bt = sapply(post.bt,median)              ## only used to determine minimum year
		#minYr  = yrs[is.element(med.bt,min(med.bt))] ## overrides GUI value or user's value
		## The following takes minimum depletion (BT/Bavg) across year(s) designated as minimum;
		##     perhaps should be the median though if only one minimum year is specified, then it makes no difference.
		#LRP.mcmc = apply(post.bt[,is.element(yrs,minYr),drop=FALSE],1,min); ## across years therefore 1000 mins (+ Burn mins)
		#USR.mcmc = 2. * LRP.mcmc
		LRP.mcmc = bLRPs
		USR.mcmc = bUSRs

		Bcurr.Bavg.mcmc = Bcurr.mcmc/Bavg.mcmc
		Bcurr.LRP.mcmc  = Bcurr.mcmc/LRP.mcmc
		Bcurr.USR.mcmc  = Bcurr.mcmc/USR.mcmc
#browser(); return()

		#post.ft   = mcmc2(A$mc.ft, start=Burn+1) #[,1:Nyears]
		#post.ht   = 1 - exp(-post.ft)                     ## harvest rate  (same as HRPs$post.ft)
		#post.ha   = post.ht[,is.element(yr,aveYr)]        ## harvest rates in years for average  (same as HRPs$post.aveut)
		#Uavg.mcmc = apply(post.ha,1,function(x){mean(x)}) ## average harvest rates (1000 MCMC samples) (same as HRPs$post.aut)
		Uavg.mcmc = post.aut                              ## (unpacked from HRPs)
		Ucurr.Uavg.mcmc = Ucurr.mcmc/Uavg.mcmc

		if (!val$useHRP) {
			# 2. Add Biomass-based values
			mcmcParamTable <- cbind(mcmcParamTable, zeros, 0.2*B0.mcmc, 0.4*B0.mcmc, B0.mcmc, Bcurr.mcmc, Bcurr.mcmc/B0.mcmc, Ucurr.mcmc)
			mpdParamVector <- c(mpdParamVector, 0, 0.2*B0.mpd, 0.4*B0.mpd, B0.mpd, Bcurr.mpd, Bcurr.mpd/B0.mpd, Ucurr.mpd)
			paramNames     <- c(paramNames, "Model-based", paste0(c("0.2","0.4",""),"$B_0$"), paste0("$B_{",A$currYr,"}",c("$","/B_0$")), paste0("$u_{",A$lastYr,"}$"))

			# 3. Add MSY-based values
			mcmcParamTable <- cbind(mcmcParamTable, zeros, 0.4*Bmsy.mcmc, 0.8*Bmsy.mcmc, Bmsy.mcmc, Bmsy.mcmc/B0.mcmc, Bcurr.mcmc/Bmsy.mcmc, msy.mcmc, Umsy.mcmc, Ucurr.mcmc/Umsy.mcmc)
			mpdParamVector <- c(mpdParamVector, 0, 0.4*Bmsy.mpd, 0.8*Bmsy.mpd, Bmsy.mpd, Bmsy.mpd/B0.mpd, Bcurr.mpd/Bmsy.mpd, msy.mpd, Umsy.mpd, Ucurr.mpd/Umsy.mpd)
			paramNames     <- c(paramNames, "MSY-based", paste0(c("0.4","0.8",""),"$B_\\mathrm{MSY}$"), paste0(c("$",paste0("$B_{",A$currYr,"}/")),"B_\\mathrm{MSY}",c("/B_0$","$")), "MSY", paste0(c("$",paste0("$u_{",A$lastYr,"}/")),"u_\\mathrm{MSY}$") )
		} else {
			# 4 Add HRP-based values
			mcmcParamTable <- cbind(mcmcParamTable, zeros, Bcurr.mcmc, Bavg.mcmc, LRP.mcmc, USR.mcmc, Bcurr.Bavg.mcmc, Bcurr.LRP.mcmc, Bcurr.USR.mcmc, Uavg.mcmc, Ucurr.Uavg.mcmc )
			mpdParamVector <- c(mpdParamVector, 0, Bcurr.mpd, Bavg.mpd, LRP.mpd, USR.mpd, Bcurr.Bavg.mpd, Bcurr.LRP.mpd, Bcurr.USR.mpd, Uavg.mpd, Ucurr.Uavg.mpd)
			paramNames     <- c(paramNames, "HRP-based", paste0("$B_{",A$currYr,"}$"), "$B_\\mathrm{avg}$", paste0(c("$\\mathrm{LRP}=","$\\mathrm{USR}=2"), "B_{", minYr, "}$"),  paste0(paste0("$B_{",A$currYr,"}/"),c("B_\\mathrm{avg}$",paste0(c("B_{","2B_{"),minYr,"}$"))), paste0(c("$",paste0("$u_{",A$lastYr,"}/")),"u_\\mathrm{avg}$") )
		}
		colnames(mcmcParamTable) <- paramNames
		table.mcmc.mpd(mcmcData  = mcmcParamTable,
			burnin    = 0,
			probs     = quantProbs,
			mpdData   = mpdParamVector,
			colLabels = paramNames,
			formatOut = "%1.6f",  # Decimal places shown regardless of what roundDec is (RH: actually, both this and roundDec need to be the same)
			roundDec  = 6,        # Number of decimal places
			tableName = "PBMQuants")

  } else {
    cat("No decision tables for age-structured model or bridging analysis\n")
  }

 # try(mb.table.all.param.est(roundDec=2,formatOut="%1.2f"), silent=silent)
 # try(table.i(), silent=silent)
 # try(table.h(mle=F,tableType="ssb"), silent=silent)
 # try(table.h(mle=F,tableType="depletion"), silent=silent)
 # try(table.h(mle=F,tableType="f40spr"), silent=silent)
 
  
  #try(table.b(), silent=silent)
  #try(table.c(), silent=silent)
  #try(table.d(), silent=silent)
  #try(table.e1(), silent=silent)
  #try(table.e2(), silent=silent)
  #try(table.f(), silent=silent)
  #try(table.h(mle=T,tableType="ssb"), silent=silent)
  #try(table.h(mle=T,tableType="depletion"), silent=silent)
  #try(table.h(mle=T,tableType="f40spr"), silent=silent)
  #try(table.1(), silent=silent)
  #try(table.2(), silent=silent)
  #try(table.3(), silent=silent)
  #try(table.4(), silent=silent)  
}

.writeAllPlots <- function(silent=F, gui=FALSE){
  # write all figures for all scenarios to disk
  # scenarioList <- as.numeric(rownames(viewHeader))  ## already defined globally
  if(exists("scenarioList",envir=.GlobalEnv, inherits=FALSE))
    scenarioList = get("scenarioList",envir=.GlobalEnv)
  else stop("PROBLEM: scenarioList not found in Global environment")
  for(scenario in scenarioList){
    assignGlobals(scenario, gui=gui)
#browser();return()
    .writePlots(silent=silent, gui=gui)
  }  
}

.writePlots <- function(silent=FALSE, gui=TRUE){
  ### write all figures for the given scenario to disk
  windows(width=7,height=7,record=TRUE); frame(); dev.box = dev.cur()
  on.exit(dev.off(dev.box))
  .setBurnThin(silent=silent, gui=gui)
  if(exists("scenarioCurr",envir=.GlobalEnv, inherits=FALSE))
    scenarioCurr = get("scenarioCurr",envir=.GlobalEnv)
  else stop("PROBLEM: scenarioCurr not found in Global environment")
  setWinVal(list(entryScenario=scenarioCurr))
  val <- getWinVal()
  assign("saveon",TRUE,envir=.GlobalEnv)

  aside.calls = c(
    "fig.a","fig.effort","fig.e1","fig.e2","fig.g","fig.h","fig.i","fig.j",
    "fig.comm.age.residuals1","fig.comm.age.residuals2","fig.comm.age.props",
    "fig.survey.age.residuals1","fig.survey.age.residuals2","fig.survey.age.props",
    "fig.comm.age.props.fit","fig.survey.age.props.fit","fig.selectivity","fig.phase",
    "fig.catchFit","fig.mcmc.density","fig.mcmc.geweke","fig.estimated.params.pairs2",
    "fig.estimated.params.pairs.key","fig.estimated.params.pairs.no.log.key","fig.variance.partitions"
    )
  simple.calls = c(
    "fig.catchFit", "fig.surveybiomass.fit", "fig.Fmcmc", "fig.weightFit",
    "fig.mcmc.trace", "fig.RecAnomMCMC", "fig.mcmc.autocor", "fig.estimated.params.pairs"
    )
  for (s in simple.calls) 
    if (exists(s)) eval(call(s))

  ### Calls that require argument specifications
  ### ------------------------------------------
  ### Spawning biomass (MPD only)
  if(exists("fig.biomass.mpd"))
    fig.biomass.mpd(ylimit=val$biomassYlim, useMaxYlim=val$maxBiomassYlim, main="Spawning biomass", xlab="Year")
  ### Spawning depletion (MPD only)
  if(exists("fig.depletion.mpd"))
    fig.depletion.mpd(ylimit=val$depletionYlim, useMaxYlim=val$maxDepletionYlim, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr)
  ### Biomass MCMC with MPD overlay
  if(exists("fig.b"))
    fig.b(includeMPD=T, ylimit=val$biomassYlim, useMaxYlim=val$maxBiomassYlim, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr, tac=val$currTAC)
  ### Depletion MCMC with MPD overlay 
  if(exists("fig.c"))
    fig.c(includeMPD=T, ylimit=val$depletionYlim, useMaxYlim=val$maxDepletionYlim, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr)
  ### Recruitment MCMC
  if(exists("fig.dmcmc"))
    fig.dmcmc(ylimit=val$recruitmentYlim, useMaxYlim=val$maxRecruitmentYlim)
  ### Recruitment MPD
  if(exists("fig.d"))
    fig.d(ylimit=val$recruitmentYlim, useMaxYlim=val$maxRecruitmentYlim)
  ### Biomass and recruitment, two-panel plot
  #try(fig.biomass.recruits(yBiomassYlim=val$biomassYlim,
  #     useMaxBiomassYlim=val$maxBiomassYlim,
  #     yRecruitmentYlim=val$recruitmentYlim,
  #     useMaxRecruitmentYlim=val$maxRecruitmentYlim), silent=silent)
  #if(exists("fig.time.varying.selectivity")) {
  #  fig.time.varying.selectivity(2)
  #  fig.time.varying.selectivity(1)
  #}
  if(exists("fig.mcmc.priors.vs.posts"))    fig.mcmc.priors.vs.posts(exFactor=1.0, showEntirePrior=T)
  #if(exists("fig.mcmc.priors.vs.posts2"))   fig.mcmc.priors.vs.posts2(exFactor=1.0)
  #if(exists("fig.mcmc.priors.vs.postskey")) fig.mcmc.priors.vs.postskey(exFactor=1.0)

  if(exists("plotChains") && isScenMcmc()){
    dmcmc = mcmc2(A$mc[,parEstimated()],start=Burn+1,thin=Thin)
    plotChains(dmcmc,pdisc=0,axes=TRUE,between=list(x=.15,y=.2),col.trace=c("green3","red","blue"),xlab="Sample",ylab="Cumulative Frequency")
  }
  if(exists("fig.Allcontrol.pts.Box") && isScenMcmc()){
    fig.Allcontrol.pts.Box(tac.use=val$currTAC, useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr)
  }
#  if(delaydiff==1){
#    #if(nyr>=2012) {
#    dd.calls = c("fig.Allcontrol.pts.Box") #,"fig.Benchmarks"), "fig.MSYcontrol.pts", "fig.Histcontrol.pts"
#    for (d in dd.calls) 
#      if (exists(d)) eval(call(d))
#  }
  if(exists("plotSnail_old") && isScenMcmc()){  ##RH: taken from PBSawatea
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
    Bmsy.mcmc = mcmc.prj[,"BMSY"]
    Umsy.mcmc = mcmc.prj[,"UMSY"]
    BoverBmsy = Bt.mcmc/Bmsy.mcmc
    UoverUmsy = Ut.mcmc/Umsy.mcmc
    plotSnail(BoverBmsy, UoverUmsy, p=c(0.1, 0.9), png=T, path=figDir, PIN=c(7,7))
  }
  if(exists("fig.snail") && isScenMcmc()){  ##RH: taken from PBSawatea
    fig.snail(useHRP=val$useHRP, p=c(0.1, 0.9), png=T, path=figDir, PIN=c(7,7))
  }
  assign("saveon",FALSE,envir=.GlobalEnv)
  return(invisible())
}

.writeSensPlots <- function(silent=F){
  # write overlay sensitivity plots
  windows(width=7,height=7,record=TRUE);frame()
  on.exit(dev.off(dev.cur()))
  assign("saveon",TRUE,envir=.GlobalEnv)
  val <- getWinVal()
  uniqueSensitivityGroups <- c()  # base must be 0
  for(scenario in 1:length(opList)){
    # count number of unique sensitivity groups
    if(!is.element(opList[[scenario]][[4]]$SensitivityGroup,uniqueSensitivityGroups) && opList[[scenario]][[4]]$SensitivityGroup != 0){
        uniqueSensitivityGroups <- c(uniqueSensitivityGroups,opList[[scenario]][[4]]$SensitivityGroup)
    }
  }
  for(sensitivityGroup in uniqueSensitivityGroups){
    try(fig.base.vs.sens(sensitivityGroup=sensitivityGroup,
                         whichPlot="biomass",
                         useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                         ylimit=val$biomassYlim,
                         useMaxYlim=val$maxBiomassYlim),silent=silent)
    try(fig.base.vs.sens(sensitivityGroup=sensitivityGroup,
                         whichPlot="depletion",
                         useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                         ylimit=val$depletionYlim,
                         useMaxYlim=val$maxDepletionYlim),silent=silent)
    try(fig.base.vs.sens(sensitivityGroup=sensitivityGroup,
                         whichPlot="recruits",
                         useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                         ylimit=val$recruitmentYlim,
                         useMaxYlim=val$maxRecruitmentYlim),silent=silent)
    try(fig.base.vs.sens(sensitivityGroup=val$entrySensitivityGroup,
                         whichPlot="refpts",
                         useHRP=val$useHRP, minYr=val$minYr, aveYr=val$aveYr,
                         ylimit=val$RefptSensYlim,
                         useMaxYlim=val$maxRefptSensYlim), silent=silent)
  }
  assign("saveon",FALSE,envir=.GlobalEnv)
}

.runCurrScenario <- function(scenario=val$entryScenario, deleteFiles=F, mleOverride=F, copyADMBExecutables=F, silent=F)
{
  # Steps:
  # 1. Copy the ADMB model executable to the current scenario's folder 
  #    or put it on the path defined on line 51 above.
  # 2. Run the scenario using a system call with either MLE or MCMC
  # 3. If MCMC, use another system call to run mceval
  
  goodToGo <- T
  val <- getWinVal()
  keepMCMC <- F
  if(val$executeType=="sMPD" || mleOverride){
    keepMCMC <- T
  }
  if(deleteFiles){
    deleteOldRun(scenario,keepMCMC=keepMCMC)
  }
  modelEXE <- exModel
	if (copyADMBExecutables)
		copyExecutableToScenarioFolder(scenario=scenario,silent=silent)
	else {
		modPath = sub("\\.;",paste0(".;",modelDir,";"),Sys.getenv()["PATH"])
		Sys.setenv(PATH=modPath)
	}

	rscriptsDir <- getwd() # Save the rscripts directory so we can get back to it
	setwd(opList[[scenario]][[1]])  # change to this scenario's directory
	if(val$delaydiff==1){  
		if(val$executeType=="sMPD" || mleOverride){
	    if(is.na(val$maxfn)){
	      modelCall <- paste(modelEXE, "-delaydiff")
	    }else{
	      modelCall <- paste(modelEXE,"-maxfn",val$maxfn, "-delaydiff")
	    }
	  }else if(val$executeType=="sMCMC"){
	    if(is.na(val$mcmc) || is.na(val$mcsave)){
	      cat("Error - check your mcmc and mcsave boxes for valid values.\n")
	      goodToGo <- F
	    }else{
	      mcevalCall <- paste(modelEXE,"-mceval", "-delaydiff")
	      if(is.na(val$maxfn)){
		modelCall <- paste(modelEXE,"-mcmc",val$mcmc,"-mcsave",val$mcsave, "-delaydiff")
	      }else{
		modelCall <- paste(modelEXE,"-mcmc",val$mcmc,"-mcsave",val$mcsave,"-maxfn",val$maxfn, "-delaydiff")
	      }
	    }
	  }
	} else { #end if delaydiff
	if(val$executeType=="sMPD" || mleOverride){
	    if(is.na(val$maxfn)){
	      modelCall <- modelEXE
	    }else{
	      modelCall <- paste(modelEXE,"-maxfn",val$maxfn)
	    }
	  }else if(val$executeType=="sMCMC"){
	    if(is.na(val$mcmc) || is.na(val$mcsave)){
	      cat("Error - check your mcmc and mcsave boxes for valid values.\n")
	      goodToGo <- F
	    }else{
	      mcevalCall <- paste(modelEXE,"-mceval")
	      if(is.na(val$maxfn)){
		modelCall <- paste(modelEXE,"-mcmc",val$mcmc,"-mcsave",val$mcsave)
	      }else{
		modelCall <- paste(modelEXE,"-mcmc",val$mcmc,"-mcsave",val$mcsave,"-maxfn",val$maxfn)
	      }
	    }
	  }
	} #end else delaydiff
#browser();return()
  if(goodToGo){
    shell(modelCall)
    if(val$executeType=="sMCMC"){                                 
      shell(mcevalCall)
    }
  }
  setwd(rscriptsDir)
  loadScenario(scenario,silent=silent)
  assignGlobals(scenario,silent=silent)
  
}

.loadPlottingLimits <- function(){
  # load the plotting limits into proper entry and checkboxes
  # BIOMASS YLIMITS
  winList <- NULL
  val <- getWinVal()
  scenario <- val$entryScenario
  if(!is.null(opList[[scenario]][[4]]$biomassYlim)){
    winList <- c(winList,biomassYlim=opList[[scenario]][[4]]$biomassYlim)
    winList <- c(winList,biomassSensYlim=opList[[scenario]][[4]]$biomassYlim)
    winList <- c(winList,biomassRetroYlim=opList[[scenario]][[4]]$biomassYlim)
  }
  if(!is.null(opList[[scenario]][[4]]$maxBiomassYlim)){
    winList <- c(winList,maxBiomassYlim=opList[[scenario]][[4]]$maxBiomassYlim)
    winList <- c(winList,maxBiomassSensYlim=opList[[scenario]][[4]]$maxBiomassYlim)
    winList <- c(winList,maxBiomassRetroYlim=opList[[scenario]][[4]]$maxBiomassYlim)
  }
  #TOTAL BIOMASS LIMITS
   if(!is.null(opList[[scenario]][[4]]$tbiomassYlim)){
      winList <- c(winList,tbiomassYlim=opList[[scenario]][[4]]$tbiomassYlim)
      winList <- c(winList,biomassSensYlim=opList[[scenario]][[4]]$tbiomassYlim)
      winList <- c(winList,biomassRetroYlim=opList[[scenario]][[4]]$tbiomassYlim)
    }
    if(!is.null(opList[[scenario]][[4]]$maxBiomassYlim)){
      winList <- c(winList,maxtBiomassYlim=opList[[scenario]][[4]]$maxtBiomassYlim)
      winList <- c(winList,maxtBiomassSensYlim=opList[[scenario]][[4]]$maxtBiomassYlim)
      winList <- c(winList,maxtBiomassRetroYlim=opList[[scenario]][[4]]$maxtBiomassYlim)
  }  
  # DEPLETION YLIMITS
  if(!is.null(opList[[scenario]][[4]]$depletionYlim)){
    winList <- c(winList,depletionYlim=opList[[scenario]][[4]]$depletionYlim)
    winList <- c(winList,depletionSensYlim=opList[[scenario]][[4]]$depletionYlim)
    winList <- c(winList,depletionRetroYlim=opList[[scenario]][[4]]$depletionYlim)
  }
  if(!is.null(opList[[scenario]][[4]]$maxDepletionYlim)){
    winList <- c(winList,maxDepletionYlim=opList[[scenario]][[4]]$maxDepletionYlim)
    winList <- c(winList,maxDepletionSensYlim=opList[[scenario]][[4]]$maxDepletionYlim)
    winList <- c(winList,maxDepletionRetroYlim=opList[[scenario]][[4]]$maxDepletionYlim)
  }
  # RECRUITMENT YLIMITS
  if(!is.null(opList[[scenario]][[4]]$recruitmentYlim)){
    winList <- c(winList,recruitmentYlim=opList[[scenario]][[4]]$recruitmentYlim)
    winList <- c(winList,recruitmentSensYlim=opList[[scenario]][[4]]$recruitmentYlim)
    winList <- c(winList,recruitmentRetroYlim=opList[[scenario]][[4]]$recruitmentYlim)
  }
  if(!is.null(opList[[scenario]][[4]]$maxRecruitmentYlim)){
    winList <- c(winList,maxRecruitmentYlim=opList[[scenario]][[4]]$maxRecruitmentYlim)
    winList <- c(winList,maxRecruitmentSensYlim=opList[[scenario]][[4]]$maxRecruitmentYlim)
    winList <- c(winList,maxRecruitmentRetroYlim=opList[[scenario]][[4]]$maxRecruitmentYlim)
  }
  try(setWinVal(winList), silent=silent)
}

.checkEntries <- function(){
  # Ensures that the entry in the Scenarios box on the GUI is within proper limits
  # Issues an alert box if they are not, and returns FALSE
  # If they are within limits, returns TRUE
  val <- getWinVal()
  scenarioList <- as.numeric(rownames(viewHeader))
  currScenario <- val$entryScenario
#browser();return()
  if(currScenario<min(scenarioList) | currScenario>max(scenarioList)){
    showAlert(paste("Your scenario must be between ",
                    min(scenarioList)," and ",
                    max(scenarioList),".\nNo plot will be drawn.",sep=""),
              title="Scenario Error",icon="warning")
    return(FALSE)
  }
  return(TRUE)
}

# .getWinName  (get the current winName)
# Purpose:     Determine which GUI is active (guiSim, guiView, guiPerf, etc.)
# Parameters:  None
# Returns:     A character containing the name of the current GUI window
# Source:      A.R. Kronlund, modified from PBSref (helper_funs.r)
.getWinName <- function(){
  win <- .PBSmod$.activeWin
  # This is only required if PBSask is used, leave it for now.
  if(win=="PBSask"){
    win <- getWinVal("win", winName="PBSask")[[1]]   # Hidden field in PBSask
    win <- gsub("\n", "", win)                       # Remove the linefeed \n
  }
  return(win)
}

.closeActWin <- function(){
  closeWin(.getWinName())
}
