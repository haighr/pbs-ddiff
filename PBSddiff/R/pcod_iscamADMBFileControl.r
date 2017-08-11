
#**********************************************************************************
# admbFileControl.r
# This file contains:
# 1. Functions for creating DAT and CTL files for the stock of your choice.  
# 2. Functions for copying these files to the execution directory
# 3. Functions for copying the outputs from the admb runs back to the scenario directories
# 4. Function to copy model executable to each Scenarios' folder
#
# The data structure used is an opList, which is a list of lists, see loadScenarios.r for
# details on this structure. This file assumes that an object called 'opList'
# exists and is a valid opList.
#
# Assumes that the tcl/tk window is open and ready to be read from/written to
# using PBSmodelling commands such as setWinVal() and getWinVal().  Also assumes
# that the variable fdScenario is in the .GlobalEnv workspace and that it points to the
# current Scenario.
# Assumes loadScenarios.r has been sourced. It instantiates the path and file name variables
#
# Author            : Chris Grandin
# Development Date  : December 2011 - January 2012
#
#**********************************************************************************

copyExecutableToScenarioFolder <- function(scenario,silent=F){
# Copies the model executable from the execution directory to
# the given scenario folder
  src <- paste(fdADMB,exModel,".exe",sep="")
  des <- paste(opList[[scenario]][[1]],exModel,".exe",sep="")
  file.copy(src,des,overwrite=T)
  if(!silent){
    cat("Copied ",src," to ",des,"\n")
  }  
}

copyExecutableToScenariosFolders <- function(silent=F){
# Copies the model executable from the execution directory to
# all scenarios folders.
  for(scenario in 1:length(dirList)){
    copyExecutableToScenarioFolder(scenario,silent)
  }
}

deleteOldRun <- function(scenario,keepMCMC=F,silent=F){
  # select all non-mcmc files
  if(keepMCMC){
    opFilesCurrDirTmp <- opFilesCurrDir[!opFilesCurrDirIsMCMC]
  }else{
    opFilesCurrDirTmp <- opFilesCurrDir
  }
  # Need to add retrospective RET* files to deletion list if they exist..
  dirScenario <- dir(opList[[scenario]][[1]])
  retLoc <- grep(paste(fnRetro,"[0-9]+",sep=""),dirScenario,T)
  if(length(retLoc)>0){
    opFilesCurrDirTmp <- c(opFilesCurrDirTmp,dirScenario[retLoc])
  }
  # End of adding retrospective files to the list

  # delete opFilesCurrDir
  des <- paste(viewHeader[scenario,],opFilesCurrDirTmp,sep="")
  unlink(des)
  # delete opFilesADMBDir
  des <- paste(viewHeader[scenario,],opFilesADMBDir,sep="")
  unlink(des)
}

