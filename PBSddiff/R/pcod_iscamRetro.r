#**********************************************************************************
# ccamRetro.r
# This file contains the code necessary to plot Retrospectives.
# 
# Assumes the opList has been built and has retrospective data in position 5,
# i.e. opList[[scenario]][[5]] contains a list of the retrospective outputs
#
# Author            : Chris Grandin
# Development Date  : December 2011 - January 2012
#
#**********************************************************************************

.runAllRetros <- function(silent=F, gui=FALSE){
  # Run retrospectives for all scenarios. Use value in entry box for years.
  val <- getWinVal()
  for(scenario in 1:length(opList)){
    .runRetros(scenario, gui=gui)
  }
  
}

.runRetros <- function(scenario=val$entryScenario, silent=FALSE, gui=TRUE){
  # Steps:
  # 1. Save the current REP file by copying using file.copy
  # 2. Run the scenario using a system call for MLE retro, 1 for each retro year.
  # 3. Rename each of these runs' REP files to RET* where * is a number.
  # 4. Use file.copy to restore the original REP file
  
  val <- getWinVal()
  retroYears <- val$entryRetro
  showOutput <- val$showRetroOutput
  
  modelEXE <- exModel
  # Save the rscripts directory so we can get back to it
  rscriptsDir <- getwd()
  # change to this scenario's directory
  setwd(opList[[scenario]][[1]])
  # Save the REP file from the main non-retro run by copying to a backup file
  file.copy("pcod_iscam.rep","pcod_iscam.backup.rep")

  for(retro in 1:retroYears){
    modelCall <- paste(modelEXE,"-retro",retro)
    system(modelCall,wait=T,show.output.on.console=showOutput) 
    file.copy("pcod_iscam.rep",paste("pcod_iscam.ret",retro,sep=""),overwrite=T)
  }
  
  # Reinstantiate the REP file from the main non-retro run
  file.copy("pcod_iscam.backup.rep","pcod_iscam.rep",overwrite=T)
  
  setwd(rscriptsDir)
  loadScenario(scenario, silent=silent)
  assignGlobals(scenario, gui=gui)
}

.writeRetroPlots <- function(silent=F){
  windows(width=7,height=7,record=TRUE);frame()
  on.exit(dev.off(dev.cur()))
  assign("saveon",TRUE,envir=.GlobalEnv)
  val <- getWinVal()
  try(fig.retro(whichPlot="biomass",
            ylimit=val$biomassYlim,
            useMaxYlim=val$maxBiomassYlim), silent=silent)
  try(fig.retro(whichPlot="depletion",
            ylimit=val$depletionYlim,
            useMaxYlim=val$maxDepletionYlim), silent=silent)
  try(fig.retro(whichPlot="recruits",
                ylimit=val$recruitmentYlim,
                useMaxYlim=val$maxRecruitmentYlim), silent=silent)
  assign("saveon",FALSE,envir=.GlobalEnv)
}

fig.retro <- function(whichPlot="biomass",ylimit=6,useMaxYlim=T,lty=1,lwd=2,pch=20){
  # plots Spawning stock biomiass, depletion, and recruitment retrospectives
  # whichPlot can be:
  # 1. "biomass"
  # 2. "depletion"
  # 3. "recruits"
  # Assumes that opList[[1]] is populated (i.e. there is at least one scenario)
	if (length(opList[[1]][[5]])==0) {
		cat("WARNING (fig.retro): No Retros generated for this scenario\n"); return(invisible("No Retro data")) }

  op	<- par(no.readonly=T)
  val <- getWinVal()
  scenario <- val$entryScenario

  baseRunName <- strsplit(opList[[scenario]][[1]],"/")[[1]][3] # Gets the run's folder name out of the reletive path
  runNames <- baseRunName
  baseRep <- opList[[scenario]][[4]]
  base <- 1
  color <- 1
  colors <- color
  # Base run is loaded, now plot it and add the retrospectives to the plot
  
  # *************************************** RECRUITMENT RETROSPECTIVE ******************************************  
  if(whichPlot == "recruits"){
    if(useMaxYlim){
      # get max of all retros first
      yUpperLimit <- max(baseRep$rt)
      for(retro in 1:length(opList[[scenario]][[5]])){
        yUpperLimit <- max(yUpperLimit,opList[[scenario]][[5]][[retro]]$rt)
      }
    }else{
      yUpperLimit <- ylimit
    }

    plot(baseRep$yr[1:(length(baseRep$yr)-1)],
         baseRep$rt,
         type="b",
         ylab="Age-1 recruits",
         ylim=c(0,yUpperLimit),
         xlab="Year",
         pch=pch,
         col=color,
         lwd=lwd,
         lty=lty)
    
    for(retro in 1:length(opList[[scenario]][[5]])){ # number of retrospectives
      color <- color + 1
      colors <- c(colors,color)
      lines(opList[[scenario]][[5]][[retro]]$yr[1:(length(opList[[scenario]][[5]][[retro]]$yr)-1)],
            opList[[scenario]][[5]][[retro]]$rt,
            lty=lty,
            col=color,
            lwd=lwd,
            ylim=c(0,yUpperLimit),
            type="b",
            xaxt="n")
      runNames <- c(runNames,paste(baseRunName," - ",retro))
    }
    legend("topright",runNames,lty=lty,col=colors,bty="n",lwd=lwd) 
    filename <- paste("fig.retro.scenario",scenario,".recruits",sep="")
    saveFig(filename)

  # *************************************** BIOMASS RETROSPECTIVE ******************************************  
  }else if(whichPlot == "biomass"){
    if(useMaxYlim){
      yUpperLimit <- max(baseRep$sbt)
      for(retro in 1:length(opList[[scenario]][[5]])){
        yUpperLimit <- max(yUpperLimit,opList[[scenario]][[5]][[retro]]$sbt)
      }
    }else{
      yUpperLimit <- ylimit
    }
    plot(baseRep$yrs,
         baseRep$sbt,
         type="l",
         col=color,
         lty=lty,
         lwd=lwd,
         ylim=c(0,yUpperLimit),
         xlab="Year",
         ylab="Spawning biomass",
         main="Spawning biomass",
         las=1)
    
    points(baseRep$yrs[1]-0.8,
           baseRep$sbo,
           col=color,
           pch=pch,
           lwd=lwd)
    for(retro in 1:length(opList[[scenario]][[5]])){ # number of retrospectives
      color <- color + 1
      colors <- c(colors,color)
      lines(opList[[scenario]][[5]][[retro]]$yrs,
            opList[[scenario]][[5]][[retro]]$sbt,
            type="l",
            col=color,
            lty=lty,
            lwd=lwd,
            ylim=c(0,yUpperLimit),
            xlab="",
            ylab="",
            las=1,
            xaxt="n")
      points(opList[[scenario]][[5]][[retro]]$yrs[1]-0.8,
             opList[[scenario]][[5]][[retro]]$sbo,
             col=color,
             pch=pch,
             lwd=lwd)
      runNames <- c(runNames,paste(baseRunName," - ",retro))
    }
    legend("topright",runNames,lty=lty,col=colors,bty="n", lwd=lwd) 
    filename <- paste("fig.retro.scenario",scenario,".biomass",sep="")
    saveFig(filename)
    
  # *************************************** DEPLETION RETROSPECTIVE ******************************************  
  }else if(whichPlot == "depletion"){
    if(useMaxYlim){
      yUpperLimit <- max(baseRep$sbt/baseRep$sbo)
      for(retro in 1:length(opList[[scenario]][[5]])){
        yUpperLimit <- max(yUpperLimit,opList[[scenario]][[5]][[retro]]$sbt/opList[[scenario]][[5]][[retro]]$sbo)
      }
    }else{
      yUpperLimit <- ylimit
    }
    plot(baseRep$yrs,
         baseRep$sbt/baseRep$sbo,
         type="l",
         col=color,
         lty=lty,
         lwd=lwd,
         ylim=c(0,yUpperLimit),
         xlab="Year",
         ylab="Spawning Depletion",
         main="Spawning Depletion",
         las=1)
    abline(h=0.40, lwd=mtLineWidth, lty=mtLineType, col=mtLineColor)     
    for(retro in 1:length(opList[[scenario]][[5]])){ # number of retrospectives
      color <- color + 1
      colors <- c(colors,color)
      lines(opList[[scenario]][[5]][[retro]]$yrs,
            opList[[scenario]][[5]][[retro]]$sbt/opList[[scenario]][[5]][[retro]]$sbo,
            type="l",
            col=color,
            lty=lty,
            lwd=lwd,
            ylim=c(0,yUpperLimit),
            xlab="",
            ylab="",
            las=1,
            xaxt="n")          
      runNames <- c(runNames,paste(baseRunName," - ",retro))
    }
    lty <- c(rep(1,(length(runNames))),mtLineType)
    lwd <- c(rep(2,(length(runNames))),mtLineWidth)
    colors <- c(colors,mtLineColor)
    legend("topright",c(runNames,"Management target"),lty=lty,col=colors,bty="n", lwd=lwd) 
    filename <- paste("fig.retro.scenario",scenario,".depletion",sep="")
    saveFig(filename)
  }
	par(op)
}
