# ccamUtils.r by Chris Grandin

destroyWorkspace <- function(){
  # remove everything from the user's workspace
  rm(list=ls(all=T,envir=.GlobalEnv),envir=.GlobalEnv)
}

#rsprintf <- function(obj,formatOut="%1.2f",roundDec=2){
#  sprintf(formatOut,round(obj,roundDec))
#}

gletter <- function(i,cex=1.0,font=1){
  #Adds letters to plots (i.e., figure numbers)
  usr <- par("usr")
  inset.x <- 0.05*(usr[2]-usr[1])
  inset.y <- 0.05*(usr[4]-usr[3])
  text(usr[1]+inset.x,usr[4]-inset.y,paste("(",letters[i],")",sep=""),cex=cex,font=font)
}

getShade <- function(color,opacity){
  # getShade()
  # returns an rgb string of the specified color and opacity:
  # - color - an R color either text or decimal number
  # - opacity - 2-decimal-digit string (00-99), i.e. "20" means 20%
  # Notes: format of returned string is #RRGGBBAA
  #        where RR=red, a 2-hexadecimal-digit string
  #        GG=green, a 2-hexadecimal-digit string
  #        BB=blue, a 2-hexadecimal-digit string
  #        AA=alpha or opacity
  
  colorDEC <- col2rgb(color)
  colorHEX <- sprintf("%X", colorDEC)
  for(i in 1:length(colorHEX)){
    if(nchar(colorHEX[i])==1){
      colorHEX[i] <- paste("0",colorHEX[i],sep="")
    }
  }
  shade <- paste("#",colorHEX[1],colorHEX[2],colorHEX[3],opacity,sep="")
  return(shade)
}

rsprintf <- function(obj,
                     formatOut = "%1.6f",
                     roundDec  = 6,
                     formatThousands = FALSE){
  # foramts the obj to the sprintf format given, rounded to roundDec
  # If formatThousands is true, the result will be a decimal integer seperated by a comma
  # in the thousands range.
  if(formatThousands){
    formatC(obj, big.mark=",", format="d")
  }else{
    sprintf(formatOut,round(obj,roundDec))
  }
}

### RH: Original function package `coda' doesn't work as desired
mcmc2 = function (data = NA, start = 1, end = numeric(0), thin = 1) 
{
	if (is.matrix(data) | is.data.frame(data)) {
		niter <- nrow(data)
		nvar <- ncol(data)
	} else {
		niter <- length(data)
		nvar <- 1
	}
	thin <- round(thin)
	if (length(start) > 1) 
		stop("Invalid start")
	if (length(end) > 1) 
		stop("Invalid end")
	if (length(thin) != 1) 
		stop("Invalid thin")
	if (missing(end)) 
		end <- start + (niter - 1) * thin
	else if (missing(start)) 
		start <- end - (niter - 1) * thin
	nobs <- floor((end - start)/thin + 1)
	if (niter < nobs) 
		stop("Start, end and thin incompatible with data")
	else {
		end <- start + thin * (nobs - 1)
		index = if (nobs<niter) 1:nobs else seq(start,nobs,thin)
		if (is.matrix(data) | is.data.frame(data))
			data <- data[index, , drop = FALSE]
		else
			data <- data[index]
	}
	#attr(data, "mcpar") <- c(start, end, thin)
	#attr(data, "class") <- "mcmc"
	data
}

## Taken from PBStools
.flush.cat = function (...)
{
    cat(...)
    flush.console()
    invisible()
}

renderVals = function(..., simplify=FALSE)
{
	dots=list(...)
	if (length(dots)==0) return(NULL)
	unpackList(dots)
	out = list()
	if (exists("minYr", inherits=FALSE)){
		if (is.numeric(minYr)) minYr = deparse(minYr)
		eval(parse(text=paste("minYr=c(",minYr,")")))
		if (!all(minYr%in%yr))
			stop(paste0("Select year(s) in ",yr[1]," to ", rev(yr)[1]," for minimum biomass"))
		out[["minYr"]] = minYr
	}
	if (exists("aveYr", inherits=FALSE)){
		if (is.numeric(aveYr)) aveYr = deparse(aveYr)
		eval(parse(text=paste("aveYr=c(",aveYr,")")))
		if (!all(aveYr%in%yr))
			stop(paste0("Select year(s) in ",yr[1]," to ", rev(yr)[1]," for average biomass"))
		out[["aveYr"]] = aveYr
	}
	if (simplify && length(out)==1) out = out[[1]]
	return(out)
}
