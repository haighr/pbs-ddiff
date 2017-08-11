### Return string w/o leading or trailing whitespace
### http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
	trim <- function (x) gsub("^\\s+|\\s+$", "", x)

### Use for rep file only
reptoRlist <- function(fn){
  #ifile <- scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  A <- list()
  ifile = readLines(fn,warn=FALSE)
  is.blank = sapply(ifile,nchar)==0
  ifile = ifile[!is.blank]
  is.comment = grepl("#",ifile)
  ifile = ifile[!is.comment]
  ifile = trim(ifile)  ## get rid of leading and trailing whitespace
  if (any(grepl("ControlFile",ifile))) {
    ii = grep("ControlFile",ifile)
    A[["ControlFile"]] =  ifile[ii+1]
    ifile = ifile[setdiff(1:length(ifile),ii+c(0,1))]
  }
  ilist = strsplit(ifile,split=" ")
  warn = options()$warn; options(warn=-1)
  ddata = sapply(ilist,as.double)
  #idx <- sapply(as.double(ifile),is.na)
  options(warn=warn)
  is.data = sapply(ddata,function(x){all(!is.na(x))})
  is.name = !is.data
  ddata[is.name] = ifile[is.name]
  vnam = ifile[is.name]
  vnam = gsub(",","_",vnam)  ## RH get rid of commas in names
  nv <- length(vnam)         ## number of objects
  idx = (1:length(ddata))[is.name]
  names(idx)=vnam
  idd = (1:length(ddata))[is.data]
  nrdat = diff(c(idx,rev(idd)[1]+1))-1
  names(nrdat) = vnam
#browser();return()
  for (i in vnam) {
    ii = nrdat[i]
    if (ii==0)  A[[i]] = NA
    else {
      idat =  ddata[idx[i]+(1:ii)]
      if (ii == 1)
        A[[i]] = idat[[1]]
      else {
        icols = sapply(idat,length)
        nc = max(icols)  ## number of columns
        ijmat = matrix(NA,nrow=ii,ncol=nc)
        for (iii in 1:ii)
          ijmat[iii,1:icols[iii]] = idat[[iii]]
        A[[i]] = ijmat
      }
    }
  }
#browser();return()
  return(A)
}
### Use for dat file only
dattoRlist <- function(fn){
  #ifile <- scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  A <- list()
#browser();return()
  ifile = readLines(fn,warn=FALSE)
  ifile = rev(rev(ifile)[4:length(ifile)])  ## get rid of the last three lines
  is.blank = sapply(ifile,nchar)==0
  ifile = ifile[!is.blank]
  is.comment = grepl("##",ifile)
  ifile = ifile[!is.comment]
  ifile = trim(ifile)  ## get rid of leading and trailing whitespace
  ifile = gsub("\\t"," ",ifile)

  dimdat = sapply(strsplit(ifile[1:5],split=" #"),trim,simplify=F)
  for (i in 1:length(dimdat))
    A[paste0("#",dimdat[[i]][2])] = as.numeric(dimdat[[i]][1])
  nocomment = sapply(strsplit(ifile[6:7],split=" "),as.numeric,simplify=F)
  A[["Selec_allo"]] = nocomment[[1]]
  A[["Catch_type"]] = nocomment[[2]]
  ifile = ifile[8:length(ifile)]

  if (any(grepl(" #",ifile))) {
    hangcomm = grep(" #",ifile)
    hangdat  = sapply(strsplit(ifile[hangcomm],split=" #"),trim,simplify=F)
    for (i in 1:length(hangcomm)){
      ii = hangcomm[i]
      ifile[ii] = hangdat[[i]][1]
    }
  }
  ifile = ifile[setdiff(1:length(ifile),grep("^#$",ifile))]  ## get rid of any standalone comment markers

  ilist = strsplit(ifile,split=" ")
  warn = options()$warn; options(warn=-1)
  ddata = sapply(ilist,as.double)
  #idx <- sapply(as.double(ifile),is.na)
  options(warn=warn)
  is.data = sapply(ddata,function(x){all(!is.na(x))})
  is.name = !is.data
  ddata[is.name] = ifile[is.name]
  vnam = ifile[is.name]
  vnam = gsub(",","_",vnam)  ## RH get rid of commas in names
  nv <- length(vnam)         ## number of objects
  idx = (1:length(ddata))[is.name]
  names(idx)=vnam
  idd = (1:length(ddata))[is.data]
  nrdat = diff(c(idx,rev(idd)[1]+1))-1
  names(nrdat) = vnam

  for (i in vnam) {
    ii = nrdat[i]
    if (ii==0)  A[[i]] = NA
    else {
      idat =  ddata[idx[i]+(1:ii)]
      if (ii == 1)
        A[[i]] = idat[[1]]
      else {
        icols = sapply(idat,length)
        nc = max(icols)  ## number of columns
        ijmat = matrix(NA,nrow=ii,ncol=nc)
        for (iii in 1:ii)
          ijmat[iii,1:icols[iii]] = idat[[iii]]
#browser();return()
        A[[i]] = ijmat
      }
    }
  }
#browser();return()
  return(A)
}
poo = function(fn) {  ### Stuff below wasn't working
  vnam <- ifile[idx] #list names
  ir <- 0
  for(i in 1:nv){
    ir <- match(vnam[i],ifile)
    if(i!=nv){
      irr <- match(vnam[i+1],ifile)
    }else{
      irr <- length(ifile)+1 #next row
    }
    dum <- NA
    if(irr-ir==2){
      dum <- as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    }
    if(irr-ir>2){
      dum <- as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    }
#browser()
    if(is.numeric(dum)) #Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]] <- dum
      }
  }
  return(A)
}

read.fit <- function(file){ 
  # Function to read a basic AD Model Builder fit. 
  # Use for instance by: 
  #   simple.fit <- read.fit('c:/admb/examples/simple') 
  # 
  # Then the object 'simple.fit' is a list containing sub
  # 'names', 'est', 'std', 'cor', and 'cov' for all model
  # parameters and sdreport quantities.
  #
  ret <- list()
  parfile <- as.numeric(scan(paste0(file,".par"),what="", n=16, quiet=TRUE)[c(6,11,16)])
  ret$nopar <- as.integer(parfile[1]) 
  ret$nlogl <- parfile[2] 
  ret$maxgrad <- parfile[3] 
  file <- paste(file,".cor", sep="") 
  lin <- readLines(file) 
  ret$npar <- length(lin)-2 
  ret$logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2]) 
  sublin <- lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
  ret$names <- unlist(lapply(sublin,function(x)x[2])) 
  ret$est <- as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
  ret$std <- as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
  ret$cor <- matrix(NA, ret$npar, ret$npar) 
  corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
  ret$cor[upper.tri(ret$cor, diag=TRUE)] <- as.numeric(corvec) 
  ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
  ret$cov <- ret$cor*(ret$std%o%ret$std) 
  return(ret) 
} 

bubble.plot <- function (x=1:dim(z)[1], y=1:dim(z)[2], z, scale = 1, log.scale = FALSE,fill=T,add=F, ...){
  zo <- z
  if (log.scale){
    zo <- log(abs(z) + 1)
  }
  n <- dim(z)[1]
  ny <- dim(z)[2]
  xo <- outer(x, rep(1, length = length(y)))
  yo <- t(outer(y, rep(1, length = length(x))))
  zo <- zo/(max(abs(zo),na.rm=T) *scale) #apply(zo, 2, "/", max(abs(zo))) * length(y) * scale
  zo[abs(zo)<=0.001] <- NA
  nt <- rowSums(z)
  #zo[zo==NA]=0
  if(!add){
    matplot(xo, yo, type = "n", ...)
  }
  for (i in 1:dim(z)[1]) {
    iclr <- rep("transparent", length = ny)
    iclr[z[i, ] <= 0] <- "salmon"
    if(fill){
      points(xo[i, 1:ny], yo[i, 1:ny], cex=abs(zo[i, ]), pch=16, col=iclr)
    }
    points(xo[i, 1:ny], yo[i, 1:ny], cex = abs(zo[i, ]), pch=1, col="black")
  }   
}

### Use for pfc file only
pfctoRlist <- function(fn){
  #ifile <- scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  A <- list()
  ifile = readLines(fn,warn=FALSE)
  is.blank = sapply(ifile,nchar)==0
  ifile = ifile[!is.blank]
  is.comment = grepl("##",ifile)
  ifile = ifile[!is.comment]
  ifile = trim(ifile)  ## get rid of leading and trailing whitespace
  if (any(grepl("loop",ifile))) {
    ii = grep("loop",ifile)
    A[["ncatpol"]] =  as.numeric(ifile[ii+1])
    ifile = ifile[setdiff(1:length(ifile),ii+c(0,1))]
  }
  if (any(grepl("tac",ifile))) {
    ii = grep("tac",ifile)
    A[["catpol"]] =  as.numeric(ifile[ii+(1:A$ncatpol)])
    ifile = ifile[setdiff(1:length(ifile),ii+c(0,1:A$ncatpol))]
  }
  return(A)
}

#out=pfctoRlist("WAPSassess04dd.pfc")