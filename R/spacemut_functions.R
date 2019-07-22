##' Infer mutational components
##' @param rate matrix of mutation rates of mutation types in genomic windows
##' @param n.comp number of components
##' @return matrices of components spectra and intensities
##' @export
extract.comp <- function (rate, n.comp){
  t <- apply(rate, 2, function(x) (x - mean(x))/sd(x))
  wh <- svd(t)

  X <- (wh$v)[, 1:n.comp]
  ic <- ica::icafast( (X), nc = n.comp)
  icM <- t(t(ic$S) + ( t(solve(t(ic$M)))%*%apply(X,2,mean))[,1])

  rownames(icM) <- colnames(t)
  icS <- (wh$u[, 1:n.comp] %*% diag(wh$d)[1:n.comp, 1:n.comp] %*% ic$M)
  moms.M <- apply(icM, 2, function(x) {
    mean((x - mean(x))^3)/sd(x)^3
  })
  moms.S <- apply(icS, 2, function(x) {
    mean((x - mean(x))^3)/sd(x)^3
  })
  moms <- moms.M + moms.S
  icM <- t(t(icM) * sign(moms))
  icS <- t(t(icS) * sign(moms))
  colnames(icS) <- paste("comp.", 1:ncol(icS), sep = "")
  colnames(icM) <- paste("comp.", 1:ncol(icM), sep = "")
  return(list(icM, icS))
}



##' Plot heatmap of refelction correlation between all pairs of components
##' @param icM matrix of spectra of mutational components
##' @param par boolean, TRUE if external par() is used. TRUE by default
##' @return heatmap
##' @export
plot.reflection.matrix <- function(icM,par=FALSE){
  mat = (cor(icM,icM[revert.context(rownames(icM)),],method="spearman"));
  if (par==FALSE) {par(mar=c(7,8,2,2),mgp=c(5,1,0))}
  image( mat^1,axes=FALSE,names=TRUE,col=adjustcolor( colorRampPalette(c("grey95","grey80","darkred"))(n = 60) ,1),
         xlab="spectrum",ylab="reverse compementary\nspectrum",cex.lab=1.5,main="Reflection matrix",font.main=1,cex.main=2)
  box(lwd=2)
  axis( 1,outer=FALSE, at= seq(0,1,length.out = nrow(mat)),labels = rownames(mat),las=2,
        col.axis="black",cex.axis=1)
  axis( 2,outer=FALSE, at= seq(0,1,length.out = nrow(mat)),labels = rownames(mat),las=1,
        col.axis="black",cex.axis=1)
}

##' Scatterplot of spectrum and reverse complementary spectrum for two selected components
##' @param icM matrix of spectra of mutational components
##' @param i index of first component
##' @param j index of second component, spectrum of which will be reverse complemented
##' @param par boolean, TRUE if external par() is used. TRUE by default
##' @return scatterplot
##' @export
reflection.scatter <- function(i,j,icM,par=FALSE){
  if (par==FALSE) {par(mar=c(6,6,3,3),mgp=c(2.5,1,0))}
  plot(icM[,i],icM[revert.context(rownames(icM)),j],cex=1,pch=19,main="",cex.main=1.7,font.main=2,
       xlab=paste("spectrum,",colnames(icM)[i]),ylab=paste("reverse complementary\nspectrum,",colnames(icM)[j]),cex.lab=1.5,cex.axis=1.5)
  abline(h=0,lty=2);abline(v=0,lty=2)
  rh <- round(cor(icM[,i],icM[revert.context(rownames(icM)),j],method="spearman"),2)
  legend("bottomright", legend=bquote(rho == .(rh)),fill=NA,border=NA,bty="n",cex=1.5)
}


##' Classify components and group in processes based on reflection property
##' @param icM matrix of spectra of mutational components
##' @param cutoff cutoff on reflection property
##' @return matrices of components classification and processes annotation in strand-independent/strand-dependent/noise/ambiguous
##' @export
reflection.test <- function(icM,cutoff=0.8){
  # reflection matrix
  mat = (cor(icM,icM[revert.context(rownames(icM)),],method="spearman"));

  # classify components
  prop <- do.call(rbind,lapply(1:nrow(mat),function(i){
    if ( max(abs(mat[i,])) < cutoff ){
      return( c(NA,"noise") )
    }
    if (mat[i,i] > cutoff) {
      return( c(i,"symmetric") )
    }else if (mat[i,i] < -cutoff){
      return( c(i,"asymmetric") )
    }else{
      j <- which.max(abs(mat[i,]))
      if (i == which.max(abs(mat[j,]))){
        return( c(j,"asymmetric") )
      }else{
        return( c(NA,"ambiguous") )
      }
    }
  }))
  prop <- as.data.frame(prop[,c(2,1)])
  colnames(prop) <- c("type","ref.comp")
  rownames(prop) <- colnames(icM)

  # annotate processes
  proc.real <- do.call(rbind,lapply( 1:nrow(prop),function(i){
    if (!is.na(prop$ref.comp[i])){
      c(ifelse(i==prop$ref.comp[i],"strand-independent","strand-dependent"), sort(c(i,as.character(prop$ref.comp[i]))) )
    }
  }))
  proc.real <- data.frame(unique(proc.real))
  colnames(proc.real) <- c("type","comp. X","comp. Y")
  rownames(proc.real) <- paste("process.",1:nrow(proc.real),sep="")

  return( list(comp = prop, proc = proc.real) )
}

##' Estimate number of components with reflected property and processes for different input extracted components
##' @param  rate matrix of mutation rates of mutation types in genomic windows
##' @param n.min min number of components to extract
##' @param n.max max number of components to extract
##' @param cutoff cutoff on reflection property
##' @param n.cores nunmber of cores to use
##' @return number of components (with reflection property) and processes per each input number of extracted components
##' @export
select.ncomponents <- function(rate,n.min=2,n.max=30,cutoff=0.8,n.cores=1){
  xx <- do.call(rbind,parallel::mclapply(n.min:n.max,function(n){
    comp.info <- extract.comp(rate,n)
    icM <- comp.info[[1]]
    plot.reflection.matrix(icM)
    ref.prop <- reflection.test(icM,cutoff)
    c( n,sum(!is.na(ref.prop[[1]]$ref.comp)),nrow(ref.prop[[2]]) )
  }))
  colnames(xx) <- c("n","comp","proc")
  cat('optimal parameters:');cat('\n')
  print(xx[which.max(xx[,2]),]);cat('\n')
  return(as.data.frame(xx))
}

##' Show number of real mutational components and processes depending on number of extracted components
##' @param stat.extract matrix output of select.ncomponents()
##' @export
show.optimal.comp <- function(stat.extract,par=FALSE){
  if (par==FALSE) {par(mfrow=c(1,1),mar=c(7,8,0.5,0.5),mgp=c(4.5,1,0))}
  plot(stat.extract$n, stat.extract$comp,ylim=c(0,max(stat.extract$comp)+1),
       type=c("l"),col="darkgreen",lwd=3,
       xlab="number of\nextracted components",ylab="number of components\npassed reflection test",cex.lab=1.5,cex.axis=2)
  points(stat.extract$n, stat.extract$comp,
         col="darkgreen",cex=0.75,lwd=3,pch=19)
  lines(stat.extract$n,stat.extract$proc,
        type=c("l"),col="red",lwd=3)
  points(stat.extract$n,stat.extract$proc,
         col="red",cex=0.75,lwd=3,pch=19)
  abline(v=stat.extract[,1][which.max(stat.extract[,2])],lwd=2,lty=2)
  legend("bottomright",c("components","processes"),fill=NA,border=NA,bty="n",col=c("darkgreen","red"),pch=19,cex=1.5)
}

##' Draw component spectrum
##' @param sign vector of mutational component spectrum
##' @param n.sign name of component to show in legend
##' @export
draw.signature <- function(sign,n.sign=NULL){
  sub.alphabet <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  sub.cols <- cbind(sub.alphabet,c("blue","black","red","grey60","darkgreen","pink")); colnames(sub.cols)=c("sub","col")
  sign.context <- parse.context(names(sign)); sign.reverted <- revert.matrix(sign.context)
  sign.r <- sign[sign.context[,3] %in% sub.alphabet]; sign.context <- sign.context[sign.context[,3] %in% sub.alphabet,]
  sign.c <- sign[sign.reverted[,3] %in% sub.alphabet]; sign.reverted <- sign.reverted[sign.reverted[,3] %in% sub.alphabet,]
  tbl <- rbind(sign.r[order(sign.context[,3],sign.context[,1])],sign.c[order(sign.reverted[,3],sign.reverted[,1])])
  #cols <- as.character( (base:::merge)(sign.context,sub.cols,by = c("sub"))[,"col"])[order(sign.context[,3],sign.context[,1])]
  cols <- as.character(sub.cols[match(sign.context[,"sub"],sub.cols[,"sub"]),"col"])[order(sign.context[,3],sign.context[,1])]
  cols1 <- rep("grey",2*length(cols))
  cols1[seq(1,2*length(cols),2)] = adjustcolor(cols,0.3);cols1[seq(2,2*length(cols),2)] = adjustcolor(cols,1)
  n.plots=3
  nf <- layout(matrix(1:n.plots ,n.plots,1, byrow=F), heights=c(0.05,0.05,rep(1,n.plots-2)))#, respect=T)
  #layout.show(nf)
  par(mar=c(0,3,0,0))
  image( t(rbind(1:6,1:6)),col="white",xaxt="na",yaxt="na",border=NA)
  for (k in 0:5) {text(k/5,0.5,sub.alphabet[k+1],cex=2,font=2)}
  par(mar=c(0,3,0,0))
  #image( t(rbind(1:6,1:6)),col=sub.cols[,2],xaxt="na",yaxt="na")
  image( t(rbind(1:6,7:12)),col=c(adjustcolor(sub.cols[,2],0.3),sub.cols[,2]),xaxt="na",yaxt="na")
  box(lwd=0.5);abline(h=0.5,lwd=0.5)
  par(mar=c(3,3,1,0),mgp=c(1.5,0.5,0),lwd=0.15)
  barplot( tbl,border= TRUE,space=rep(0.0,2*ncol(tbl)),
           xaxs = "i", beside=T,ylab="loading",cex.lab = 1.5,cex.axis = 1.5, legend.text=F,names.arg=sign.context[,1],cex.names = 0.6,las=2,col=cols1)
  if (!is.null(n.sign)) legend("topright",paste("Signature",n.sign),cex=1,text.font=2,box.lwd=0)
}

##' Draw component intensities along the genome
##' @param track vector of component intensities
##' @param info data.frame containing chr and start/end information about intensities
##' @param chr chromosome to show
##' @param span.wind number of windows to use for loess smoothing
##' @param start start of the region to show
##' @param end end of the region to show
##' @export
plot.intensities <- function(track, info, chr=1,span.wind=20, start=NULL, end=NULL){
  start.chr <- min(info$start[info$chr==chr]); end.chr <- max(info$end[info$chr==chr])
  if (!is.null(start)){ start.chr <- max(start,start.chr) }
  if (!is.null(end)){ end.chr <- min(end,end.chr) }
  ind <- info$chr == chr & info$start >= start.chr & info$end <= end.chr
  ft <- loess(track[ind] ~ info$start[ind],span=span.wind/sum(ind))
  par(mfrow=c(1,1),mar=c(4,4,1,1),mgp=c(2,0.5,0))
  plot(info$start[ind],track[ind],pch=19,cex=0.3,col="grey50",
       ylab="intensity",xlab="position",cex.lab=2)
  lines(info$start[ind],ft$fitted,col="blue",lwd=1,lty=1)
  abline(h=0,lty=2,col="black")
}


##' Estimate robustness of components spectra using bootstrap
##' @param  icM matrix of components spectra
##' @param rate matrix of mutation rates of mutation types in genomic windows
##' @param n.comp number of extract components in bootstrap. By default ,the same as in icM
##' @param n.boot number of bootstraps
##' @param n.cores number of cores to use
##' @return matrix of Spearman correlations between original and bootstrapped spectra of components
##' @export
spectra.bootstrap <- function(icM,rate,n.comp=ncol(icM),n.boot=100,n.cores=1){
  boot.stat <- do.call(rbind,parallel::mclapply(1:n.boot,function(k){
    rgns <- sample(1:nrow(rate),replace = TRUE)
    comp.info.boot <- extract.comp(rate[rgns,],n.comp)
    icM.boot <- comp.info.boot[[1]]
    icS.boot <- comp.info.boot[[2]]
    mat = (cor(icM[,],icM.boot[,],method="spearman")); #mat[mat < 0] = 0
    apply(mat,1,function(x) max(abs(x)) )
  },mc.cores = n.cores))
  rownames(boot.stat) <- paste("run_",1:n.boot,sep="")
  return(boot.stat)
}

##' Estimate robustness of components spectra using bootstrap
##' @param  boot matrix of components spectra
##' @param  icM matrix of components spectra. If supplied than reflection correlations are also shown.
##' @param par boolean. If TRUE, external plot par() configuration is used. FALSE by default
##' @return boxplot of bootstrap correlations.
##' @export
visualize.bootstrap <- function(boot,icM=NULL,par=FALSE){
  if (par==FALSE){par(mar=c(5,7,0.5,0.5),mgp=c(4,1,0))}
  # draw boxplot of bootstrap
  boxplot(boot[],
          ylim=c(0.6,1),xlab="component",ylab="spectrum similarity,\nSpearman correlation",names=1:13,cex.lab=1.5,cex.axis=1.5,
          range=1e-6,outline=FALSE,boxwex=0.6,staplewex=0,outwex=0,lwd=1,col=adjustcolor("blue",0.1),border="blue" )#,
  # add bootstrap points
  yy <- unlist(lapply(1:ncol(boot),function(i) boot[,i]))
  xx <- unlist(lapply(1:ncol(boot),function(i) rep(i,nrow(boot)) ))
  stripchart(yy~xx, vertical = TRUE,method = "jitter", add = TRUE, pch = 19,cex=0.3, col = adjustcolor('blue',0.1))
  # add points from reflection test
  if (!is.null(icM)){
    xx <- 1:13
    mat = (cor(icM[,],icM[revert.context(rownames(icM)),],method="spearman")); #mat[mat < 0] = 0
    cor.refl <- apply(mat,1,function(x) max(abs(x)) )
    stripchart(cor.refl~ xx, vertical = TRUE,method = "jitter", add = TRUE, pch = 19,cex=1.0, col = 'red')
    legend("bottomright",c("reflection","bootstrap"),pch=19,col=c("red","blue"),bty="n",cex=1.5)
  }
}


##' Return reverse complementary mutation types
##' @param  sign.names vector of mutation types
##' @return vector of reverse complementary mutation types
##' @export
revert.context <- function(sign.names){
  parse.matrix <- parse.context(sign.names)
  revert.entry <- function(entry){
    x1=paste(com(substring(entry[1],3,3)),com(substring(entry[1],2,2)),com(substring(entry[1],1,1)),sep="")
    #c(x1,com(entry[4]),paste( com(substring(entry[1],2,2)),">",com( substring(entry[3],3,3)),sep=""),com(entry[2]))
    paste(x1,com(substring(entry[3],3,3)),sep="_")
  }
  if (is.null(ncol(parse.matrix))){
    revert.entry(parse.matrix)
  }else{
    apply(parse.matrix,1,function(entry) revert.entry(entry))}
}


rearrange.components <- function(comp.info,process.order){
  if ( min(process.order) < 1 | max(process.order) > ncol(comp.info[[1]]) ) {stop("wrong order indices")}
  icM <- comp.info[[1]][,process.order]
  icS <- comp.info[[2]][,process.order]

  colnames(icS) <- paste("comp.", 1:ncol(icS), sep = "")
  colnames(icM) <- paste("comp.", 1:ncol(icM), sep = "")
  return(list(icM, icS))

}

##'  Parsing coded mutation type
##'  E.g. converts 'TCT_A' to a vector c('TCT','T','C>A','T')
##' @param  sign.names vector of mutation types
##' @return matrix of parsed mutation types
parse.context <- function(sign.names){
  info = do.call(rbind,lapply(strsplit(sign.names,"_"),function(pattern){
    c(pattern[1],substring(pattern[1],1,1),paste(substring(pattern[1],2,2),">",pattern[2],sep=""),substring(pattern[1],3,3))
  }))
  colnames(info) = c("context","left","sub","right")
  return(info)
}


##'  Reverse parsed matrix of mutation type
##' @param  parse.matrix matrix of parsed vector of mutation types
##' @return matrix of reversed parsed mutation types
revert.matrix <- function(parse.matrix){
  revert.entry <- function(entry){
    x1=paste(com(substring(entry[1],3,3)),com(substring(entry[1],2,2)),com(substring(entry[1],1,1)),sep="")
    c(x1,com(entry[4]),paste( com(substring(entry[1],2,2)),">",com( substring(entry[3],3,3)),sep=""),com(entry[2]))
  }
  if (is.null(ncol(parse.matrix))){
    revert.entry(parse.matrix)
  }else{
    t(apply(parse.matrix,1,function(entry) revert.entry(entry)))}
}


##' Return complementary nucleotide
##' @param  x nucleotide (A, T, G or C)
##' @return complementary nucleotide
com = function(x){
  if (x=="A") {"T"
  }else if (x=="T"){"A"
  }else if (x=="G"){"C"
  }else if (x=="C"){"G"}
}

