library(parallel)
library(gtools)

N.mutations <- sum(mt.all[,-1])

simcase.list <- mclapply(1:5000, function(i){
  simcase <- simul.comps(n.proc = 9, n.asymmetric = 5, win.size = 1e+4, 
                         rates = rate.rec, N.mutations = N.mutations, contexts = sites[,-1])
  return(simcase)
}, mc.cores = 30, mc.preschedule = TRUE)

saveRDS(simcase.list,"/d0/home/solrust/mutations/paper_analysis/simcase.rds")


vr.list <- mclapply(simcase.list, function(simcase){
  vol <- simcase$vol
  vrx <- vrnmf::volnmf_main(vol, n.comp = 14, wvol = 1e-3,
                            iter.nmf = 3e+3, n.iter = 1e+3, 
                            #accelerate = TRUE, acc.C = 1, acc.R = 1,
                            verbose = FALSE)
  
  C <- vrx$C
  rownames(C) <- rownames(simcase$spec)
  return(C)
}, mc.cores = 30, mc.preschedule = TRUE)

saveRDS(vr.list,"/d0/home/solrust/mutations/paper_analysis/vrsim.rds")

n.context <- apply(sites[,-1],2,mean)
n.context <- n.context[ unlist(lapply(strsplit(names(rate.rec),"_"),function(x)x[1]))]

simul.stat <- do.call(rbind,lapply(1:length(simcase.list), function(i){
  spec1 <- simcase.list[[i]]$spec/simcase.list[[i]]$vol$col.factors
  affinity <- apply(cor(spec1, vr.list[[i]], method="pearson"),1,max)
  cbind(recovery = affinity, hl = simcase.list[[i]]$hl, load = simcase.list[[i]]$load/14, 
        sparsity = simcase.list[[i]]$sparsity)#, coef = coefs )
}))
simul.stat <- as.data.frame(simul.stat)

ind <- simul.stat$sparsity < 0.5
par(mfrow=c(1,1))
plot(simul.stat$hl[ind],simul.stat$load[ind],col = ifelse(simul.stat$recovery[ind] > 0.8,"red","black"),
     pch=19,cex=0.02,log="xy")


### Dependence of components recovery of spectrum degeneracy
pdf("/d0/home/solrust/mutations/paper_analysis/figures/simulations_spectum_degeneracy.pdf",width = 4.5,height = 4.5)
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1))
smpl <-  1:nrow(simul.stat) #sample(nrow(simul.stat),20000) # sample(nrow(simul.stat),20000)
plot(simul.stat$sparsity[smpl],simul.stat$recovery[smpl],xlab="process spectrum degeneracy",ylab="recovery",
     pch=19,cex=0.5,col=adjustcolor("black",0.02),cex.axis=1.5,cex.lab=1.5,log="x")
abline(h=0.8,lty=1,lwd=2,col="darkgreen")
abline(v=1.0,col="red",lwd=2)
dev.off()


### Draw a heatmap of scale vs loading portrait
ind <- simul.stat$sparsity < 0.5
bins <- 100
xrng <- seq( min(log(simul.stat$hl[ind])),max(log(simul.stat$hl[ind])),length.out=bins)
yrng <- seq( min(  log(simul.stat$load[ind])),max( log(simul.stat$load[ind])),length.out=bins)
#yrng <- seq( min(  log(5e-04)),max( log(5e-02)),length.out=bins)
bin.neibor <- 200 
res.mat <- do.call(rbind,mclapply( xrng,function(i){
  unlist(lapply(yrng,function(j){
    dst <- sqrt( ( log(simul.stat$load[ind])-j)^2/max(abs(log(simul.stat$load[ind])))^2 +
                   ( log(simul.stat$hl[ind])-i)^2/max(abs( log(simul.stat$hl[ind])))^2)
    mean(simul.stat$recovery[ind][rank(dst) < bin.neibor])
  }))
},mc.cores = 30))

# examples of local neighbor
dst <- sqrt( ( log(simul.stat$load[ind])-j)^2/max(abs(log(simul.stat$load[ind])))^2 +
               ( log(simul.stat$hl[ind])-i)^2/max(abs( log(simul.stat$hl[ind])))^2)
inds <- rank(dst) < bin.neibor
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(simul.stat$hl[ind],simul.stat$load[ind],col="grey80",pch=19,cex=0.02,log="xy")
points(simul.stat$hl[ind][inds],simul.stat$load[ind][inds],col="black",pch=19,cex=0.02)

pdf("/d0/home/solrust/mutations/paper_analysis/figures/simulations_heatmap.pdf",width = 5,height = 5)
par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),mgp=c(3,1.0,0))
image( (res.mat),axes=FALSE,names=TRUE,col=adjustcolor(colorRampPalette(c("darkblue","grey90","darkred"))(n = 60),0.75),
       xlab="spatial scale, kb",ylab="spatial loading, % mutations",cex.lab=2)
lblx <- c( 1e+2, 1e+3, 1e+4, 1e+5, 1e+6 )
axis( 1,outer=FALSE, at= (log(lblx)-min(xrng))/(max(xrng)-min(xrng)),labels= lblx/1e+3 ,
      col.axis="black",cex.axis=2)
lbly <- c(2e-4,1e-3,4e-3,1e-2,2e-2)
axis( 2,outer=FALSE, at= (log(lbly)-min(yrng))/(max(yrng)-min(yrng)),labels= lbly*100 ,
      col.axis="black",cex.axis=2)
box(lwd=2)
abline(v=(log(1e+4)-min(xrng))/(max(xrng)-min(xrng)),lty=2,lwd=2)
abline(h=(log(4e-3)-min(yrng))/(max(yrng)-min(yrng)),lty=2,lwd=2)
legend("topright","recovered processes",bty="n",cex=1.5, text.font = 1)
legend("bottomleft","non-recovered\nprocesses",bty="n",cex=1.5,x.intersp = 0, text.font = 1)
dev.off()



###################################### functions #########################################

### wrapper 
simul.comps <-function(n.proc = 9, n.asymmetric = 5, win.size = 1, rates, N.mutations = 4e+8, contexts){
  # set up half-life and spatial variation
  hl.ou <- runif(n.proc, log(1e+2 / win.size), log(1e+6 / win.size)) 
  hl.ou <- exp(hl.ou)
  
  # define stationary variance
  var.ou <-  runif(n.proc, log(1e-4), log(16e-2)) 
  var.ou <- exp(var.ou)
  
  ### simulate intensities
  intensity <- intensity.generate(n.proc, n.asymmetric, win.n = nrow(rates),
                                  mean.ou = rep(1, n.proc), hl.ou = hl.ou, var.ou = var.ou, n.cores=1)
  
  ### generate spectra 
  # set up dirichlet parameters reflecting spectrum sparsity
  alpha.dir <- exp(runif(n.proc, log(1e-2), log(1e+1)))
  # simulate specta
  mname <- colnames(rates)
  mname.compl <- spacemut::revert.context(mname)
  spec <- spectra.generate(n.proc = n.proc, n.asymmetric = n.asymmetric, alpha.dir = alpha.dir, mname, mname.compl)
  
  # re-normalize mutation spectra to obtain overall dataset spectrum rate.mean
  rate.mean <- colMeans(rates)
  im <- apply(intensity,2,mean)
  ct <- unlist(lapply(1:192, function(t) {
    rate.mean[t] / sum(im * spec[t,])
  }))
  spec <- spec * ct
  
  ### simulate mutation rates
  
  # estimate average number of contexts 
  n.context <- apply(contexts,2,mean)
  n.context <- n.context[ unlist(lapply(strsplit(names(rates),"_"),function(x)x[1]))]
  
  mut.est <- mutations.generate(spec, intensity, 
                                N.mutations = N.mutations, n.context = n.context,
                                loads = rep(1, 14)/14 )
  
  require(vrnmf)
  vol <- vrnmf::vol_preprocess(mut.est)
  
  hl <- c(hl.ou,hl.ou[(n.proc-n.asymmetric+1):n.proc])*win.size
  sparsity <- c(alpha.dir,alpha.dir[(n.proc-n.asymmetric+1):n.proc])
  
  #fracs <-  unlist(lapply(1:14,function(i) (t(spec)%*%n.context)[i,]))  
  spatial.load <- sqrt(c(var.ou, var.ou[(n.proc-n.asymmetric+1):n.proc]))*(1/14)
  
  vol$X.process <- NULL
  return(list( vol = vol, spec = spec, hl = hl, sparsity = sparsity, load = spatial.load))
}


#### function to generate spectra
spectra.generate <- function(n.proc, n.asymmetric = 0, alpha.dir = 1, mname, mname.compl){
  if (length(alpha.dir)==1){alpha.dir <- rep(alpha.dir,n.proc)}
  #spec <- 
  spec <- do.call(rbind,lapply(1:n.proc,function(i) rdirichlet(1, rep(alpha.dir[i],192) )))
  colnames(spec) <- mname
  spec[1:(n.proc-n.asymmetric),] <- (spec[1:(n.proc-n.asymmetric),]+spec[1:(n.proc-n.asymmetric),mname.compl])/2
  if (n.asymmetric > 0){
    spec <- rbind(spec,spec[(n.proc-n.asymmetric+1):n.proc,][,mname.compl] )
    colnames(spec) <- mname
  }
  return(t(spec))
}


#### function to generate intensities
intensity.generate <- function(n.proc, n.asymmetric = 0, win.n, 
                               mean.ou = rep(1,n.proc), hl.ou = rep(10, n.proc), var.ou = 0.3, 
                               n.cores = 1 ){
  
  if (n.asymmetric > 0) {
    hl.ou <- c(hl.ou, hl.ou[(n.proc-n.asymmetric+1):n.proc])
    mean.ou <- c(mean.ou, mean.ou[(n.proc-n.asymmetric+1):n.proc])
    if (length(var.ou)==1) {
      var.ou <- rep(var.ou, n.proc+n.asymmetric)
    }else{
      var.ou <- c(var.ou, var.ou[(n.proc-n.asymmetric+1):n.proc])
    }
  }
  
  # define parameters of Ornstein-Uhlenbeck process
  lambdas <-  log(2)/hl.ou
  sigmas <- sqrt(2*var.ou*lambdas)
  # define parameters of autoregressive process resulting from window-averaging of OU
  delta <- 1
  var.auto <- sigmas^2 * ((delta - (1 - exp(-lambdas * delta)) / lambdas) / lambdas^2)
  cov.neigbor <- sigmas^2 * ((1 - exp((-lambdas * delta)))^2 / (2 * lambdas^3))
  
  do.call(cbind,mclapply(1:(n.proc+n.asymmetric),function(k){
    #lambda <- lambdas[k] 
    #sigma <- sigmas[k] 
    x <- mean.ou[k]
    v.ar <- var.auto[k]
    cov.ar <- cov.neigbor[k]
    dt <- 1
    for (i in 2:win.n){
      #Ext <- x[i-1]*exp(-lambda*dt) + mean.ou[k]*(1-exp(-lambda*dt))
      #Vxt <- sigma^2/(2*lambda)*(1-exp(-2*lambda*dt))
      Ext <- mean.ou[k] + (cov.ar / v.ar) * (x[i-1] - mean.ou[k])
      Vxt <- (1 - (cov.ar / v.ar)^2) * v.ar
      x[i] <- rnorm(1, Ext, sqrt(Vxt))
      if (x[i] < 0) {x[i] <- 0}
    }
    x <- x / sum(x)
    x
  },mc.cores=n.cores))
}


#### function to generate observed mutations from processes
#spectra <- spec.proc
#intensities <- intensity.proc
#N.mutations <- 400e+6
mutations.generate <- function(spectra, intensities, N.mutations, n.context, loads){
  loads <- loads/sum(loads)
  intensities.n <- t(t(intensities)*loads)
  m.est <- t(t(intensities.n%*%t(spectra))*n.context)
  m.est <- m.est*N.mutations/sum(m.est)
  yy <- rpois(length(as.numeric(m.est)),lambda=as.numeric(m.est) )
  m.est.sample <- matrix(yy,nrow=nrow(m.est),ncol=ncol(m.est))
  rate.est.sum <- t(t(m.est.sample)/n.context)
  return(rate.est.sum)
}
