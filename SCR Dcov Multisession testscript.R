#simulation scenario set up to vary M, J, and K by session to make sure that functionality works

library(nimble)
library(coda)
source("sim.SCR.Dcov.R") #single session data simulator
source("sim.SCR.Dcov.Multi.R") #multisession wrapper that pastes together session-level data
source("init.SCR.Dcov.Multi.R") #formats data, constants, etc for nimble
source("NimbleModel SCR Dcov Multisession.R") #model file
source("NimbleFunctions SCR Dcov Multisession.R") #custom nimble functions and updates
source("sSampler Dcov Multisession.R") #custom activity center update
source("mask.check.R") #function to check the habitat mask

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

#simulate some data
N.session <- 3
#baseline detection probability
#one per session
p0 <- rep(0.25,N.session)
sigma <- rep(0.50,N.session) #detection spatial scale
K <- c(10,11,9) #number of occasions per session
buff <- rep(3,N.session) #state space buffer. Should be at least 3 sigma (generally).

#make an SCR trapping array. Making the trapping array size vary by session
X <- vector("list",N.session)
X[[1]] <- as.matrix(expand.grid(1:9,1:9))
X[[2]] <- as.matrix(expand.grid(1:8,1:8))
X[[3]] <- as.matrix(expand.grid(1:10,1:10))

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
for(g in 1:N.session){
  xlim[g,] <- range(X[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X[[g]][,1] <- X[[g]][,1]- x.shift
  X[[g]][,2] <- X[[g]][,2]- y.shift
}

res <- rep(0.25,N.session) #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- y.vals <- dSS <- cells <- vector("list",N.session)
n.cells <- n.cells.x <- n.cells.y <- rep(NA,N.session)
for(g in 1:N.session){
  x.vals[[g]] <- seq(xlim[g,1]+res[g]/2,xlim[g,2]-res[g]/2,res[g]) #x cell centroids
  y.vals[[g]] <- seq(ylim[g,1]+res[g]/2,ylim[g,2]-res[g]/2,res[g]) #y cell centroids
  dSS[[g]] <- as.matrix(cbind(expand.grid(x.vals[[g]],y.vals[[g]])))
  cells[[g]] <- matrix(1:nrow(dSS[[g]]),nrow=length(x.vals[[g]]),ncol=length(y.vals[[g]]))
  n.cells[g] <- nrow(dSS[[g]])
  n.cells.x[g] <- length(x.vals[[g]])
  n.cells.y[g] <- length(y.vals[[g]])
}

#create a density covariate - one for each session
library(geoR)
D.cov <- vector("list",N.session)
#need a simulated landscape with individuals living around traps to be captured
#these are pretty good
D.seeds <- c(13216,13216,13218)
for(g in 1:N.session){
  set.seed(D.seeds[g])
  D.cov.tmp <- grf(n.cells[g],grid=dSS[[g]],cov.pars=c(1000,1000),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
  D.cov.tmp <- as.numeric(scale(D.cov.tmp)) #scale
  par(mfrow=c(1,1),ask=FALSE)
  D.cov[[g]] <- D.cov.tmp
  image(x.vals[[g]],y.vals[[g]],matrix(D.cov[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," D.cov"),xlab="X",ylab="Y",col=cols1)
}

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
InSS <- vector("list",N.session)
for(g in 1:N.session){
  dSS.tmp <- dSS[[g]] - res[g]/2 #convert back to grid locs
  InSS[[g]] <- rep(1,length(D.cov[[g]]))
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Habitat"))
}

#Density covariates
D.beta0 <- rep(-1,N.session)
D.beta1 <- rep(0.5,N.session)
#what is implied expected N in state space?
for(g in 1:N.session){
  lambda.cell <- InSS[[g]]*exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  print(sum(lambda.cell)) #expected N in state space
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Expected Density"))
  points(X[[g]],pch=4,cex=0.75) #SCR traps
}

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(13435342) #change seed for new data set
data <- sim.SCR.Dcov.Multi(N.session=N.session,D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                  p0=p0,sigma=sigma,K=K,X=X,xlim=xlim,ylim=ylim,res=res)

#simulated N per session
unlist(lapply(data,function(x){x$N}))
#SCR-detected individuals by session. 
unlist(lapply(data,function(x){x$n}))

#Visualize activity centers
for(g in 1:N.session){
  lambda.cell <- InSS[[g]]*exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Expected Density"))
  points(X[[g]],pch=4,cex=0.75)
  points(data[[g]]$s,pch=16)
}


for(g in 1:N.session){
  #function to test for errors in mask set up. 
  mask.check(dSS=data[[g]]$dSS,cells=data[[g]]$cells,n.cells=data[[g]]$n.cells,n.cells.x=data[[g]]$n.cells.x,
             n.cells.y=data[[g]]$n.cells.y,res=data[[g]]$res,xlim=data[[g]]$xlim,ylim=data[[g]]$ylim,
             x.vals=data[[g]]$x.vals,y.vals=data[[g]]$y.vals)
}


#Data augmentation level
M <- c(175,175,200) #one per session

#initialize multisession data structures
nimbuild <- init.SCR.Dcov.Multi(data,M)

#plot to make sure initialized activity centers are in habitat
for(g in 1:N.session){
  image(data[[g]]$x.vals,data[[g]]$y.vals,matrix(data[[g]]$InSS,data[[g]]$n.cells.x,data[[g]]$n.cells.y))
  points(nimbuild$s.init[g,,],pch=16)
}

#inits for nimble - MUST use z init and N init for data augmentation scheme to work. should use s.init, too.
Niminits <- list(z=nimbuild$z.init,N=rowSums(nimbuild$z.init),s=nimbuild$s.init,
                 p0=rep(0.5,N.session),sigma=rep(1,N.session),
                 D0=rowSums(nimbuild$z.init)/(rowSums(nimbuild$InSS)*nimbuild$cellArea),D.beta1=rep(0,N.session))

#constants for Nimble
constants <- list(N.session=N.session,M=M,J=nimbuild$J,K1D=nimbuild$K1D,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                  D.cov=nimbuild$D.cov,cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells,
                  res=nimbuild$res)

#supply data to nimble
Nimdata <- list(y=nimbuild$y,z=nimbuild$z.data,X=nimbuild$X,
                dummy.data=nimbuild$dummy.data,cells=nimbuild$cells,InSS=nimbuild$InSS)

# set parameters to monitor
parameters <- c('N','lambda.N','p0','sigma','D0',"D.beta1")
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 25 #thinning rate for parameters2

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0','sigma','D0',"D.beta1")
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2,thin2=nt2,nodes=config.nodes,
                      useConjugacy=FALSE)


# ###*required* sampler replacements
z.ups <- round(M*0.25) # how many z proposals per iteration??? 25% of M generally seems good, but no idea what is optimal
for(g in 1:N.session){
  conf$addSampler(target = paste("N[",g,"]"),
                  type = 'zSampler',control = list(g=g,inds.detected=1:nimbuild$n[g],z.ups=z.ups[g],
                                                   J=nimbuild$J[g],K=nimbuild$K[g],M=M[g]),
                  silent = TRUE)
}

for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste0("s[",g,",",i,", 1:2]"),
                    type = 'sSampler',control=list(g=g,i=i,J=nimbuild$J[g],K=nimbuild$K[g],
                                                   res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                   n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,]),silent = TRUE)
  }
}

#often better to block these 2 sets of parameters
for(g in 1:N.session){
  #AF_slice mixes better, but runs more slowly, with reduced runtime a function of how
  #costly it is to evaluate the likelihoods involved.
  #AF_slice is fast for D covs
  conf$addSampler(target = c(paste0("D0[",g,"]"),paste0("D.beta1[",g,"]")),
                  type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
  #AF_slice causes larger slowdown in run time. May still be worth it. But RW_block is faster.
  #Not sure which is most efficient in terms of effective sample size/unit time
  conf$addSampler(target = c(paste0("p0[",g,"]"),paste0("sigma[",g,"]")),
                  type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
  
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
# can ignore warning about pi.cell having NAs
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[burnin:nrow(mvSamples),]))

#truth
unlist(lapply(data,function(x){x$N}))
unlist(lapply(data,function(x){x$lambda}))


#plot to make sure final iteration activity centers are in habitat
for(g in 1:N.session){
  image(data[[g]]$x.vals,data[[g]]$y.vals,matrix(data[[g]]$InSS,data[[g]]$n.cells.x,data[[g]]$n.cells.y))
  points(Cmodel$s[g,,],pch=16)
}


#Look at cell-level expected density estimates, compare to truth
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
n.cells.max <- max(n.cells)
lambda.cell.idx <- matrix(lambda.cell.idx,N.session,n.cells.max)
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10

#compare expected D plot to truth
#image will show posterior means
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- array(NA,dim=c(N.session,n.cells.max,length(n.iter.use)))
lambda.cell <- lambda.cell.ests <- array(NA,dim=c(N.session,n.cells.max))
lambda.cell.HPDs <- array(NA,dim=c(N.session,n.cells.max,2))
for(g in 1:N.session){
  lambda.cell[g,1:n.cells[g]] <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  lambda.cell.post[g,1:n.cells[g],] <- t(cellArea[g]*mvSamples2[n.iter.use,D0.idx[g]]*
    mvSamples2[n.iter.use,lambda.cell.idx[g,1:n.cells[g]]])
  lambda.cell.ests[g,1:n.cells[g]] <- rowMeans(lambda.cell.post[g,1:n.cells[g],])
  lambda.cell.HPDs[g,1:n.cells[g],] <- HPDinterval(mcmc(t(lambda.cell.post[g,1:n.cells[g],])))
  #remove nonhabitat (or not, comment out)
  lambda.cell[g,InSS[[g]]==0] <- NA
  lambda.cell.ests[g,InSS[[g]]==0] <- NA
}

par(mfrow=c(1,1),ask=FALSE)
for(g in 1:N.session){
  zlim <- range(c(lambda.cell[g,],lambda.cell.ests[g,]),na.rm=TRUE) #use same zlim for plots below
  #truth
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"True Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
  #estimate, posterior means
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell.ests[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"Est Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
}
#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
for(g in 1:N.session){
  idx <- order(lambda.cell[g,1:n.cells[g]])
  plot(lambda.cell.ests[g,1:n.cells[g]][idx]~lambda.cell[g,1:n.cells[g]][idx],type="l",lwd=2,
       main=paste("Session",g,"True vs. Estimated Density"),ylim=range(lambda.cell.HPDs[g,1:n.cells[g],]))
  lines(lambda.cell.HPDs[g,1:n.cells[g],1][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  lines(lambda.cell.HPDs[g,1:n.cells[g],2][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  abline(0,1,col="darkred",lwd=2) #1:1 expectation
}