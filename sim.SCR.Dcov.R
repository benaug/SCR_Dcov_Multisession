e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.Dcov <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,lam0=NA,
           p0=NA,sigma=NA,K=NA,X=NA,xlim=NA,ylim=NA,res=NA,seed=NA){
    if(!is.na(seed)){
      set.seed(seed)
    }
    #get expected N
    cellArea <- res^2
    lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
    lambda.N <- sum(lambda.cell)
    #simulate realized N
    N <- rpois(1,lambda.N)
    
    #recreate some Dcov things so we can pass fewer arguments into this function
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    # simulate a population of activity centers
    pi.cell <- lambda.cell/sum(lambda.cell)
    s.cell <- sample(1:n.cells,N,prob=pi.cell,replace=TRUE)
    #distribute activity centers uniformly inside cells
    s <- matrix(NA,nrow=N,ncol=2)
    for(i in 1:N){
      tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
    }
    D <- e2dist(s,X)
    J <- nrow(X)
    
    # Capture individuals - SCR
    y <- array(0,dim=c(N,J,K))
    y.state <- array(0,dim=c(N,J,K))
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          y[i,j,k] <- rbinom(1,1,pd[i,j])
        }
      }
    }
    
    #Get observed data
    caught <- which(apply(y,c(1),sum)>0)
    n <- length(caught)
    # y.SCR.true <- y
    y <- y[caught,,]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
    }
    
    #create trap operation. assumed to be complete in data simulator
    K2D <- matrix(1,J,K)
    
    out <- list(y=y,X=X,K=K,K2D=K2D,s=s,n=n,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N,seed=seed)
    return(out)
  }
