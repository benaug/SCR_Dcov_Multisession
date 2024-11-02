init.SCR.Dcov.Multi <-
  function(data=data,M=M){
    N.session <- length(data)
    if(length(M)!=N.session)stop("M and data must be of length 'N.session'")
    M.max <- max(M)
    J <- unlist(lapply(data,function(x){nrow(x$X)}))
    J.max <- max(J)
    K <- unlist(lapply(data,function(x){x$K}))
    K.max <- max(K)
    n.cells <- unlist(lapply(data,function(x){x$n.cells}))
    n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
    n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
    n.cells.max <- max(n.cells)
    n.cells.x.max <- max(n.cells.x)
    n.cells.y.max <- max(n.cells.y)
    #Structure data for nimble
    y <- array(0,dim=c(N.session,M.max,J.max)) #maximal augmentation across sessions
    X <- array(0,dim=c(N.session,J.max,2))
    xlim <- ylim <- matrix(0,N.session,2)
    n <- rep(NA,N.session)
    K1D <- array(0,dim=c(N.session,J.max))
    res <- unlist(lapply(data,function(x){x$res}))
    cellArea <- res^2
    x.vals <- matrix(NA,N.session,n.cells.x.max)
    y.vals <- matrix(NA,N.session,n.cells.y.max)
    dSS <- array(NA,dim=c(N.session,n.cells.max,2))
    InSS <- array(0,dim=c(N.session,n.cells.max))
    D.cov <- array(NA,dim=c(N.session,n.cells.max))
    cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
    for(g in 1:N.session){
      y[g,1:data[[g]]$n,1:J[g]] <- apply(data[[g]]$y,c(1,2),sum)
      X[g,1:J[g],1:2] <- data[[g]]$X
      xlim[g,] <- data[[g]]$xlim
      ylim[g,] <- data[[g]]$ylim
      n[g] <- data[[g]]$n
      K1D[g,1:J[g]] <- rowSums(data[[g]]$K2D)
      x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
      y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
      dSS[g,1:n.cells[g],] <- data[[g]]$dSS
      InSS[g,1:n.cells[g]] <- data[[g]]$InSS
      D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
      cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
    }
    
    s.init <- array(NA,dim=c(N.session,M.max,2))
    y.2D <- apply(y,c(1,2),sum)
    for(g in 1:N.session){
      s.init[g,1:M[g],] <- cbind(runif(M[g],xlim[g,1],xlim[g,2]), runif(M[g],ylim[g,1],ylim[g,2])) #assign random locations
      idx <- which(y.2D[g,]>0) #switch for those actually caught
      for(i in idx){
        trps <- matrix(X[g,y[g,i,]>0,1:2],ncol=2,byrow=FALSE)
        if(nrow(trps)>1){
          s.init[g,i,] <- c(mean(trps[,1]),mean(trps[,2]))
        }else{
          s.init[g,i,] <- trps
        }
      }
    }
    #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
    e2dist  <-  function (x, y){
      i <- sort(rep(1:nrow(y), nrow(x)))
      dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
      matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
    }
    getCell  <-  function(s,res,cells){
      cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
    }
    for(g in 1:N.session){
      alldists <- e2dist(s.init[g,,],data[[g]]$dSS)
      alldists[,data[[g]]$InSS==0] <- Inf
      for(i in 1:M[g]){
        this.cell <- data[[g]]$cells[trunc(s.init[g,i,1]/data[[g]]$res)+1,trunc(s.init[g,i,2]/data[[g]]$res)+1]
        if(data[[g]]$InSS[this.cell]==0){
          cands <- alldists[i,]
          new.cell <- which(alldists[i,]==min(alldists[i,]))
          s.init[g,i,] <- data[[g]]$dSS[new.cell,]
        }
      }
    }
    
    z.init <- 1*(apply(y,c(1,2),sum)>0)
    dummy.data <- matrix(0,N.session,M.max) #dummy data not used, doesn't really matter what the values are
    
    z.data <- matrix(NA,N.session,M.max)
    for(g in 1:N.session){
      z.data[g,1:n[g]] <- 1
    }
    
    return(list(y=y,X=X,xlim=xlim,ylim=ylim,K=K,J=J,
                n=n,K1D=K1D,res=res,cellArea=cellArea,x.vals=x.vals,
                y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
                n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data,
                s.init=s.init,z.init=z.init,z.data=z.data))
  }
