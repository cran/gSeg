### main functions
# single change-point
gseg1 = function(n, E, statistics=c("all","o","w","g","m"), n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100){
  r1 = list()
  n0 = ceiling(n0)
  n1 = floor(n1)
  Ebynode = vector("list", n)
  
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  
  n0_us = n0
  n1_us = n1
  
  if(n0<2){
  	n0=2
  }
  if(n1>(n-2)){
  	n1=n-2
  }

  r1$scanZ = gcp1bynode(n,Ebynode,statistics,n0,n1)
    
  if (pval.appr==TRUE){
    mypval1 = pval1(n,E,Ebynode,r1$scanZ,statistics, skew.corr,n0,n1)
    r1$pval.appr = mypval1
  }
  if (pval.perm==TRUE){
    mypval2 = permpval1(n,Ebynode,r1$scanZ,statistics,B,n0,n1)
    r1$pval.perm = mypval2
  }

  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
  	if(n0_us<=1){
      cat("  Note: Starting index has been set to n0 = 2 as the original edge-count test statistic is not well-defined for t<2. \n")
    }
    if(n1_us>=n-1){
      cat("  Note: Ending index has been set to n1 =", n-2, " as the original edge-count test statistic is not well-defined for t>",n-2,". \n")
    }
    cat("Original edge-count scan statistic: \n")
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    cat("Weighted edge-count statistic: \n")
    if(n0_us<=1){
      cat("  Note: Starting index has been set to n0 = 2 as the weighted edge-count test statistic is not well-defined for t<2. \n")
    }
    if(n1_us>=n-1){
      cat("  Note: Ending index has been set to n1 =", n-2, " as the weighted edge-count test statistic is not well-defined for t>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    cat("Generalized edge-count statistic: \n")
    if(n0_us<=1){
      cat("  Note: Starting index has been set to n0 = 2 as the generalized edge-count test statistic is not well-defined for t<2. \n")
    }
    if(n1_us>=n-1){
      cat("  Note: Ending index has been set to n1 =", n-2, " as the generalized edge-count test statistic is not well-defined for t>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    cat("Max-type edge-count statistic: \n")
    if(n0_us<=1){
      cat("  Note: Starting index has been set to n0 = 2 as the max-type edge-count test statistic is not well-defined for t<2. \n")
    }
    if(n1_us>=n-1){
      cat("  Note: Ending index has been set to n1 =", n-2, " as the max-type edge-count test statistic is not well-defined for t>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type$pval, "\n")
    }
  }

  return(r1)
}

# changed interval
gseg2 = function(n, E, statistics=c("all", "o", "w", "g", "m"), l0=0.05*n, l1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100){
  l0 = ceiling(l0)
  l1 = floor(l1)
  Ebynode = vector("list", n)
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }

	l0_us = l0
	l1_us = l1
	
  if(l0<=1){
     l0 = 2
   }
  if(l1>=(n-1)){
     l1=n-2
   }
  r1 = list()
  r1$scanZ = gcp2bynode(n,Ebynode,statistics,l0,l1)

  if (pval.appr==TRUE){
    mypval1 = pval2(n,E,Ebynode,r1$scanZ,statistics, skew.corr,l0,l1)
    r1$pval.appr = mypval1
  }
  if (pval.perm==TRUE){
    mypval2 = permpval2(n,Ebynode,r1$scanZ,statistics,B,l0,l1)
    r1$pval.perm = mypval2
  }

  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    cat("Original edge-count scan statistic: \n")
    if(l0_us<=1){
      cat("  Note: Minimum interval length has been set to l0 = 2 as the original edge-count test statistic is not well-defined for interval length<2. \n")
    }
    if(l1_us>=n-1){
      cat("  Note: Maximum interval length has been set to l1 =", n-2, " as the original edge-count test statistic is not well-defined for interval length>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    cat("Weighted edge-count statistic: \n")
    if(l0_us<=1){
      cat("  Note: Minimum interval length has been set to l0 = 2 as the weighted edge-count test statistic is not well-defined for interval length<2. \n")
    }
    if(l1_us>=n-1){
      cat("  Note: Maximum interval length has been set to l1 =", n-2, " as the weighted edge-count test statistic is not well-defined for interval length>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    cat("Generalized edge-count statistic: \n")
    if(l0_us<=1){
      cat("  Note: Minimum interval length has been set to l0 = 2 as the generalized edge-count test statistic is not well-defined for interval length< 2. \n")
    }
    if(l1_us>=n-1){
      cat("  Note: Maximum interval length has been set to l1 =", n-2, " as the generalized edge-count test statistic is not well-defined for interval length>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    cat("Max-type edge-count statistic: \n")
    if(l0_us<=1){
      cat("  Note: Minimum interval length has been set to l0 = 2 as the max-type edge-count test statistic is not well-defined for interval length< 2. \n")
    }
    if(l1_us>=n-1){
      cat("  Note: Maximum interval length has been set to l1 =", n-2, " as the max-type edge-count test statistic is not well-defined for interval length>",n-2,". \n")
    }
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$Zmax, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type$pval, "\n")
    }
  }

  return(r1)
}


# the Nu function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

# single change-point
gcp1bynode = function(n, Ebynode, statistics="all", n0=ceiling(0.05*n), n1=floor(0.95*n)){
  # "n" is the total number of nodes.
  # Ebynode[[i]] is the list of nodes that are connect to i by an edge.
  # The nodes are numbered by their order in the sequence.
  # To estimate the change-point, we find the maximum of Z(t), the standardized
  # version of R(t), between n1 and n2.

  nodedeg = rep(0,n)
  for(i in 1:n) nodedeg[i] = length(Ebynode[[i]])
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2

  g = rep(1,n)
  R = rep(0,n)
  R1 = rep(0,n)
  R2 = rep(0,n)

  for(i in 1:(n-1)){
    g[i] = 0  # update g
    links = Ebynode[[i]]

    if(i==1){
      if(length(links)>0){
        R[i] = sum(rep(g[i],length(links)) != g[links])
      } else {
        R[i] = 0
      }
      R1[i]=0
      R2[i]=nE-length(links)
    } else {
      if(length(links)>0){
        add = sum(rep(g[i],length(links)) != g[links])
        subtract = length(links)-add
        R[i] = R[i-1]+add-subtract
        R1[i] = R1[i-1]+subtract
      } else {
        R[i] = R[i-1]
        R1[i]=R1[i-1]
      }
    }

    R2[i] = nE-R[i]-R1[i]
  }

  tt = 1:n
  temp=n0:n1

  scanZ = list()
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    mu.t = nE* 2*tt*(n-tt)/(n*(n-1))
    p1.tt = 2*tt*(n-tt)/(n*(n-1))
    p2.tt = tt*(n-tt)*(n-2)/(n*(n-1)*(n-2))
    p3.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
    A.tt = (p1.tt-2*p2.tt+p3.tt)*nE+(p2.tt-p3.tt)*sumEisq+p3.tt*nE^2
    Z = (mu.t-R)/sqrt(A.tt-mu.t^2)
    Z[n] = 0

    tauhat = temp[which.max(Z[n0:n1])]
    ori = list(tauhat=tauhat, Zmax=Z[tauhat], Z=Z, R=R)
    scanZ$ori = ori
  }
  if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    Rw = ((n-tt-1)*R1+(tt-1)*R2)/(n-2)
    mu.Rw = nE*((n-tt-1)*tt*(tt-1)+(tt-1)*(n-tt)*(n-tt-1))/(n*(n-1)*(n-2))

    mu.R1 = nE*tt*(tt-1)/(n*(n-1))
    mu.R2 = nE*(n-tt)*(n-tt-1)/(n*(n-1))
    v11 = mu.R1*(1-mu.R1) + 2*(0.5*sumEisq-nE)*(tt*(tt-1)*(tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*(tt*(tt-1)*(tt-2)*(tt-3))/(n*(n-1)*(n-2)*(n-3))
    v22 = mu.R2*(1-mu.R2) + 2*(0.5*sumEisq-nE)*((n-tt)*(n-tt-1)*(n-tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*((n-tt)*(n-tt-1)*(n-tt-2)*(n-tt-3))/(n*(n-1)*(n-2)*(n-3))

    v12 = (nE*(nE-1)-2*(0.5*sumEisq-nE))*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3)) - mu.R1*mu.R2

    var.Rw=((n-tt-1)/(n-2))^2*v11 + 2*((n-tt-1)/(n-2))*((tt-1)/(n-2))*v12+((tt-1)/(n-2))^2*v22
    Zw = -(mu.Rw-Rw)/sqrt(apply(cbind(var.Rw,rep(0,n)),1,max))
    
    if (length(which(!is.na(match(c("w","weighted","m","all"),statistics))))>0){
      tauhat = temp[which.max(Zw[n0:n1])]
      weighted = list(tauhat=tauhat, Zmax=Zw[tauhat], Zw=Zw, Rw=Rw)
      scanZ$weighted = weighted
      
    }
    if (length(which(!is.na(match(c("m","max","g","generalized","all"),statistics))))>0){
      Rd = R1-R2
      Zd = (Rd-(mu.R1-mu.R2))/sqrt(apply(cbind(v11+v22-2*v12,rep(0,n)),1,max))
      
      if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
        M = apply(cbind(abs(Zd),Zw),1,max)
        tauhat = temp[which.max(M[n0:n1])]
        max.type = list(tauhat=tauhat, Zmax=M[tauhat], M=M)
        scanZ$max.type = max.type
      }
      if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
        Z = Zw^2 + Zd^2
        tauhat = temp[which.max(Z[n0:n1])]
        generalized = list(tauhat=tauhat, Zmax=Z[tauhat], S=Z)
        scanZ$generalized = generalized
      }
    }
  }
  return(scanZ)
}


# changed interval
gcp2bynode = function(n, Ebynode, statistics="all", l0=ceiling(0.05*n), l1=floor(0.95*n)){

  nodedeg = rep(0,n)
  for(i in 1:n) nodedeg[i] = length(Ebynode[[i]])
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2

  Rtmp = R1 = R2 = Rw = matrix(0,n,n)
 
  
  for (i in 1:(n-1)){
    g = rep(0,n)
    for (j in (i+1):n){
      g[j] = 1 # update g
      links = Ebynode[[j]]
      if (j == (i+1)){
        if (length(links)>0){
          Rtmp[i,j] = sum(rep(g[j],length(links)) != g[links])
          R1[i,j]=0
          R2[i,j]=nE-length(links)
        }
      }else{
        if (length(links)>0){
          add = sum(rep(g[j],length(links)) != g[links])
          subtract = length(links)-add
          Rtmp[i,j] = Rtmp[i,j-1]+add-subtract
          R1[i,j] = R1[i,j-1]+subtract
        }else{
          Rtmp[i,j] = Rtmp[i,j-1]
          R1[i,j] = R1[i,j-1]
        }
      }
      R2[i,j] = nE-Rtmp[i,j]-R1[i,j]
      Rw[i,j] = ((n-j+i-1)/(n-2))*R1[i,j] + (j-i-1)/(n-2)*R2[i,j]
    }
  }

  dif = matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      dif[i,j] = j-i
    }
  }
  difv = as.vector(t(dif))
  ids = which(difv>0)
  ids2 = which((difv>=l0) & (difv<=l1))

  scanZ = list()

  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    tt = 1:n
    mu.t = nE* 2*tt*(n-tt)/(n*(n-1))
    p1.tt = 2*tt*(n-tt)/(n*(n-1))
    p2.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
    V.tt = p2.tt*nE+(p1.tt/2-p2.tt)*sumEisq+(p2.tt-p1.tt^2)*nE^2

    Rv = as.vector(t(Rtmp))
    Zv = rep(0,n*n)
    Zv[ids] = (mu.t[difv[ids]]-Rv[ids])/sqrt(V.tt[difv[ids]])
    Zmax = max(Zv[ids2])
    tauhat0 = which(Zv == Zmax)
    tauhat = c(floor(tauhat0/n)+1, (tauhat0-1)%%n+1)

    ori = list(tauhat=tauhat, Zmax=Zmax, Z=matrix(Zv,n,byrow=TRUE), R=Rtmp, Zv=Zv)
    scanZ$ori = ori
  }
  if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"), statistics))))>0){
    if(l0<=1){
      l0 = 2
    }
    if(l1>=(n-1)){
      l1=n-2
    }
    ids2 = which((difv>=l0) & (difv<=l1))
    
    tt = 1:n

    mu.r1 = nE*tt*(tt-1)/(n*(n-1))
    mu.r2 = nE*(n-tt)*(n-tt-1)/(n*(n-1))

    sig11 = mu.r1*(1-mu.r1) + 2*(0.5*sumEisq-nE)*(tt*(tt-1)*(tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*(tt*(tt-1)*(tt-2)*(tt-3))/(n*(n-1)*(n-2)*(n-3))
    sig22 = mu.r2*(1-mu.r2) + 2*(0.5*sumEisq-nE)*((n-tt)*(n-tt-1)*(n-tt-2))/(n*(n-1)*(n-2)) + (nE*(nE-1)-2*(0.5*sumEisq-nE))*((n-tt)*(n-tt-1)*(n-tt-2)*(n-tt-3))/(n*(n-1)*(n-2)*(n-3))

    sig12=sig21 = (nE*(nE-1)-2*(0.5*sumEisq-nE))*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3)) - mu.r1*mu.r2

    p = (tt-1)/(n-2)
    q = 1-p
    muRw.tt=q*mu.r1+p*mu.r2
    sigRw=q^2*sig11 + p^2*sig22+2*p*q*sig12

    Rw_v = as.vector(t(Rw))
    Zv2 = rep(0,n*n)
    Zv2[ids] = -(muRw.tt[difv[ids]]-Rw_v[ids])/sqrt(sigRw[difv[ids]])

    if (length(which(!is.na(match(c("w","weighted","m","all"), statistics))))>0){
      Zmax = max(Zv2[ids2])
      tauhat0 = which(Zv2 == Zmax)
      tauhat = c(floor(tauhat0/n)+1, (tauhat0-1)%%n+1)
      weighted = list(tauhat=tauhat, Zmax=Zmax, Zw=matrix(Zv2,n,byrow=TRUE), Rw=Rw_v, Zw_v=Zv2)
      scanZ$weighted = weighted
    }
    if (length(which(!is.na(match(c("m","max","g","generalized","all"), statistics))))>0){
      Rsub = R1-R2
      Rsub_v = as.vector(t(Rsub))
      mu1.tt=mu.r1-mu.r2
      sig1 = sig11+sig22-2*sig12
      Zv1 = rep(0,n*n)
      Zv1[ids] = -(mu1.tt[difv[ids]]-Rsub_v[ids])/sqrt(sig1[difv[ids]])
      if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
        # for(i in 1:length(Zv2)){
        #   if(Zv2[i]=='NaN'){
        #     Zv2[i]=0
        #   }
        # }
        M = rep(0,n*n)
        M[ids] =apply(cbind(abs(Zv1[ids]), Zv2[ids]),1,max)
        Zmax = max(M[ids2])
        tauhat0 = which(M == Zmax)
        tauhat = c(floor(tauhat0/n)+1, (tauhat0-1)%%n+1)
        max.type = list(tauhat=tauhat, Zmax=Zmax, M=matrix(M,n,byrow=TRUE), Mv=M)
        scanZ$max.type = max.type
      }
      if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
        Zv = rep(0,n*n)
        Zv[ids] = (Zv1[ids])^2 + (Zv2[ids])^2
        Zmax = max(Zv[ids2])
        tauhat0 = which(Zv == Zmax)
        tauhat = c(floor(tauhat0/n)+1, (tauhat0-1)%%n+1)
        generalized = list(tauhat=tauhat, Zmax=Zmax, S=matrix(Zv,n,byrow=TRUE), Sv=Zv)
        scanZ$generalized = generalized
      }
    }
  }
  return(scanZ)
}



# rho_one = n h_G
rho_one = function(n, s, sumE, sumEisq){
  f1 = 4*(n-1)*(2*s*(n-s)-n)
  f2 = ((n+1)*(n-2*s)^2-2*n*(n-1))
  f3 = 4*((n-2*s)^2-n)
  f4 = 4*n*(s-1)*(n-1)*(n-s-1)
  f5 = n*(n-1)*((n-2*s)^2-(n-2))
  f6 = 4*((n-2)*(n-2*s)^2-2*s*(n-s)+n)
  n*(n-1)*(f1*sumE + f2*sumEisq - f3*sumE^2)/(2*s*(n-s)*(f4*sumE + f5*sumEisq - f6*sumE^2))
}

rho_one_Rw = function(n, t){
  -((2*t^2 - 2*n*t + n)*(n^2 - 3*n + 2)^4)/(2*t*(n - 1)^3*(n - 2)^4*(t - 1)*(n^2 - 2*n*t - n + t^2 + t))
}


# p value approximation for single change-point
pval1 = function(n, E, Ebynode, scanZ, statistics="all", skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)){
  output = list()
  deg = rep(0,n)
  for(i in 1:n) deg[i] = length(Ebynode[[i]])
  sumE = sum(deg)/2
  sumEisq = sum(deg^2)

  if (skew.corr==FALSE){
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    		  # if(lower<2){
       		# lower = 2
    	   	  # }
    		  # if(upper>(n-2)){
       		# upper = n-2
     	  # }
      b = scanZ$ori$Zmax
      if (b>0){
        integrandO = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          x*Nu(sqrt(2*b^2*x))
        }
        pval.ori = dnorm(b)*b*integrate(integrandO, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.ori = 1
      }
      output$ori = min(pval.ori,1)
    }
    
    # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
     # if(lower<2){
       # lower = 2
     # }
     # if(upper>(n-2)){
       # upper = n-2
     # }
    # }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b = scanZ$weighted$Zmax
      if (b>0){
        integrandW = function(t){
          x = rho_one_Rw(n,t)
          x*Nu(sqrt(2*b^2*x))
        }
        pval.weighted = dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.weighted = 1
      }
      output$weighted = min(pval.weighted,1)
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax
      if (b>0){
        integrand1 = function(t){
          x1 = n/(2*t*(n - t))
          x1*Nu(sqrt(2*b^2*x1))
        }
        integrand2 = function(t){
          x2 = rho_one_Rw(n,t)
          x2*Nu(sqrt(2*b^2*x2))
        }
        pval_u1 = 2*dnorm(b)*b*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        pval_u2 = dnorm(b)*b*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
      }else{
        pval.max.type = 1
      }
      output$max.type = pval.max.type
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      b = scanZ$generalized$Zmax
      if (b>0){
        integrandG = function(t,w){
          x1 = n/(2*t*(n - t))
          x2 = rho_one_Rw(n,t)
          2*(x1*cos(w)^2+x2*sin(w)^2)*b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
        }
        integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
        pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.generalized = 1
      }
      output$generalized = min(pval.generalized,1)
    }
    return(output)
  }

  x1 = sum(deg*(deg-1))
  x2 = sum(deg*(deg-1)*(deg-2))
  x3 = 0
  for (i in 1:nrow(E)){
    x3 = x3 + (deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  }
  x4 = sum(deg*(deg-1)*(sumE-deg))
  x5 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    x5 = x5 + length(which(!is.na(match(Ebynode[[j]], Ebynode[[k]]))))
  }

  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    b = scanZ$ori$Zmax
    if (b>0){
      s = 1:n
      x = rho_one(n,s,sumE,sumEisq)
      p1 = 2*s*(n-s)/(n*(n-1))
      p2 = 4*s*(s-1)*(n-s)*(n-s-1)/(n*(n-1)*(n-2)*(n-3))
      p3 = s*(n-s)*((n-s-1)*(n-s-2) + (s-1)*(s-2))/(n*(n-1)*(n-2)*(n-3))
      p4 = 8*s*(s-1)*(s-2)*(n-s)*(n-s-1)*(n-s-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
      mu = p1*sumE
      sig = sqrt(apply(cbind(p2*sumE + (p1/2-p2)*sumEisq + (p2-p1^2)*sumE^2, rep(0,n)), 1, max))  # sigma
      ER3 = p1*sumE + p1/2*3*x1 + p2*(3*sumE*(sumE-1)-3*x1) + p3*x2 + p2/2*(3*x4-6*x3) + p4*(sumE*(sumE-1)*(sumE-2)-x2-3*x4+6*x3)- 2*p4*x5
      r = (mu^3 + 3*mu*sig^2 - ER3)/sig^3
      theta_b = rep(0,n)
      pos = which(1+2*r*b>0)
      theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
      ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
      a = x*Nu(sqrt(2*b^2*x)) * ratio
      nn = n-length(pos)
      if (nn>0.75*n){
        cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
        integrand = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          x*Nu(sqrt(2*b^2*x))
        }
        pval.ori = dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori = min(pval.ori,1)
      }else{
        if (nn>=(lower-1)+(n-upper)){
          neg = which(1+2*r*b<=0)
          dif = neg[2:nn]-neg[1:(nn-1)]
          id1 = which.max(dif)
          id2 = id1 + ceiling(0.03*n)
          id3 = id2 + ceiling(0.09*n)
          inc = (a[id3]-a[id2])/(id3-id2)
          a[id2:1] = a[id2+1]-inc*(1:id2)
          a[(n/2+1):n] = a[(n/2):1]
          neg2 = which(a<0)
          a[neg2] = 0
        }
        integrand = function(s){
          a[s]
        }
        
        result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
        if (is.numeric(result)){
          output$ori = min(result,1)
        }else{
          cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
          b = scanZ$ori$Zmax
          integrand = function(s){
            x = rho_one(n,s,sumE,sumEisq)
            x*Nu(sqrt(2*b^2*x))
          }
          pval.ori = dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$ori = min(pval.ori,1)
        }
      }
    }else{
      output$ori = 1
    }
  }

  if (length(which(!is.na(match(c("w","weighted","m","max","all"), statistics))))>0){
    t = 1:(n-1)
    A1 = sumE*t*(t-1)/(n*(n-1)) + 3*x1*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*(t*(t-1)*(t-2)*(t-3))/(n*(n-1)*(n-2)*(n-3)) +
    2*x5*(t*(t-1)*(t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*t*(t-1)*(t-2)*(t-3)*(t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    B1 = (sumE*(sumE-1)-x1)*(t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(t*(t-1)*(t-2)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    C1 = (sumE*(sumE-1)-x1)*(n-t)*(n-t-1)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(n-t)*(n-t-1)*(n-t-2)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    D1 =  sumE*(n-t)*(n-t-1)/(n*(n-1)) + 3*x1*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*((n-t)*(n-t-1)*(n-t-2)*(n-t-3))/(n*(n-1)*(n-2)*(n-3)) +
    2*x5*((n-t)*(n-t-1)*(n-t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    r1=sumE*(t*(t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))

    r2=sumE*((n-t)*(n-t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))

    r12= (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))

    t = 1:(n-1)
    x = rho_one_Rw(n,t)
    # for(i in 1:length(x)){
    #   if ((abs(x[i]))=="Inf"){
    #     x[i]=0
    #   }
    # }
    q=(n-t-1)/(n-2)
    p=(t-1)/(n-2)

    mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
    sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
    sig = sqrt(sig1)  # sigma
    ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
    r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
  
    b = scanZ$weighted$Zmax
    result.u2 = pval1_sub_2(n,b,r,x,lower,upper)
    r.Rw = r
    x.Rw = x

    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      if (is.numeric(result.u2) && result.u2 > 0){
        output$weighted = min(result.u2,1)
      }else{
        if (result.u2 ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Weighted edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$weighted$Zmax
        if (b>0){
          integrandW = function(t){
            x = rho_one_Rw(n,t)
            x*Nu(sqrt(2*b^2*x))
          }
          pval.weighted = dnorm(b)*b*integrate(integrandW, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        }else{
          pval.weighted = 1
        }
        output$weighted = min(pval.weighted,1)
      }
    }

    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax

      t = 1:(n-1)
      x = n/(2*t*(n - t))
      
      q=1
      p=-1
      mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
      sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
      sig = sqrt(apply(cbind(sig1, rep(0,n-1)), 1, max))  # sigma
      ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
      r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
      
      
      result.u1 = pval1_sub_1(n,b,r,x,lower,upper)
      result.u2 = pval1_sub_2(n,b,r.Rw,x.Rw,lower,upper)

      if (!is.numeric(result.u1) || !is.numeric(result.u2) || result.u1 ==0 || result.u2 ==0){
        if(result.u1 ==0 || result.u2 == 0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Max-type edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$max.type$Zmax
        if (b>0){
          integrand1 = function(t){
            x1 = n/(2*t*(n - t))
            x1*Nu(sqrt(2*b^2*x1))
          }
          integrand2 = function(t){
            x2 = rho_one_Rw(n,t)
            x2*Nu(sqrt(2*b^2*x2))
          }
          pval_u1 = 2*dnorm(b)*b*integrate(integrand1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_u2 = dnorm(b)*b*integrate(integrand2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
        }else{
          pval.max.type = 1
        }
        output$max.type = pval.max.type
      }else{
        output$max.type = 1-(1-min(result.u1,1))*(1-min(result.u2,1))
      }
    }
  }

  # for generalized edge-count test, the approximated p-value without skewness correction is reported
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    b = scanZ$generalized$Zmax
    if (b>0){
      integrandG = function(t,w){
        x1 = n/(2*t*(n - t))
        x2 = rho_one_Rw(n,t)
        2*(x1*cos(w)^2+x2*sin(w)^2)*b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
      }
      integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.generalized = 1
    }
    output$generalized = min(pval.generalized,1)
  }


  #if (length(!is.na(match(c("g","generalized","all"), statistics)))>0){
  # cat("Generalized edge-count statistic: p-value approximation without skewness correction is reported.\n")
  #}

  return(output)
}

# p-value approximation for single change-point, sub functions

pval1_sub_1 = function(n,b,r,x,lower,upper){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    return(0)
  }
  if (nn.l>=lower){
    neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  if (nn.r>=(n-upper)){
    neg = which(1+2*r[ceiling(n/2):(n-1)]*b<=0 )
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    a[(n/2):(n-1)] = (x*Nu(sqrt(2*b^2*x)) * ratio)[(n/2):(n-1)] # update a after extrapolation 
  }
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]
  }
  result = try(2*dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)

}

pval1_sub_2 = function(n,b,r,x,lower,upper){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=(lower-1)+(n-upper)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]
  }
  result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}

# p value approximation for changed interval
pval2 = function(n, E, Ebynode, scanZ, statistics="all", skew.corr=TRUE, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  output = list()
  deg = rep(0,n)
  for(i in 1:n) deg[i] = length(Ebynode[[i]])
  sumE = sum(deg)/2
  sumEisq = sum(deg^2)

  if (skew.corr==FALSE){
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      b = scanZ$ori$Zmax
      if (b>0){
        integrand = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-s)
        }
        pval.ori = dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.ori = 1
      }
      output$ori = min(pval.ori,1)
    }
    # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
      # if(l0<=1){
        # l0 = 2
      # }
      # if(l1>=(n-1)){
        # l1=n-2
      # }
    # }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b = scanZ$weighted$Zmax
      if (b>0){
        integrandW = function(t){
          x = rho_one_Rw(n,t)
          (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-t)
        }
        pval.weighted = try(dnorm(b)/b*integrate(integrandW, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
      }else{
        pval.weighted = 1
      }
      output$weighted = min(pval.weighted,1)
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax
      if (b>0){
        integrand1 = function(t){
          x1 = n/(2*t*(n - t))
          (b^2*x1*Nu(sqrt(2*b^2*x1)))^2*(n-t)
        }
        integrand2 = function(t){
          x2 = rho_one_Rw(n,t)
          (b^2*x2*Nu(sqrt(2*b^2*x2)))^2*(n-t)
        }
        pval_u1 = try(2*dnorm(b)/b*integrate(integrand1, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,silent=T)
        pval_u2 = try(dnorm(b)/b*integrate(integrand2, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,silent=T)
        pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
      }else{
        pval.max.type = 1
      }
      output$max.type = pval.max.type
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      b = scanZ$generalized$Zmax
      if (b>0){
        integrandG = function(t,w){
          x1 = n/(2*t*(n - t))
          x2 = rho_one_Rw(n,t)
          (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2) *b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
        }
        integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
        pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0), l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
      }else{
        pval.generazlied = 1
      }
      output$generalized = min(pval.generalized,1)
    }
    return(output)
  }

  x1 = sum(deg*(deg-1))
  x2 = sum(deg*(deg-1)*(deg-2))
  x3 = 0
  for (i in 1:nrow(E)){
    x3 = x3 + (deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  }
  x4 = sum(deg*(deg-1)*(sumE-deg))
  x5 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    x5 = x5 + length(which(!is.na(match(Ebynode[[j]], Ebynode[[k]]))))
  }

  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    b = scanZ$ori$Zmax
    if (b>0){
      s = 1:n
      x = rho_one(n,s,sumE,sumEisq)
      p1 = 2*s*(n-s)/(n*(n-1))
      p2 = 4*s*(s-1)*(n-s)*(n-s-1)/(n*(n-1)*(n-2)*(n-3))
      p3 = s*(n-s)*((n-s-1)*(n-s-2) + (s-1)*(s-2))/(n*(n-1)*(n-2)*(n-3))
      p4 = 8*s*(s-1)*(s-2)*(n-s)*(n-s-1)*(n-s-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
      mu = p1*sumE
      sig = sqrt(p2*sumE + (p1/2-p2)*sumEisq + (p2-p1^2)*sumE^2)  # sigma
      ER3 = p1*sumE + p1/2*3*x1 + p2*(3*sumE*(sumE-1)-3*x1) + p3*x2 + p2/2*(3*x4-6*x3) + p4*(sumE*(sumE-1)*(sumE-2)-x2-3*x4+6*x3) - 2*p4*x5
      r = (mu^3 + 3*mu*sig^2 - ER3)/sig^3
      theta_b = rep(0,n)
      pos = which(1+2*r*b>0)
      theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
      ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
      a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio
      nn = n-length(pos)
      if (nn>0.75*n){
        #output$ori = 0
        cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$ori$Zmax
        integrand = function(s){
          x = rho_one(n,s,sumE,sumEisq)
          (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-s)
        }
        pval.ori = dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori = min(pval.ori,1)
      } else {
        if (nn>=2*l0-1){
          neg = which(1+2*r*b<=0)
          dif = neg[2:nn]-neg[1:(nn-1)]
          id1 = which.max(dif)
          id2 = id1 + ceiling(0.03*n)
          id3 = id2 + ceiling(0.09*n)
          inc = (a[id3]-a[id2])/(id3-id2)
          a[id2:1] = a[id2+1]-inc*(1:id2)
          a[(n/2+1):n] = a[(n/2):1]
          neg2 = which(a<0)
          a[neg2] = 0
        }
        integrand = function(s){
          a[s]*(n-s)
        }
        result = try(dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
        if (is.numeric(result)){
          output$ori = min(result,1)
        }else{
          cat("Original edge-count statistic: p-value approximation without skewness correction is reported.\n")
          b = scanZ$ori$Zmax
          integrand = function(s){
            x = rho_one(n,s,sumE,sumEisq)
            (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-s)
          }
          pval.ori = dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
          output$ori = min(pval.ori,1)
        }
      }
    }else{
      output$ori = 1
    }
  }

  if (length(which(!is.na(match(c("w","weighted","m","max","all"), statistics))))>0){
    t = 1:(n-1)
    A1 = sumE*t*(t-1)/(n*(n-1)) + 3*x1*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*(t*(t-1)*(t-2)*(t-3))/(n*(n-1)*(n-2)*(n-3)) +
    2*x5*(t*(t-1)*(t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*t*(t-1)*(t-2)*(t-3)*(t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    B1 = (sumE*(sumE-1)-x1)*(t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(t*(t-1)*(t-2)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    C1 = (sumE*(sumE-1)-x1)*(n-t)*(n-t-1)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)) + (x4+2*x5-4*x3)*(n-t)*(n-t-1)*(n-t-2)*t*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*t*(t-1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    D1 =  sumE*(n-t)*(n-t-1)/(n*(n-1)) + 3*x1*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (3*sumE*(sumE-1)-3*x1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))  + x2*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3)) + (6*x3 - 6*x5)*((n-t)*(n-t-1)*(n-t-2)*(n-t-3))/(n*(n-1)*(n-2)*(n-3)) +
    2*x5*((n-t)*(n-t-1)*(n-t-2))/(n*(n-1)*(n-2)) + (3*x4+6*x5-12*x3)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4)) +
    (sumE*(sumE-1)*(sumE-2)+6*x3-2*x5-x2-3*x4)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

    r1=sumE*(t*(t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*t*(t-1)*(t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))

    r2=sumE*((n-t)*(n-t-1)/(n*(n-1))) + 2*(0.5*sumEisq-sumE)*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)) + (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))

    r12= (sumE*(sumE-1)-(2*(0.5*sumEisq-sumE)))*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))

    t = 1:(n-1)
    x = rho_one_Rw(n,t)
    
    q=(n-t-1)/(n-2)
    p=(t-1)/(n-2)

    mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
    sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
    sig = sqrt(sig1)  # sigma
    ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
    r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
   
    b = scanZ$weighted$Zmax
    # if(l0<=1){
      # l0 = 2
    # }
    # if(l1>=(n-1)){
      # l1=n-2
    # }
    
    result.u2 = pval2_sub_2(n,b,r,x,l0,l1)
    r.Rw = r
    x.Rw = x

    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      if (is.numeric(result.u2) && result.u2 >0){
        output$weighted = min(result.u2,1)
      }else{
        if(result.u2 ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Weighted edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$weighted$Zmax
        if (b>0){
          integrandW = function(t){
            x = rho_one_Rw(n,t)
            (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-t)
          }
          pval.weighted = try(dnorm(b)/b*integrate(integrandW, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
        }else{
          pval.weighted = 1
        }
        output$weighted = min(pval.weighted,1)
      }
    }

    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b = scanZ$max.type$Zmax

      t = 1:(n-1)
      x = n/(2*t*(n - t))
     
      q=1
      p=-1
      mu = sumE*(q*t*(t-1)+p*(n-t)*(n-t-1))/(n*(n-1))
      sig1= q^2*r1 + 2*q*p*r12 + p^2*r2 - mu^2
      sig = sqrt(apply(cbind(sig1, rep(0,n-1)), 1, max))  # sigma
      ER3 = q^3*A1 + 3*q^2*p*B1 + 3*q*p^2*C1 + p^3*D1
      r =  (ER3- 3*mu*sig^2 - mu^3)/sig^3
      # for(i in 1:length(r)){
      #   if (is.na(r[i])==TRUE | abs(r[i])=='Inf'){
      #     r[i]=0 }
      # }
      # if (r[n/2]==0) {r[n/2]=r[n/2+1]}

      result.u1 = pval2_sub_1(n,b,r,x,l0,l1)
      result.u2 = pval2_sub_2(n,b,r.Rw,x.Rw,l0,l1)

      if (!is.numeric(result.u1) || !is.numeric(result.u2) || result.u1 ==0 || result.u2 ==0){
        if(result.u1 ==0 || result.u2 ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Max-type edge-count statistic: p-value approximation without skewness correction is reported.\n")
        b = scanZ$max.type$Zmax
        if (b>0){
          integrand1 = function(t){
            x1 = n/(2*t*(n - t))
            (b^2*x1*Nu(sqrt(2*b^2*x1)))^2*(n-t)
          }
          integrand2 = function(t){
            x2 = rho_one_Rw(n,t)
            (b^2*x2*Nu(sqrt(2*b^2*x2)))^2*(n-t)
          }
          pval_u1 = try(2*dnorm(b)/b*integrate(integrand1, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,silent=T)
          pval_u2 = try(dnorm(b)/b*integrate(integrand2, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,silent=T)
          pval.max.type = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
        }else{
          pval.max.type = 1
        }
        output$max.type = pval.max.type
      }else{
        output$max.type = 1-(1-min(result.u1,1))*(1-min(result.u2,1))
      }
    }
  }

  # for generalized edge-count test, the approximated p-value without skewness correction is reported
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    b = scanZ$generalized$Zmax
    if (b>0){
      integrandG = function(t,w){
        x1 = n/(2*t*(n - t))
        x2 = rho_one_Rw(n,t)
        (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2) *b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
      }
      integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0), l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
    }else{
      pval.generazlied = 1
    }
    output$generalized = min(pval.generalized,1)
  }

  return(output)

}

# p-value approximation for changed interval, sub functions
pval2_sub_1 = function(n,b,r,x,l0,l1){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point alternative
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    return(0)
  }
  neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
  if (nn.l>=l0){  ## this part also differs from single change-point
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  neg = which(1+2*r[ceiling(n/2):n]*b<=0)
  if (nn.r>=(n-l1)){ ## this part also differs from single change-point
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point alternative
  }
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]*(n-s)   ## this is different from single change-point alternative
  }
  result = try(2*dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)   ## this is different from single change-point alternative
  return(result)

}


pval2_sub_2 = function(n,b,r,x,l0,l1){
  if (b<0){
    return(1)
  }
  theta_b = rep(0,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)

  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  neg = which(1+2*r*b<=0)
  if (nn>=(l0-1)+(n-l1)){
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE )
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]*(n-s)
  }
  result = try(dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}


# p value from permutation for single change point
permpval1 = function(n, Ebynode, scanZ, statistics="all", B=100, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  # Computes the pvalue P(max_{n1<=t<=n2} Z(t)>b) by permuting the nodes in the graph.
  Z.ori = Z.weighted = Z.max.type = Z.generalized = matrix(0,B,n)
  for(b in 1:B){
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    perm = sample(n)
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]
    }
    gcpstar=gcp1bynode(n,Ebnstar,statistics,n0,n1)

    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      Z.ori[b,] = gcpstar$ori$Z
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      Z.weighted[b,] = gcpstar$weighted$Zw
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      Z.max.type[b,] = gcpstar$max.type$M
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      Z.generalized[b,] = gcpstar$generalized$Z
    }
  }

  output = list()
  p=1-(0:(B-1))/B

  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
  	# if((n0<=1 & n1>=(n-2)) | (n0<=2 & n1>=(n-1))){
      # n0 = 2
      # n1 = n-2
    # }
    maxZ = apply(Z.ori[,n0:n1],1,max)
    maxZs = sort(maxZ)
    output$ori = list(pval=length(which(maxZs>=scanZ$ori$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.ori)
  }
  

  # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    # if((n0<=1 & n1>=(n-2)) | (n0<=2 & n1>=(n-1))){
      # n0 = 2
      # n1 = n-2
    # }
  # }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    
    maxZ = apply(Z.weighted[,n0:n1],1,max)
    maxZs = sort(maxZ)
    output$weighted = list(pval=length(which(maxZs>=scanZ$weighted$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.weighted)
  }
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    maxZ = apply(Z.max.type[,n0:n1],1,max)
    maxZs = sort(maxZ)
    output$max.type = list(pval=length(which(maxZs>=scanZ$max.type$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.max.type)
  }
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    maxZ = apply(Z.generalized[,n0:n1],1,max)
    maxZs = sort(maxZ)
    output$generalized = list(pval=length(which(maxZs>=scanZ$generalized$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z.generalized)
  }

  return(output)
}

# p value from permutation for changed interval
permpval2 = function(n,Ebynode,scanZ,statistics="all", B=100,l0=ceiling(0.05*n),l1=floor(0.95*n)){
  # Computes the pvalue for changed interval by permuting the nodes in the graph.
  Zmax.ori = Zmax.weighted = Zmax.max.type = Zmax.generalized = rep(0,n)
  for(b in 1:B){
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    perm = sample(n)
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]
    }
    gcpstar=gcp2bynode(n,Ebnstar,statistics,l0,l1)

    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      Zmax.ori[b] = gcpstar$ori$Zmax
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      Zmax.weighted[b] = gcpstar$weighted$Zmax
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      Zmax.max.type[b] = gcpstar$max.type$Zmax
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      Zmax.generalized[b] = gcpstar$generalized$Zmax
    }
  }

  output = list()
  p=1-(0:(B-1))/B

  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    maxZ = max(Zmax.ori)
    maxZs = sort(maxZ)
    output$ori = list(pval=length(which(maxZs>=scanZ$ori$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Zmax=Zmax.ori)
  }
  # if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    # if(l0<=1){
      # l0 = 2
    # }
    # if(l1>=(n-1)){
      # l1=n-2
    # }
  # }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    maxZ = max(Zmax.weighted)
    maxZs = sort(maxZ)
    output$weighted = list(pval=length(which(maxZs>=scanZ$weighted$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Zmax=Zmax.weighted)
  }
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    maxZ = max(Zmax.max.type)
    maxZs = sort(maxZ)
    output$max.type = list(pval=length(which(maxZs>=scanZ$max.type$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Zmax=Zmax.max.type)
  }
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    maxZ = max(Zmax.generalized)
    maxZs = sort(maxZ)
    output$generalized = list(pval=length(which(maxZs>=scanZ$generalized$Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Zmax=Zmax.generalized)
  }

  return(output)
}




