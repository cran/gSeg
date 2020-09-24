### main functions
# single change-point
gseg1_discrete = function(n, E, id, statistics=c("all","o","w","g","m"), n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100) {
  r1 = list()
  n0 = ceiling(n0)
  n1 = floor(n1)
  
  r1$scanZ = gcp1bynode_discrete(n, E, id, statistics, n0, n1)
  if (pval.appr==TRUE){
    mypval1_discrete = pval1_discrete(n, E, id, r1$scanZ, statistics, skew.corr, n0, n1)
    r1$pval.appr = mypval1_discrete
  }
  if (pval.perm==TRUE){
    mypval2_discrete = permpval1_discrete(n, E, id, r1$scanZ, statistics, B, n0, n1)
    r1$pval.perm = mypval2_discrete
  }
  
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    cat("Original edge-count scan statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zo_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori_a$pval, "\n")
    }
    cat("Original edge-count scan statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zo_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    cat("Weighted edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zw_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted_a$pval, "\n")
    }
    cat("Weighted edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zw_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    cat("Generalized edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$S_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized_a$pval, "\n")
    }
    cat("Generalized edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$S_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    cat("Max-type edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$M_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type_a$pval, "\n")
    }
    cat("Max-type edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$M_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type_u$pval, "\n")
    }
  }
  
  return(r1)
}


# Changed Interval
gseg2_discrete = function(n, E, id, statistics=c("all","o","w","g","m"), l0=0.05*n, l1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100) {
  r1 = list()
  l0 = ceiling(l0)
  l1 = floor(l1)
  
  r1$scanZ = gcp2bynode_discrete(n, E, id, statistics, l0, l1)
  if (pval.appr==TRUE){
    mypval1_discrete = pval2_discrete(n, E, id, r1$scanZ, statistics, skew.corr, l0, l1)
    r1$pval.appr = mypval1_discrete
  }
  if (pval.perm==TRUE){
    mypval2_discrete = permpval2_discrete(n, E, id, r1$scanZ, statistics, B, l0, l1)
    r1$pval.perm = mypval2_discrete
  }
  
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    cat("Original edge-count scan statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zo_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori_a$pval, "\n")
    }
    cat("Original edge-count scan statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$ori$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$ori$Zo_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$ori_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$ori_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    cat("Weighted edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zw_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted_a$pval, "\n")
    }
    cat("Weighted edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$weighted$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$weighted$Zw_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$weighted_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$weighted_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    cat("Generalized edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$S_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized_a$pval, "\n")
    }
    cat("Generalized edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$generalized$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$generalized$S_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$generalized_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$generalized_u$pval, "\n")
    }
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    cat("Max-type edge-count statistic (a) : \n")
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat_a, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$M_a_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type_a, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type_a$pval, "\n")
    }
    cat("Max-type edge-count statistic (u) : \n")
    cat("  Estimated change-point location:", r1$scanZ$max.type$tauhat_u, "\n")
    cat("  Test statistic:", r1$scanZ$max.type$M_u_max, "\n")
    if (pval.appr==TRUE){
      cat("  Approximated p-value:", r1$pval.appr$max.type_u, "\n")
    }
    if (pval.perm==TRUE){
      cat("  p-value from", B, "permutations:", r1$pval.perm$max.type_u$pval, "\n")
    }
  }
  
  return(r1)
}


# single change-point
gcp1bynode_discrete = function(n, E, id, statistics="all", n0=ceiling(0.05*n), n1=floor(0.95*n)) {
  # "n" is the total number of nodes (=obs).
  # E is the similarity graph (MST, NNL, etc).
  # The nodes are numbered by their order in the sequence.
  # id is the index of observations (the order of observations)
  # To estimate the change-point, we find the maximum of Z(t), the standardized
  # version of R(t), between n1 and n2.
  
  temp = getMV_discrete1(E,id)
  muo_a = temp$muo_a
  muo_u = temp$muo_u
  varo_a = temp$varo_a
  varo_u = temp$varo_u
  mu1_a = temp$mu1_a
  mu2_a = temp$mu2_a
  var1_a = temp$var1_a
  var2_a = temp$var2_a
  var12_a = temp$var12_a
  mu1_u = temp$mu1_u
  mu2_u = temp$mu2_u
  var1_u = temp$var1_u
  var2_u = temp$var2_u
  var12_u = temp$var12_u
  
  p_hat = rep(0,n)
  t = 1:n
  p_hat = (t-1)/(n-2)
  q_hat = 1-p_hat
  
  # mu and variance of the weighted edge-count test (average method)
  muw_a = q_hat*mu1_a + p_hat*mu2_a 
  varw_a = q_hat^2*var1_a + p_hat^2*var2_a + 2*q_hat*p_hat*var12_a
  # mu and variance of the difference of two with-in group edge-counts (average method)
  mud_a = mu1_a - mu2_a 
  vard_a = apply(cbind(var1_a+var2_a-2*var12_a,rep(0,n)),1,max)
  
  # mu and variance of the weighted edge-count test (union method)
  muw_u = q_hat*mu1_u + p_hat*mu2_u 
  varw_u = q_hat^2*var1_u + p_hat^2*var2_u + 2*q_hat*p_hat*var12_u
  # mu and variance of the difference of two with-in group edge-counts (union method)
  mud_u = mu1_u - mu2_u 
  vard_u = apply(cbind(var1_u+var2_u-2*var12_u,rep(0,n)),1,max)
  
  temp = getR1R2_discrete1(E,id)
  R1_a = temp$R1_a
  R2_a = temp$R2_a
  R1_u = temp$R1_u
  R2_u = temp$R2_u
  Ro_a = temp$Ro_a
  Ro_u = temp$Ro_u
  
  Rw_a = q_hat*R1_a + p_hat*R2_a
  Rd_a = R1_a - R2_a	
  Zw_a = (Rw_a-muw_a)/sqrt(varw_a)
  Zd_a = (Rd_a-mud_a)/sqrt(vard_a)	
  
  Rw_u = q_hat*R1_u + p_hat*R2_u
  Rd_u = R1_u - R2_u	
  Zw_u = (Rw_u-muw_u)/sqrt(varw_u)
  Zd_u = (Rd_u-mud_u)/sqrt(vard_u)
  
  temp = n0:n1
  
  scanZ = list()
  
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    Zo_a = -(Ro_a-muo_a)/sqrt(varo_a)
    Zo_u = -(Ro_u-muo_u)/sqrt(varo_u)
    tauhat_a = temp[which.max(Zo_a[n0:n1])]
    tauhat_u = temp[which.max(Zo_u[n0:n1])]
    ori = list(tauhat_a=tauhat_a, Zo_a_max=Zo_a[tauhat_a], Zo_a=Zo_a, Ro_a=Ro_a, 
               tauhat_u=tauhat_u, Zo_u_max=Zo_u[tauhat_u], Zo_u=Zo_u, Ro_u=Ro_u)
    scanZ$ori = ori
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    tauhat_a = temp[which.max(Zw_a[n0:n1])]
    tauhat_u = temp[which.max(Zw_u[n0:n1])]
    weighted = list(tauhat_a=tauhat_a, Zw_a_max=Zw_a[tauhat_a], Zw_a=Zw_a, Rw_a=Rw_a, 
                    tauhat_u=tauhat_u, Zw_u_max=Zw_u[tauhat_u], Zw_u=Zw_u, Rw_u=Rw_u)
    scanZ$weighted = weighted
  }
  if (length(which(!is.na(match
                          (c("g","generalized","all"),statistics))))>0){
    S_a = Zw_a^2 + Zd_a^2
    S_u = Zw_u^2 + Zd_u^2
    tauhat_a = temp[which.max(S_a[n0:n1])]
    tauhat_u = temp[which.max(S_u[n0:n1])]
    generalized = list(tauhat_a=tauhat_a, S_a_max=S_a[tauhat_a], S_a=S_a, 
                       tauhat_u=tauhat_u, S_u_max=S_u[tauhat_u], S_u=S_u)
    scanZ$generalized = generalized
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    M_a = apply(cbind(Zw_a,abs(Zd_a)), 1, max)
    M_u = apply(cbind(Zw_u,abs(Zd_u)), 1, max)
    tauhat_a = temp[which.max(M_a[n0:n1])]
    tauhat_u = temp[which.max(M_u[n0:n1])]
    max.type = list(tauhat_a=tauhat_a, M_a_max=M_a[tauhat_a], M_a=M_a, 
                    tauhat_u=tauhat_u, M_u_max=M_u[tauhat_u], M_u=M_u)
    scanZ$max.type = max.type
  }
  
  return(scanZ)
}


# supporting function for gcp1bynode_discrete
# give the test statistics such as R1_a(t), R2_u(t), Ro_a(t)
getR1R2_discrete1 = function(E, id) {
  n = length(id) # number of obs
  K = length(unique(id)) # number of distinct obs(= number of category)
  
  nE = nrow(E)
  mk = as.numeric(table(id))
  temp1 = sum(mk^2-mk)/2
  temp2 = 0
  for (i in 1:nE){
    e1 = E[i,1]
    e2 = E[i,2]
    temp2 = temp2 + mk[e1]*mk[e2]
  }
  
  R1_a = R2_a = R1_u = R2_u = rep(0,n)
  Ro_a = Ro_u = rep(0,n)
  
  for (i in 1:n) {
    if (i == 1) {
      V = matrix(0, 2, K)
      V[2,] = mk
      V[1,id[i]] = 1
      V[2,id[i]] = V[2,id[i]] - 1
      So_a1 = sum(V[1,]*V[2,]*2/mk)
      So_u1 = sum(V[1,]*V[2,])
      So_a2 = 0
      So_u2 = 0
      for (k in i:nE) {
        So_a2 = So_a2 + (V[1,E[k,1]]*V[2,E[k,2]]+V[1,E[k,2]]*V[2,E[k,1]])/mk[E[k,1]]/mk[E[k,2]]
        So_u2 = So_u2 + (V[1,E[k,1]]*V[2,E[k,2]]+V[1,E[k,2]]*V[2,E[k,1]])
      }
      R1_a[i] = 0
      Ro_a[i] = So_a1 + So_a2
      R2_a[i] = n - K + nE - R1_a[i] - Ro_a[i]
      R1_u[i] = 0
      Ro_u[i] = So_u1 + So_u2
      R2_u[i] = temp1 + temp2 - R1_u[i] - Ro_u[i]
    }
    else {
      V[1,id[i]] = V[1,id[i]] + 1
      V[2,id[i]] = V[2,id[i]] - 1
      A = which(E==id[i], arr.ind = T)
      A[,2] = 3 - A[,2]
      R1_a[i] = R1_a[i-1] + 2*(V[1,id[i]]-1)/mk[id[i]] + sum(V[1,E[A]]/mk[E[A]])/mk[id[i]]
      Ro_a[i] = Ro_a[i-1] + ( 2 - (4*V[1,id[i]]-2)/mk[id[i]] ) + sum((V[2,E[A]]-V[1,E[A]])/mk[E[A]])/mk[id[i]]
      R1_u[i] = R1_u[i-1] + (V[1,id[i]]-1) + sum(V[1,E[A]])
      Ro_u[i] = Ro_u[i-1] + (mk[id[i]] - 2*V[1,id[i]] + 1) + sum((V[2,E[A]]-V[1,E[A]]))
    }
    R2_a[i] = n - K + nE - R1_a[i] - Ro_a[i]
    R2_u[i] = temp1 + temp2 - R1_u[i] - Ro_u[i]
  }
  
  return( list(R1_a = R1_a, R2_a = R2_a, R1_u = R1_u, R2_u = R2_u, Ro_a = Ro_a, Ro_u = Ro_u) )
}


# supporting function for gcp1bynode_discrete
# give the mu, var or cov of test statistics such as E(R1_a(t))
getMV_discrete1 = function(E, id) {
  n = length(id) 
  K = length(unique(id)) 
  nE = dim(E)[1] # the number of edge
  mu1_a=mu2_a=var1_a=var2_a=var12_a=mu1_u=mu2_u=var1_u=var2_u=var12_u= rep()
  muo_a = muo_u = varo_a = varo_u = rep()
  
  t = 1:n
  p1 = t*(t-1)/n/(n-1)
  p2 = t*(t-1)*(t-2)/n/(n-1)/(n-2)
  p3 = t*(t-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = (n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)
  q3 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)/n/(n-1)/(n-2)/(n-3)
  f1 = t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  p0 = t*(n-t)/n/(n-1)
  
  muo_a = (n-K+nE)*2*p0
  mu1_a = (n-K+nE)*p1
  mu2_a = (n-K+nE)*q1
  
  nodedeg_a = nodedeg_u = rep(0,K)
  quan = quan_u = 0
  mk = as.numeric(table(id))
  
  for (i in 1:nE){
    e1 = E[i,1]
    e2 = E[i,2]
    nodedeg_a[e1] = nodedeg_a[e1]+1
    nodedeg_a[e2] = nodedeg_a[e2]+1
    nodedeg_u[e1] = nodedeg_u[e1] + mk[e2]
    nodedeg_u[e2] = nodedeg_u[e2] + mk[e1]
    quan = quan + 1/mk[e1]/mk[e2]
    quan_u = quan_u + mk[e1]*mk[e2]
  }
  
  temp1 = n-K+2*nE+sum(nodedeg_a^2/4/mk)-sum(nodedeg_a/mk)
  temp2 = K-sum(1/mk)
  temp3 = quan
  temp4 = (n-K+nE)^2
  
  varo_a = 4*(p0-4*f1)*temp1+(24*f1-4*p0)*temp2+4*f1*temp3+(4*f1-4*p0^2)*temp4
  var1_a = 4*(p2-p3)*temp1+2*(p1-4*p2+3*p3)*temp2+(p1-2*p2+p3)*temp3+(p3-p1^2)*temp4
  var2_a = 4*(q2-q3)*temp1+2*(q1-4*q2+3*q3)*temp2+(q1-2*q2+q3)*temp3+(q3-q1^2)*temp4
  var12_a = f1*(-4*temp1+6*temp2+temp3)+(f1-p1*q1)*temp4
  
  G = sum(mk*(mk-1))/2 + quan_u
  nodedeg_u = mk-1+nodedeg_u
  
  muo_u = G*2*p0
  mu1_u = G*p1
  mu2_u = G*q1
  
  varo_u = (2*p0-4*f1)*G + (p0-4*f1)*sum(mk*nodedeg_u*(nodedeg_u-1)) + 4*(f1-p0^2)*G^2
  var1_u = (p1-p3)*G + (p2-p3)*sum(mk*nodedeg_u*(nodedeg_u-1)) + (p3-p1^2)*G^2
  var2_u = (q1-q3)*G + (q2-q3)*sum(mk*nodedeg_u*(nodedeg_u-1)) + (q3-q1^2)*G^2
  var12_u = f1*(G^2-G-sum(mk*nodedeg_u*(nodedeg_u-1))) - p1*q1*G^2
  
  return(list(mu1_a=mu1_a,mu2_a=mu2_a,var1_a=var1_a,var2_a=var2_a,var12_a=var12_a,
              mu1_u=mu1_u,mu2_u=mu2_u,var1_u=var1_u,var2_u=var2_u,var12_u=var12_u,
              muo_a=muo_a,muo_u=muo_u,varo_a=varo_a,varo_u=varo_u))
}


# Changed Interval
gcp2bynode_discrete = function(n, E, id, statistics="all", l0=ceiling(0.05*n), l1=floor(0.95*n)) {
  # "n" is the total number of nodes (=obs).
  # E is the similarity graph (MST, NNL, etc).
  # The nodes are numbered by their order in the sequence.
  # id is the index of the data (the order of data)
  # To estimate the change-point, we find the maximum of Z(t), the standardized
  # version of R(t), between n1 and n2.
  
  # use the same function in gcp1bynode_discrete
  temp = getMV_discrete1(E,id)
  muo_a = temp$muo_a
  muo_u = temp$muo_u
  varo_a = temp$varo_a
  varo_u = temp$varo_u
  mu1_a = temp$mu1_a
  mu2_a = temp$mu2_a
  var1_a = temp$var1_a
  var2_a = temp$var2_a
  var12_a = temp$var12_a
  mu1_u = temp$mu1_u
  mu2_u = temp$mu2_u
  var1_u = temp$var1_u
  var2_u = temp$var2_u
  var12_u = temp$var12_u
  
  p_hat = rep(0,n)
  t = 1:n
  p_hat = (t-1)/(n-2)
  q_hat = 1-p_hat
  
  # mu and variance of the weighted edge-count test (average method)
  muw_a = q_hat*mu1_a + p_hat*mu2_a 
  varw_a = q_hat^2*var1_a + p_hat^2*var2_a + 2*q_hat*p_hat*var12_a
  # mu and variance of the difference of two with-in group edge-counts (average method)
  mud_a = mu1_a - mu2_a 
  vard_a = var1_a + var2_a - 2*var12_a
  
  # mu and variance of the weighted edge-count test (union method)
  muw_u = q_hat*mu1_u + p_hat*mu2_u 
  varw_u = q_hat^2*var1_u + p_hat^2*var2_u + 2*q_hat*p_hat*var12_u
  # mu and variance of the difference of two with-in group edge-counts (union method)
  mud_u = mu1_u - mu2_u 
  vard_u = var1_u + var2_u - 2*var12_u
  
  temp = getR1R2_discrete2(E,id)
  R1_a = temp$R1_a
  R2_a = temp$R2_a
  R1_u = temp$R1_u
  R2_u = temp$R2_u
  Ro_a = temp$Ro_a
  Ro_u = temp$Ro_u
  
  Rw_a = temp$Rw_a
  Rw_u = temp$Rw_u
  Rd_a = R1_a - R2_a	
  Rd_u = R1_u - R2_u
  
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
    Ro_av = as.vector(t(Ro_a))
    Zv_a = rep(0,n*n)
    Zv_a[ids] = -(Ro_av[ids]-muo_a[difv[ids]])/sqrt(varo_a[difv[ids]])
    Zo_a_max = max(Zv_a[ids2])
    tauhat_a0 = which(Zv_a == Zo_a_max)
    tauhat_a = c(floor(tauhat_a0/n)+1, (tauhat_a0-1)%%n+1)
    Ro_uv = as.vector(t(Ro_u))
    Zv_u = rep(0,n*n)
    Zv_u[ids] = -(Ro_uv[ids]-muo_u[difv[ids]])/sqrt(varo_u[difv[ids]])
    Zo_u_max = max(Zv_u[ids2])
    tauhat_u0 = which(Zv_u == Zo_u_max)
    tauhat_u = c(floor(tauhat_u0/n)+1, (tauhat_u0-1)%%n+1)
    
    ori = list(tauhat_a=tauhat_a, Zo_a_max=Zo_a_max, Zo_a=matrix(Zv_a,n,byrow=T), Ro_a=Ro_a,
               tauhat_u=tauhat_u, Zo_u_max=Zo_u_max, Zo_u=matrix(Zv_u,n,byrow=T), Ro_u=Ro_u)
    scanZ$ori = ori
  }
  
  if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"), statistics))))>0){
    if (l0 <= 1) {
      l0 = 2
    }
    if (l1 >= (n-1)) {
      l1 = n-2
    }
    ids2 = which((difv>=l0) & (difv<=l1))
    
    Rw_av = as.vector(t(Rw_a))
    Zwv_a = rep(0,n*n)
    Zwv_a[ids] = (Rw_av[ids]-muw_a[difv[ids]])/sqrt(varw_a[difv[ids]])
    Rw_uv = as.vector(t(Rw_u))
    Zwv_u = rep(0,n*n)
    Zwv_u[ids] = (Rw_uv[ids]-muw_u[difv[ids]])/sqrt(varw_u[difv[ids]])
    
    
    if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0) {
      Zw_a_max = max(Zwv_a[ids2])
      tauhat_a0 = which(Zwv_a == Zw_a_max)
      tauhat_a = c(floor(tauhat_a0/n)+1, (tauhat_a0-1)%%n+1)
      Zw_u_max = max(Zwv_u[ids2])
      tauhat_u0 = which(Zwv_u == Zw_u_max)
      tauhat_u = c(floor(tauhat_u0/n)+1, (tauhat_u0-1)%%n+1)
      
      weighted = list(tauhat_a=tauhat_a, Zw_a_max=Zw_a_max, Zw_a=matrix(Zwv_a,n,byrow=T), Rw_a=Rw_a,
                      tauhat_u=tauhat_u, Zw_u_max=Zw_u_max, Zw_u=matrix(Zwv_u,n,byrow=T), Rw_u=Rw_u)
      scanZ$weighted = weighted
    }
    
    if (length(which(!is.na(match(c("m","max","g","generalized","all"), statistics))))>0) {
      Rd_av = as.vector(t(Rd_a))
      Zdv_a = rep(0,n*n)
      Zdv_a[ids] = (Rd_av[ids]-mud_a[difv[ids]])/sqrt(vard_a[difv[ids]])
      Rd_uv = as.vector(t(Rd_u))
      Zdv_u = rep(0,n*n)
      Zdv_u[ids] = (Rd_uv[ids]-mud_u[difv[ids]])/sqrt(vard_u[difv[ids]])
      
      if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0) {
        Sv_a = rep(0, n*n)
        Sv_a[ids] = (Zwv_a[ids])^2 + (Zdv_a[ids])^2
        S_a_max = max(Sv_a[ids2])
        tauhat_a0 = which(Sv_a == S_a_max)
        tauhat_a = c(floor(tauhat_a0/n)+1, (tauhat_a0-1)%%n+1)
        Sv_u = rep(0, n*n)
        Sv_u[ids] = (Zwv_u[ids])^2 + (Zdv_u[ids])^2
        S_u_max = max(Sv_u[ids2])
        tauhat_u0 = which(Sv_u == S_u_max)
        tauhat_u = c(floor(tauhat_u0/n)+1, (tauhat_u0-1)%%n+1)
        
        generalized = list(tauhat_a=tauhat_a, S_a_max=S_a_max, S_a=matrix(Sv_a,n,byrow=T),
                           tauhat_u=tauhat_u, S_u_max=S_u_max, S_u=matrix(Sv_u,n,byrow=T))
        scanZ$generalized = generalized
      }
      
      if (length(which(!is.na(match(c("m","max","all"), statistics))))>0) {
        for(i in 1:length(Zwv_a)){
          if(Zwv_a[i]=='NaN'){
            Zwv_a[i]=0
          }
        }
        M_a = rep(0,n*n)
        M_a[ids] = apply(cbind(abs(Zdv_a[ids]),Zwv_a[ids]), 1, max)
        M_a_max = max(M_a[ids2])
        tauhat_a0 = which(M_a == M_a_max)
        tauhat_a = c(floor(tauhat_a0/n)+1, (tauhat_a0-1)%%n+1)
        for(i in 1:length(Zwv_u)){
          if(Zwv_u[i]=='NaN'){
            Zwv_u[i]=0
          }
        }
        M_u = rep(0,n*n)
        M_u[ids] = apply(cbind(abs(Zdv_u[ids]),Zwv_u[ids]), 1, max)
        M_u_max = max(M_u[ids2])
        tauhat_u0 = which(M_u == M_u_max)
        tauhat_u = c(floor(tauhat_u0/n)+1, (tauhat_u0-1)%%n+1)
        
        max.type = list(tauhat_a=tauhat_a, M_a_max=M_a_max, M_a=matrix(M_a,n,byrow=TRUE), 
                        tauhat_u=tauhat_u, M_u_max=M_u_max, M_u=matrix(M_u,n,byrow=TRUE))
        scanZ$max.type = max.type
      }
    }
  }
  return(scanZ)
}

# supporting function for gcp2bynode_discrete
# give the test statistics such as R1_a(t1,t2), R2_u(t1,t2), Ro_a(t1,t2), Rw_a(t1,t2)
getR1R2_discrete2 = function(E, id) {
  n = length(id) # number of obs
  K = length(unique(id)) # number of distinct obs(= number of category)
  
  nE = nrow(E)
  mk = as.numeric(table(id))
  temp1 = sum(mk^2-mk)/2
  temp2 = 0
  for (i in 1:nE){
    e1 = E[i,1]
    e2 = E[i,2]
    temp2 = temp2 + mk[e1]*mk[e2]
  }
  
  R1_a = R2_a = R1_u = R2_u = matrix(0,n,n)
  Ro_a = Ro_u = Rw_a = Rw_u = matrix(0,n,n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (j == (i+1)) {
        V = matrix(0, 2, K)
        V[2,] = mk
        V[1,id[j]] = 1
        V[2,id[j]] = V[2,id[j]] - 1
        So_a1 = sum(V[1,]*V[2,]*2/mk)
        So_u1 = sum(V[1,]*V[2,])
        So_a2 = 0
        So_u2 = 0
        for (k in 1:nE) {
          So_a2 = So_a2 + (V[1,E[k,1]]*V[2,E[k,2]]+V[1,E[k,2]]*V[2,E[k,1]])/mk[E[k,1]]/mk[E[k,2]]
          So_u2 = So_u2 + (V[1,E[k,1]]*V[2,E[k,2]]+V[1,E[k,2]]*V[2,E[k,1]])
        }
        
        R1_a[i,j] = 0
        Ro_a[i,j] = So_a1 + So_a2
        R2_a[i,j] = n - K + nE - R1_a[i,j] - Ro_a[i,j]
        R1_u[i,j] = 0
        Ro_u[i,j] = So_u1 + So_u2
        R2_u[i,j] = temp1 + temp2 - R1_u[i,j] - Ro_u[i,j]
      }
      else {
        V[1,id[j]] = V[1,id[j]] + 1
        V[2,id[j]] = V[2,id[j]] - 1
        probmat = which(E==id[j], arr.ind = T)
        probmat[,2] = 3 - probmat[,2]
        R1_a[i,j] = R1_a[i,j-1] + 2*(V[1,id[j]]-1)/mk[id[j]] + sum(V[1,E[probmat]]/mk[E[probmat]])/mk[id[j]]
        Ro_a[i,j] = Ro_a[i,j-1] + ( 2 - (4*V[1,id[j]]-2)/mk[id[j]] ) + sum((V[2,E[probmat]]-V[1,E[probmat]])/mk[E[probmat]])/mk[id[j]]
        R1_u[i,j] = R1_u[i,j-1] + (V[1,id[j]]-1) + sum(V[1,E[probmat]])
        Ro_u[i,j] = Ro_u[i,j-1] + (mk[id[j]] - 2*V[1,id[j]] + 1) + sum((V[2,E[probmat]]-V[1,E[probmat]]))
      }
      R2_a[i,j] = n - K + nE - R1_a[i,j] - Ro_a[i,j]
      R2_u[i,j] = temp1 + temp2 - R1_u[i,j] - Ro_u[i,j]
      Rw_a[i,j] = ((n-j+i-1)/(n-2))*R1_a[i,j] + (j-i-1)/(n-2)*R2_a[i,j]
      Rw_u[i,j] = ((n-j+i-1)/(n-2))*R1_u[i,j] + (j-i-1)/(n-2)*R2_u[i,j]
    }
  }
  
  return( list(R1_a = R1_a, R2_a = R2_a, R1_u = R1_u, R2_u = R2_u, Ro_a = Ro_a, Ro_u = Ro_u, Rw_a = Rw_a, Rw_u = Rw_u) )
}






# p value from permutation for single change point (permutation p-value)
permpval1_discrete = function(n, E, id, scanZ, statistics="all", B=100, n0=ceiling(0.05*n), n1=floor(0.95*n)) {
  # Computes the pvalue P(max_{n1<=t<=n2}Z(t) > b) by permuting the nodes(obs) in the graph.
  Z.ori_a=Z.ori_u=Z.weighted_a=Z.weighted_u=Z.max.type_a=Z.max.type_u = matrix(0,B,n)
  Z.generalized_a=Z.generalized_u = matrix(0,B,n)
  
  for(b in 1:B) {
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    id0 = sample(id, replace = FALSE)    # permute the data
    gcpstar = gcp1bynode_discrete(n, E, id0, statistics, n0, n1)
    
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      Z.ori_a[b,] = gcpstar$ori$Zo_a
      Z.ori_u[b,] = gcpstar$ori$Zo_u
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      Z.weighted_a[b,] = gcpstar$weighted$Zw_a
      Z.weighted_u[b,] = gcpstar$weighted$Zw_u
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      Z.max.type_a[b,] = gcpstar$max.type$M_a
      Z.max.type_u[b,] = gcpstar$max.type$M_u
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      Z.generalized_a[b,] = gcpstar$generalized$S_a
      Z.generalized_u[b,] = gcpstar$generalized$S_u
    }
  }
  
  output = list()
  p=1-(0:(B-1))/B
  
  if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    if((n0<=1 & n1>=(n-2)) | (n0<=2 & n1>=(n-1))){
      n0 = 2
      n1 = n-2
    }
  }
  
  # pval : permutation p-value, curve : distribution of B max(Z(t))
  # maxZs : B max(Z(t)) after calculation by B permutation
  # Z : B Z(t) (B x n matrix)
  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    maxZ_a = apply(Z.ori_a[,n0:n1],1,max)
    maxZs_a = sort(maxZ_a)
    output$ori_a = list(pval=length(which(maxZs_a>=scanZ$ori$Zo_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Z.ori_a)
    maxZ_u = apply(Z.ori_u[,n0:n1],1,max)
    maxZs_u = sort(maxZ_u)
    output$ori_u = list(pval=length(which(maxZs_u>=scanZ$ori$Zo_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Z.ori_u)
  }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    maxZ_a = apply(Z.weighted_a[,n0:n1],1,max)
    maxZs_a = sort(maxZ_a)
    output$weighted_a = list(pval=length(which(maxZs_a>=scanZ$weighted$Zw_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Z.weighted_a) 
    maxZ_u = apply(Z.weighted_u[,n0:n1],1,max)
    maxZs_u = sort(maxZ_u)
    output$weighted_u = list(pval=length(which(maxZs_u>=scanZ$weighted$Zw_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Z.weighted_u)
  }
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    maxZ_a = apply(Z.max.type_a[,n0:n1],1,max)
    maxZs_a = sort(maxZ_a)
    output$max.type_a = list(pval=length(which(maxZs_a>=scanZ$max.type$M_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Z.max.type_a)
    maxZ_u = apply(Z.max.type_u[,n0:n1],1,max)
    maxZs_u = sort(maxZ_u)
    output$max.type_u = list(pval=length(which(maxZs_u>=scanZ$max.type$M_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Z.max.type_u)
  }
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    maxZ_a = apply(Z.generalized_a[,n0:n1],1,max)
    maxZs_a = sort(maxZ_a)
    output$generalized_a = list(pval=length(which(maxZs_a>=scanZ$generalized$S_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Z.generalized_a)
    maxZ_u = apply(Z.generalized_u[,n0:n1],1,max)
    maxZs_u = sort(maxZ_u)
    output$generalized_u = list(pval=length(which(maxZs_u>=scanZ$generalized$S_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Z.generalized_u)
  }
  
  return(output)
}
  

# p value from permutation for changed interval (permutation p-value)
permpval2_discrete = function(n, E, id, scanZ, statistics="all", B=100, l0=ceiling(0.05*n), l1=floor(0.95*n)) {
  # Computes the pvalue P(max_{n1<=t<=n2}Z(t) > b) by permuting the nodes(obs) in the graph.
  Zmax.ori_a=Zmax.ori_u=Zmax.weighted_a=Zmax.weighted_u=Zmax.max.type_a=Zmax.max.type_u = rep(0,n)
  Zmax.generalized_a=Zmax.generalized_u = rep(0,n)
  
  for(b in 1:B) {
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    id0 = sample(id, replace = FALSE)    # permute the data
    gcpstar = gcp2bynode_discrete(n, E, id0, statistics, l0, l1)
    
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      Zmax.ori_a[b] = gcpstar$ori$Zo_a_max
      Zmax.ori_u[b] = gcpstar$ori$Zo_u_max
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      Zmax.weighted_a[b] = gcpstar$weighted$Zw_a_max
      Zmax.weighted_u[b] = gcpstar$weighted$Zw_u_max
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      Zmax.max.type_a[b] = gcpstar$max.type$M_a_max
      Zmax.max.type_u[b] = gcpstar$max.type$M_u_max
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      Zmax.generalized_a[b] = gcpstar$generalized$S_a_max
      Zmax.generalized_u[b] = gcpstar$generalized$S_u_max
    }
  }
  
  output = list()
  p=1-(0:(B-1))/B
  
  # pval : permutation p-value, curve : distribution of B max(Z(t))
  # maxZs : B max(Z(t)) after calculation by B permutation
  # Z : B Z(t) (B x n matrix)
  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    maxZ_a = Zmax.ori_a
    maxZs_a = sort(maxZ_a)
    output$ori_a = list(pval=length(which(maxZs_a>=scanZ$ori$Zo_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Zmax.ori_a)
    maxZ_u = Zmax.ori_u
    maxZs_u = sort(maxZ_u)
    output$ori_u = list(pval=length(which(maxZs_u>=scanZ$ori$Zo_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Zmax.ori_u)
  }
  if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
    if(l0<=1){
      l0 = 2
    }
    if(l1>=(n-1)){
      l1=n-2
    }
  }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    maxZ_a = Zmax.weighted_a
    maxZs_a = sort(maxZ_a)
    output$weighted_a = list(pval=length(which(maxZs_a>=scanZ$weighted$Zw_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Zmax.weighted_a) 
    maxZ_u = Zmax.weighted_u
    maxZs_u = sort(maxZ_u)
    output$weighted_u = list(pval=length(which(maxZs_u>=scanZ$weighted$Zw_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Zmax.weighted_u)
  }
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    maxZ_a = Zmax.max.type_a
    maxZs_a = sort(maxZ_a)
    output$max.type_a = list(pval=length(which(maxZs_a>=scanZ$max.type$M_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Zmax.max.type_a)
    maxZ_u = Zmax.max.type_u
    maxZs_u = sort(maxZ_u)
    output$max.type_u = list(pval=length(which(maxZs_u>=scanZ$max.type$M_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Zmax.max.type_u)
  }
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    maxZ_a = Zmax.generalized_a
    maxZs_a = sort(maxZ_a)
    output$generalized_a = list(pval=length(which(maxZs_a>=scanZ$generalized$S_a_max))/B, curve=cbind(maxZs_a,p), maxZs_a=maxZs_a, Z=Zmax.generalized_a)
    maxZ_u = Zmax.generalized_u
    maxZs_u = sort(maxZ_u)
    output$generalized_u = list(pval=length(which(maxZs_u>=scanZ$generalized$S_u_max))/B, curve=cbind(maxZs_u,p), maxZs_u=maxZs_u, Z=Zmax.generalized_u)
  }
  
  return(output)
}



## supporting function for analytical p-value (pval1_discrete)
# the v(x) function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}
# c_w(t) = rho_w(t) for Zw(t)
rho_w = function(n, t) {
  ((2*t-1)*(n-t)*(n-t-1) - t*(t-1)*(2*t-2*n+1))/2/t/(t-1)/(n-t)/(n-t-1)
}
# c_d(t) = rho_d(t) for Zd(t)
rho_d = function(n, t) {
  n/2/t/(n-t)
}
# rho_one_discrete = n h_G
rho_one_discrete = function(n, t, one, two, three){
  f1 = 4*(n-1)*(2*t*(n-t)-n)
  f2 = ((n+1)*(n-2*t)^2-2*n*(n-1))
  f3 = 4*((n-2*t)^2-n)
  f4 = 4*n*(t-1)*(n-1)*(n-t-1)
  f5 = n*(n-1)*((n-2*t)^2-(n-2))
  f6 = 4*((n-2)*(n-2*t)^2-2*t*(n-t)+n)
  n*(n-1)*(f1*one + f2*two - f3*three)/(2*t*(n-t)*(f4*one + f5*two - f6*three))
}


# p value approximation for single change-point (analytical p-value)
pval1_discrete = function(n, E, id, scanZ, statistics="all", skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)) {
  output = list()
  t = 1:n
  t = as.numeric(t)
  k = length(unique(id))
  nE = nrow(E)
  mk = as.numeric(table(id))
  
  nodedeg_a = nodedeg_u = rep(0,k)
  quan = quan_u = 0
  for (i in 1:nE) {
    e1 = E[i,1]
    e2 = E[i,2]
    nodedeg_a[e1] = nodedeg_a[e1] + 1
    nodedeg_a[e2] = nodedeg_a[e2] + 1
    nodedeg_u[e1] = nodedeg_u[e1] + mk[e2]
    nodedeg_u[e2] = nodedeg_u[e2] + mk[e1]
    quan = quan + 1/mk[e1]/mk[e2]
    quan_u = quan_u + mk[e1]*mk[e2]
  }
  temp3 = sum(1/mk)
  temp1 = sum(nodedeg_a/mk)
  temp8 = sum(nodedeg_a^2/mk)
  
  G = sum(mk*(mk-1))/2 + quan_u
  nodedeg_u = nodedeg_u + mk - 1
  
  one_a = 2*(2*k - 2*temp3 + quan)
  two_a = 4*(n-2*k+2*nE+temp8/4-temp1+temp3)
  three_a = (n-k+nE)^2
  one_u = G
  two_u = sum(mk*nodedeg_u^2)
  three_u = G^2

  
  # approximated p-value without skewness correction
  if (skew.corr==FALSE) {
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      b_a = scanZ$ori$Zo_a_max
      b_u = scanZ$ori$Zo_u_max
      integrando_a = function(t){
        c1 = rho_one_discrete(n,t,one_a,two_a,three_a)
        c1*Nu(b_a*sqrt(2*c1))
      }
      integrando_u = function(t){
        c1 = rho_one_discrete(n,t,one_u,two_u,three_u)
        c1*Nu(b_u*sqrt(2*c1))
      }
      pval.ori_a = b_a*dnorm(b_a)*integrate(integrando_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.ori_u = b_u*dnorm(b_u)*integrate(integrando_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      output$ori_a = min(pval.ori_a,1)
      output$ori_u = min(pval.ori_u,1)
    }
    if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
      if(lower<2){
        lower = 2
      }
      if(upper>(n-2)){
        upper = n-2
      }
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b_a = scanZ$weighted$Zw_a_max
      b_u = scanZ$weighted$Zw_u_max
      integrandW_a = function(t){
        c1 = rho_w(n, t)
        c1*Nu(b_a*sqrt(2*c1))
      }
      integrandW_u = function(t){
        c1 = rho_w(n, t)
        c1*Nu(b_u*sqrt(2*c1))
      }
      pval.weighted_a = b_a*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.weighted_u = b_u*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      output$weighted_a = min(pval.weighted_a,1)
      output$weighted_u = min(pval.weighted_u,1)
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b_a = scanZ$max.type$M_a_max
      b_u = scanZ$max.type$M_u_max
      integrandD_a = function(t){
        c1 = rho_d(n, t)
        c1*Nu(b_a*sqrt(2*c1))
      }
      integrandD_u = function(t){
        c1 = rho_d(n, t)
        c1*Nu(b_u*sqrt(2*c1))
      }
      integrandW_a = function(t){
        c1 = rho_w(n, t)
        c1*Nu(b_a*sqrt(2*c1))
      }
      integrandW_u = function(t){
        c1 = rho_w(n, t)
        c1*Nu(b_u*sqrt(2*c1))
      }
      pval_a1 = 2*b_a*dnorm(b_a)*integrate(integrandD_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_a2 = b_a*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.max.type_a = as.numeric(1-(1-min(pval_a1,1))*(1-min(pval_a2,1)))
      pval_u1 = 2*b_u*dnorm(b_u)*integrate(integrandD_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_u2 = b_u*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.max.type_u = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
      output$max.type_a = pval.max.type_a
      output$max.type_u = pval.max.type_u
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      b_a = scanZ$generalized$S_a_max
      b_u = scanZ$generalized$S_u_max
      integrandG_a = function(t,w){
        x1 = rho_d(n,t)
        x2 = rho_w(n,t)
        2*(x1*cos(w)^2+x2*sin(w)^2)*b_a*Nu(sqrt(2*b_a*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
      }
      integrand0_a = function(t) {integrate(integrandG_a,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      integrandG_u = function(t,w){
        x1 = rho_d(n,t)
        x2 = rho_w(n,t)
        2*(x1*cos(w)^2+x2*sin(w)^2)*b_u*Nu(sqrt(2*b_u*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
      }
      integrand0_u = function(t) {integrate(integrandG_u,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      pval.generalized_a = dchisq(b_a,2)*integrate(Vectorize(integrand0_a),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      pval.generalized_u = dchisq(b_u,2)*integrate(Vectorize(integrand0_u),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      output$generalized_a = min(pval.generalized_a,1)
      output$generalized_u = min(pval.generalized_u,1)
    }
     
    return(output)
  }
  
  
  ##########
  # approximated p-value with skewness correction
  if (length(which(!is.na(match(c("o","ori","original","w","weighted","m","max","all"), statistics))))>0) {
    
    temp = skewcorr(E, id) # give statistics for third moment of E(Rw_a^3) etc
    ER3w_a = temp$ER3w_a
    ER3w_u = temp$ER3w_u
    ER3d_a = temp$ER3d_a
    ER3d_u = temp$ER3d_u
    ER3_a = temp$ER3_a
    ER3_u = temp$ER3_u
    
    temp = getMV_discrete1(E, id)
    muo_a = temp$muo_a
    muo_u = temp$muo_u
    varo_a = temp$varo_a
    varo_u = temp$varo_u
    mu1_a = temp$mu1_a
    mu2_a = temp$mu2_a
    var1_a = temp$var1_a
    var2_a = temp$var2_a
    var12_a = temp$var12_a
    mu1_u = temp$mu1_u
    mu2_u = temp$mu2_u
    var1_u = temp$var1_u
    var2_u = temp$var2_u
    var12_u = temp$var12_u
    
    p_hat = rep(0,n)
    p_hat = (t-1)/(n-2)
    q_hat = 1-p_hat
    
    # mu and variance of the weighted edge-count test (average method)
    muw_a = q_hat*mu1_a + p_hat*mu2_a 
    varw_a = q_hat^2*var1_a + p_hat^2*var2_a + 2*q_hat*p_hat*var12_a
    muw_u = q_hat*mu1_u + p_hat*mu2_u 
    varw_u = q_hat^2*var1_u + p_hat^2*var2_u + 2*q_hat*p_hat*var12_u
    # mu and variance of the difference of two with-in group edge-counts (average method)
    mud_a = mu1_a - mu2_a 
    vard_a = var1_a + var2_a - 2*var12_a
    mud_u = mu1_u - mu2_u
    vard_u = var1_u + var2_u - 2*var12_u
    
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0) {
      A = n-k+nE
      ER3o_a = A^3-3*A^2*muo_a+3*A*(varo_a+muo_a^2) - ER3_a
      ER3o_u = G^3-3*G^2*muo_u+3*G*(varo_u+muo_u^2) - ER3_u
      
      ro_a = (-ER3o_a + 3*muo_a*varo_a + muo_a^3)/(varo_a^(3/2)) # E(Zo_a^3)
      for(i in 1:length(ro_a)){
        if (is.na(ro_a[i])==TRUE){
          ro_a[i]=0
        }
      }
      ro_u = (-ER3o_u + 3*muo_u*varo_u + muo_u^3)/(varo_u^(3/2)) # E(Zo_u^3)
      for(i in 1:length(ro_u)){
        if (is.na(ro_u[i])==TRUE){
          ro_u[i]=0
        }
      }
      c1_o_a = rho_one_discrete(n, t, one_a, two_a, three_a)
      for(i in 1:length(c1_o_a)){
        if ((abs(c1_o_a[i]))=="Inf"){
          c1_o_a[i]=0
        }
      }
      c1_o_u = rho_one_discrete(n, t, one_u, two_u, three_u)
      for(i in 1:length(c1_o_u)){
        if ((abs(c1_o_u[i]))=="Inf"){
          c1_o_u[i]=0
        }
      }
      
      b_a = scanZ$ori$Zo_a_max
      b_u = scanZ$ori$Zo_u_max
      result_o_a = pval1_discrete_sub_2(n, b_a, ro_a, c1_o_a, lower, upper)
      result_o_u = pval1_discrete_sub_2(n, b_u, ro_u, c1_o_u, lower, upper)
      # average
      if (is.numeric(result_o_a) && result_o_a > 0){
        output$ori_a = min(result_o_a, 1)
      }else{
        if (result_o_a ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Original edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
        b_a = scanZ$ori$Zo_a_max
        integrando_a = function(t){
          c1 = rho_one_discrete(n,t,one_a,two_a,three_a)
          c1*Nu(b_a*sqrt(2*c1))
        }
        pval.ori_a = b_a*dnorm(b_a)*integrate(integrando_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori_a = min(pval.ori_a,1)
      }
      # union
      if (is.numeric(result_o_u) && result_o_u > 0){
        output$ori_u = min(result_o_u, 1)
      }else{
        if (result_o_u ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Original edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
        b_u = scanZ$ori$Zo_u_max
        integrando_u = function(t){
          c1 = rho_one_discrete(n,t,one_u,two_u,three_u)
          c1*Nu(b_u*sqrt(2*c1))
        }
        pval.ori_u = b_u*dnorm(b_u)*integrate(integrando_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori_u = min(pval.ori_u,1)
      }
    }
    
    if (length(which(!is.na(match(c("w","weighted","m","max","all"), statistics))))>0) {
      rw_a =  (ER3w_a- 3*muw_a*varw_a - muw_a^3)/(varw_a^(3/2)) # E(Zw_a^3)
      for(i in 1:length(rw_a)){
        if (is.na(rw_a[i])==TRUE){
          rw_a[i]=0
        }
      }
      rw_u =  (ER3w_u- 3*muw_u*varw_u - muw_u^3)/(varw_u^(3/2)) # E(Zw_u^3)
      for(i in 1:length(rw_u)){
        if (is.na(rw_u[i])==TRUE){
          rw_u[i]=0
        }
      }
      c1_w = rho_w(n, t)
      for(i in 1:length(c1_w)){
        if ((abs(c1_w[i]))=="Inf"){
          c1_w[i]=0
        }
      }
      if(lower<2){
        lower = 2
        #print(lower)
      }
      if(upper>(n-2)){
        upper = n-2
        #print(upper)
      }
      # weighted case
      if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
        b_a = scanZ$weighted$Zw_a_max
        b_u = scanZ$weighted$Zw_u_max
        result_w_a = pval1_discrete_sub_2(n, b_a, rw_a, c1_w, lower, upper)
        result_w_u = pval1_discrete_sub_2(n, b_u, rw_u, c1_w, lower, upper)
        # average
        if (is.numeric(result_w_a) && result_w_a > 0){
          output$weighted_a = min(result_w_a, 1)
        }else{
          if (result_w_a ==0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Weighted edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
          b_a = scanZ$weighted$Zw_a_max
          integrandW_a = function(t){
            c1 = rho_w(n, t)
            c1*Nu(b_a*sqrt(2*c1))
          }
          pval.weighted_a = b_a*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$weighted_a = min(pval.weighted_a,1)
        }
        # union
        if (is.numeric(result_w_u) && result_w_u > 0){
          output$weighted_u = min(result_w_u, 1)
        }else{
          if (result_w_u ==0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Weighted edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
          b_u = scanZ$weighted$Zw_u_max
          integrandW_u = function(t){
            c1 = rho_w(n, t)
            c1*Nu(b_u*sqrt(2*c1))
          }
          pval.weighted_u = b_u*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$weighted_u = min(pval.weighted_u,1)
        }
      }
      # max-type case
      if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
        c1_d = rho_d(n, t)
        for(i in 1:length(c1_d)){
          if ((abs(c1_d[i]))=="Inf"){
            c1_d[i]=0
          }
        }
        rd_a =  (ER3d_a - 3*mud_a*vard_a - mud_a^3)/(vard_a^(3/2))
        for(i in 1:length(rd_a)){
          if (is.na(rd_a[i])==TRUE){
            rd_a[i]=0
          }
        }
        if (rd_a[n/2]==0) {rd_a[n/2]=rd_a[n/2+1]}
        
        rd_u =  (ER3d_u - 3*mud_u*vard_u - mud_u^3)/(vard_u^(3/2))
        for(i in 1:length(rd_u)){
          if (is.na(rd_u[i])==TRUE){
            rd_u[i]=0
          }
        }
        if (rd_u[n/2]==0) {rd_u[n/2]=rd_u[n/2+1]}
        
        b_a = scanZ$max.type$M_a_max
        b_u = scanZ$max.type$M_u_max
        result_d_a = pval1_discrete_sub_1(n, b_a, rd_a, c1_d, lower, upper)
        result_w_a = pval1_discrete_sub_2(n, b_a, rw_a, c1_w, lower, upper)
        result_d_u = pval1_discrete_sub_1(n, b_u, rd_u, c1_d, lower, upper)
        result_w_u = pval1_discrete_sub_2(n, b_u, rw_u, c1_w, lower, upper)
        
        # average
        if (!is.numeric(result_d_a) || !is.numeric(result_w_a) || result_d_a ==0 || result_w_a ==0){
          if(result_d_a ==0 || result_w_a == 0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Max-type edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
          b_a = scanZ$max.type$M_a_max
          integrandD_a = function(t){
            c1 = rho_d(n, t)
            c1*Nu(b_a*sqrt(2*c1))
          }
          integrandW_a = function(t){
            c1 = rho_w(n, t)
            c1*Nu(b_a*sqrt(2*c1))
          }
          pval_a1 = 2*b_a*dnorm(b_a)*integrate(integrandD_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_a2 = b_a*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type_a = as.numeric(1-(1-min(pval_a1,1))*(1-min(pval_a2,1)))
          output$max.type_a = pval.max.type_a
        }else{
          output$max.type_a = 1-(1-min(result_d_a,1))*(1-min(result_w_a,1))
        }
        # union
        if (!is.numeric(result_d_u) || !is.numeric(result_w_u) || result_d_u ==0 || result_w_u ==0){
          if(result_d_u ==0 || result_w_u == 0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Max-type edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
          b_u = scanZ$max.type$M_u_max
          integrandD_u = function(t){
            c1 = rho_d(n, t)
            c1*Nu(b_u*sqrt(2*c1))
          }
          integrandW_u = function(t){
            c1 = rho_w(n, t)
            c1*Nu(b_u*sqrt(2*c1))
          }
          pval_u1 = 2*b_u*dnorm(b_u)*integrate(integrandD_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_u2 = b_u*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type_u = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
          output$max.type_u = pval.max.type_u
        }else{
          output$max.type_u = 1-(1-min(result_d_u,1))*(1-min(result_w_u,1))
        }
      }
    }
  }
  
  # for generalized edge-count test, the approximated p-value without skewness correction is reported
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0) {
    b_a = scanZ$generalized$S_a_max
    b_u = scanZ$generalized$S_u_max
    integrandG_a = function(t,w){
      x1 = rho_d(n,t)
      x2 = rho_w(n,t)
      2*(x1*cos(w)^2+x2*sin(w)^2)*b_a*Nu(sqrt(2*b_a*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
    }
    integrand0_a = function(t) {integrate(integrandG_a,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    integrandG_u = function(t,w){
      x1 = rho_d(n,t)
      x2 = rho_w(n,t)
      2*(x1*cos(w)^2+x2*sin(w)^2)*b_u*Nu(sqrt(2*b_u*(x1*cos(w)^2+x2*sin(w)^2)))/(2*pi)
    }
    integrand0_u = function(t) {integrate(integrandG_u,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    pval.generalized_a = dchisq(b_a,2)*integrate(Vectorize(integrand0_a),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    pval.generalized_u = dchisq(b_u,2)*integrate(Vectorize(integrand0_u),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    output$generalized_a = min(pval.generalized_a,1)
    output$generalized_u = min(pval.generalized_u,1)
  }
    
    return(output)  
}
  
# Support function for approximated p-value (pval1_discrete, pval2_discrete)
# This provides E(Rw_a^3),E(Rw_u^3), E(Rd_a^3) and E(Rd_u^3)
skewcorr = function(E, id) {
  n = length(id) # number of obs
  k = length(unique(id)) # number of distinct obs(= number of category)
  mk = as.numeric(table(id))
  t = 1:n
  
  p_hat = rep(0,n)
  p_hat = (t-1)/(n-2)
  q_hat = 1-p_hat
  
  p1 = t*(t-1)/n/(n-1)
  p2 = t*(t-1)*(t-2)/n/(n-1)/(n-2)
  p3 = t*(t-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
  p4 = t*(t-1)*(t-2)*(t-3)*(t-4)/n/(n-1)/(n-2)/(n-3)/(n-4)
  p5 = t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  p6 = t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  p7 = t*(t-1)*(t-2)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)
  p8 = t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = (n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)
  q3 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)/n/(n-1)/(n-2)/(n-3)
  q4 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/n/(n-1)/(n-2)/(n-3)/(n-4)
  q5 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  q7 = t*(t-1)*(n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)/(n-3)/(n-4)
  q8 = t*(t-1)*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  nE = nrow(E)
  
  Ebynode = vector("list", k)
  for(i in 1:k) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  
  temp28 = temp29 = temp30 = temp31 = temp41 = temp45 = 0
  nodedeg_a = nodedeg_u = rep(0,k)
  quan = quan1 = quan2 = quan3 = quan4 = quan5 = quan_u = 0
  for (i in 1:nE){
    e1 = E[i,1]
    e2 = E[i,2]
    nodedeg_a[e1] = nodedeg_a[e1] + 1
    nodedeg_a[e2] = nodedeg_a[e2] + 1
    nodedeg_u[e1] = nodedeg_u[e1] + mk[e2]
    nodedeg_u[e2] = nodedeg_u[e2] + mk[e1]
    quan = quan + 1/mk[e1]/mk[e2]
    quan_u = quan_u + mk[e1]*mk[e2]
    quan1 = quan1 + (mk[e1]+mk[e2])/mk[e1]^2/mk[e2]^2
    quan2 = quan2 + (mk[e1]+mk[e2])^2/mk[e1]/mk[e2]
    quan3 = quan3 + (mk[e1]^2+mk[e2]^2)/mk[e1]/mk[e2]
    quan4 = quan4 + 1/mk[e1]/mk[e1]/mk[e2]/mk[e2]
    quan5 = quan5 + (mk[e1]^2+mk[e2]^2)/mk[e1]/mk[e1]/mk[e2]/mk[e2]
    
    e3 = intersect(Ebynode[[e1]], Ebynode[[e2]])
    e4 = length(which(!is.na(match(Ebynode[[e1]], Ebynode[[e2]]))))
    temp28 = temp28 + e4
    temp29 = temp29 + sum( 1/mk[e3]/mk[e1]/mk[e2] )
    temp30 = temp30 + sum( (mk[e3]+mk[e1]+mk[e2])/mk[e3]/mk[e1]/mk[e2] )
    temp31 = temp31 + sum( (mk[e3]*mk[e1]+mk[e3]*mk[e2]+mk[e1]*mk[e2])/mk[e3]/mk[e1]/mk[e2] )
    temp41 = temp41 + (1/mk[e1]+1/mk[e2])*e4
    temp45 = temp45 + (mk[e1]+mk[e2]-1)*e4/mk[e1]/mk[e2]
  }
  
  G = sum(mk*(mk-1))/2 + quan_u
  nodedeg_u = nodedeg_u + mk - 1
  
  temp47 = sum(nodedeg_u*(nodedeg_u-1)*mk)
  temp49 = sum(nodedeg_u*(nodedeg_u-1)*(G-nodedeg_u)*mk)
  
  temp34 = temp35 = temp36 = temp37 = temp42 = 0
  temp51 = temp52 = 0
  for (i in 1:nE) {
    e1 = E[i,1]
    e2 = E[i,2]
    e3 = intersect(Ebynode[[e1]], Ebynode[[e2]])
    e4 = length(which(!is.na(match(Ebynode[[e1]], Ebynode[[e2]]))))
    temp34 = (nodedeg_a[e1]-1)*(nodedeg_a[e2]-1)
    temp35 = temp35 + 6*(temp34-e4)/mk[e1]/mk[e2]
    temp36 = temp36 + 6*(temp34-e4)*((mk[e1]+mk[e2]-2)/mk[e1]/mk[e2])
    temp37 = temp37 + 6*(temp34-e4)*((mk[e1]*mk[e2]-mk[e1]-mk[e2]+1)/mk[e1]/mk[e2])
    temp42 = temp42 + temp34
    temp51 = temp51 + sum(mk[e1]*mk[e2]*mk[e3]) + 3*mk[e1]*(mk[e1]-1)*mk[e2]/2 + 3*mk[e2]*(mk[e2]-1)*mk[e1]/2
    temp52 = temp52 + mk[e1]*mk[e2]*(nodedeg_u[e1]-1)*(nodedeg_u[e2]-1)
  }
  temp51 = temp51 + 3*sum(mk*(mk-1)*(mk-2)/6)
  temp52 = temp52 + sum(mk*(mk-1)*(nodedeg_u-1)^2/2)
  
  temp1 = sum(nodedeg_a/mk)
  temp2 = sum(nodedeg_a/(mk^2))
  temp3 = sum(1/mk)
  temp4 = sum(mk*nodedeg_a)
  temp5 = sum(mk^2)
  temp6 = sum(1/(mk^2))
  temp7 = sum(nodedeg_a^2)
  temp8 = sum(nodedeg_a^2/mk)
  temp9 = sum(nodedeg_a^2*mk)
  temp10 = sum(nodedeg_a^2*mk^2)
  temp12 = sum(nodedeg_a^2/mk/mk)
  temp13 = sum(nodedeg_a^3)
  temp15 = sum(nodedeg_a^3/mk)
  temp33 = sum(nodedeg_a^3/mk/mk)
  temp43 = sum(nodedeg_a*(nodedeg_a-1)*(3*nE-2*nodedeg_a-2))
  
  temp16 = temp19 = temp20 = temp23 = temp38 = rep()
  for (i in 1:k) {
    temp16 = c(temp16, sum(mk[Ebynode[[i]]]))
    temp19 = c(temp19, sum((nodedeg_a[Ebynode[[i]]]-1)*(nodedeg_a[Ebynode[[i]]]-2)))
    temp20 = c(temp20, sum(nodedeg_a[setdiff(c(1:k)[-i], Ebynode[[i]])]*(nodedeg_a[setdiff(c(1:k)[-i], Ebynode[[i]])]-1)))
    temp23 = c(temp23, sum(1/mk[Ebynode[[i]]]))
    temp38 = c(temp38, sum(nodedeg_a[Ebynode[[i]]]))
  }
  temp17 = sum( (nodedeg_a-1)*(1-1/mk)*((n-mk)*nodedeg_a-2*temp16) ) 
  temp18 = sum( (nodedeg_a-1)*((n-mk)*nodedeg_a-2*temp16)/mk ) 
  temp21 = sum( (mk-1)*((nE-nodedeg_a)*(nE-nodedeg_a-1)-temp19-temp20) )
  temp25 = sum( (nodedeg_a-1)*temp23 )
  temp26 = sum( (nodedeg_a-1)*temp23/mk )
  temp27 = sum( (nodedeg_a-1)*temp23/mk^2 )
  temp39 = sum( (nodedeg_a-1)*temp38 )
  temp40 = sum( (nodedeg_a-1)*temp38/mk )
  
  temp22 = temp32 = temp46 = temp44 = w = w1 = w2 = w3 = 0
  for (i in 1:k) {
    for (j in Ebynode[[i]]) {
      temp22 = temp22 + (nodedeg_a[j]-1)*((4*p3-10*p4+6*p5)/mk[j]+(2*p4-2*p5)*mk[i]/mk[j]-(4*p3-8*p4+4*p5)/mk[i]/mk[j])
      temp32 = temp32 + (nE-nodedeg_a[i]-nodedeg_a[j]+1)*(p5+(p3-2*p4+p5)/mk[i]/mk[j]+(p4-p5)*(1/mk[i]+1/mk[j]))
      temp44 = temp44 + (nodedeg_a[j]-1)*((18*p8-10*p7)/mk[j]+mk[i]*(2*p7-6*p8)/mk[j]+(8*p7-12*p8)/mk[i]/mk[j])
      temp46 = temp46 + (nE-nodedeg_a[i]-nodedeg_a[j]+1)*(3*p8+(p6-2*p7+3*p8)/mk[i]/mk[j]+(p7-3*p8)*(1/mk[i]+1/mk[j]))
      w = w + (nodedeg_a[j]-1)*((4*q3-10*q4+6*q5)/mk[j]+(2*q4-2*q5)*mk[i]/mk[j]-(4*q3-8*q4+4*q5)/mk[i]/mk[j])
      w1 = w1 + (nE-nodedeg_a[i]-nodedeg_a[j]+1)*(q5+(q3-2*q4+q5)/mk[i]/mk[j]+(q4-q5)*(1/mk[i]+1/mk[j]))
      w2 = w2 + (nodedeg_a[j]-1)*((18*q8-10*q7)/mk[j]+mk[i]*(2*q7-6*q8)/mk[j]+(8*q7-12*q8)/mk[i]/mk[j])
      w3 = w3 + (nE-nodedeg_a[i]-nodedeg_a[j]+1)*(3*q8+(p6-2*q7+3*q8)/mk[i]/mk[j]+(q7-3*q8)*(1/mk[i]+1/mk[j]))
    }
  }
  
  R1_a3 = p1*(4*temp3-4*temp6+quan4) + p2*(32*k-96*temp3+64*temp6+12*temp1-12*temp2+24*quan-9*quan1-6*quan4+3*temp27+2*temp29) + p3*(6*(n-k)*k+32*n-216*k+(6*k+96)*nE-(6*n-6*k+6*nE-412)*temp3-228*temp6+temp33-12*temp12+83*temp2+12*temp8-132*temp1+(3*n-3*k-75)*quan+13*quan4+30*quan1+quan5-6*temp29+2*temp30+9*temp26-12*temp27+temp35) + p4*(12*(n-k)*(n-3*k)-72*(n-5*k)+24*nE^2+(36*n-60*k-216)*nE-3*temp33+3*temp7-141*temp2+3*(nE-k-5)*temp8+30*temp12-(15*nE+9*n-12*k-240)*temp1+24*(n-k+nE-24)*temp3+288*temp6-6*(n-k-13)*quan-12*quan4-33*quan1+3*quan3-3*quan2-3*quan5+6*temp29-4*temp30+2*temp31+3*temp18+3*temp25-18*temp26+15*temp27+temp36-6*temp40+3*temp41) + p5*((n-k)^3-3*(n-k)*(4*n-10*k)+40*n-176*k-42*nE^2+(3*(n-k)^2-33*n+57*k+95)*nE-2*temp13+3*(nE-k+6)*temp7+2*temp33+70*temp2-3*(nE-k-1)*temp8-18*temp12+(9*n-12*k+15*nE-120)*temp1-(18*n-18*k+18*nE-256)*temp3-120*temp6-3*temp9+(6*nE-3)*temp4+3*(n-k-9)*quan+4*quan4+12*quan1-3*quan3+2*quan5+3*quan2+3*temp17+3*temp21-3*temp25+9*temp26-6*temp27-2*temp29+2*temp30-2*temp31+temp37-6*temp39+6*temp40+6*temp28-3*temp41+6*temp42-temp43+nE*(nE-1)*(nE-2)) + 3*temp22 + 1.5*temp32 
  R1_a2R2_a = p6*(2*(n-k)*k-8*k+2*k*nE-2*(nE+n-k-10)*temp3-12*temp6-4*temp1+4*temp2+(n-k-3)*quan+quan1+quan4+temp26-temp27) + p7*(8*nE^2+(12*n-20*k-72)*nE+4*(n-k)*(n-3*k)-24*n+120*k+8*(nE+n-k-24)*temp3+96*temp6+temp7-temp33-(5*nE+3*n-4*k-80)*temp1-47*temp2+(nE-k-5)*temp8+10*temp12-(2*n-2*k-26)*quan+quan3-11*quan1-quan2-quan5-4*quan4+2*temp29-temp30+temp45+temp18+temp25-6*temp26+5*temp27+temp36/3-2*temp40+temp41) + p8*((n-k)^3-3*(n-k)*(4*n-10*k)+40*n-176*k-42*nE^2+(3*(n-k)^2-33*n+57*k+95)*nE+2*temp33-2*temp13-3*temp9+(6*nE-3)*temp4-(18*n-18*k+18*nE-256)*temp3-120*temp6+3*(nE-k+6)*temp7+(15*nE+9*n-12*k-120)*temp1+70*temp2+(3*k-3*nE+3)*temp8-18*temp12+3*(n-k-9)*quan+12*quan1-3*quan3+3*quan2+2*quan5+4*quan4+3*temp17+3*temp21-3*temp25+9*temp26-6*temp27-2*temp29+2*temp30-2*temp31+temp37-6*temp39+6*temp40+6*temp28-3*temp41+nE*(nE-1)*(nE-2)+6*temp42-temp43) + temp44 + 0.5*temp46
  R1_aR2_a2 = p6*(2*(n-k)*k-8*k+2*k*nE-2*(nE+n-k-10)*temp3-12*temp6-4*temp1+4*temp2+(n-k-3)*quan+quan1+quan4+temp26-temp27) + q7*(8*nE^2+(12*n-20*k-72)*nE+4*(n-k)*(n-3*k)-24*n+120*k+8*(nE+n-k-24)*temp3+96*temp6+temp7-temp33-(5*nE+3*n-4*k-80)*temp1-47*temp2+(nE-k-5)*temp8+10*temp12-(2*n-2*k-26)*quan+quan3-11*quan1-quan2-quan5-4*quan4+2*temp29-temp30+temp45+temp18+temp25-6*temp26+5*temp27+temp36/3-2*temp40+temp41) + q8*((n-k)^3-3*(n-k)*(4*n-10*k)+40*n-176*k-42*nE^2+(3*(n-k)^2-33*n+57*k+95)*nE+2*temp33-2*temp13-3*temp9+(6*nE-3)*temp4-(18*n-18*k+18*nE-256)*temp3-120*temp6+3*(nE-k+6)*temp7+(15*nE+9*n-12*k-120)*temp1+70*temp2+(3*k-3*nE+3)*temp8-18*temp12+3*(n-k-9)*quan+12*quan1-3*quan3+3*quan2+2*quan5+4*quan4+3*temp17+3*temp21-3*temp25+9*temp26-6*temp27-2*temp29+2*temp30-2*temp31+temp37-6*temp39+6*temp40+6*temp28-3*temp41+nE*(nE-1)*(nE-2)+6*temp42-temp43) + w2 + 0.5*w3
  R2_a3 = q1*(4*temp3-4*temp6+quan4) + q2*(32*k-96*temp3+64*temp6+12*temp1-12*temp2+24*quan-9*quan1-6*quan4+3*temp27+2*temp29) + q3*(6*(n-k)*k+32*n-216*k+(6*k+96)*nE-(6*n-6*k+6*nE-412)*temp3-228*temp6+temp33-12*temp12+83*temp2+12*temp8-132*temp1+(3*n-3*k-75)*quan+13*quan4+30*quan1+quan5-6*temp29+2*temp30+9*temp26-12*temp27+temp35) + q4*(12*(n-k)*(n-3*k)-72*(n-5*k)+24*nE^2+(36*n-60*k-216)*nE-3*temp33+3*temp7-141*temp2+3*(nE-k-5)*temp8+30*temp12-(15*nE+9*n-12*k-240)*temp1+24*(n-k+nE-24)*temp3+288*temp6-6*(n-k-13)*quan-12*quan4-33*quan1+3*quan3-3*quan2-3*quan5+6*temp29-4*temp30+2*temp31+3*temp18+3*temp25-18*temp26+15*temp27+temp36-6*temp40+3*temp41) + q5*((n-k)^3-3*(n-k)*(4*n-10*k)+40*n-176*k-42*nE^2+(3*(n-k)^2-33*n+57*k+95)*nE-2*temp13+3*(nE-k+6)*temp7+2*temp33+70*temp2-3*(nE-k-1)*temp8-18*temp12+(9*n-12*k+15*nE-120)*temp1-(18*n-18*k+18*nE-256)*temp3-120*temp6-3*temp9+(6*nE-3)*temp4+3*(n-k-9)*quan+4*quan4+12*quan1-3*quan3+2*quan5+3*quan2+3*temp17+3*temp21-3*temp25+9*temp26-6*temp27-2*temp29+2*temp30-2*temp31+temp37-6*temp39+6*temp40+6*temp28-3*temp41+6*temp42-temp43+nE*(nE-1)*(nE-2)) + 3*w + 1.5*w1 
  R1_u3 = p1*G + p2*(3*temp47+2*temp51) + p3*(3*G*(G-1)-3*temp47+sum(nodedeg_u*(nodedeg_u-1)*(nodedeg_u-2)*mk)+6*temp52-6*temp51) + p4*(3*temp49+6*temp51-12*temp52) + p5*(G*(G-1)*(G-2)+6*temp52-2*temp51-sum(nodedeg_u*(nodedeg_u-1)*(3*G-2*nodedeg_u-2)*mk))
  R1_u2R2_u = p6*(G*(G-1)-temp47) + p7*(temp49+2*temp51-4*temp52) + p8*(G*(G-1)*(G-2)+6*temp52-2*temp51-sum(nodedeg_u*(nodedeg_u-1)*(nodedeg_u-2)*mk)-3*temp49)
  R1_uR2_u2 = p6*(G*(G-1)-temp47) + q7*(temp49+2*temp51-4*temp52) + q8*(G*(G-1)*(G-2)+6*temp52-2*temp51-sum(nodedeg_u*(nodedeg_u-1)*(nodedeg_u-2)*mk)-3*temp49)
  R2_u3 = q1*G + q2*(3*temp47+2*temp51) + q3*(3*G*(G-1)-3*temp47+sum(nodedeg_u*(nodedeg_u-1)*(nodedeg_u-2)*mk)+6*temp52-6*temp51) + q4*(3*temp49+6*temp51-12*temp52) + q5*(G*(G-1)*(G-2)+6*temp52-2*temp51-sum(nodedeg_u*(nodedeg_u-1)*(3*G-2*nodedeg_u-2)*mk))
  
  p_hat = rep(0,n)
  p_hat = (t-1)/(n-2)
  q_hat = 1-p_hat
  
  ER3w_a = q_hat^3*R1_a3 + 3*q_hat^2*p_hat*R1_a2R2_a + 3*q_hat*p_hat^2*R1_aR2_a2 + p_hat^3*R2_a3 # E(Rw_a^3)
  ER3w_u = q_hat^3*R1_u3 + 3*q_hat^2*p_hat*R1_u2R2_u + 3*q_hat*p_hat^2*R1_uR2_u2 + p_hat^3*R2_u3 # E(Rw_u^3)
  
  q_hat = 1
  p_hat = -1
  
  ER3d_a = q_hat^3*R1_a3 + 3*q_hat^2*p_hat*R1_a2R2_a + 3*q_hat*p_hat^2*R1_aR2_a2 + p_hat^3*R2_a3 # E(Rd_a^3)
  ER3d_u = q_hat^3*R1_u3 + 3*q_hat^2*p_hat*R1_u2R2_u + 3*q_hat*p_hat^2*R1_uR2_u2 + p_hat^3*R2_u3 # E(Rd_u^3)
  
  q_hat = 1
  p_hat = 1
  
  ER3_a = q_hat^3*R1_a3 + 3*q_hat^2*p_hat*R1_a2R2_a + 3*q_hat*p_hat^2*R1_aR2_a2 + p_hat^3*R2_a3 
  ER3_u = q_hat^3*R1_u3 + 3*q_hat^2*p_hat*R1_u2R2_u + 3*q_hat*p_hat^2*R1_uR2_u2 + p_hat^3*R2_u3 
  
  
  return( list(ER3w_a = ER3w_a, ER3w_u = ER3w_u, ER3d_a = ER3d_a, ER3d_u = ER3d_u, ER3_a = ER3_a, ER3_u = ER3_u) )
}



# Support function(skewness correction) for approximated p-value in single change point setting
# approximated p-value with extrapolation for max-count statistic
pval1_discrete_sub_1 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
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
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):n]*b>0))
  
  if (nn.l>0.35*n){
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
    neg = which(1+2*r[ceiling(n/2):n]*b<=0)
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:n] = ratio[id2-1]+inc*((id2:n)-id2)
    ratio[ratio<0]=0
    a[(n/2):n] = (x*Nu(sqrt(2*b^2*x))*ratio)[(n/2):n] # update a after extrapolation 
  }
 
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]
  }
  
  result = try(2*dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}

# approximated p-value(skewness correction) with extrapolation for weighted statistic
pval1_discrete_sub_2 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b) # S(t) in integrand
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-length(pos)
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




# p value approximation for changed interval (analytical p-value)
pval2_discrete = function(n, E, id, scanZ, statistics="all", skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)) {
  output = list()
  t = 1:n
  t = as.numeric(t)
  k = length(unique(id))
  nE = nrow(E)
  mk = as.numeric(table(id))
  
  nodedeg_a = nodedeg_u = rep(0,k)
  quan = quan_u = 0
  for (i in 1:nE) {
    e1 = E[i,1]
    e2 = E[i,2]
    nodedeg_a[e1] = nodedeg_a[e1] + 1
    nodedeg_a[e2] = nodedeg_a[e2] + 1
    nodedeg_u[e1] = nodedeg_u[e1] + mk[e2]
    nodedeg_u[e2] = nodedeg_u[e2] + mk[e1]
    quan = quan + 1/mk[e1]/mk[e2]
    quan_u = quan_u + mk[e1]*mk[e2]
  }
  temp3 = sum(1/mk)
  temp1 = sum(nodedeg_a/mk)
  temp8 = sum(nodedeg_a^2/mk)
  
  G = sum(mk*(mk-1))/2 + quan_u
  nodedeg_u = nodedeg_u + mk - 1
  
  one_a = 2*k - 2*temp3 + quan
  two_a = 4*(n-2*k+2*nE+temp8/4-temp1+temp3)
  three_a = (n-k+nE)^2
  one_u = G
  two_u = sum(mk*nodedeg_u^2)
  three_u = G^2
  
  # approximated p-value without skewness correction
  if (skew.corr==FALSE) {
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      b_a = scanZ$ori$Zo_a_max
      b_u = scanZ$ori$Zo_u_max
      integrando_a = function(t){
        c1 = rho_one_discrete(n,t,one_a,two_a,three_a)
        (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
      }
      integrando_u = function(t){
        c1 = rho_one_discrete(n,t,one_u,two_u,three_u)
        (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
      }
      pval.ori_a = b_a^3*dnorm(b_a)*integrate(integrando_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.ori_u = b_u^3*dnorm(b_u)*integrate(integrando_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      output$ori_a = min(pval.ori_a,1)
      output$ori_u = min(pval.ori_u,1)
    }
    if (length(which(!is.na(match(c("w","weighted","m","max","g","generalized","all"),statistics))))>0){
      if(lower<2){
        lower = 2
      }
      if(upper>(n-2)){
        upper = n-2
      }
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b_a = scanZ$weighted$Zw_a_max
      b_u = scanZ$weighted$Zw_u_max
      integrandW_a = function(t){
        c1 = rho_w(n, t)
        (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
      }
      integrandW_u = function(t){
        c1 = rho_w(n, t)
        (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
      }
      pval.weighted_a = b_a^3*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.weighted_u = b_u^3*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      output$weighted_a = min(pval.weighted_a,1)
      output$weighted_u = min(pval.weighted_u,1)
    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      b_a = scanZ$max.type$M_a_max
      b_u = scanZ$max.type$M_u_max
      integrandD_a = function(t){
        c1 = rho_d(n, t)
        (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
      }
      integrandD_u = function(t){
        c1 = rho_d(n, t)
        (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
      }
      integrandW_a = function(t){
        c1 = rho_w(n, t)
        (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
      }
      integrandW_u = function(t){
        c1 = rho_w(n, t)
        (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
      }
      pval_a1 = 2*b_a^3*dnorm(b_a)*integrate(integrandD_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_a2 = b_a^3*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.max.type_a = as.numeric(1-(1-min(pval_a1,1))*(1-min(pval_a2,1)))
      pval_u1 = 2*b_u^3*dnorm(b_u)*integrate(integrandD_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval_u2 = b_u^3*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
      pval.max.type_u = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
      output$max.type_a = pval.max.type_a
      output$max.type_u = pval.max.type_u
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      b_a = scanZ$generalized$S_a_max
      b_u = scanZ$generalized$S_u_max
      integrandG_a = function(t,w){
        x1 = rho_d(n,t)
        x2 = rho_w(n,t)
        (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2)*b_a*Nu(sqrt(2*b_a*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
      }
      integrand0_a = function(t) {integrate(integrandG_a,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      integrandG_u = function(t,w){
        x1 = rho_d(n,t)
        x2 = rho_w(n,t)
        (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2)*b_u*Nu(sqrt(2*b_u*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
      }
      integrand0_u = function(t) {integrate(integrandG_u,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
      pval.generalized_a = dchisq(b_a,2)*integrate(Vectorize(integrand0_a),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      pval.generalized_u = dchisq(b_u,2)*integrate(Vectorize(integrand0_u),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
      output$generalized_a = min(pval.generalized_a,1)
      output$generalized_u = min(pval.generalized_u,1)
    }
    
    return(output)
  }
  
  ##
  # approximated p-value with skewness correction
  if (length(which(!is.na(match(c("o","ori","original","w","weighted","m","max","all"), statistics))))>0) {
    
    temp = skewcorr(E, id) # give statistics for third moment of E(Rw_a^3) etc
    ER3w_a = temp$ER3w_a
    ER3w_u = temp$ER3w_u
    ER3d_a = temp$ER3d_a
    ER3d_u = temp$ER3d_u
    ER3_a = temp$ER3_a
    ER3_u = temp$ER3_u
    
    temp = getMV_discrete1(E, id)
    muo_a = temp$muo_a
    muo_u = temp$muo_u
    varo_a = temp$varo_a
    varo_u = temp$varo_u
    mu1_a = temp$mu1_a
    mu2_a = temp$mu2_a
    var1_a = temp$var1_a
    var2_a = temp$var2_a
    var12_a = temp$var12_a
    mu1_u = temp$mu1_u
    mu2_u = temp$mu2_u
    var1_u = temp$var1_u
    var2_u = temp$var2_u
    var12_u = temp$var12_u
    
    p_hat = rep(0,n)
    p_hat = (t-1)/(n-2)
    q_hat = 1-p_hat
    
    # mu and variance of the weighted edge-count test (average method)
    muw_a = q_hat*mu1_a + p_hat*mu2_a 
    varw_a = q_hat^2*var1_a + p_hat^2*var2_a + 2*q_hat*p_hat*var12_a
    muw_u = q_hat*mu1_u + p_hat*mu2_u 
    varw_u = q_hat^2*var1_u + p_hat^2*var2_u + 2*q_hat*p_hat*var12_u
    # mu and variance of the difference of two with-in group edge-counts (average method)
    mud_a = mu1_a - mu2_a 
    vard_a = var1_a + var2_a - 2*var12_a
    mud_u = mu1_u - mu2_u
    vard_u = var1_u + var2_u - 2*var12_u
    
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0) {
      A = n-k+nE
      ER3o_a = A^3-3*A^2*muo_a+3*A*(varo_a+muo_a^2) - ER3_a
      ER3o_u = G^3-3*G^2*muo_u+3*G*(varo_u+muo_u^2) - ER3_u
      
      ro_a = (-ER3o_a + 3*muo_a*varo_a + muo_a^3)/(varo_a^(3/2)) # E(Zo_a^3)
      for(i in 1:length(ro_a)){
        if (is.na(ro_a[i])==TRUE){
          ro_a[i]=0
        }
      }
      ro_u = (-ER3o_u + 3*muo_u*varo_u + muo_u^3)/(varo_u^(3/2)) # E(Zo_u^3)
      for(i in 1:length(ro_u)){
        if (is.na(ro_u[i])==TRUE){
          ro_u[i]=0
        }
      }
      c1_o_a = rho_one_discrete(n, t, one_a, two_a, three_a)
      for(i in 1:length(c1_o_a)){
        if ((abs(c1_o_a[i]))=="Inf"){
          c1_o_a[i]=0
        }
      }
      c1_o_u = rho_one_discrete(n, t, one_u, two_u, three_u)
      for(i in 1:length(c1_o_u)){
        if ((abs(c1_o_u[i]))=="Inf"){
          c1_o_u[i]=0
        }
      }
      
      b_a = scanZ$ori$Zo_a_max
      b_u = scanZ$ori$Zo_u_max
      result_o_a = pval2_discrete_sub_2(n, b_a, ro_a, c1_o_a, lower, upper)
      result_o_u = pval2_discrete_sub_2(n, b_u, ro_u, c1_o_u, lower, upper)
      # average
      if (is.numeric(result_o_a) && result_o_a > 0){
        output$ori_a = min(result_o_a, 1)
      }else{
        if (result_o_a ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Original edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
        b_a = scanZ$ori$Zo_a_max
        integrando_a = function(t){
          c1 = rho_one_discrete(n,t,one_a,two_a,three_a)
          (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
        }
        pval.ori_a = b_a^3*dnorm(b_a)*integrate(integrando_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori_a = min(pval.ori_a,1)
      }
      # union
      if (is.numeric(result_o_u) && result_o_u > 0){
        output$ori_u = min(result_o_u, 1)
      }else{
        if (result_o_u ==0){
          cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
        }
        cat("Original edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
        b_u = scanZ$ori$Zo_u_max
        integrando_u = function(t){
          c1 = rho_one_discrete(n,t,one_u,two_u,three_u)
          (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
        }
        pval.ori_u = b_u^3*dnorm(b_u)*integrate(integrando_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
        output$ori_u = min(pval.ori_u,1)
      }
    }
    
    if (length(which(!is.na(match(c("w","weighted","m","max","all"), statistics))))>0) {
      
      rw_a =  (ER3w_a- 3*muw_a*varw_a - muw_a^3)/(varw_a^(3/2)) # E(Zw_a^3)
      for(i in 1:length(rw_a)){
        if (is.na(rw_a[i])==TRUE){
          rw_a[i]=0
        }
      }
      rw_u =  (ER3w_u- 3*muw_u*varw_u - muw_u^3)/(varw_u^(3/2)) # E(Zw_u^3)
      for(i in 1:length(rw_u)){
        if (is.na(rw_u[i])==TRUE){
          rw_u[i]=0
        }
      }
      
      c1_w = rho_w(n, t)
      for(i in 1:length(c1_w)){
        if ((abs(c1_w[i]))=="Inf"){
          c1_w[i]=0
        }
      }
      
      if(lower<2){
        lower = 2
        #print(lower)
      }
      if(upper>(n-2)){
        upper = n-2
        #print(upper)
      }
      
      # weighted case
      if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
        b_a = scanZ$weighted$Zw_a_max
        b_u = scanZ$weighted$Zw_u_max
        result_w_a = pval2_discrete_sub_2(n, b_a, rw_a, c1_w, lower, upper)
        result_w_u = pval2_discrete_sub_2(n, b_u, rw_u, c1_w, lower, upper)
        
        # average
        if (is.numeric(result_w_a) && result_w_a > 0){
          output$weighted_a = min(result_w_a, 1)
        }else{
          if (result_w_a ==0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Weighted edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
          b_a = scanZ$weighted$Zw_a_max
          integrandW_a = function(t){
            c1 = rho_w(n, t)
            (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
          }
          pval.weighted_a = b_a^3*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$weighted_a = min(pval.weighted_a,1)
        }
        # union
        if (is.numeric(result_w_u) && result_w_u > 0){
          output$weighted_u = min(result_w_u, 1)
        }else{
          if (result_w_u ==0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Weighted edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
          b_u = scanZ$weighted$Zw_u_max
          integrandW_u = function(t){
            c1 = rho_w(n, t)
            (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
          }
          pval.weighted_u = b_u^3*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          output$weighted_u = min(pval.weighted_u,1)
        }
      }
      
      # max-type case
      if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
        c1_d = rho_d(n, t)
        for(i in 1:length(c1_d)){
          if ((abs(c1_d[i]))=="Inf"){
            c1_d[i]=0
          }
        }
        
        rd_a =  (ER3d_a - 3*mud_a*vard_a - mud_a^3)/(vard_a^(3/2))
        for(i in 1:length(rd_a)){
          if (is.na(rd_a[i])==TRUE){
            rd_a[i]=0
          }
        }
        if (rd_a[n/2]==0) {rd_a[n/2]=rd_a[n/2+1]}
        
        rd_u =  (ER3d_u - 3*mud_u*vard_u - mud_u^3)/(vard_u^(3/2))
        for(i in 1:length(rd_u)){
          if (is.na(rd_u[i])==TRUE){
            rd_u[i]=0
          }
        }
        if (rd_u[n/2]==0) {rd_u[n/2]=rd_u[n/2+1]}
        
        b_a = scanZ$max.type$M_a_max
        b_u = scanZ$max.type$M_u_max
        result_d_a = pval2_discrete_sub_1(n, b_a, rd_a, c1_d, lower, upper)
        result_w_a = pval2_discrete_sub_2(n, b_a, rw_a, c1_w, lower, upper)
        result_d_u = pval2_discrete_sub_1(n, b_u, rd_u, c1_d, lower, upper)
        result_w_u = pval2_discrete_sub_2(n, b_u, rw_u, c1_w, lower, upper)
        
        # average
        if (!is.numeric(result_d_a) || !is.numeric(result_w_a) || result_d_a ==0 || result_w_a ==0){
          if(result_d_a ==0 || result_w_a == 0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Max-type edge-count statistic (a): p-value approximation without skewness correction is reported.\n")
          b_a = scanZ$max.type$M_a_max
          integrandD_a = function(t){
            c1 = rho_d(n, t)
            (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
          }
          integrandW_a = function(t){
            c1 = rho_w(n, t)
            (c1*Nu(b_a*sqrt(2*c1)))^2*(n-t)
          }
          pval_a1 = 2*b_a^3*dnorm(b_a)*integrate(integrandD_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_a2 = b_a^3*dnorm(b_a)*integrate(integrandW_a, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type_a = as.numeric(1-(1-min(pval_a1,1))*(1-min(pval_a2,1)))
          output$max.type_a = pval.max.type_a
        }else{
          output$max.type_a = 1-(1-min(result_d_a,1))*(1-min(result_w_a,1))
        }
        # union
        if (!is.numeric(result_d_u) || !is.numeric(result_w_u) || result_d_u ==0 || result_w_u ==0){
          if(result_d_u ==0 || result_w_u == 0){
            cat("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
          }
          cat("Max-type edge-count statistic (u): p-value approximation without skewness correction is reported.\n")
          b_u = scanZ$max.type$M_u_max
          integrandD_u = function(t){
            c1 = rho_d(n, t)
            (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
          }
          integrandW_u = function(t){
            c1 = rho_w(n, t)
            (c1*Nu(b_u*sqrt(2*c1)))^2*(n-t)
          }
          pval_u1 = 2*b_u^3*dnorm(b_u)*integrate(integrandD_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval_u2 = b_u^3*dnorm(b_u)*integrate(integrandW_u, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
          pval.max.type_u = as.numeric(1-(1-min(pval_u1,1))*(1-min(pval_u2,1)))
          output$max.type_u = pval.max.type_u
        }else{
          output$max.type_u = 1-(1-min(result_d_u,1))*(1-min(result_w_u,1))
        }
      }
    }
  }
    
    
  # for generalized edge-count test, the approximated p-value without skewness correction is reported
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0) {
    b_a = scanZ$generalized$S_a_max
    b_u = scanZ$generalized$S_u_max
    integrandG_a = function(t,w){
      x1 = rho_d(n,t)
      x2 = rho_w(n,t)
      (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2)*b_a*Nu(sqrt(2*b_a*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
    }
    integrand0_a = function(t) {integrate(integrandG_a,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    integrandG_u = function(t,w){
      x1 = rho_d(n,t)
      x2 = rho_w(n,t)
      (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2)*b_u*Nu(sqrt(2*b_u*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
    }
    integrand0_u = function(t) {integrate(integrandG_u,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    pval.generalized_a = dchisq(b_a,2)*integrate(Vectorize(integrand0_a),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    pval.generalized_u = dchisq(b_u,2)*integrate(Vectorize(integrand0_u),lower,upper,subdivisions=3000, stop.on.error=FALSE)$value
    output$generalized_a = min(pval.generalized_a,1)
    output$generalized_u = min(pval.generalized_u,1)
  }
  
  return(output)  
}

# Support function(skewness correction) for approximated p-value in changed interval setting
# approximated p-value with extrapolation for max-count statistic
pval2_discrete_sub_1 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2*ratio  ## this is different from single change-point alternative
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):n]*b>0))
  
  if (nn.l>0.35*n){
    return(0)
  }
  
  neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
  if (nn.l>=lower){      ## this part also differs from single change-point
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  neg = which(1+2*r[ceiling(n/2):n]*b<=0) ## this part is different from gSeg
  if (nn.r>=(n-upper)){        ## this part also differs from single change-point
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:n] = ratio[id2-1]+inc*((id2:n)-id2)
    ratio[ratio<0]=0
    a[(n/2):n] = ((b^2*x*Nu(sqrt(2*b^2*x)))^2*ratio)[(n/2):n] ########### differs from gSeg
  }
  
  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]*(n-s)       # different from single change point
  }
  
  result = try(2*dnorm(b)/b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T) # differs from single
  return(result)
}

# approximated p-value(skewness correction) with extrapolation for weighted statistic
pval2_discrete_sub_2 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b) # S(t) in integrand
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  neg = which(1+2*r*b<=0)  ## this is different from single change-point
  if (nn>=(lower-1)+(n-upper)){
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
    a[s]*(n-s)    # different from single change point
  }
  
  result = try(dnorm(b)/b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T) ## this is different from single 
  return(result)
}





