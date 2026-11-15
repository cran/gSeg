# Compatibility boundary for the native MST replacement of ade4::mstree.
# Unusual conversion methods retain the original R sequence and neig class.
.gseg_mst_conversion_ok <- function(d) {
  tryCatch({
    if (!.gseg_dist_conversion_ok() ||
        !identical(attr(d,"class",exact=TRUE),"dist") ||
        !(typeof(d) %in% c("integer","double"))) return(FALSE)
    n<-attr(d,"Size",exact=TRUE);labels<-attr(d,"Labels",exact=TRUE)
    if(is.object(n)||length(n)!=1L||!(typeof(n)%in%c("integer","double")))return(FALSE)
    if(!is.null(labels)&&(is.object(labels)||!is.character(labels)||length(labels)!=n))return(FALSE)
    b<-asNamespace("base");st<-asNamespace("stats")
    if(bindingIsActive("as.double",b)||!identical(get("as.double",b),.Primitive("as.double")))return(FALSE)
    for(name in c("t","t.default"))if(bindingIsActive(name,b))return(FALSE)
    tg<-get("t",b);td<-get("t.default",b)
    gt<-function(x)UseMethod("t");dt<-".Internal(t.default(x))"
    if(!identical(environment(tg),b)||!identical(environment(td),b)||
       !identical(body(tg),body(gt))||!identical(formals(tg),formals(gt))||
       !identical(paste(deparse(body(td)),collapse="\n"),dt)||!identical(formals(td),formals(gt))||
       isdebugged(tg)||isdebugged(td))return(FALSE)
    scopes<-list(.GlobalEnv,st)
    extra<-c(paste0("t.",c("matrix","array","double","integer","numeric","neig")),"[.neig","as.matrix.neig")
    for(scope in scopes){
      for(name in extra)if(exists(name,scope,inherits=TRUE))return(FALSE)
      env<-scope
      repeat{
        if(exists("t.default",env,inherits=FALSE)){
          if(bindingIsActive("t.default",env)||!identical(get("t.default",env,inherits=FALSE),td))return(FALSE)
          break
        }
        if(identical(env,emptyenv()))return(FALSE)
        env<-parent.env(env)
      }
    }
    table<-get(".__S3MethodsTable__.",b)
    for(name in extra)if(exists(name,table,inherits=FALSE))return(FALSE)
    if(exists("t.default",table,inherits=FALSE)&&
       (bindingIsActive("t.default",table)||!identical(get("t.default",table,inherits=FALSE),td)))return(FALSE)
    if(bindingIsActive("as.matrix.dist",st))return(FALSE)
    method<-get("as.matrix.dist",st)
    current<-quote({
      size<-attr(x,"Size");df<-matrix(0,size,size)
      if(size>1L){
        n..1<-(size-1L):1L;s.1<-size+1L
        up<-sequence.default(n..1,from=seq.int(s.1,length(df),s.1),by=size)
        lo<-sequence.default(n..1,from=seq.int(2L,length(df)+1L,s.1))
        df[up]<-df[lo]<-x
      }
      labels<-attr(x,"Labels") %||% seq_len(size)
      dimnames(df)<-list(labels,labels);df
    })
    older<-quote({
      size<-attr(x,"Size");df<-matrix(0,size,size)
      lower<-row(df)>col(df);df[lower]<-x;df<-t(df);df[lower]<-x
      labels<-attr(x,"Labels")
      dimnames(df)<-if(is.null(labels))list(seq_len(size),seq_len(size))else list(labels,labels)
      df
    })
    newer<-quote({
      size<-attr(x,"Size");df<-matrix(0,size,size)
      if(size>1L){
        n..1<-(size-1L):1L;s.1<-size+1L;s2<-size*size
        up<-sequence.default(n..1,seq.int(s.1,s2,s.1),by=size)
        lo<-sequence.default(n..1,seq.int(2L,s2,s.1))
        df[up]<-df[lo]<-x
      }
      labels<-attr(x,"Labels") %||% seq_len(size)
      dimnames(df)<-list(labels,labels);df
    })
    if(!identical(environment(method),st)||isdebugged(method)||
       !identical(formals(method),formals(function(x,...)NULL))||
       !(identical(body(method),current)||identical(body(method),older)||identical(body(method),newer)))return(FALSE)
    for(scope in scopes){
      env<-scope
      repeat{
        if(exists("as.matrix.dist",env,inherits=FALSE)){
          if(bindingIsActive("as.matrix.dist",env)||!identical(get("as.matrix.dist",env,inherits=FALSE),method))return(FALSE)
          break
        }
        if(identical(env,emptyenv()))break
        env<-parent.env(env)
      }
    }
    if(!exists("as.matrix.dist",table,inherits=FALSE)||bindingIsActive("as.matrix.dist",table)||
       !identical(get("as.matrix.dist",table,inherits=FALSE),method))return(FALSE)
    TRUE
  },error=function(e)FALSE)
}
.gseg_mst_legacy_conversion <- function(d,k) {
  if(!inherits(d,"dist"))stop("Object of class 'dist' expected")
  D<-as.matrix(d)
  nlig<-nrow(D)
  values<-as.double(D)
  if(k<=1)k<-1
  if(k>=nlig)k<-1
  k<-as.integer(k)
  mat01<-.Call(C_gseg_mst_matrix,values,nlig,k)
  # GPL>=2: ade4 neig(mat01=...) / neig.util.GtoL, retaining conversion order.
  G<-as.matrix(mat01);n<-nrow(G)
  if(ncol(G)!=n)stop("Square matrix expected")
  if(any(t(G)!=G))stop("Symetric matrix expected")
  if(sum(G==0|G==1)!=n*n)stop("0-1 values expected")
  if(sum(diag(G)!=0))stop("Null diagonal expected")
  G<-G*(row(G)<col(G))
  G<-(row(G)+0+(0+(0+1i))*col(G))*G
  G<-as.vector(G);G<-G[G!=0];output<-cbind(Re(G),Im(G))
  degrees<-as.vector(apply(mat01,1,sum))
  if(!is.null(rownames(mat01)))names(degrees)<-rownames(mat01)
  else if(!is.null(colnames(mat01)))names(degrees)<-colnames(mat01)
  else names(degrees)<-as.character(1:nrow(mat01))
  attr(output,"degrees")<-degrees
  attr(output,"call")<-quote(neig(mat01=mst))
  class(output)<-"neig"
  output
}
.gseg_reference_adjacency <- function(E,n) {
  out<-vector("list",n)
  for(i in 1:n)out[[i]]<-rep(0,0)
  for(i in 1:nrow(E)){
    out[[E[i,1]]]<-c(out[[E[i,1]]],E[i,2])
    out[[E[i,2]]]<-c(out[[E[i,2]]],E[i,1])
  }
  out
}
