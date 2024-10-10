###############################################################################
#### Functions for the optimal location of the cut-off points based on GAM ####
###############################################################################

# x <- los datos de la variable que queremos categoriza
# y <- prediccion modelo GAM con terms
# w <- se.fit^(-2) (inversa de la varianza)
# pred <- prediccion de la funci?n a trozos.
# ymed la media de las y en cada nivel del factor

# Esta funci?n calcula el error asociado (comparaci?n GAM y funci?n a trozos) para cada punto de corte
CortesGAM=function(x,y,w=rep(1,length(y)),cortes,nmin=0.02*length(x)){
  
  cortes=unique(sort(c(min(x), cortes, max(x))))
  pred=numeric(length(y))
  
  
  x2=cut(x, include.lowest=TRUE,breaks=cortes)
  niveles=levels(x2); nn=length(niveles)
  minimo=min(table(x2)) # n?mero m?nimo de valores en cortes
  
  error=Inf
  if (minimo >= nmin) {
    ymed=tapply(y,x2,mean); 
    for (j in 1:nn) pred[x2==niveles[j]]=ymed[j]
    error=sum(w*(y-pred)**2)}
  list(x=x,y=y,w=w,cortes=cortes,pred=pred,err=error)
}

# Busqueda BackAddFor para los puntos de corte ?ptimos. 
# El modelo lo ajusto fuera de la funci?n.
# nc= n?mero de puntos de corte que quiero buscar.

Opt.cpoints.GAM <- function(x, y, w, nfino=100, nc, eps=0.001,repmax=5 ){
  
  cfino=seq(min(x),max(x),length=nfino) # grid donde va a buscar
  
  mat.aux <- NULL
  
  # Nodos Iniciales
  c0=seq(quantile(x,0.025),quantile(x,0.975),length=nc+2)[-c(1,nc+2)] # elimino el primero y el ?ltimo
  err0=CortesGAM(x = x , y = y , w=w , c=c0)$err
  
  # B?squeda Sequencial
  cfino=seq(min(x),max(x),length=nfino+2)[-c(1,nfino+2)]
  
  seguir=TRUE
  nrep=0
  err_old=err0
  
  while (seguir) {
    nrep=nrep+1 
    err_old=err0
    
    for ( j in 1:nc) { 
      for (ifino in 1:nfino) {
        c1=c0
        c1[j]=cfino[ifino]
        err1=CortesGAM(x = x , y = y , w=w , c=c1)$err
        mat.aux = rbind(mat.aux,c(c1,err1))
        if (err1<err0) {err0=err1; c0=c1}
      }
    }
    
    if (nrep>1 & abs((err0-err_old)/err_old)<eps) {seguir=FALSE}
    if (nrep>repmax) {seguir=FALSE} 
    err_old=err0  
  }
  
  list(res=CortesGAM(x,y,w,c0), errores = mat.aux)
}


###############################################################################
#### Function to simulate the theoretical logistic model                   ####
###############################################################################

# Docle búsqueda multivariante
SimLOG.multi=function(cortes1=c(0.2,0.5),cortes2=c(0.2,0.5),
                      fx1=c(0.5,1.5,2),fx2=c(0.5,1.5,2),
                      x1=runif(200), x2=runif(200), z=runif(200)) {
  n=length(x1); K1=length(cortes1) ; K2=length(cortes2) ; 
  pred1=numeric(n);pred2=numeric(n); pred=numeric(n);
  pred1[x1<=cortes1[1]]=fx1[1]
  pred1[x1>cortes1[K1]]=fx1[K1+1]
  pred2[x2<=cortes2[1]]=fx2[1]
  pred2[x2>cortes2[K2]]=fx2[K2+1]
  for (k in 1:(K1-1)) pred1[cortes1[k]<x1 & x1<=cortes1[k+1]]=fx1[k+1]
  for (j in 1:(K2-1)) pred2[cortes2[j]<x2 & x2<=cortes2[j+1]]=fx2[j+1]
  
  pred = pred1 + pred2 + 0.1*z
  
  pred = exp(pred)/(1+exp(pred))
  
  list(x1=x1,x2=x2, z=z,pred=pred)}