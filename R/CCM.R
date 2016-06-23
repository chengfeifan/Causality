#' The original Convergent Cross Mapping
#' @param x A vector represent the time series
#' @param y A vector represent the time series
#' @return A list of SugiX, SugiY, OriginX, OriginY, CorX, CorY
#' @examples x<-c(1:100)
#' y<-c(100:1)
#' result<-CCM(x,y,tau=1,E=2)
#' @details This method is proposed by the Sugihara, which is based on the Taken's Theorem.
#' @export
CCM<-function(x,y,tau,E,LMN,...){
  if(nargs()<4){
    stop('bad input')
  }
  else if(nargs()==4){
    LMN=E+1
  }
  else {

  }

  L=length(x)
  Ti=1+(E-1)*tau
  Xm=matrix(0,nrow = L-Ti+1,ncol=E)
  Ym=matrix(0,nrow = L-Ti+1,ncol=E)
  SugiN=E+1
  N=L-Ti+1

  for(t in 1:N){
    Xm[t,]<-x[seq((Ti+t-1),(Ti+t-1-(E-1)*tau),-tau)]
    Ym[t,]<-y[seq((Ti+t-1),(Ti+t-1-(E-1)*tau),-tau)]
  }

  LMj<-array(0,dim=c(2,2,N))

  library(parallel)
  cl<-makeCluster(getOption("cl.cores",8))
  clusterEvalQ(cl,library(FNN))
  clusterExport(cl,c('SugiX','SugiY'))
  dat<-floor((L-Ti+1)/2)
  j=(dat+1):(L-Ti+1)
  Sugi<-parLapply(cl,j,function(ii,Xm,Ym){
    Xknn<-get.knnx(Xm[(ii-dat):(ii-1),],t(as.matrix(Xm[ii,])),algorithm="kd_tree",SugiN)
    Yknn<-get.knnx(Ym[(ii-dat):(ii-1),],t(as.matrix(Ym[ii,])),algorithm="kd_tree",SugiN)

    uls=exp(Xknn$nn.dist/Xknn$nn.dist[,1])
    wls=uls/sum(uls)
    SugiY<<-wls%*%y[(Xknn$nn.index+Ti-1+ii-(dat+1))]

    uls=exp(Yknn$nn.dist/Yknn$nn.dist[,1])
    wls=uls/sum(uls)
    SugiX<<-wls%*%x[(Yknn$nn.index+Ti-1+ii-(dat+1))]

    return(list(SugiX,SugiY))
  },Xm,Ym)
  Sugi<-unlist(Sugi)
  SugiX<-Sugi[seq(1,length(Sugi),2)]
  SugiY<-Sugi[seq(2,length(Sugi),2)]

  OriginY<-y[Ti:(length(y))]
  OriginY<-OriginY[(dat+1):N]
  OriginX<-x[Ti:(length(x))]
  OriginX<-OriginX[(dat+1):N]

  CorX<-cor(SugiX,OriginX)
  CorY<-cor(SugiY,OriginY)

  return(list("SugiX"=SugiX,"SugiY"=SugiY,"OriginX"=OriginX,"OriginY"=OriginY,
               "CorX"=CorX,"CorY"=CorY))
}
