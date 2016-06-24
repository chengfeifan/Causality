#' Run the Convergent Cross Mapping
#' @param x A vector represent the time series
#' @param y A vector represent the time series
#' @param E The embedding dimension
#' @param tau The time interval of sample
#' @param L Vector of computing points
#' @return $L Vector of computing points
#' @return $CorX Vector of correlation X
#' @return $CorY Vector of correlation Y
#' @details If x causes y, the correlation X will converge into a certain number.
#' @examples
#' set.seed(43)
#'x<-runif(100)
#'set.seed(50)
#'y<-runif(100)
#'L<-seq(50,100,10)
#'tau<-1
#'E=2
#'dataCCM<-runCCM(x,y,tau = tau,E=E,L=L)
#'plot(dataCCM)
#' @export
runCCM<-function(x, y, tau, E, L){
  Len<-length(L)
  CorX<-c(1:Len)
  CorY<-c(1:Len)
  for(i in 1:Len){
    CCMresult<-CCM(x[1:L[i]],y[1:L[i]],tau,E,CALL=match.call())
    CorX[i]<-CCMresult$CorX
    CorY[i]<-CCMresult$CorY
  }
  result<-list("L" = L,"CorX" = CorX,"CorY" = CorY)
  class(result)<-"CCM"
  return(result)
}

print.CCM<-function(x,...){
  CCMframe<-data.frame(x$L,x$CorX,x$CorY)
  colnames(CCMframe)<-c("L","CorX","CorY")
  print(CCMframe)
}

plot.CCM<-function(x,...){
  plotylimits<-range(c(x$CorX,x$CorY,na.rm=TRUE))
  plot(x$L,x$CorX,type='l',lwd=2,ylim=c(plotylimits[1],plotylimits[2])
       ,col=1,lty=1,xlab = 'Data Length',ylab=expression(rho))
  lines(x$L,x$CorY,type='l',lwd=2,col=2,lty=2)
  legend("topleft",c("X causes Y","Y causes X"),lty=c(1,2),col=c(1,2),
         lwd=2,bty = 'n')
}
