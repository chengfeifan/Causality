#' Cross Correlation Function
#' @param x A vector represents a time sequence
#' @param y A vector represents another time sequence
#' @param lag.max A numeric determines the max time lag
#'
#' @return $lag The sequence of time lag
#' @return $CCF The value of cross correlation function
#' @details The method of cross correlation function can be used to detect
#' the causality in the linear system.
#' @examples
#' set.seed(42)
#' x<-runif(100)
#' y<-runif(100)
#' result<-CCF(x,y,lag.max=50)
#' plot(result)
#' @export


CCF<-function(x,y,lag.max,...){
  oneTagCCF<-function(x,y,pho){
    if(pho>=0){
      varLen<-length(x)
      meanX<-mean(x)
      meanY<-mean(y)
      num<-seq(1,varLen-pho)
      numerator<-sum(unlist(lapply(num,function(num,x,y){
        (x[num]-meanX)*(y[num+pho]-meanY)
      },x,y)))
      denominator1<-sum(unlist(lapply(num,function(num,x){
        (x[num]-meanX)^2
      },x)))
      denominator2<-sum(unlist(lapply(num,function(num,y){
        (y[num+pho]-meanY)^2
      },y)))
      denominator<-sqrt(denominator1*denominator2)
      cPho<-numerator/denominator*(varLen-pho)/varLen
      return(cPho)
    }
    else{
      temp<-x
      x<-y
      y<-temp
      pho=-pho
      varLen<-length(x)
      meanX<-mean(x)
      meanY<-mean(y)
      num<-seq(1,varLen-pho)
      numerator<-sum(unlist(lapply(num,function(num,x,y){
        (x[num]-meanX)*(y[num+pho]-meanY)
      },x,y)))
      denominator1<-sum(unlist(lapply(num,function(num,x){
        (x[num]-meanX)^2
      },x)))
      denominator2<-sum(unlist(lapply(num,function(num,y){
        (y[num+pho]-meanY)^2
      },y)))
      denominator<-sqrt(denominator1*denominator2)
      cPho<-numerator/denominator*(varLen-pho)/varLen
      return(cPho)
    }
  }

  sequence<-seq(-lag.max,lag.max,1)
  values<-unlist(lapply(sequence,function(s,x,y){
    oneTagCCF(x,y,s)
  },x,y))
  result<-list("lag"=sequence,"CCF"=values)
  class(result)<-"CCF"
  return(result)
}

plot.CCF<-function(x,...){
  plot(x$lag,x$CCF,type='l',lwd=2,xlab='lag',ylab='ccf',lty=1
       ,main = 'CCF' )
}

print.CCF<-function(x,...){
  result<-data.frame(x$lag,y$CCF)
  colnames(result)<-c("lag","CCF")
  print(result)
}
