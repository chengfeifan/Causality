#' Generate the data structure for CCM
#'
#' Takes in two sequences, time lag, and row tag, E
#' @param x A vector represents the time sequence
#' @param y A vector represents another time sequence
#' @param lag Time lag between x and y
#' @param tag A vector, which selects the dots of x and y to calculate the CCM
#' @param E+1 is the embedded dimension of the manifold
#'
#' @return $dataX A matrix, which is manifold of X.
#' $dataY A vector, which is the original Y
#' @examples
#' x<-c(1:100)
#' y<-c(100:1)
#' data<-sampleCCM(x,y,l=2,t = c(1:50),E=2)
#' @export
sampleCCM<-function(x,y,lag,tag,E,...){
  xLen<-length(x)
  yLen<-length(y)
  if(xLen!=yLen){
    stop('The length does not match')
  }
  else{
    tag<-tag[which(tag+lag+E<xLen+1 & tag+lag>0)]
    dataRange<-unlist(lapply(tag,function(tagdot){
      xdot<-x[seq(tagdot+lag,tagdot+lag+E,1)]
    }))
    dataX<-matrix(dataRange,ncol=E+1,byrow = T)
    dataY<-unlist(lapply(tag,function(tagdot){
      ydot<-y[tagdot]
    }))
    return(list("dataX"=dataX,"dataY"=dataY))
  }
}


#' Generate the data structure for CCM
#'
#' Takes in two sequences, time lag, and row tag, E
#' @param x A vector represents the time sequence
#' @param y A vector represents another time sequence
#' @param lag Time lag between x and y
#' @param tag Select the dots of x and y to calculate the CCM
#' @param E+1 is the embedded dimension of the manifold
#' @examples
#' x<-c(1:100)
#' y<-c(100:1)
#' data<-sampleCCM(x,y,l=2,t = c(1:50),E=2)
#' @return $dataX A matrix, dataX is manifold of X.
#' $DataY is the manifold of y, dataY is a matrix
#' @export
multiSample<-function(x,y,lag,tag,E,...){
  xLen<-length(x)
  yLen<-length(y)
  if(xLen!=yLen){
    stop('The length does not match')
  }
  else{
    tag<-tag[which(tag+lag+E<xLen+1 & tag+E<xLen+1 & tag+lag>0)]
    dataRange<-unlist(lapply(tag,function(tagdot){
      xdot<-x[seq(tagdot+lag,tagdot+lag+E,1)]
    }))
    dataX<-matrix(dataRange,ncol=E+1,byrow = T)
    dataY<-unlist(lapply(tag,function(tagdot){
      ydot<-y[seq(tagdot,tagdot+E,1)]
    }))
    dataY<-matrix(dataY,ncol=E+1,byrow = T)
    return(list("dataX"=dataX,"dataY"=dataY))
  }
}
