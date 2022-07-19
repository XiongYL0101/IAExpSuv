#' IApre2
#'
#' IApre2 is used for survival prediction of stage IA lung adenocarcinoma based on a gene expression of a single patient
#' @param x a data.frame containing GeneID and gene expression of one or more patients
#' @example IApre2(x)
#'
#'


IApre2<-function(x){
  IAdata1<-IAdata1[order(IAdata1$GeneID),]
  IAdata2<-IAdata2[order(IAdata2$GeneID),]
  x<-x[order(x$GeneID),]
  GeneID1<-x$GeneID
  Patients<-colnames(x)[-1]
  x<-as.data.frame(apply(x[,-1],2,function(y)(y/quantile(y,0.75))))
  x$GeneID<- GeneID1
  x<-x[,c(ncol(x),c(1:(ncol(x)-1)))]
  x<-x[order(x$GeneID),]
  x<-subset(x,x$GeneID %in% IAdata2$GeneID,select=c(1:ncol(x)))
  x<-x[order(x$GeneID),]
  GeneID2<-x$GeneID
  x<-as.data.frame(sapply(x[,-1],function(z)z<-as.data.frame(t(scale(t(cbind(z,IAdata2)))))[,1]))
  x$GeneID<- GeneID2
  x<-x[,c(ncol(x),c(1:(ncol(x)-1)))]
  x<-subset(x,x$GeneID %in% IAdata1$GeneID,select=c(1:ncol(x)))
  x<-x[order(x$GeneID),]
  rownames(x)<-x$GeneID
  x<-as.data.frame(sapply(x[,-1],function(y)sum(y*IAdata1$coef)))
  x$Patients<-Patients
  colnames(x)<-c("LAS","Patients")
  fiveyears<-sapply(x[,1],function(y)(1/(1+exp(-(-0.8857 + (-0.1072*y))))))
  threeyears<-sapply(x[,1],function(y)(1/(1+exp(-(-1.6411 + (-0.0818*y))))))
  oneyear<-sapply(x[,1],function(y)(1/(1+exp(-(-3.3685 + (-0.0446*y))))))
  x<-data.frame(x[,2],fiveyears,threeyears,oneyear)
  colnames(x)<-c("Patients","Fiveyears","threeyears","oneyear")
  x}
