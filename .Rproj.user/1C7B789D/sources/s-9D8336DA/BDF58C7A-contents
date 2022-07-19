#' IApre1
#'
#' IApre1 is used for survival prediction of stage IA lung adenocarcinoma based on a z-scale of a gene matrix
#' @param x a data.frame containing GeneID and gene expression of one or more patients
#' @example IApre1(x)
#'
#'


IApre1<-function(x){
  IAdata1<-IAdata1[order(IAdata1$GeneID),]
  IAdata2<-IAdata2[order(IAdata2$GeneID),]
  x<-subset(x,x$GeneID %in% IAdata1$GeneID,select=c(1:ncol(x)))
  x<-x[order(x$GeneID),]
  Patients<-colnames(x)[-1]
  rownames(x)<-x$GeneID
  x<-as.data.frame(t(scale(t(x[,-1]))))
  x<-as.data.frame(sapply(x,function(y)sum(y*IAdata1$coef)))
  x$Patients<-Patients
  colnames(x)<-c("LAS","Patients")
  fiveyears<-sapply(x[,1],function(y)(1/(1+exp(-(-0.6271 + (-2.2610*y))))))
  threeyears<-sapply(x[,1],function(y)(1/(1+exp(-(-1.5455 + (-1.1097*y))))))
  oneyear<-sapply(x[,1],function(y)(1/(1+exp(-(-3.4167 + (-0.7781*y))))))
  x<-data.frame(x[,2],fiveyears,threeyears,oneyear)
  colnames(x)<-c("Patients","Fiveyears","threeyears","oneyear")
  x}

