#' IApres
#'
#' IApres is used for predicting the five-year risk of death for a single stage IA lung adenocarcinoma patient
#' @param x a data.frame containing GeneID and gene expression of a single patient
#'
#'
#'


IApres<-function(x){
             x<-x[order(x$GeneID),]
             GeneID1<-x$GeneID
             Patients<-colnames(x)[-1]
             x1<-as.data.frame(apply(x[,-1],2,function(y)(y/quantile(y,0.75))))
             x2<-as.data.frame(apply(x[,-1],2,function(y)(y-median(y))/mad(y)))
             x1$GeneID<- GeneID1
             x1<-x1[,c(ncol(x1),c(1:(ncol(x1)-1)))]
             x1<-x1[order(x1$GeneID),]
             x1<-subset(x1,x1$GeneID %in% IAdata2$GeneID,select=c(1:ncol(x1)))
             x1<-x1[order(x1$GeneID),]
             GeneID2<-x1$GeneID
             x1<-as.data.frame(sapply(x1[,-1],function(z)z<-as.data.frame(t(scale(t(cbind(z,IAdata2)))))[,1]))
             x1$GeneID<- GeneID2
             x1<-x1[,c(ncol(x1),c(1:(ncol(x1)-1)))]
             x1<-subset(x1,x1$GeneID %in% IAdata1$GeneID,select=c(1:ncol(x1)))
             x1<-x1[order(x1$GeneID),]
             rownames(x1)<-x1$GeneID
             x1<-as.data.frame(sapply(x1[,-1],function(y)sum(y*IAdata1$coef)))
             x1$Patients<-Patients
             colnames(x1)<-c("LAS","Patients")
             UQ<-sapply(x1[,1],function(y)(1/(1+exp(-(-0.3628 + (-1.8313*y))))))
             x1<-data.frame(x1[,2],UQ)
             colnames(x1)<-c("Patients","UQ")
             x2$GeneID<- GeneID1
             x2<-x2[,c(ncol(x2),c(1:(ncol(x2)-1)))]
             x2<-x2[order(x2$GeneID),]
             x2<-subset(x2,x2$GeneID %in% IAdata2$GeneID,select=c(1:ncol(x2)))
             x2<-x2[order(x2$GeneID),]
             GeneID2<-x2$GeneID
             x2<-as.data.frame(sapply(x2[,-1],function(z)z<-as.data.frame(t(scale(t(cbind(z,IAdata2)))))[,1]))
             x2$GeneID<- GeneID2
             x2<-x2[,c(ncol(x2),c(1:(ncol(x2)-1)))]
             x2<-subset(x2,x2$GeneID %in% IAdata1$GeneID,select=c(1:ncol(x2)))
             x2<-x2[order(x2$GeneID),]
             rownames(x2)<-x2$GeneID
             x2<-as.data.frame(sapply(x2[,-1],function(y)sum(y*IAdata1$coef)))
             x2$Patients<-Patients
             colnames(x2)<-c("LAS","Patients")
             MAD<-sapply(x2[,1],function(y)(1/(1+exp(-(-0.3628 + (-1.8313*y))))))
             x2<-data.frame(x2[,2],MAD)
             colnames(x2)<-c("Patients","MAD")
             P<-merge(x1,x2,by=c("Patients"))
             P$Survivial.rate.of.five.years<-((P$UQ+P$MAD)/2)
             S<-P[,c(1,4)]
             colnames(S)<-c("Patients","Risk Probability")
             S}
