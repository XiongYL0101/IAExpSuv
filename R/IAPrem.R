#' IAprem
#'
#' IAprem is used for predicting the five-year risk of death for a group of stage IA lung adenocarcinoma patients
#' @param x a data.frame containing GeneID and gene expression of a group of patients
#'
#'
#'


IAprem<-function(x){
        x<-subset(x,x$GeneID %in% IAdata1$GeneID,select=c(1:ncol(x)))
        x<-x[order(x$GeneID),]
        Patients<-colnames(x)[-1]
        rownames(x)<-x$GeneID
        x<-as.data.frame(t(scale(t(x[,-1]))))
        x<-as.data.frame(sapply(x,function(y)sum(y*IAdata1$coef)))
        x$Patients<-Patients
        colnames(x)<-c("LAS","Patients")
        fiveyears<-sapply(x[,1],function(y)(1/(1+exp(-(-0.3628 + (-1.8313*y))))))
        x<-data.frame(x[,2],fiveyears)
        colnames(x)<-c("Patients","Risk Probability")
        x}
