
#' load and clean DAISY iggdata
#' 
#' This function loads the serum anti-carbohydrate
#' antibody data from the DAISY cohort along
#' with the DAISY phenotype data
#' It keeps only the earliest measured sample
#' for each individual and performs background
#' thresholding
#' 
#' @param dataloc Path to DAISY Iggdata
#' @param phenoloc Path to DAISY phenotype data
#' @param igg_cols a vector containing names of the molecules or glycans
#' @param sub_id_col name of the column where subject ID can be located
#' @param removeoutliers Logical indicating whether or not to remove outliers in the iggdata
#' @param backgroundThesholdval Value for background thresholding of Iggdata
#' @param log2 Logical indicating whether or not to log2 transform the iggdata
#' 
#' @return a list containing two data frames, igg and pheno
#' @export
load_and_clean_DAISY_data<-function(dataloc,
                              phenoloc,
                              igg_cols,
                              mols_to_remove,
                              sub_id_col,
                              removeoutliers=F, 
                              backgroundThesholdval=1,
                              log2=F){
  
  
  #load glycan data
  iggdata<-read.csv(dataloc)
  
  #remove GLYRDC_078
  iggdata<-iggdata[,-which(colnames(iggdata)==mols_to_remove)]
  
  #keep earliest sample
  iggdata<-keep_earliest_samples(iggdata)
  rownames(iggdata)<-iggdata[,sub_id_col]
  
  #keep only igg data
  iggdata<-iggdata[,igg_cols]
  
  
  #load updated pheno data
  if (!is.null(phenoloc) && phenoloc != "") {
    phenodata <- xlsx::read.xlsx(phenoloc, 1)[-c(1:5), ]
    rownames(phenodata) <- phenodata[, sub_id_col]
  }
  
  #remove outliers, 7734 outliers
  if (removeoutliers==T){
  iggdata<-outliers_to_na(iggdata)
  }
  #remove glycans with more than 70 na's 
  #GLYPW_057, GLYPW_079, GLYPW_095, GLYPW_107, GLYRDC_022, GLYRDC_057,
  #GLYRDC_068, GLYRDC_080, GLYRDC_081, GLYRDC_082, GLYRDC_098, GLYRDC_099, 
  #GLYRDC_0742
  #45 gly
  #gly_torm<-names(which(apply(iggdata,2,function(x)sum(is.na(x)))>30))
  #iggdata<-iggdata[,-which(colnames(iggdata)%in%gly_torm)]
  
  #correct all negative values to 0.1
  if (!is.null(backgroundThesholdval) && backgroundThesholdval != "") {
  iggdata[iggdata<backgroundThesholdval]<-backgroundThesholdval+0.1
  }
  
  # #merge with updated phenotype data
  newdatafile<-merge.data.frame(iggdata,phenodata,by="row.names")
  rownames(newdatafile)<-newdatafile$Row.names
  newdatafile$Age.of.the.first<-as.numeric(newdatafile$Age.of.the.first)
  colnames(newdatafile)
  cols_torm<-c("Row.names","Subject.ID","age.at.last","For.T1D.case.","For.loose.case..age.on.event.visit.date..For.non.cases..age.on.last.visit.date..4AB.")
  newdatafile<-newdatafile[,!colnames(newdatafile)%in%cols_torm]
  newdatafile[,c((dim(iggdata)[1]+1):dim(newdatafile)[1])]
  
  daisy=list('igg'=newdatafile[,c(1:dim(iggdata)[2])],'pheno'=newdatafile[,c((dim(iggdata)[2]+1):dim(newdatafile)[2])])
  
  colnames(daisy$pheno)<-c("Group","Sex","FDR","Draw_Age","HLA_risk")
  
  if (log2==T){
    daisy$igg=log2(daisy$igg)
  }
  daisy$pheno$Group<-factor(daisy$pheno$Group,
                            levels=c("Progressor",
                                     "Non-progressor",
                                     "Control"))
  
  daisy
}

#' load and clean PAGODA iggdata
#' 
#' This function loads the serum anti-carbohydrate
#' antibody data from the PAGODA cohort along
#' with the PAGODA phenotype data
#' It keeps only the earliest measured sample
#' for each individual and performs background
#' thresholding
#' 
#' @param fileloc Path to PAGODA data
#' @param removeoutliers Logical indicating whether or not to remove outliers in the iggdata
#' @param backgroundThesholdval Value for background thresholding of Iggdata
#' @param log2 Logical indicating whether or not to log2 transform the iggdata
#' 
#' @return a list containing two data frames, igg and pheno
#' @export
load_and_clean_PAGODA_data<-function(fileloc,
                                     removeoutliers = F, 
                                     backgroundThesholdval = 1,
                                     log2=F){
  iggdata<-read.csv(fileloc,row.names = "Person_ID",na.strings = "")
  pagoda<-list("igg"=iggdata[,which(colnames(iggdata)%in%colnames(daisy$igg))],
               "pheno"=iggdata[,c("Group","Sex","FDR_T1D","Dr_Age","Genetic_Risk")])
  rownames(pagoda$igg)<-pagoda$pheno$Person_ID
  colnames(pagoda$pheno)<-c("Group","Sex","FDR","Draw_Age","HLA_risk")
  #remove outliers
  if (removeoutliers==T){
    pagoda$igg<-outliers_to_na(pagoda$igg)
  }
  
  #correct all negative values to 0.1
  pagoda$igg[pagoda$igg<backgroundThesholdval]<-backgroundThesholdval+0.1
  
  if (log2==T){
    pagoda$igg=log2(pagoda$igg)
  }
  
  
  pagoda
  
}


#' Load glycan classes information
#' 
#' This function loads the classification information for glycans
#' 
#' @param fileloc Path to location of glycan class data file
#' @return data frame for classification of glycans
#' @export
load_glycan_classes<-function(fileloc="/Users/paultran/Downloads/glycan_classes.csv"){
  glycan_classes<-read.csv(fileloc,na.strings = "")
  glycan_classes<-glycan_classes[-which(duplicated(glycan_classes$Glycan.ID)),]
  glycan_classes<-glycan_classes[-which(is.na(glycan_classes$Glycan.ID)),]
  rownames(glycan_classes)<-glycan_classes$Glycan.ID
  glycan_classes
}

#' Keep earliest sample only for multiple samples
#' 
#' 
#' @param iggdata data frame with repeat samples
#' @return data frame  with only earliest samples
#' @export
keep_earliest_samples<-function(iggdata){
  subjects<-as.character(unique(iggdata$Subject.ID))
  iggdata$Date.of.Sample.Collection<-as.Date(iggdata$Date.of.Sample.Collection,format="%m/%d/%Y")
  keep<-c()
  
  for (i in 1:length(subjects)){
    samps<-which(iggdata$Subject.ID==subjects[i])
    dd<-cbind.data.frame(dates=iggdata$Date.of.Sample.Collection[samps],samps)[order(iggdata$Date.of.Sample.Collection[samps]),]
    index<-dd$samps[1]
    keep<-c(keep,index)
    rm(dd,i,index,samps)
  }
  
  iggdata<-iggdata[keep,]
  iggdata
}
