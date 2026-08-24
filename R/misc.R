
#' entropy/shannon diversity index function
#' 
#' This function calculates the entropy of shannon diversity index for a given vector
#' 
#' @param pop vector
#' 
#' @return entropy value
#' @export
entropy<-function(pop){
  n = pop
  N = sum(pop, na.rm = T)
  p = n/N
  H = -sum(p*log(p),na.rm = T)
  H
}


#' outliers to na
#' 
#' This function takes a data frame and converts all column-wise outliers to na
#' An outlier is a value less than the 1st quartile minus the interquartile range
#' or a value greater than the 3rd quartile plus the interquartile range
#' 
#' @param iggdata numeric matrix
#' 
#' @return numeric matrix with outliers converted to "NA"
#' @export
outliers_to_na<-function(iggdata){
  sum_igg<-apply(iggdata,2,quantile,na.rm=T)
  igg_iqr<-sum_igg[4,]-sum_igg[2,]
  igg_lower<-sum_igg[2]-igg_iqr
  igg_higher<-sum_igg[4]+igg_iqr
  
  for (i in 1:dim(iggdata)[2]){
    iggdata[which(iggdata[,i]<igg_lower[i]),i]<-NA
    iggdata[which(iggdata[,i]>igg_higher[i]),i]<-NA
  }
  iggdata
}


#' rader/spider chart per glycan class
#' 
#' This function plots the significance of association for each anti-carbohydrate
#' antibody in a class against islet autoimmunity and progression to type 1 diabetes
#' accounting for covariates sex, first-degree relative, HLA risk, and age. This uses the 
#' ??univariate_models function
#' 
#' @param daisy igg data
#' @param gly_class glycan class to assess. Default = "Aminoglycoside antibiotics"
#' @param class_col glycan grouping columns for your proteins, genes, glycans eg class_col = "Glycan.Class" or "group"
#' @param id_col name of molecules stored in  a separate file . eg  default is id_col = "Glycan.ID"
#' @param group_col grouping variable. Default is group_col = "Group" values are (c("Control", "Non-progressor", "Progressor")
#' @param group_levels order of levels in your grouping column default is group_levels = c("Control", "Non-progressor", "Progressor")
#' @param scale values are TRUE or FALSE, if you supplied scaled data then change it to FALSE. default is scale=TRUE
#' @param legend_labels comming from your grouping column default is group_levels = c("Control", "Non-progressor", "Progressor")

#' 
#' @return spider plot
#' @export
univariate_radar_chart <- function(daisy = daisy, 
                                   covariates = c("Group", "Sex",  "Draw_Age"),
                                   gly_class = "Aminoglycoside antibiotics", 
                                   label = TRUE,
                                   class_col = "Glycan.Class", 
                                   id_col = "Glycan.ID",
                                   glygroups = glycan_classes,
                                   group_col = "Group",
                                   group_levels = c("Control", "Non-progressor", "Progressor"),
                                  scale=TRUE,
                                  legend_labels = c("Control", "Non-progressor", "Progressor")) {
  
  # Dynamically extract classes and IDs using custom column names and trim whitespace
  classes_vec <- trimws(as.character(glygroups[[class_col]]))
  ids_vec <- trimws(as.character(glygroups[[id_col]]))
  
  # Filter IDs matching the requested gly_class
  selected_ids <- ids_vec[classes_vec == trimws(gly_class)]
  selected_ids <- selected_ids[!is.na(selected_ids)]
  selected_ids <- intersect(selected_ids, colnames(daisy$igg))
  
  # Safety check: Stop and print available classes if none match
  if (length(selected_ids) == 0) {
    stop(paste("No matching molecules found for group:", gly_class, 
               "\nAvailable groups are:", paste(unique(classes_vec), collapse = ", ")))
  }
  
  if (scale == TRUE) {
    daisy$combined <- cbind.data.frame(
      apply(daisy$igg[, selected_ids, drop = FALSE], 2, scale), 
      daisy$pheno
    )
  } else {
    daisy$combined <- cbind.data.frame(
      daisy$igg[, selected_ids, drop = FALSE], 
      daisy$pheno
    )
  }
  
  # 4. Create/overwrite the Group column using your specific phenotype column (e.g., ptype)
  # Replace 'ptype' with whatever column name exists in daisy$pheno if it changes
  daisy$combined$Group <- factor(
    daisy$pheno[[group_col]], 
    levels = group_levels
  )
  
  daisy_aminoglyc_models <- univariate_models(
    daisy$combined,
    covariates, group_col = group_col,
    depVarList=selected_ids
      )
  
  radar_df <- rbind.data.frame(
    1.5, 0,
    -log10(daisy_aminoglyc_models$small_table$prog_p),
    -log10(0.05),
    -log10(daisy_aminoglyc_models$small_table$nonprog_p)
  )
  
  colnames(radar_df) <- rownames(daisy_aminoglyc_models$small_table)
  # 2. Dynamically calculate the maximum value in your data and add 10% padding
  max_val <- ceiling(max(radar_df, na.rm = TRUE) * 1.1)
  
  # Ensure the max ceiling is at least 1.5 just in case all p-values are very high/low
  max_val <- max(max_val, 1.5) 
  
  # 3. Prepend the Max (row 1) and Min (row 2) rows required by fmsb::radarchart
  radar_df <- rbind(
    rep(max_val, ncol(radar_df)), # Row 1: Dynamic Maximum Ceiling
    rep(0, ncol(radar_df)),       # Row 2: Minimum (0)
    radar_df                      # Rows 3+: Your actual data
  )
  
  areas <- c(rgb(0, 0, 1, 0.15), rgb(0, 0, 0, 0.1), rgb(1, 0, 0, 0.15))
  
  if (label == FALSE) {
    fmsb::radarchart(
      radar_df, pfcol = areas, pty = c(16, 32, 16),
      pcol = c("blue", "black", "red"), axistype = 0, # Changed axistype to 1 so you can see scale rings
      calcex = 1, vlabels = NA, axislabcol = "black",
      cglcol = "gray",
    )
  } else {
    fmsb::radarchart(
      radar_df, pfcol = areas, pty = c(16, 32, 16),
      pcol = c("blue", "black", "red"), axistype = 0, # Changed axistype to 1
      calcex = 1, axislabcol = "black",
      cglcol = "gray",
    )
  }
  # Add the Legend
  legend(
    x = "topright",                          # Position of the legend (e.g., "topright", "bottomleft", etc.)
    legend = legend_labels, # Labels matching your groups
    col = c("blue", "black", "red"),         # Line colors matching pcol
    lty = 1,                                 # Line type
    lwd = 2,                                 # Line width
    bty = "n"                                # Box type ("n" removes the border box)
  )
}


#' Calculate the AUC of binomial glmnet model from random number of glycans
#' 
#' @param daisy_prog igg data subseted two only include two factors for comparison. In this case, control subjects and progressors.
#' @param no_glycans Number of glycans to randomly sample
#' @return AUC value
#' @export
glmnet_auc<-function(daisy_prog=daisy_prog,no_glycans){
  subigg<-daisy_prog$igg[,sample(1:dim(daisy_prog$igg)[2],no_glycans)]
  x<-data.matrix(cbind(daisy_prog$pheno[,-1],subigg))
  cv_fit <- glmnet::cv.glmnet(y=as.vector(daisy_prog$pheno$Group),x=x , alpha = 0,family="binomial")
  mypred<-predict(cv_fit$glmnet.fit,newx=as.matrix(x),type = "response",s=cv_fit$lambda.min)
  
  suppressMessages(pROC::roc(fastDummies::dummy_cols(daisy_prog$pheno$Group)[,2],as.numeric(mypred))$auc)
}

