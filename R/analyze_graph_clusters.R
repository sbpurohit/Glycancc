
#' Construct igraph based on pairwise spearman correlations
#' 
#' This function is not recommended since a K-nearest neighbors approach
#' with louvain clustering is a more generalizable approach to clustering.
#' See the ??make_umap and ??knn_network_from_umap functions for more details. 
#' This function calculates the pairwise spearman correlaton
#' for all anti-carbohydrate antibodies, then constructs a graph based
#' on the spearman's rho and prunes weak edges based on p-value and 
#' spearmna's rho cutoffs.
#' 
#' @param iggdata data frame of iggdata
#' @param pcutoff p-value cutoff for edge removal from graph. Default = 0.05
#' @param rcutoff Spearman's rho cutoff for edge removal from graph. Default = 0.5
#' @param n_neighbor keep edges based on number of nearest neighbors. Default = 2
#' 
#' @return igraph object
#' @export
make_correlation_network<-function(iggdata,
                                   pcutoff=0.05,
                                   rcutoff=0.5,
                                   n_neighbor=2){
  my_cor_data<-Hmisc::rcorr(as.matrix(as.data.frame(iggdata)),type = "spearman")
  
  diag(my_cor_data$r)<-0
  
  maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
  k_ind<-apply(my_cor_data$r,2,maxn(n_neighbor))#x[which(x!=max(x))]<0)
  
  for (i in 1:dim(daisy$igg[,aminoglycosides])[2]){
    #which(my_cor_data$P[,i]<0.05)
    my_cor_data$r[which(my_cor_data$P[,i]>=pcutoff),i]<-0
    my_cor_data$r[which(is.na(my_cor_data$r[,i])),i]<-0
    my_cor_data$r[which(my_cor_data$r[,i]<my_cor_data$r[k_ind[i],i]),i]<-0
    
    my_cor_data$r[which(my_cor_data$r[,i]<rcutoff),i]<-0
    
    # my_cor_data$r[which(my_cor_data$r[,i]<0),i]<-0
    
    
    
  }
  
  #generate undirected graph
  mygraph<-igraph::graph_from_adjacency_matrix(my_cor_data$r,mode="undirected",weighted = T)
  mygraph
}


#' glycan class enrichment analysis
#' 
#' This function calculates the significance of enrichment of glycans for
#' each glycan class in each glycan cluster using fisher's exact test
#' 
#' @param co_cluster data frame of iggdata
#' @param glygroups p-value cutoff for edge removal from graph. Default = 0.05
#' 
#' @return list including a full table and partial table of significantly enriched classes
#' @export
glycan_class_enrichment<-function(co_cluster=co_cluster,
                                  glygroups=glycan_classes){
  
  results2<-data.frame("index"=1:length(unique(glycan_classes$Glycan.Class)))
  for (j in 1:max(co_cluster)){
    # j=1
    gly_clust<-names(which(co_cluster==j))
    group_dummyvars<-fastDummies::dummy_cols(glygroups$Glycan.Class)
    
    gly_clust_vec<-rep(0,dim(glygroups)[1])
    gly_clust_vec[which(glygroups$Glycan.ID%in%gly_clust)]=1
    results<-data.frame("Class","p-val")[-1,]
    for (i in 2:dim(group_dummyvars)[2]){
      # i=2
      tab<-table(gly_clust_vec,group_dummyvars[,i])
      myvec<-c(paste(colnames(group_dummyvars)[i]),
               fisher.test(tab,alternative="greater")$p.val)
      results<-rbind.data.frame(results,myvec)
    }
    
    results2<-cbind.data.frame(results2,results)
  }
  results2<-results2[,c(2,seq(3,dim(results2)[2],2))]
  colnames(results2)<-c("Glycan.Class",paste0("Cluster",1:max(co_cluster)))
  glyclass<-results2$Glycan.Class
  glyclass<-unlist(lapply(strsplit(glyclass,"_"),'[[',2))
  
  results2<-results2[,-1]
  results2<-apply(results2,2,function(x)as.numeric(as.vector(x)))
  
  rows_to_keep<-which(apply(results2,1,function(x)min(x)<0.05))
  table3<-results2[rows_to_keep,]
  
  rownames(results2)<-glyclass
  rownames(table3)<-glyclass[rows_to_keep]
  
  list("fulltab"=results2,"smalltab"=table3)  
}

#' plot result of glycan enrichment analysis as bubble chart on ggplot2
#' 
#' This function plots the results from ??glycan_class_enrichment as a bubble chart
#' 
#' @param cluster_class_enrichment output from ??glycan_class_enrichment
#' 
#' @return bubble chart
#' @export
plot_significant_cluster_classes<-function(cluster_class_enrichment){
  
  sig_class_cluster<-cbind.data.frame("Glycan_Class"=rownames(cluster_class_enrichment$smalltab),
                                      tidyr::gather(data.frame(cluster_class_enrichment$smalltab),
                                                    "cluster","p.val",colnames(cluster_class_enrichment$smalltab)))
  sig_class_cluster$p.val<- -log10(sig_class_cluster$p.val)
  sig_class_cluster$cluster<- factor(sig_class_cluster$cluster,
                                     levels = paste0("Cluster",1:dim(cluster_class_enrichment$smalltab)[2]))
  sig_class_cluster$Glycan_Class<- factor(sig_class_cluster$Glycan_Class,
                                          levels = sort(unique(sig_class_cluster$Glycan_Class),decreasing = T))
  
  g<-ggplot(data=sig_class_cluster,
            aes(y=Glycan_Class,x=cluster,size=p.val,color=p.val))+
    geom_point()+ 
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90))
  g
}


#' regression models
#' 
#' This function performs linear regression for each column in a given data frame
#' except for those listed as a dependent variable
#' @param newdatafile data frame with all independent and dependent variables to be modeled
#' @param covariates a vector of all dependent variables. Example: groups<-c("Group", "Sex", "FDR", "HLA_risk","Draw_Age")
#' @param depVarList a vector of all y variables. Example: depVarList<-c("Cluster 1", "Cluster 2", "Cluster 3", ....) or c("CRP","SAA","sgp130","sIL6R" )
#' 
#' @return list object including abbreviated summary table of significant associations, summary table of all associations, and all models
#' @export
univariate_models <- function(newdatafile, 
                              covariates, 
                              depVarList,
                              group_col = "Group") {
  
  # If you pass custom variables, it uses them. Otherwise, 
  # it falls back to automatic detection using group_col:
  if (missing(depVarList) || is.null(depVarList)) {
    glys_torm <- suppressWarnings(names(which(is.na(apply(newdatafile, 2, sd, na.rm = TRUE)))))
    depVarList <- setdiff(colnames(newdatafile), group_col)
    depVarList <- setdiff(depVarList, intersect(depVarList, glys_torm))
  }
  
  # Apply over them and create model for each safely
  allModels <- lapply(depVarList, function(x) {
    
    # Use dynamic group_col instead of hardcoded "Group"
    form <- as.formula(paste0("`", x, "` ~ `", group_col, "` + ", paste(setdiff(covariates, group_col), collapse = " + ")))
    
    # Check if target 'x' or covariates contain Inf/NaN for this specific iteration
    vars_needed <- c(x, group_col, setdiff(covariates, group_col))
    sub_df <- newdatafile[, vars_needed, drop = FALSE]
    
    # Keep only rows where values are strictly finite (drops NA, NaN, Inf)
    valid_rows <- apply(sub_df, 1, function(row) {
      all(is.finite(suppressWarnings(as.numeric(as.factor(row))))) || all(is.finite(suppressWarnings(as.numeric(row))))
    })
    
    # Safely fit lm, returns NULL if an error occurs
    tryCatch({
      lm(formula = form, data = newdatafile, subset = valid_rows)
    }, error = function(e) {
      NULL 
    })
  })
  
  # Name the list and filter out any NULLs if models failed completely
  names(allModels) <- depVarList
  allModels <- purrr::compact(allModels) 
  
  univ_mod_sum <- lapply(allModels, broom::tidy)
  
  estimates <- lapply(univ_mod_sum, '[[', 2)
  pvals <- lapply(univ_mod_sum, '[[', 5)
  
  univ_mod_sum_tab <- cbind.data.frame(
    nonprog_est = unlist(lapply(estimates, '[[', 2)),
    nonprog_p   = unlist(lapply(pvals, '[[', 2)),
    prog_est    = unlist(lapply(estimates, '[[', 3)),
    prog_p      = unlist(lapply(pvals, '[[', 3))
  )
  
  list(small_table = univ_mod_sum_tab, fulltable = univ_mod_sum, models = allModels)
}

#convert model stats to ggplot format
convert_to_ggplotformat<-function(daisy_models=daisy_eigen_models$fulltable){
  #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/es-calc.html
  fulltab<-data.frame("effect_size","lowerCI","upperCI","Cluster", "Variable")[-1,]
  
  for (j in 2:dim(daisy_models$Cluster1)[1]){
    # rowid=3
    # j=2
    rowid=j
    resultstab<-data.frame(row.names = c(
      "effect_size",
      "lowerCI",
      "upperCI"))
    for (i in 1:length(daisy_models)){
      # i=11
      
      # Then, we use the weights to calculate the pooled effect
      pooled_effect <- daisy_models[[i]][rowid,2]
      lowerCI<-pooled_effect-1.96*daisy_models[[i]][rowid,3]
      upperCI<-pooled_effect+1.96*daisy_models[[i]][rowid,3]
      resultstab<-cbind.data.frame(resultstab,unlist(c(pooled_effect,
                                                       lowerCI,
                                                       upperCI)))
      
    }
    
    colnames(resultstab)<-paste0("Cluster",1:length(daisy_models))
    resultstab<-data.frame(t(resultstab))
    resultstab$Cluster<-rownames(resultstab)
    resultstab$Variable<-daisy_models$Cluster1$term[rowid]
    fulltab<-rbind.data.frame(fulltab,resultstab)
    
  }
  fulltab$Cluster<-factor(fulltab$Cluster,levels=paste0("Cluster",length(daisy_models):1))
  
  fulltab
}
#' convert model stats to ggplot format
#' 
#' This function converts the output from ??univariate_models to ggplot2 format
#' to prepare for using ??forest_plot_lm_model function
#' 
#' @param daisy_models output from ??univariate_models
#' 
#' @return table of linear model results, including effect size and 95% CI
#' @export
convert_to_ggplotformat<-function(daisy_models=daisy_eigen_models$fulltable){
  #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/es-calc.html
  fulltab<-data.frame("effect_size","lowerCI","upperCI","Cluster", "Variable")[-1,]
  
  for (j in 2:dim(daisy_models$Cluster1)[1]){
    # rowid=3
    # j=2
    rowid=j
    resultstab<-data.frame(row.names = c(
      "effect_size",
      "lowerCI",
      "upperCI"))
    for (i in 1:length(daisy_models)){
      # i=11
      
      # Then, we use the weights to calculate the pooled effect
      pooled_effect <- daisy_models[[i]][rowid,2]
      lowerCI<-pooled_effect-1.96*daisy_models[[i]][rowid,3]
      upperCI<-pooled_effect+1.96*daisy_models[[i]][rowid,3]
      resultstab<-cbind.data.frame(resultstab,unlist(c(pooled_effect,
                                                       lowerCI,
                                                       upperCI)))
      
    }
    
    colnames(resultstab)<-paste0("Cluster",1:length(daisy_models))
    resultstab<-data.frame(t(resultstab))
    resultstab$Cluster<-rownames(resultstab)
    resultstab$Variable<-daisy_models$Cluster1$term[rowid]
    fulltab<-rbind.data.frame(fulltab,resultstab)
    
  }
  fulltab$Cluster<-factor(fulltab$Cluster,levels=paste0("Cluster",length(daisy_models):1))
  
  fulltab
}


#' plot ggplot converted lm model as forest plot
#' 
#' This function plots the results from ??convert_to_ggplotformat as a forest plot
#' 
#' @param daisy_model output from ??convert_to_ggplotformat
#' 
#' @return forest plot as a ggplot2 object
#' @export
forest_plot_lm_model<-function(daisy_model){
  clrid<-rep('black',dim(daisy_model)[1])
  clrid[which(daisy_model$upperCI<0)] <- "red"
  clrid[which(daisy_model$lowerCI>0)] <- "red"
  
  g<-ggplot2::ggplot(data=daisy_model,
            ggplot2::aes(x=Cluster,y=effect_size))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=0))+
    ggplot2::geom_point(col=clrid)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lowerCI,ymax=upperCI),col=clrid,width=0.3)+
    ggplot2::facet_grid(cols=ggplot2::vars(Variable))+
    ggplot2::coord_flip()+
    ggplot2::theme_classic()
  g
}
