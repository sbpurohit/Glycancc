
#' run umap algorithm on igg data
#' 
#' 
#' @param datamat data matrix containing igg data
#' @param n_neighbors number of neighbors to use in umap hyperparameters. Default = 15
#' @param min_dist hyperparameter minimum distance for umap. Default = 0.1
#' @param metric hyperparameter distance metric to use for umap. Default = "cosine"
#' @param random_state umap variable for reproducible results. Default = 123
#' 
#' @return umap object including k-nearest neighbors results and 2D projection
#' @export
make_umap<-function(datamat,
                    n_neighbors=15,
                    min_dist=0.1,
                    metric="cosine",
                    random_state=123)
{
  custom.config<-umap::umap.defaults
  custom.config$n_neighbors=n_neighbors
  custom.config$min_dist=min_dist
  custom.config$metric=metric
  custom.config$random_state=random_state
  
  myumap<-umap::umap(datamat,config = custom.config)
  myumap
}

#knn from umap
#' Construct K-nearest neighbors graph from UMAP implementation of KNN
#' 
#' 
#' @param myumap output from make_umap function
#' @param metric hyperparameter distance metric used for umap function. Default = "cosine"
#' 
#' @return kNN graph as an igraph object
#' @export
knn_network_from_umap<-function(myumap,metric="cosine"){
  
  edge_list<-as.matrix(dplyr::left_join(
    cbind.data.frame("gly1"=rownames(myumap$knn$indexes),
                     "index1"=myumap$knn$indexes[,1],
                     "index2"=as.vector(myumap$knn$indexes)),
    cbind.data.frame("index2"=myumap$knn$indexes[,1],"gly2"=rownames(myumap$knn$indexes))
  )[,c(1,4)])
  
  #generate undirected graph object
  knn_graph<-igraph::graph_from_edgelist(edge_list,
                                 directed = F)
  if (metric=="euclidean"){
    #weights are equal to the inverse square distance
    igraph::E(knn_graph)$weight<-1/(as.vector(myumap$knn$distances)+0.1)^2
  }
  else igraph::E(knn_graph)$weight<-as.vector(myumap$knn$distances)
  #remove multiple edges and self loops
  knn_graph<-igraph::simplify(knn_graph,remove.multiple = T,remove.loops = T)
  knn_graph
}


#' combine two igraph objects
#' 
#' @param graph1 First igraph object
#' @param graph2 Second igraph object
#' @return igraph object
#' @export
join_graphs<-function(graph1=daisy_knn_graph,
                      graph2=pagoda_knn_graph){
  
  #extract edge lists
  edgelist_1<-cbind(data.frame(igraph::get.edgelist(graph1)),
                    igraph::edge_attr(daisy_knn_graph))
  edgelist_2<-cbind(data.frame(igraph::get.edgelist(graph2)),
                    igraph::edge_attr(pagoda_knn_graph))
  
  #merge edge lists
  merged_edgelist<-merge.data.frame(edgelist_1,edgelist_2,by=c("X1","X2"),all = T)
  
  #calculate mutual information of multiple edges
  merged_edgelist$sum<-rowSums(merged_edgelist[,3:4],na.rm = T)
  merged_edgelist$prod<-apply(merged_edgelist[,3:4],1,prod,na.rm=F)
  merged_edgelist$weight<-rowSums(cbind(merged_edgelist$sum,-merged_edgelist$prod),na.rm = T)
  
  #create merged graph
  merged_graph<-igraph::graph_from_edgelist(as.matrix(merged_edgelist[,c(1,2)]),directed = F)
  igraph::E(merged_graph)$weight <- merged_edgelist$weight
  
  merged_graph
}

#' run umap algorithm on igg data
#' 
#' @param knn_graph igraph object
#' @param cluster_method graph cluster method. Default = "louvain"
#' @param vertex.label provide a vector of names of genes, proteins, glycans
#' @return igraph cluster object
#' @export
cluster_graph<-function(knn_graph,cluster_method="prop"){
  
  #fast optimal
  # t_optimal_clusters<-igraph::cluster_optimal(knn_graph)
  
  #fast leiden
  # t_leiden_clusters<-igraph::cluster_leiden(knn_graph)
  
  #fast louvain
  if (cluster_method=="louvain"){
    t_clusters<-igraph::cluster_louvain(knn_graph)
  }
  #fast leading eigen
  # t_eigen_clusters<-igraph::cluster_leading_eigen(knn_graph)
  
  #fast greedy
  if (cluster_method=="fast_greedy"){
    t_clusters<-igraph::cluster_fast_greedy(knn_graph)
  }
  
  #edge
  if (cluster_method=="edge"){
    t_clusters<-igraph::cluster_edge_betweenness(knn_graph)
  }
  #prop
  if (cluster_method=="prop"){
    t_clusters<-cluster_label_prop(knn_graph)
  }
  
  t_clusters
}


plot_cluster_graph<-function(knn_graph,
                             myumap,
                             layout="umap",
                             cluster_method="prop")
{
  
  
  if (layout=="forcedirected"){
    l=layout_with_drl(knn_graph)
  }
  if (layout=="umap"){
    l=myumap$layout
  }
  # l=layout_with_fr(knn_graph)
  # l=layout_with_gem(knn_graph)
  # l=layout_with_graphopt(knn_graph)
  # l=layout_with_kk(knn_graph)
  # l=layout_with_lgl(knn_graph)
  # l=layout_with_mds(knn_graph)
  
  #fast optimal
  # t_optimal_clusters<-igraph::cluster_optimal(knn_graph)
  
  #fast leiden
  # t_leiden_clusters<-igraph::cluster_leiden(knn_graph)
  
  #fast louvain
  if (cluster_method=="louvain"){
    t_clusters<-igraph::cluster_louvain(knn_graph)
  }
  #fast leading eigen
  # t_eigen_clusters<-igraph::cluster_leading_eigen(knn_graph)
  
  #fast greedy
  if (cluster_method=="fast_greedy"){
    t_clusters<-igraph::cluster_fast_greedy(knn_graph)
  }
  
  #edge
  if (cluster_method=="edge"){
    t_clusters<-igraph::cluster_edge_betweenness(knn_graph)
  }
  #prop
  if (cluster_method=="prop"){
    t_clusters<-cluster_label_prop(knn_graph)
  }
  
  
  
  plot(t_clusters,knn_graph,layout=l)
  
}

#compare clusters alluvial
#' Alluvial plot to compare igraph cluster objects
#' 
#' @param cluster1 igraph cluster object 1
#' @param cluster2 igraph cluster object 2
#' 
#' @return alluvial plot object and table of membership
#' @export
compare_clusters_alluvial<-function(cluster1=daisy_clusters,
                                    cluster2=pagoda_clusters){
  daisy_mem<-membership(cluster1)
  pagoda_mem<-membership(cluster2)
  t1d_mem<-cbind.data.frame("row.names"=names(pagoda_mem),"daisy_mem"=as.vector(daisy_mem[which(names(daisy_mem)%in%names(pagoda_mem))]),
                            "pagoda_mem"=as.vector(pagoda_mem))
  
  clinfactors_forsankey<-data.frame(
    table(
      t1d_mem$daisy_mem,
      t1d_mem$pagoda_mem
    )
  )
  
  colnames(clinfactors_forsankey) <-
    c(
      "DAISY",
      "PAGODA",
      "Freq"
    )
  
  g<-alluvial::alluvial(clinfactors_forsankey[,1:2], freq=clinfactors_forsankey$Freq,
                        hide = clinfactors_forsankey$Freq < 5,
                        #col = ifelse( clinfactors_forsankey$TCGA_Histology == "LUAD", "#1B9E77", "#D95F02"),
                        blocks=T,cw=0.25,cex=1
  )
  
  list(g,"t1d_mem"=t1d_mem)
}

#make eigen matrix
#' Make eigen matrix
#' 
#' for each cluster of variables, calculate the first
#' principle component the cluster
#' 
#' @param mydata igg data
#' @param myclusters vector indicating cluster membership for each variable
#' 
#' @return data frame with rows as samples and columns as PC1 of each cluster
#' @export
make_eigen_matrix<-function(mydata,myclusters){
  # mydata<-daisy$igg
  eigenmatrix<-data.frame(row.names = rownames(mydata))
  if (class(myclusters)=="communities"){
    for (i in 1:max(membership(myclusters))){
      # i=2
      gly_names<-names(which(membership(myclusters)==i))
      subdata<-mydata[,colnames(mydata)%in%gly_names]
      mypca<-gmodels::fast.prcomp(subdata)
      eigenmatrix<-cbind.data.frame(eigenmatrix,mypca$x[,1])
      
    }
    colnames(eigenmatrix)<-paste0("Cluster",1:max(membership(myclusters)))
  }
  else {
    for (i in 1:max(myclusters)){
      # i=2
      gly_names<-names(which(myclusters==i))
      subdata<-mydata[,colnames(mydata)%in%gly_names]
      mypca<-gmodels::fast.prcomp(subdata)
      eigenmatrix<-cbind.data.frame(eigenmatrix,mypca$x[,1])
      
    }
    colnames(eigenmatrix)<-paste0("Cluster",1:max(myclusters))  
  }
  eigenmatrix
}


