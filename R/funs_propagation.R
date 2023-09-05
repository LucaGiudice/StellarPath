#' Harmonizes the matrix of sample's profiles with respect the network
#'
#' Remove all the features/rows that don't have their name matching with a node inside the network
#'
#' @param counts Numeric matrix with rows matching nodes of the network
#' @param net Adjacency matrix of a network s.t. nodes have same names of the rows in the matrix of profiles
#' @return The matrix of sample's profiles without the rows not having a match with a node of the network
#' @export
#'
harm_countsWnet = function(counts,net){
  message(">>>harmonizing count matrix with a biological network")
  #Recover patient genes and network genes
  pat_genes=rownames(counts)
  network_genes=rownames(net)

  #Find genes which exist in the profiles but do no exist in the network
  network_genes_2rem=setdiff(pat_genes,network_genes)
  #Remove
  if(length(network_genes_2rem)!=0){
    message(">>>the following number of sample's rownames are missing as nodes: ",length(network_genes_2rem))
    counts_tmp=as.matrix(counts[-match(network_genes_2rem,pat_genes),])
    colnames(counts_tmp)=colnames(counts)
    counts=counts_tmp
    pat_genes=rownames(counts)
  }

  #Find genes which exist in the network but do no exist in the profiles
  network_genes_2add=setdiff(network_genes,pat_genes)
  #Normalize
  if(length(network_genes_2add)!=0){
    m2bind=matrix(0,length(network_genes_2add),length(colnames(counts)))
    counts=rbind(counts,m2bind)
    rownames(counts)=c(pat_genes,network_genes_2add)
  }
  #Finish
  return(counts)
}

#' Network based propagation based on Random Walk with Restart
#'
#' This function performs propagation of the sample's profiles passed as columns of a matrix
#' Features of the profiles are mapped to the nodes of the network
#' Individually for each sample, the values of the features are projected on the corresponding nodes
#' Individually for each sample, the values of the features are propagated through the network due to the node's connectivity
#'
#' @param net Adjacency matrix or Edge list of a non-directed non-weighted network/graph
#' @param counts Numeric matrix of sample's profiles
#' @param norm Default row, options: row or column or laplacian, it indicates the graph normalization performed on the network
#' @param n_cores Default 2, an integer value greater than 0 indicating the number of cores to use for parallel computing
#' @param r Default 0.8, double value lower than 1 indicating the percentage of information that a gene keeps (0.8 is 80 percentage)
#' @param stop_step Default 200, integer value greater than 0 indicating the number of iterations of the propagation
#' @param stop_delta Default 1e-06, double value lower than 0 indicating the threshold under which all imputed propagation values are set 0
#' @param keep_zero_profiles Default FALSE, Boolean: TRUE if to keep profiles without information
#' @param keep_no_nodes Default FALSE, Boolean: TRUE if to keep profile's rows without node in the network
#' @return Matrix of propagated sample's profiles
#' @importFrom  igraph simplify
#' @importFrom  igraph degree
#' @importFrom  igraph delete.vertices
#' @importFrom  igraph get.adjacency
#' @importFrom  igraph graph_from_edgelist
#' @importFrom  igraph as_adjacency_matrix
#' @importFrom Rfast rowsums
#' @import data.table
#' @import Matrix
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import stats
#' @export
#'
get_propagated=function(net, counts, norm = "row", n_cores=2, r = 0.8, stop_step=200, stop_delta = 1e-06,
                 keep_zero_profiles=FALSE, keep_no_nodes=FALSE){

  message(">>network-based propagation")
  #Convert network in adjacency matrix ----
  if(class(net)[1]=="igraph"){
    #message(">converting igraph network in adjacency matrix")
    net=igraph::simplify(net)
    Isolated=which(igraph::degree(net)==0)
    net = igraph::delete.vertices(net, Isolated)
    adjM=igraph::get.adjacency(net, type = "both", attr = NULL, edges = F,
                               names = T, sparse = getIgraphOpt("sparsematrices"))

  }
  if(class(net)[1]=="matrix" | class(net)[1]=="data.frame"){
    if(ncol(net)==2){
      #message(">converting edge list network in adjacency matrix")
      net=igraph::graph_from_edgelist(net, directed = F)
      net=igraph::simplify(net)
      Isolated = which(igraph::degree(net)==0)
      net = igraph::delete.vertices(net, Isolated)
      adjM=igraph::as_adjacency_matrix(net, type = "both", attr = NULL, edges = F,
                                       names = T, sparse = igraph::getIgraphOpt("sparsematrices"))
    }else{
      #message(">network is an adjacency matrix")
      adjM=net
    }
  }

  #Prepare the profile matrix to the propagation -----
  #Remove gene rows that don't have a match withing the network as node
  data=as.matrix(harm_countsWnet(counts,adjM))
  #Divide profiles with starting information from one without
  zero_data=data[,colSums(data)==0]
  #Keep only profiles with information to propagate
  data1=as.matrix(data[,colSums(data)!=0])
  colnames(data1)=colnames(data)[colSums(data)!=0]
  data=data1;rm(data1)
  cnames=colnames(data)
  rnames=rownames(data)

  #Apply network normalization ----
  A <- adjM != 0
  if(norm == "row"){
    #message(">row normalization of the network")
    D <- Matrix::Diagonal(x = (Matrix::rowSums(A))^(-1))
    nadjM <- adjM %*% D
  }else if(norm == "column"){
    #message(">column normalization of the network")
    D <- Matrix::Diagonal(x = (Matrix::colSums(A))^(-1))
    nadjM <- D %*% adjM
  }else if(norm == "laplacian"){
    #message(">laplacian normalization of the network")
    D <- Matrix::Diagonal(x = (Matrix::colSums(A))^(-0.5))
    nadjM <- D %*% adjM %*% D
  }else{
    nadjM <- adjM
  }

  ind <- match(rownames(data), rownames(adjM))
  P0matrix <- matrix(0, nrow = nrow(nadjM), ncol = ncol(data))
  P0matrix[ind[!is.na(ind)], ] <- as.matrix(data[!is.na(ind), ])

  #Prepare matrix for the propagation ----
  rownames(P0matrix) <- rownames(adjM)
  colnames(P0matrix) <- cnames
  P0matrix <- Matrix::Matrix(P0matrix, sparse = T)

  #Set parallel computing -----
  if(.Platform$OS.type == "unix") {
    message(">>>parallel for Linux and Mac")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    message(">>>parallel for Windows")
    cl <- makeCluster(n_cores);
  }
  doParallel::registerDoParallel(cl)

  #Divide the profile matrix in chuncks for the propagation made in parallel -----
  ind=seq(1,ncol(P0matrix))
  sets=split(ind, ceiling(seq_along(ind)/round(length(ind)/n_cores)))
  P0_l=lapply(sets,function(x){Matrix::Matrix(as.matrix(P0matrix[,x]), sparse = T)})

  #Random walk with restart -----
  PTmatrix=foreach(k = 1:length(P0_l),.inorder = T,.combine = c("cbind"),.packages=c("Matrix")) %dopar% {

    P0matrix1=P0_l[[k]]
    for(j in 1:ncol(P0matrix1)){
      P0 <-P0matrix1[, j]; step <- 0; delta <- 1; PT <- P0;
      while(step <= stop_step){
        PX <- (1 - r) * nadjM %*% PT + r * P0
        delta <- sum(abs(PX - PT))
        PT <- PX
        step <- step + 1
      }
      PT=as.matrix(PT)
      P0matrix1[,j]=PT
    }
    return(P0matrix1)

  }

  #Finish the matrix, close parallel
  PTmatrix[PTmatrix < stop_delta]= 0
  colnames(PTmatrix)=colnames(P0matrix)
  parallel::stopCluster(cl)

  #Add profiles with zero startign information if we want to keep them
  if(keep_zero_profiles == TRUE & ncol(zero_data)!=0){
    message(">>>keeping samples without information")
    tmp=merge(PTmatrix,zero_data,by="row.names")
    rownames(tmp)=tmp[,1];tmp=tmp[,-1]
    PTmatrix=tmp
  }

  if(keep_no_nodes == TRUE){
    message(">>>keeping samples's molecules not represented as nodes")
    genes_exis=rownames(PTmatrix)
    genes_miss=setdiff(rownames(counts),genes_exis)
    genes_miss_df=counts[match(genes_miss,rownames(counts)),]
    PTmatrix=rbind(PTmatrix,genes_miss_df)
    rownames(PTmatrix)=c(genes_exis,genes_miss)
  }

  #End the script
  PTmatrix=as.matrix(PTmatrix)
  PTmatrix=PTmatrix[Rfast::rowsums(as.matrix(PTmatrix))!=0,]
  message(">>end propagation")
  return(PTmatrix)
}
