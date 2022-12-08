#' Human gene-gene interaction network
#'
#' @format list
#' \describe{
#'   \item{NETWORK}{gene-gene interaction network in adjacent and edge list format}
#' }
#' @source{
#'   \url{https://www.ebi.ac.uk/intact/home}
#' }
"dbs_hsa_net_l"

#' Mouse gene-gene interaction network
#'
#' @format list
#' \describe{
#'   \item{NETWORK}{gene-gene interaction network in adjacent and edge list format}
#' }
#' @source{
#'   \url{https://www.ebi.ac.uk/intact/home}
#' }
"dbs_mmu_net_l"

#' Mouse specific list of pathways about different omics
#'
#' @format The list has two levels, first level indicates the specie, second level indicates the set of database data:
#' \describe{
#'   \item{MIRWALK}{list where for each miRNA contains the set of its target genes}
#'   \item{RNARNA}{list that contains RNA-RNA interactions (e.g. lncrna_gene_targ_l contains the target genes for each lncRNA)}
#'   \item{PATHWAYS}{list that contains gene canonical pathways}
#' }
#' @source{
#'   \url{http://mirwalk.umm.uni-heidelberg.de/}
#'   \url{http://www.rnainter.org/}
#'   \url{http://baderlab.org/GeneSets}
#' }
"dbs_hsa_paths_l"

#' Human specific list of pathways about different omics
#'
#' @format The list has two levels, first level indicates the specie, second level indicates the set of database data:
#' \describe{
#'   \item{MIRWALK}{list where for each miRNA contains the set of its target genes}
#'   \item{RNARNA}{list that contains RNA-RNA interactions (e.g. lncrna_gene_targ_l contains the target genes for each lncRNA)}
#'   \item{PATHWAYS}{list that contains gene canonical pathways}
#' }
#' @source{
#'   \url{http://mirwalk.umm.uni-heidelberg.de/}
#'   \url{http://www.rnainter.org/}
#'   \url{http://baderlab.org/GeneSets}
#' }
"dbs_hsa_paths_l"

#' Vector of words banned to appear in the name of a pathway
#'
#' GO, KEGG, MSIGDB include pathways which are not informative or very context specific that would be
#' strange to use in most of the analysis, this list is removed from the database
#'
#' @format Vector of words banned to appear in the name of a pathway
#' \describe{
#'   \item{ban_list}{Vector of words banned to appear in the name of a pathway}
#' }
#' @source{
#'   \url{experience}
#' }
"ban_list"

#' OGD WT vs Normoxia WT
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{dataframe containing the mRNA count matrix}
#'   \item{miRNA}{dataframe containing the miRNA count matrixs}
#' }
#' @source{
#'   \url{}
#' }
"ogd_l"

#' Preview of a matrix like object
#'
#' Allows to quickly see the first rows and columns of a matrix like object
#'
#' @param m A matrix like object
#' @param r Default 5, integer indicating number of rows
#' @param c Default 5, integer indicating number of columns
#' @return Preview of a matrix like object
#' @export
#'
see = function(m, r=5, c=5){
  if(nrow(m)<5){
    r=nrow(m)
  }
  if(ncol(m)<5){
    c=ncol(m)
  }
  m[1:r,1:c]
}

#' Checks vectors of elements
#'
#' This function checks if two vectors are exactly the same both in the elements that contain and their order
#'
#' @param x A vector
#' @param y A vector
#' @return The result of the check. 0 not match. 1 match.
#' @export
#'
check_vectors = function(x, y){
  #Check if the vectors have the same length
  if(length(x)==length(y)){

    #Check if all elements in y are in x
    indx_ord=match(y,x)
    if(sum(!is.na(indx_ord))==0){
      stop("ERROR: not maching vectors")
      return(0)
    }

    #Check if all elements in x are in y
    indx_ord=match(x,y)
    if(sum(!is.na(indx_ord))==0){
      stop("ERROR: not maching vectors")
      return(0)
    }

    #Check if all elements are also in the same position
    corr_ord=seq(1,length(y))
    err=0
    for(k in corr_ord){
      check=indx_ord[k]==k
      if(is.na(check)){
        err=c(err,k)
      }else{
        if(!check){
          err=c(err,k)
        }
      }
    }
    err=err[-1]

    if(length(err)==0){return(1)}
    if(length(err)!=0){stop("ERROR: not maching vectors");return(0)}
  }else{
    stop("ERROR: not maching vectors");return(0)
  }
}

#' Compute standard deviation per group
#'
#' Compute standard deviation per column group of each feature row
#'
#' @param m numeric matrix
#' @param groups character or factor vector s.t. one element is the label indicating the group of one sample column
#' @param mad boolean, default is FALSE, set to true if you want to use median absolute deviations
#' @return numeric matrix wit same column of m, number of columns equal to the number of groups, a cell's value is
#' equal to the standard deviation of the values of the row feature in the original columns of the same group
#' @importFrom matrixStats colSds
#' @importFrom Rfast colMads
#' @export
#'
get_sdXgr = function(m, groups, mad=FALSE) {
  s <- split(as.data.frame(t(m)), groups)
  if(!mad){
    sdXgr = as.matrix(sapply(s, function(x){matrixStats::colSds(as.matrix(x))}))
  }else{
    sdXgr = as.matrix(sapply(s, function(x){Rfast::colMads(as.matrix(x))}))
  }
  rownames(sdXgr) = rownames(m)
  colnames(sdXgr) = names(s)
  indxs=match(unique(groups),colnames(sdXgr))
  sdXgr=sdXgr[,indxs]
  return(sdXgr)
}

#' Compute mean per group
#'
#' Compute mean per column group of each feature row
#'
#' @param m numeric matrix
#' @param groups character or factor vector s.t. one element is the label indicating the group of one sample column
#' @param median boolean, default is FALSE, set to true if you want to use the median instead of the mean
#' @return numeric matrix wit same column of m, number of columns equal to the number of groups, a cell's value is
#' equal to the mean of the values of the row feature in the original columns of the same group
#' @importFrom Rfast colmeans
#' @importFrom Rfast colMedians
#' @export
#'
get_meanXgr = function(m, groups, median=FALSE) {
  s <- split(as.data.frame(t(m)), groups)
  if(!median){
    meanXgr = as.matrix(sapply(s, function(x){Rfast::colmeans(as.matrix(x))}))
  }else{
    meanXgr = as.matrix(sapply(s, function(x){Rfast::colMedians(as.matrix(x))}))
  }
  rownames(meanXgr) = rownames(m)
  colnames(meanXgr) = names(s)
  indxs=match(unique(groups),colnames(meanXgr))
  meanXgr=meanXgr[,indxs]
  return(meanXgr)
}

#' Matthews correlation coefficient (MCC) score
#'
#' Calculate the Matthews correlation coefficient (MCC) score
#'
#' @param act actual values (vector), 1 (positive), or 0 (negative)
#' @param pred predict values (vector), 1 (positive), or 0 (negative)
#' @return MCC score
#' @export
get_mcc = function(real, predicted, one_group){
  predicted_binary = ifelse(predicted == one_group, 1, 0)
  real_binary = ifelse(real == one_group, 1, 0)

  TP <- sum(real_binary == 1 & predicted_binary == 1)
  TN <- sum(real_binary == 0 & predicted_binary == 0)
  FP <- sum(real_binary == 0 & predicted_binary == 1)
  FN <- sum(real_binary == 1 & predicted_binary == 0)

  denom <- as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (any((TP+FP) == 0, (TP+FN) == 0, (TN+FP) == 0, (TN+FN) == 0)) denom <- 1
  mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)
  return(mcc)
}

#' 0-1 standardization
#'
#' Vectorize the matrix, standardize based on 0 and 1, return the standardized matrix
#'
#' @param m numeric matrix
#' @param min Default 0, number to set the low limit of the range standardization
#' @param max Default 1,  number to set the high limit of the range standardization
#' @param capping Default FALSE, do not perform capping before standardize
#' @return range standardized input matrix
#' @importFrom scales rescale
#' @export
do_standard=function(m,min=0,max=1,capping=FALSE){
  if(capping==TRUE){
    m=apply(m,2,function(x){
      qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
      caps <- quantile(x, probs=c(.1, .90), na.rm = T)
      H <- 1.5 * IQR(x, na.rm = TRUE)
      x[x < (qnt[1] - H)] <- caps[1]
      x[x > (qnt[2] + H)] <- caps[2]
      return(x)
    })
  }

  #Convert matrix to vector
  geno_p_v=as.vector(m)
  #Scale the vector
  geno_p_v_sc=scales::rescale(geno_p_v, to=c(min,max))
  #Convert vector to matrix with same original dimensions
  geno_p2=matrix(geno_p_v_sc,nrow(m),ncol(m))
  #Apply original labels
  colnames(geno_p2)=colnames(m)
  rownames(geno_p2)=rownames(m)
  #Return
  return(geno_p2)
}
