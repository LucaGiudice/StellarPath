#' Normalization and Rank transformation of RNAseq data
#'
#' Library normalize count values and then apply rank transformation
#'
#' @param counts Numeric matrix of raw count values from RNA sequencing (e.g. genes x samples)
#' @param groups Character vector, each element indicates the group of the sample in the same column position in m
#' @param winz Default FALSE, not adjusting outliers by Winsorizing
#' @param min Default 0.1, double value indicating the minimum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and max value
#' @param max Default 0.9, double value indicating the maximum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and greater than min value
#' @param type Default total_rank, character string which defines the direction of ranking, options in the help of the get_rank_01 function
#' @return Ratio transformed and library normalized count values
#' @export
#'
prepare_RNAseq = function(counts, groups=NULL, winz=FALSE, min=0.1, max=0.9, type="total_rank"){
  #Library Normalize and rank transform
  lcpm=suppressMessages(get_lcpm(counts=counts, groups=groups, winz=winz))
  rlcpm=get_rank_01(m=lcpm, type=type)
  return(rlcpm)
}

#' Library normalize a raw count matrix
#'
#' This function filters, library normalizes and log transform count values with edgeR
#'
#' @param counts Numeric matrix of raw count values (e.g. molecules x samples)
#' @param groups Character vector, each element indicates the group of the sample in the same column position in counts
#' @param filtering Default TRUE, removing row molecules with low counts on the overall samples
#' @param winz Default FALSE, not adjusting outliers by Winsorizing
#' @param min Default 0.1, double value indicating the minimum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and max value
#' @param max Default 0.9, double value indicating the maximum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and greater than min value
#' @return library normalized count matrix in lcpm
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR cpm
#' @importFrom countTransformers l2Transformer
#' @export
#'
get_lcpm = function(counts, groups=NULL, filtering=TRUE, winz=FALSE, min=0.1, max=0.9){
  #Create edgelist object and compute L and M variables
  options(warn=-1)
  x=edgeR::DGEList(counts = counts, remove.zeros = TRUE)
  if(filtering){
    if(!is.null(groups)){
      #Remove genes which are low expressed and non informative of the difference between groups
      keep.exprs <- edgeR::filterByExpr(x, group=groups)
      x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    }else{
      #Remove genes which are low expressed and non informative of the difference between groups
      keep.exprs <- edgeR::filterByExpr(x)
      x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    }
  }
  #Normalize data by library size
  x <- edgeR::calcNormFactors(x, method = "TMM")
  norm_counts <- edgeR::cpm(x, log=F);

  #Adjust outliers
  if(!is.null(groups) & winz){
    norm_counts=do_winz(m=norm_counts, groups=groups, min=min ,max=max)
  }

  #Log normalization
  norm_counts[norm_counts<=1]=1.1
  lcpm=countTransformers::l2Transformer(mat=norm_counts)$mat2
  options(warn=0)
  return(lcpm)
}

#' Apply 0-1 standardized rank transformation
#'
#' This function transforms the sample's counts in 0-1 standardized ranks.
#'
#' By default is type "sample_rank": For each sample, a molecule is ranked based on its count value and then the ranks are 0-1 standardized.
#' The most expressed molecule of a sample has the highest rank (e.g. 20000) divided by the number of molecules get 1
#' The least expressed molecule get a value equal to 0.
#'
#' If type is "total_rank": a molecule is ranked based on its count value with respect all the count values in the matrix and then
#' the ranks are 0-1 standardized.
#'
#' If type is "expression": no transformation is applied, lcpm count values are used
#'
#' @param m Numeric matrix of count values (e.g. molecules x samples)
#' @param type Default sample, character string which defines the direction of ranking, it can be "sample_rank", "total_rank", "expression"
#' @return 0-1 standardized ranks of the counts for each sample (default) or for the total matrix (type: total_rank) or no transformation (type: expression)
#' @export
#'
get_rank_01 = function(m, type="total_rank"){
  if(type=="sample_rank"){
    m1=apply(m,2,function(x){
      rr=rank(x)
      ratio=rr/length(x)
      return(ratio)
    })
  }
  if(type=="total_rank"){
    vm=as.vector(m)
    rr=rank(vm)
    ratio=rr/length(rr)
    m1=matrix(ratio,nrow(m),ncol(m))
    colnames(m1)=colnames(m)
    rownames(m1)=rownames(m)
  }
  if(type=="expression"){
    m1=m
  }
  return(m1)
}

#' Winsorizing
#'
#' This function performs Winsorizing of the row values based on column groups
#'
#' @param m numeric matrix of samples divided in groups at the columns
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @return Input winzorised numeric matrix
#' @importFrom DescTools Winsorize
#' @export
#'
do_winz = function(m, groups, min=0.05, max=0.95){
  groups=as.character(groups)
  for(gr in unique(groups)){
    indxs=which(gr==groups)
    for(ge in 1:nrow(m)){
      m[ge,indxs]=DescTools::Winsorize(m[ge,indxs], probs = c(min, max))
    }
  }
  return(m)
}

#' Find the specific differential role/expressed molecules between two groups
#'
#' This function detects the different molecules between two groups based on the differential in standard deviation
#' It computes the average of the values (ranks or expression) associated to a molecule in each group
#' It computes the standard deviation of the values (ranks or expression) associated to a molecule in each group
#' It ranks the molecules based on the difference in average and standard deviation between the two groups
#' It selects the top 10% which the highest average and standard deviation because it means that:
#' the molecule's values differ between the two groups
#' the molecule's values are stable in one group but vary in the other
#'
#' The meaning of the selected molecules depends by how the count matrix has been normalized in "prepare_data" due to the
#' variable "type".
#'
#' By default is type "sample_rank": For each sample, a molecule is ranked based on its count value and then the ranks are 0-1 standardized.
#' The most expressed molecule of a sample has the highest rank (e.g. 20000) divided by the number of molecules get 1
#' The least expressed molecule get a value equal to 0. In this case, the selected molecules differ in how much they have
#' been activated or silenced in each sample of the group with respect any other sample of the opposite group.
#'
#' If type is "total_rank": a molecule is ranked based on its count value with respect all the count values in the matrix and then
#' the ranks are 0-1 standardized. In this case, the selected molecules differ in how much they have been activated or
#' silenced with respect all the other molecules in study.
#'
#' If type is "expression": no transformation and standardization is applied. In this case, the
#' result reflects the typical differential expression analysis. A selected molecule is up or down regulated with respect
#' only the same molecule in the opposite group.
#'
#' @param counts count matrix of samples at the columns and molecules at the row (e.g. gene expression matrix)
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param top_molecules_th Default 0.1, double value of the molecules that are selected as best ones
#' @param force_advanced Default FALSE, the best molecules are not labelled and divided based on their fold change.
#' This means that a pathway/target set will be enriched for a set of molecules which fold change is not coordinate in one direction.
#' Keep FALSE when the number of molecules composing the count matrix is small or when the direction of fold change does not matter.
#' @return Dataframe of descriptive statistics related to the row molecules in a count matrix having two groups of samples in comparison
#' logFoldChange measures the change in count values between first vs second group defined in the groups variable (unique(groups)[1] vs unique(groups)[2])
#' Stability measures how much is differently stable in first vs second group (higher positive and more is stable for first group)
#' Activity_1 shows the mean of the count values (rank or expression) in the first group
#' Activity_2 shows the mean of the count values (rank or expression) in the second group
#' Stable_for indicates in which sample group the molecule is stable for
#' Significant indicates if the molecule has been tested significant
#' Barcode is an internal variable used by the software to understand which test the molecule has passed to be significant
#' @export
#'
find_SDR = function(counts, groups, top_molecules_th=0.05, force_directional=FALSE){
  if(nrow(counts)<=4000 & !force_directional){
    DR_res = find_SDR_undirectional(counts=counts, groups=groups, top_molecules_th=top_molecules_th)
  }else{
    DR_res = find_SDR_directional(counts=counts, groups=groups, top_molecules_th=top_molecules_th)
  }
  return(DR_res)
}

#' Find undirectional specific differential role/expressed molecules between two groups
#'
#' Internal function of find_SDR
#'
#' @param counts count matrix of samples at the columns and molecules at the row (e.g. gene expression matrix)
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param top_molecules_th Default 0.1, double value of the molecules that are selected as best ones
#' @return Dataframe of descriptive statistics related to the row molecules in a count matrix having two groups of samples in comparison
#' logFoldChange measures the change in count values between first vs second group defined in the groups variable (unique(groups)[1] vs unique(groups)[2])
#' Stability measures how much is differently stable in first vs second group (higher positive and more is stable for first group)
#' Activity_1 shows the mean of the count values (rank or expression) in the first group
#' Activity_2 shows the mean of the count values (rank or expression) in the second group
#' Stable_for indicates in which sample group the molecule is stable for
#' Significant indicates if the molecule has been tested significant
#' Barcode is an internal variable used by the software to understand which test the molecule has passed to be significant
#'
find_SDR_undirectional = function(counts, groups, top_molecules_th=0.05){
  #Determine how many molecules to retrieve signficant as the prob parameter in Limma eBayes
  top_molecules_th=nrow(counts)*top_molecules_th

  #Find how much a molecule is different in standard deviation and (rank or expression) between the groups
  #independently by the group in which is weaker or stronger
  meanXgrs=get_meanXgr(m=counts, groups=groups, median=TRUE)
  sdXgrs=get_sdXgr(m=counts, groups=groups, mad=TRUE)
  diffsXgr=data.frame(row.names = rownames(meanXgrs),
                      logFoldChange=log2(meanXgrs[,1]/meanXgrs[,2]),
                      Stability=(log2(sdXgrs[,1]/sdXgrs[,2])*-1),
                      #*-1 in order to have the greatest positive Stability for the best molecule for the first class
                      Activity_1=meanXgrs[,1],Activity_2=meanXgrs[,2],
                      abs_lFC=abs(log2(meanXgrs[,1]/meanXgrs[,2])),
                      abs_Stability=abs(log2(sdXgrs[,1]/sdXgrs[,2])*-1))

  colnames(diffsXgr)=gsub("1",unique(groups)[1],colnames(diffsXgr))
  colnames(diffsXgr)=gsub("2",unique(groups)[2],colnames(diffsXgr))
  diffsXgr$Stable_for=unique(groups)[1]
  diffsXgr$Stable_for[diffsXgr$Stability<0]=unique(groups)[2]

  #Set the combinations of quantiles to iterate over in order to gradually include more and more significant molecules until
  #threshold is reached, the difference in average between the groups is the criterium that gets weak faster
  ths=cbind(rep(seq(0.85,0.25,-0.1),each=5),seq(0.7,0.1,-0.1))

  #Select the best molecules, high difference in sd and rank/expression
  i=1;n_top_molecules=0;
  while(n_top_molecules<=top_molecules_th & i<=nrow(ths)){
    mean_th=ths[i,1]
    sd_th=ths[i,2]

    top_molecules=rownames(diffsXgr)[diffsXgr$abs_lFC>=mean_th &
                                       diffsXgr$abs_Stability>=sd_th]

    n_top_molecules=length(top_molecules)
    i=i+1
  }

  #Set a variable to discern significant and not molecules
  diffsXgr$Significant=FALSE
  diffsXgr$Significant[rownames(diffsXgr) %in% top_molecules]=TRUE
  #Set a variable to discern significant and not molecules detected in an undirectional manner (9)
  diffsXgr$Barcode=0
  diffsXgr$Barcode[rownames(diffsXgr) %in% top_molecules]=9
  diffsXgr=diffsXgr[,-grep("abs_",colnames(diffsXgr))]
  return(diffsXgr)
}

#' Find directional specific differential role/expressed molecules between two groups
#'
#' Internal function of find_SDR
#'
#' @param counts count matrix of samples at the columns and molecules at the row (e.g. gene expression matrix)
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param top_molecules_th Default 0.1, double value of the molecules that are selected as best ones
#' @return Dataframe of descriptive statistics related to the row molecules in a count matrix having two groups of samples in comparison
#' logFoldChange measures the change in count values between first vs second group defined in the groups variable (unique(groups)[1] vs unique(groups)[2])
#' Stability measures how much is differently stable in first vs second group (higher positive and more is stable for first group)
#' Activity_1 shows the mean of the count values (rank or expression) in the first group
#' Activity_2 shows the mean of the count values (rank or expression) in the second group
#' Stable_for indicates for which sample group the molecule is stable for
#' Significant indicates if the molecule has been tested significant
#' Barcode is an internal variable used by the software to understand which test the molecule has passed to be significant
#'
find_SDR_directional = function(counts, groups, top_molecules_th=0.05){
  #Determine how many features to retrieve signficant as the prob in Limma eBayes
  top_molecules_th=nrow(counts)*top_molecules_th

  #Find how much a molecule is different in standard deviation and (rank or expression) between the groups
  #independently by the group in which is weaker or stronger
  meanXgrs=get_meanXgr(m=counts, groups=groups, median=TRUE)
  sdXgrs=get_sdXgr(m=counts, groups=groups, mad=TRUE)
  meanXgrs[meanXgrs==0]=0.1
  sdXgrs[sdXgrs==0]=0.1
  diffsXgr=data.frame(row.names = rownames(meanXgrs),
                      logFoldChange=log2(meanXgrs[,1]/meanXgrs[,2]),
                      Stability=(log2(sdXgrs[,1]/sdXgrs[,2])*-1),
                      Activity_1=meanXgrs[,1],Activity_2=meanXgrs[,2],
                      abs_lFC=abs(log2(meanXgrs[,1]/meanXgrs[,2])),
                      abs_Stability=abs(log2(sdXgrs[,1]/sdXgrs[,2])*-1))
  diffsXgr=na.omit(diffsXgr)

  colnames(diffsXgr)=gsub("1",unique(groups)[1],colnames(diffsXgr))
  colnames(diffsXgr)=gsub("2",unique(groups)[2],colnames(diffsXgr))
  diffsXgr$Stable_for=unique(groups)[1]
  diffsXgr$Stable_for[diffsXgr$Stability<0]=unique(groups)[2]

  #Subset the matrix to find 4 different differentially activated and disease-specific molecules
  #Positive lfc and low variation for class 1
  diffsXgr1_high=diffsXgr[diffsXgr$logFoldChange>0 & diffsXgr$Stability>0,]
  #Positive lfc and low variation for class 2
  diffsXgr2_high=diffsXgr[diffsXgr$logFoldChange<0 & diffsXgr$Stability<0,]
  #Negative lfc and low variation for class 1
  diffsXgr1_low=diffsXgr[diffsXgr$logFoldChange<0 & diffsXgr$Stability>0,]
  #Negative lfc and low variation for class 2
  diffsXgr2_low=diffsXgr[diffsXgr$logFoldChange>0 & diffsXgr$Stability<0,]
  #Compact in a list that it will be iterated over
  diffs_dfs=list(diffsXgr1_high=diffsXgr1_high,diffsXgr2_high=diffsXgr2_high,
                 diffsXgr1_low=diffsXgr1_low,diffsXgr2_low=diffsXgr2_low)

  #For each coordinated set of molecules, select the most significant
  select_significant = function(diffsXgr, top_molecules_th=300){
    ths=cbind(rep(seq(0.85,0.25,-0.1),each=5),seq(0.7,0.1,-0.1))
    i=1;n_top_molecules=0;
    while(n_top_molecules<=top_molecules_th & i<=nrow(ths)){
      mean_th=ths[i,1]
      sd_th=ths[i,2]

      top_molecules = rownames(diffsXgr)[diffsXgr$abs_lFC>=mean_th &
                                           diffsXgr$abs_Stability>=sd_th]

      n_top_molecules = length(top_molecules)
      i=i+1

    }
    return(top_molecules)
  }

  top_molecules_l=lapply(diffs_dfs,function(x){
    select_significant(x,top_molecules_th)
  })

  #Combine the table of the stats with the signficant molecule list but label each molecule due to which criterium satisfied
  diffsXgr$Significant=FALSE
  diffsXgr$Significant[rownames(diffsXgr) %in% unique(unlist(top_molecules_l))]=TRUE
  diffsXgr$Barcode=0
  #Positive lfc and low variation for class 1
  diffsXgr$Barcode[rownames(diffsXgr) %in% top_molecules_l$diffsXgr1_high]=1
  #Positive lfc and low variation for class 2
  diffsXgr$Barcode[rownames(diffsXgr) %in% top_molecules_l$diffsXgr2_high]=2
  #Negative lfc and low variation for class 1
  diffsXgr$Barcode[rownames(diffsXgr) %in% top_molecules_l$diffsXgr1_low]=-1
  #Negative lfc and low variation for class 2
  diffsXgr$Barcode[rownames(diffsXgr) %in% top_molecules_l$diffsXgr2_low]=-2
  #Remove absolute values
  diffsXgr=diffsXgr[,-grep("abs_",colnames(diffsXgr))]
  return(diffsXgr)
}
