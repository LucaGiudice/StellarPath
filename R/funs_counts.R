#' Normalization and Standardization of RNAseq data
#'
#' Library normalization of count values and standardization with rank and ratio
#'
#' @param counts Numeric matrix of raw count values from RNA sequencing (e.g. genes x samples)
#' @param groups Character vector, each element indicates the group of the sample in the same column position in m
#' @param winz Default FALSE, boolean to indicate if to adjuste outliers with winsorizing
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

#' Normalization of RNAseq data
#'
#' Library normalization and log trasformation of count values with edgeR
#' We followed the guide https://f1000research.com/articles/5-1408
#'
#' @param counts Numeric matrix of raw count values (e.g. molecules x samples)
#' @param groups Character vector, each element indicates the group of the sample in the same column position in counts
#' @param filtering Default TRUE, removing row molecules with low counts on the overall samples
#' @param winz Default FALSE, boolean to indicate if to adjuste outliers with winsorizing
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
  options(warn=-1)
  if(filtering & nrow(counts)>5000){
    counts=counts[rowSums(counts) > 10,]
    #counts=rem_extreme_sd(counts)
  }

  #Create edgelist object and compute L and M variables
  x=edgeR::DGEList(counts = counts, remove.zeros = TRUE)

  #Normalize data by library size
  x <- edgeR::calcNormFactors(x, method = "TMM")
  norm_counts <- edgeR::cpm(x, log=F);
  #Log normalization
  norm_counts[norm_counts<=1]=1.1
  norm_counts=countTransformers::l2Transformer(mat=norm_counts)$mat2

  # ids <- as.factor(colnames(counts))
  # dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts),
  #                                       colData = as.data.frame(ids),
  #                                       design = ~ ids)
  # # Perform the variance-stabilizing transformation (vst)
  # out_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  # # Access the expression values after vst transformation
  # norm_counts <- out_vst@assays@data@listData[[1]]

  #Adjust outliers
  if(!is.null(groups) & winz){
    norm_counts=do_winz(m=norm_counts, groups=groups, min=min ,max=max)
  }

  options(warn=0)
  return(norm_counts)
}

#' Rank and 0-1 standardize a normalized count matrix
#'
#' This function transforms the sample's counts in 0-1 standardized ranks.
#'
#' By default is type "total_rank": a molecule is ranked based on its count value with respect all the count values in the matrix and then
#' the ranks are 0-1 standardized.
#'
#' If type "sample_rank": For each sample, a molecule is ranked based on its count value and then the ranks are 0-1 standardized.
#' The most expressed molecule of a sample has the highest rank (e.g. 20000) divided by the number of molecules get 1
#' The least expressed molecule get a value equal to 0.
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
#' Winsorizing of the row values based on column groups
#'
#' https://www.r-bloggers.com/2011/06/winsorization/
#'
#' @param m numeric matrix of samples divided in groups at the columns
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param min Default 0.1, double value indicating the minimum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and max value
#' @param max Default 0.9, double value indicating the maximum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and greater than min value
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

#' Find significantly different molecules between groups
#'
#' Find significantly different molecules between groups evaluating both expression and stability.
#' Label the significantly different molecules in four categories:
#'
#' - Positive lfc and low variation for class 1
#'
#' - Positive lfc and low variation for class 2
#'
#' - Negative lfc and low variation for class 1
#'
#' - Negative lfc and low variation for class 2
#'
#' @param counts numeric count matrix of samples at the columns and molecules at the row (e.g. gene expression matrix)
#' @param groups character vector of groups' label s.t. each label indicates the belonging of a sample to one group
#' @param top_molecules_th Default 0.05, double value of the molecules that are selected as best ones
#' @return Dataframe of descriptive statistics related to the molecules in a count matrix having two groups of samples in comparison
#'
#' - logFoldChange measures the change in count values between first vs second group defined in the groups variable (unique(groups)[1]
#' vs unique(groups)[2])
#'
#' - stability measures how much is differently stable in first vs second group (higher positive and more is stable for first group)
#'
#' - activity_1 shows the mean of the count values (rank or expression) in the first group
#'
#' - activity_2 shows the mean of the count values (rank or expression) in the second group
#'
#' - stable_for indicates for which sample group the molecule is stable for
#'
#' - significant indicates if the molecule has been tested significant
#'
#' - sign_type identifies which test the molecule has passed to be significant.
#' 1)Positive lfc and low variation for class 1.
#' 2)Positive lfc and low variation for class 2.
#' -1)Negative lfc and low variation for class 1.
#' -2)Negative lfc and low variation for class 2.
#'
#' @export
find_SDR = function(counts, groups, top_molecules_th=0.05){
  #Determine how many features to retrieve signficant as the prob in Limma eBayes
  top_molecules_th=nrow(counts)*top_molecules_th
  if(top_molecules_th>300){
    top_molecules_th=300
  }

  #Find how much a molecule is different in standard deviation and (rank or expression) between the groups
  #independently by the group in which is weaker or stronger
  meanXgrs=get_meanXgr(m=counts, groups=groups, median=TRUE)
  sdXgrs=get_cvXgr(m=counts, groups=groups, mad=TRUE)
  meanXgrs[meanXgrs==0]=0.05
  sdXgrs[sdXgrs==0]=0.05
  diffsXgr=data.frame(row.names = rownames(meanXgrs),
                      logFoldChange=log2(meanXgrs[,1]/meanXgrs[,2]),
                      stability=(log2(sdXgrs[,1]/sdXgrs[,2])*-1),
                      activity_1=meanXgrs[,1],activity_2=meanXgrs[,2],
                      abs_lFC=abs(log2(meanXgrs[,1]/meanXgrs[,2])),
                      abs_stability=abs(log2(sdXgrs[,1]/sdXgrs[,2])*-1))
  diffsXgr=na.omit(diffsXgr)

  colnames(diffsXgr)=gsub("1",unique(groups)[1],colnames(diffsXgr))
  colnames(diffsXgr)=gsub("2",unique(groups)[2],colnames(diffsXgr))
  diffsXgr$stable_for=unique(groups)[1]
  diffsXgr$stable_for[diffsXgr$stability<0]=unique(groups)[2]

  #Subset the matrix to find 4 different differentially activated and disease-specific molecules
  #Positive lfc and low variation for class 1
  diffsXgr1_high=diffsXgr[diffsXgr$logFoldChange>0 & diffsXgr$stability>0,]
  #Positive lfc and low variation for class 2
  diffsXgr2_high=diffsXgr[diffsXgr$logFoldChange<0 & diffsXgr$stability<0,]
  #Negative lfc and low variation for class 1
  diffsXgr1_low=diffsXgr[diffsXgr$logFoldChange<0 & diffsXgr$stability>0,]
  #Negative lfc and low variation for class 2
  diffsXgr2_low=diffsXgr[diffsXgr$logFoldChange>0 & diffsXgr$stability<0,]
  #Compact in a list that it will be iterated over
  diffs_dfs=list(diffsXgr1_high=diffsXgr1_high,diffsXgr2_high=diffsXgr2_high,
                 diffsXgr1_low=diffsXgr1_low,diffsXgr2_low=diffsXgr2_low)

  #For each coordinated set of molecules, select the most significant
  select_significant = function(diffsXgr, top_molecules_th=300){
    ths=cbind(rep(seq(0.85,0.25,-0.1),each=7),seq(0.7,0.1,-0.1))
    i=1;n_top_molecules=0;
    while(n_top_molecules<=top_molecules_th & i<=nrow(ths)){
      mean_th=ths[i,1]
      sd_th=ths[i,2]

      top_molecules = rownames(diffsXgr)[diffsXgr$abs_lFC>=mean_th & diffsXgr$abs_stability>=sd_th]

      n_top_molecules = length(top_molecules)
      i=i+1

    }
    return(top_molecules)
  }

  top_molecules_l=lapply(diffs_dfs,function(x){
    select_significant(x,top_molecules_th)
  })

  #Combine the table of the stats with the signficant molecule list but label each molecule due to which criterium satisfied
  diffsXgr$significant=FALSE
  diffsXgr$significant[rownames(diffsXgr) %in% unique(unlist(top_molecules_l))]=TRUE
  diffsXgr$sign_type=0
  #Positive lfc and low variation for class 1
  diffsXgr$sign_type[rownames(diffsXgr) %in% top_molecules_l$diffsXgr1_high]=1
  #Positive lfc and low variation for class 2
  diffsXgr$sign_type[rownames(diffsXgr) %in% top_molecules_l$diffsXgr2_high]=2
  #Negative lfc and low variation for class 1
  diffsXgr$sign_type[rownames(diffsXgr) %in% top_molecules_l$diffsXgr1_low]=-1
  #Negative lfc and low variation for class 2
  diffsXgr$sign_type[rownames(diffsXgr) %in% top_molecules_l$diffsXgr2_low]=-2
  #Remove absolute values
  diffsXgr=diffsXgr[,-grep("abs_",colnames(diffsXgr))]
  return(diffsXgr)
}
