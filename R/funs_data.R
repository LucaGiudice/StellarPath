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
"dbs_mmu_paths_l"

#' Vector of words used to ban pathways that are difficult to interpret
#'
#' GO, KEGG, MSIGDB include pathways which are not informative or are very context specific.
#' Ban list includes terms like "SOUND", "STABILIZATION", "ODONTOGENESIS", "PREGNANCY"
#' By default any pathway including such terms is not used in the analysis.
#' The ban can also be avoided personalizing the function setting.
#'
#' @format Vector of words used to ban pathways that are difficult to interpret
#' \describe{
#'   \item{ban_list}{Vector of words used to ban pathways that are difficult to interpret}
#' }
#' @source{
#'   \url{experience of using the pathway analysis}
#' }
"ban_list"

#' List of data describing oxygen-glucose deprived mouse samples vs non-treated ones
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#' }
#' @source{
#'   \url{}
#' }
"ogd_l"

#' List of data describing cll_l samples vs non-treated ones
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#' }
#' @source{
#'   \url{}
#' }
"cll_l"

#' List of data describing multiple small tcga datasets (part 1)
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#'   \item{mutation}{matrix containing the binary information of mutations mapped to genes}
#' }
#' @source{
#'   \url{}
#' }
"tcga_p1_l"

#' List of data describing multiple medium-small tcga datasets (part 2)
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#'   \item{mutation}{matrix containing the binary information of mutations mapped to genes}
#' }
#' @source{
#'   \url{}
#' }
"tcga_p2_l"

#' List of data describing multiple medium-large tcga datasets (part 3)
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#'   \item{mutation}{matrix containing the binary information of mutations mapped to genes}
#' }
#' @source{
#'   \url{}
#' }
"tcga_p3_l"

#' List of data describing multiple large tcga datasets (part 4)
#'
#' @format The list has three elements
#' \describe{
#'   \item{info}{dataframe containing information about the mouse samples}
#'   \item{mRNA}{matrix containing the bulk RNA counts (both genes and lncRNAs)}
#'   \item{miRNA}{matrix containing the miRNA counts}
#'   \item{mutation}{matrix containing the binary information of mutations mapped to genes}
#' }
#' @source{
#'   \url{}
#' }
"tcga_p4_l"

#' Preview of a matrix like object
#'
#' Allows to quickly see the first rows and columns of a matrix-like object
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

#' Compute coefficient of variation per group for each feature
#'
#' Compute coefficient of variation per column group of each feature row
#'
#' @param m numeric matrix
#' @param groups character or factor vector s.t. one element is the label indicating the group of one sample column
#' @param mad boolean, default is FALSE, set to true if you want to use median absolute deviations
#' @return numeric matrix with same columns as m, number of columns equal to the number of groups, a cell's value is
#' equal to the coefficient of variation of the values of the row feature in the original columns of the same group
#' @importFrom matrixStats colSds
#' @importFrom Rfast colmeans
#' @importFrom Rfast colMads
#' @importFrom Rfast colMedians
#' @export
get_cvXgr <- function(m, groups, mad=FALSE){
  s <- split(as.data.frame(t(m)), groups)
  if (!mad) {
    sdXgr <- as.matrix(sapply(s, function(x) { matrixStats::colSds(as.matrix(x)) }))
    meanXgr <- as.matrix(sapply(s, function(x) { Rfast::colmeans(as.matrix(x)) }))
  } else {
    sdXgr <- as.matrix(sapply(s, function(x) { Rfast::colMads(as.matrix(x)) }))
    meanXgr <- as.matrix(sapply(s, function(x) { Rfast::colMedians(as.matrix(x)) }))
  }

  cvXgr <- sdXgr / meanXgr
  rownames(cvXgr) <- rownames(m)
  colnames(cvXgr) <- names(s)

  indxs <- match(unique(groups), colnames(cvXgr))
  cvXgr <- cvXgr[, indxs]

  return(cvXgr)
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

#' Calculates the classification metrics given the true and predicted labels.
#'
#' This function takes in a vector of true binary labels and a vector of predicted binary labels.
#' Returns the classification metrics (Matthews correlation coefficient, accuracy, sensitivity, and specificity).
#'
#' @param true_labels A vector of true binary labels
#' @param predicted_labels A vector of predicted binary labels
#' @return A dataframe containing the calculated metrics: Matthews correlation coefficient, accuracy, sensitivity, and specificity.
#' @examples
#' true_labels=c("LATE","LATE","LATE","EARLY","EARLY","EARLY")
#' predicted_labels=c("LATE","LATE","EARLY","EARLY","EARLY","LATE")
#' get_perf_metrics(true_labels, predicted_labels)
get_perf_metrics <- function(true_labels, predicted_labels) {

  # Count the number of true positives (TP), true negatives (TN),
  # false positives (FP), and false negatives (FN)
  A=unique(true_labels)[1]
  B=unique(true_labels)[2]

  TP <- sum(true_labels == A & predicted_labels == A)
  TN <- sum(true_labels == B & predicted_labels == B)
  FP <- sum(true_labels == B & predicted_labels == A)
  FN <- sum(true_labels == A & predicted_labels == B)

  # Calculate accuracy, sensitivity, and specificity
  accuracy <- (TP + TN) / length(true_labels)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)

  denom <- as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (any((TP+FP) == 0, (TP+FN) == 0, (TN+FP) == 0, (TN+FN) == 0)) denom <- 1
  mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)

  # Return a list containing the calculated metrics
  df_perf=data.frame(mcc = mcc, accuracy = accuracy, sensitivity = sensitivity, specificity = specificity)
  return(df_perf)
}

#' Calculates the Area Under the Receiver Operating Characteristic (ROC) Curve (AUC).
#'
#' This function takes in a vector of true labels and a vector of predicted scores/probabilities, and returns the AUC using the roc() function from the pROC package.
#'
#' @param true_labels A vector of true binary labels
#' @param predicted_scores A vector of predicted scores or probabilities for the positive class.
#' @return A numeric value representing the AUC.
#' @importFrom pROC roc
#' @importFrom pROC auc
#' @examples
#' true_labels=c("LATE","LATE","LATE","EARLY","EARLY","EARLY")
#' predicted_scores=c(0.8,0.9,0.2,0.4,0.3,0.9)
#' get_auc(true_labels, predicted_scores)
get_auc = function(true_labels, predicted_scores){
  roc_obj <- pROC::roc(true_labels, predicted_scores)
  auc <- pROC::auc(roc_obj)
  res_l=list(roc_obj=roc_obj,auc=auc)
  return(res_l)
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

#' Create a dataframe from a vector of strings containing attribute-value pairs
#'
#' This function takes a vector of strings, where each string represents multiple attribute-value pairs separated by '|'. The attribute-value pairs are in the format 'attribute;;value'. The function extracts the attribute names and values and creates a dataframe where each attribute becomes a column and the values are assigned to the corresponding rows.
#'
#' @param strings A vector of strings containing attribute-value pairs
#'
#' @return A dataframe with attribute names as columns and attribute values assigned to the respective rows
#' @examples
#' strings <- c("age at death;;93|braak score;;IV|diagnosis;;Alzheimer's disease",
#'              "age at death;;85|braak score;;II|diagnosis;;Parkinson's disease",
#'              "age at death;;78|braak score;;III|diagnosis;;Alzheimer's disease")
#' df <- create_df4attr(strings)
#' print(df)
#' @export
create_df4attr <- function(strings) {
  # Initialize an empty list to store attribute-value pairs
  attributes <- list()

  # Extract unique attribute names
  attribute_names <- unique(gsub(";;.*", "", unlist(strsplit(strings[1], "\\|", fixed = FALSE))))

  # Iterate over each string in the vector
  for (i in 1:length(strings)) {
    # Split the string by "|"
    string_parts <- strsplit(strings[i], "\\|", fixed = FALSE)[[1]]

    # Split each part by ";;" to extract attribute names and values
    for (j in 1:length(string_parts)) {
      attribute_parts <- strsplit(string_parts[j], ";;", fixed = FALSE)[[1]]

      # Extract the attribute name and value
      attribute_name <- attribute_parts[1]
      attribute_value <- attribute_parts[2]

      # Append the attribute-value pair to the list
      attributes[[attribute_name]] <- c(attributes[[attribute_name]], attribute_value)
    }
  }

  check_len=sapply(attributes,length)
  attributes=attributes[check_len==median(check_len)]

  # Create a dataframe from the attribute-value pairs
  df <- as.data.frame(attributes)

  return(df)
}


#' Remove the ending part from a word
#'
#' This function takes a word and removes the ending part if it matches the pattern '.<number>' (dot followed by a number). It returns the modified word without the ending part.
#'
#' @param word The word to remove the ending from
#'
#' @return The modified word without the ending part
#' @examples
#' word <- "example.1"
#' modified_word <- remove_ending(word)
#' print(modified_word)
#' @export
remove_ending <- function(word) {
  # Find the position of the last dot
  dot_position <- regexpr("\\.[0-9]+$", word)

  # Check if a dot and number were found at the end of the word
  if (dot_position > 0) {
    # Remove the ending part from the word
    word <- substr(word, 1, dot_position - 1)
  }

  return(word)
}

#' Remove duplicated rows based on row names and keep the highest values
#'
#' This function takes a matrix and a vector of row names. It removes duplicated rows based on the row names, keeping only the row with the highest values. The function returns a new matrix with the duplicated rows removed.
#'
#' @param matrix The input matrix with duplicated rows
#' @param row_names Default NULL, vector of row names corresponding to the matrix
#' @param avg Default TRUE, compute the column-wise average of the duplicated rows
#'
#' @return A new matrix with duplicated rows removed, keeping only the row with the highest values
#' @examples
#' matrix <- matrix(seq(1,12), ncol = 2)
#' row_names <- c("A", "B", "A", "C", "B", "C")
#' result <- remove_duplicates(matrix, row_names)
#' print(result)
#' @export
remove_duplicates <- function(matrix, row_names=NULL, avg=TRUE) {
  if(is.null(row_names)){
    row_names=rownames(matrix)
  }

  # Check if the number of row names matches the number of rows in the matrix
  if(length(row_names) != nrow(matrix)) {
    stop("Number of row names does not match the number of rows in the matrix")
  }

  if(!avg){
    # Combine the row names and matrix into a data frame
    df <- data.frame(row_names, matrix)
    # Convert row names to factors for efficient indexing
    df$row_names <- as.factor(df$row_names)
    # Get the indices of the rows to keep (highest values)
    indices <- tapply(seq_len(nrow(df)), df$row_names, function(x) x[which.max(as.matrix(rowSums(df[x, -1])))])
    # Select the rows with the highest values
    result <- df[indices, ]
    # Remove the row names column
    result$row_names <- NULL
    result=as.matrix(result)
  }else{
    counts=aggregate(matrix, by=list(row.names(matrix)), FUN=median)
    rownames(counts)=counts[,1];counts=counts[,-1]
    result=as.matrix(counts)
  }

  return(result)
}

#' Remove rows with the lowest and highest standard deviation from a numeric matrix.
#'
#' This function takes a numeric matrix as input and removes the rows with the
#' lowest and highest standard deviation based on the specified percentile thresholds.
#'
#' @param numeric_matrix A numeric matrix with rows representing samples and
#' columns representing features (variables).
#' @param low_percentile The percentile threshold (default is 5) to determine the
#' rows to remove with the lowest standard deviation. It should be a value between 0 and 100.
#' @param high_percentile The percentile threshold (default is 95) to determine the
#' rows to remove with the highest standard deviation. It should be a value between 0 and 100.
#'
#' @return A new numeric matrix with rows having the lowest and highest standard deviation
#' removed.
#'
#' @examples
#' # Sample numeric matrix
#' numeric_matrix <- matrix(
#'   c(1, 2, 1, 3, 1,
#'     5, 5, 4, 4, 6,
#'     9, 8, 8, 8, 10,
#'     12, 11, 12, 10, 13),
#'   nrow = 4, byrow = TRUE
#' )
#'
#' # Remove rows with the lowest 5% and highest 5% standard deviation
#' filtered_matrix <- remove_lowest_highest_sd_rows(numeric_matrix)
#' print(filtered_matrix)
#'
#' @export
rem_extreme_sd <- function(numeric_matrix, low_percentile = 5, high_percentile = 95) {
  # Calculate the standard deviation for each row
  row_sd <- apply(numeric_matrix, 1, sd)

  # Find the thresholds for the lowest and highest standard deviation
  low_threshold <- quantile(row_sd, low_percentile / 100)
  high_threshold <- quantile(row_sd, high_percentile / 100)

  # Identify the rows to remove (those with standard deviation below the low_threshold or above the high_threshold)
  rows_to_remove <- which(row_sd < low_threshold | row_sd > high_threshold)

  # Remove the rows from the numeric matrix
  filtered_numeric_matrix <- numeric_matrix[-rows_to_remove, ]

  return(filtered_numeric_matrix)
}

