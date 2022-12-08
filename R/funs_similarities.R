#' Create Similarity Network based on Adjusted Weighted Jaccard similarity
#'
#' Determine which is the similarity between column profiles based on Adjusted Weighted Jaccard similarity
#' Higher the variability of the same features between two profiles and lower is the similarity
#' Smaller the variability of the same features between two profiles and higher is the similarity
#'
#' @param m Numeric matrix of column profiles and row features (genes per samples)
#' @return Similarity Network based on Adjusted Weighted Jaccard function
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import stats
#' @export
build_adjWJ_PSN <- function(m, comb2=NULL) {
  if(is.null(comb2)){
    samples <- 1:ncol(m)
    comb <- data.table::CJ(samples, samples)
    comb[, i := .I]
    comb <- data.table::melt(comb, 'i')
    data.table::setorder(comb, value)
    v2 <- paste0("V", 1:2)
    comb[, variable2 := v2 , keyby = i]
    comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')
    combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes
  }else{
    combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes
  }

  XXcomb=combUnique
  XX <- apply(combUnique[, -'i'], 1, function(x) {
    x2 <- rowRanges(m, cols = x)
    x2 = scales::rescale(x2,c(0,1))
    s <- colSums2(x2)
    res=s[1]/s[2]
    return(res)
  })

  data.table::set(XXcomb, j = 'xx', value = XX)
  rez2 <- merge(comb2, XXcomb[, -'i'], by = c('V1', 'V2'), all.x = T)
  data.table::setorder(rez2, i)
  rez2 <- array(rez2$xx, dim = rep(ncol(m), 2))
  rownames(rez2) <- colnames(m)
  colnames(rez2) <- colnames(m)
  rez2=round(rez2,4)
  return(rez2)
}

#' Create Similarity Network based on Average of Pair Weights
#'
#' For each pair of sample's profiles, it determines the average of the feature's values and then
#' it determines their average to assign the score of similarity to the pair of samples.
#'
#' @param m Numeric matrix of column profiles and row features (genes per samples)
#' @return Similarity Network based on Average of Pair Weights
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import stats
#' @importFrom Rfast colmeans
#' @export
build_avg_PWN=function(m, comb2=NULL){
  if(is.null(comb2)){
    samples <- 1:ncol(m)
    comb <- data.table::CJ(samples, samples)
    comb[, i := .I]
    comb <- data.table::melt(comb, 'i')
    data.table::setorder(comb, value)
    v2 <- paste0("V", 1:2)
    comb[, variable2 := v2 , keyby = i]
    comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')
    combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes
  }else{
    combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes
  }

  YY <- apply(combUnique, 1, function(x){
    y <- cbind(m[,x[2]],m[,x[3]])
    r_avgs <- Rfast::colmeans(y)
    res=mean(r_avgs)
    return(res)
  })

  data.table::set(combUnique, j = 'yy', value = YY)
  rez3 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
  data.table::setorder(rez3, i)
  rez3 <- array(rez3$yy, dim = rep(ncol(m), 2))
  rownames(rez3) <- colnames(m)
  colnames(rez3) <- colnames(m)
  rez3=round(rez3,4)
  return(rez3)
}

#' Triangle area between three points
#'
#' One is fixed (0,0), the second is (b_value,0), the third is (0,c_value)
#' Triangle area increases more when two points increase
#' Triangle area penalizes when only one point changes
#'
#' @param b_value Value
#' @param c_value Value
#' @return Triangle area
#' @export
get_AUtri=function(b_value, c_value){
  a=c(0,0)
  b=c(b_value,0)
  c=c(0,c_value)
  res=((a[1]*(b[2]-c[2]))+(b[1]*(c[2]-a[2]))+(c[1]*(a[2]-b[2])))/2
  res=abs(res)
  return(res)
}

#' Quadrilateral area between four points
#'
#' One is fixed (0,0), the second is (b_value,0), the third is (0,c_value), the fourth is (b_value,c_value)
#' Quadrilateral area increases more when three points increase
#' Quadrilateral area penalizes when only one point changes
#'
#' @param b_value Value
#' @param c_value Value
#' @param d_value Value
#' @return Triangle area
#' @export
get_AUquad = function(b_value, c_value, d_value){
  a=c(0,0)
  b=c(b_value,0)
  c=c(0,c_value)
  d=c(b_value,d_value)

  m=rbind(a,b,c,d)
  dists=as.matrix(dist(m))+1
  diag=sqrt(((dists["a","b"])^2)+((dists["a","c"])^2))

  AreaABC=dists["a","b"]*dists["a","c"]*0.5
  s=(diag + dists["c","b"] + dists["c","d"])/2
  AreaBCD=sqrt(s*(s-diag)*(s-dists["c","b"])*(s-dists["c","d"]))
  AreaABCD=(AreaABC+AreaBCD)
  return(AreaABCD)
}

#' Create a Similarity Network based on Adjusted Weighted Jaccard similarity and Average of Pair Weights
#'
#' It creates the SN based on the Adjusted Weighted Jaccard similarity
#' It creates the SN based on the Average of Pair Weights
#' It combines the two values that each pair of samples has
#' The combination is performed computing the triangular area of the two values
#'
#' @param m Numeric matrix s.t. columns are subjects, rows are features and a cell contains a value
#' @return Similarity Network based on Adjusted Weighted Jaccard similarity and Average of Pair Weights
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import stats
#' @export
build_PSNs=function(m){
  samples <- 1:ncol(m)
  comb <- data.table::CJ(samples, samples)
  comb[, i := .I]
  comb <- data.table::melt(comb, 'i')
  data.table::setorder(comb, value)
  v2 <- paste0("V", 1:2)
  comb[, variable2 := v2 , keyby = i]
  comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')

  adjWJ_PSN=build_adjWJ_PSN(m,comb2)
  avg_PWN=build_avg_PWN(m,comb2)

  adjWJ_PSN=scales::rescale(adjWJ_PSN,to = c(0,1))
  avg_PWN=scales::rescale(avg_PWN,to = c(0,1))

  score_PSN_high=adjWJ_PSN
  for(row in 1:nrow(adjWJ_PSN)) {
    for(col in 1:ncol(adjWJ_PSN)) {
      b_value=adjWJ_PSN[row, col]
      c_value=avg_PWN[row, col]
      score_PSN_high[row,col]=get_AUtri(b_value,c_value)
    }
  }
  score_PSN_high = scales::rescale(score_PSN_high, to=c(0,1))
  diag(score_PSN_high)=1

  score_PSN_low=adjWJ_PSN
  for(row in 1:nrow(adjWJ_PSN)) {
    for(col in 1:ncol(adjWJ_PSN)) {
      b_value=adjWJ_PSN[row, col]
      c_value=(1-avg_PWN[row, col])
      score_PSN_low[row,col]=get_AUtri(b_value,c_value)
    }
  }
  score_PSN_low = scales::rescale(score_PSN_low, to=c(0,1))
  diag(score_PSN_low)=1

  PSNs_l=list(PSN_high=score_PSN_high,PSN_low=score_PSN_low)
  return(PSNs_l)
}

