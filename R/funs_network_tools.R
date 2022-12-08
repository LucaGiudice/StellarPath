#' Determines the Power of a Patient Similarity Network
#'
#' Determines if the two groups of nodes in the network separate and how much
#'
#' More a group has the intragroup similarities greater than the opposite intragroup similarities and
#' intergroup similarities and more the power is close to 10
#'
#' In order to understand the magnitude, a low percentile (e.g. 0.2) of the similarities of a group
#' is tested if is greater than a high percentile (e.g. 0.8) of the other two sets of similarities.
#'
#' Lower is the low percentile and higher is the high percentile while still keeping a separability of
#' the similarities and more the power is high.
#'
#' @param PSN Patient Similarity Network in Adjacent matrix format
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @return Dataframe indicating the group which nodes are more similar than the opposite, the power or magnitude measuring
#' how much the strongest group is more similar than the others, the smallest similarity of the strongest group which is greater
#' than the biggest similarity of the opposite group and of the similarities between the groups. Dataframe is NA if
#' the power is 0 and the two groups do not show significant different trends of similarities.
#' @export
#'
get_power=function(PSN,groups){
  min_th=0.350
  max_th=0.650
  min_lvs=seq(min_th,0.1,-0.025)
  max_lvs=seq(max_th,0.9,0.025)
  power_lvs=cbind(min_lvs,max_lvs)
  power_lvs=power_lvs[-c(10),]
  colnames(power_lvs)=c("min","max");rownames(power_lvs)=seq(1,nrow(power_lvs))

  info_df=data.frame(group="NA",power=0,min_intra_STRONG=0,max_intra_WEAK=0,max_inter_similarities=0,stringsAsFactors = F)
  power=best_distr=0;fin_label=label="NA";
  for(lv in 1:nrow(power_lvs)){
    v_info=test_power(PSN=PSN, groups=groups, min_perc=power_lvs[lv,1], max_perc=power_lvs[lv,2])
    label=v_info[1]
    e_distr=v_info[c(2,3,4)]
    if(label!="NA"){
      power=lv
      fin_label=label
      best_distr=as.numeric(e_distr)
    }
    label="NA"
  }
  info=c(fin_label,power,best_distr)
  if(power!=0){
    info_df[1,]=info
    info_df$power=as.numeric(info_df$power)
    info_df$min_intra_STRONG=as.numeric(info_df$min_intra_STRONG)
    info_df$max_intra_WEAK=as.numeric(info_df$max_intra_WEAK)
    info_df$max_inter_similarities=as.numeric(info_df$max_inter_similarities)
  }
  return(info_df)
}

#' Test if a power holds in a Patient Similarity Network
#'
#' Given a low and a high percentile, it tests if the low percentile of the distribution of the intragroup
#' similarities of one group is greater than the high percentile of the distributions of:
#' - the intragroup of similarities of the opposite group
#' - the intergroup similarities
#'
#' @param PSN Patient Similarity Network in Adjacent Matrix format
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param min_perc Double value between 0 and 1 indicating the low percentile
#' @param max_perc Double value between 0 and 1 indicating the high percentile
#' @return Vector with the label of the strongest group, the value of the low percentile of the intragroup similarities of the strongest group,
#' the value of the high percentile of the intragroup similarities of the weak group and of the intergroup similarities
#' @export
#'
test_power=function(PSN, groups, min_perc=0.3, max_perc=0.7){
  #Get indexes of the two groups in the PSN
  A_label=unique(groups)[1]
  B_label=unique(groups)[2]
  indxsAs=grep(A_label,groups)
  indxsBs=grep(B_label,groups)

  #Get similarities
  AA_PSN=PSN[indxsAs,indxsAs]
  BB_PSN=PSN[indxsBs,indxsBs]
  AB_PSN=PSN[indxsAs,indxsBs]

  #Remove similarities of non-existing edges
  sims_AA=AA_PSN[upper.tri(AA_PSN)]
  sims_AA=sims_AA[sims_AA!=0]
  sims_BB=BB_PSN[upper.tri(BB_PSN)]
  sims_BB=sims_BB[sims_BB!=0]
  sims_AB=AB_PSN
  sims_AB=sims_AB[sims_AB!=0]

  if(length(sims_AA)==0){sims_AA=0}
  if(length(sims_BB)==0){sims_BB=0}
  if(length(sims_AB)==0){sims_AB=0}

  #Measure stats
  AAe=quantile(sims_AA,probs = c(min_perc,max_perc))
  BBe=quantile(sims_BB,probs = c(min_perc,max_perc))
  ABe=quantile(sims_AB,probs = c(min_perc,max_perc))

  #Determine if one group is signature
  label="NA"
  if(AAe[1]>BBe[2]){
    if(AAe[1]>ABe[2]){
      label=c(A_label,AAe[1],BBe[2],ABe[2])
    }
  }

  if(BBe[1]>AAe[2]){
    if(BBe[1]>ABe[2]){
      label=c(B_label,BBe[1],AAe[2],ABe[2])
    }
  }
  return(label)
}

#' Adjust intragroup similarities based on intergroup ones
#'
#' Considers the similarity between two patients belonging to the same group, it increases this similarity if
#' they have low interconnections with members of other groups, it lowers down this similarity if the opposite
#'
#' @param PSN Patient Similarity Network in Adjacent Matrix format
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @return Return the patient similarity with the weights between patients adjusted
#' @export
#'
adjust_intra_sims=function(PSN, groups){
  adj_PSN=PSN
  for(gri in 1:2){
    if(gri==1){
      indxs_gr1=which(groups==unique(groups)[1])
      indxs_gr2=which(groups==unique(groups)[2])
    }else{
      indxs_gr2=which(groups==unique(groups)[1])
      indxs_gr1=which(groups==unique(groups)[2])
    }

    for(ci in 1:length(indxs_gr1)){
      p1=indxs_gr1[ci]
      indxs_ci_intra=setdiff(indxs_gr1,ci)
      for(cii in 1:length(indxs_ci_intra)){
        p2=indxs_ci_intra[cii]
        cell=PSN[p1,p2]
        p1_inter=mean(PSN[p1,indxs_gr2])
        p2_inter=mean(PSN[p2,indxs_gr2])
        score=get_AUquad(cell,1-p1_inter,1-p2_inter)
        adj_PSN[p1,p2]=adj_PSN[p2,p1]=score
      }
    }
  }
  return(adj_PSN)
}

#' Removes noise in a Patient Similarity Network
#'
#' For each group of nodes, it removes the cliques having the lowest average of weights
#'
#' @param PSN Patient Similarity Network in Adjacent Matrix format
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param n_comb Default 3, integer defining the size of the clique
#' @param dir_sign Default 1, internal parameter of the software
#' @param perc_filter Default 0.05, double between 0 and 1 defining the amount of cliques to remove because outliers
#' @return Return a sparse patient similarity network with the similarities of the lowest cliques equal to 0
#' @importFrom Rfast comb_n
#' @importFrom utils combn
#' @export
#'
remove_outliers = function(PSN, groups, n_comb=3, dir_sign=1, perc_filter=0.05){
  #Adjust intragroup similarities by intergroup similarities
  aPSN=adjust_intra_sims(PSN=PSN, groups=groups)
  v=as.vector(aPSN)
  n=nrow(aPSN)
  m=ncol(aPSN)

  #Determine the percentage of clique to filter for the two groups
  if(dir_sign==1){
    perc_filter1=perc_filter
    perc_filter2=1-perc_filter
  }else{
    perc_filter2=perc_filter
    perc_filter1=1-perc_filter
  }

  #Find indexes of each group
  indxs_gr1=which(groups==unique(groups)[1])
  indxs_gr2=which(groups==unique(groups)[2])

  #Create clique of n_comb grade
  combs_gr1=Rfast::comb_n(indxs_gr1,n_comb,simplify = T)
  combs_gr2=Rfast::comb_n(indxs_gr2,n_comb,simplify = T)

  #Compute clique weights
  v_score=apply(combs_gr1,2,function(comb){
    #Determine the similarities of a clique
    edges=combn(comb,2)

    avg_edges=mean(apply(edges,2,function(x){
      i=x[1]
      j=x[2]
      v_value=v[(n * (j-1)) + i]
      return(v_value)
    }))

    return(avg_edges)
  })

  combs_gr1=rbind(combs_gr1,v_score)

  #Compute clique weights
  v_score=apply(combs_gr2,2,function(comb){
    #Determine the similarities of a clique
    edges=combn(comb,2)

    avg_edges=mean(apply(edges,2,function(x){
      i=x[1]
      j=x[2]
      v_value=v[(n * (j-1)) + i]
      return(v_value)
    }))

    return(avg_edges)
  })

  combs_gr2=rbind(combs_gr2,v_score)

  #Select the clique to filter
  if(dir_sign==1){
    combs_gr1[4,combs_gr1[4,]<=quantile(combs_gr1[4,],perc_filter1)]=0
    combs_gr2[4,combs_gr2[4,]>=quantile(combs_gr2[4,],perc_filter2)]=0
  }else{
    combs_gr1[4,combs_gr1[4,]>=quantile(combs_gr1[4,],perc_filter1)]=0
    combs_gr2[4,combs_gr2[4,]<=quantile(combs_gr2[4,],perc_filter2)]=0
  }

  #Apply the filter to the original PSN
  fPSN=PSN
  combs=cbind(combs_gr1,combs_gr2)
  combs=combs[,combs[4,]==0]
  for(ri in 1:ncol(combs)){
    comb=combs[,ri]
    edges=t(combn(comb[1:n_comb],2))
    for(ei in 1:nrow(edges)){
      fPSN[edges[ei,1],edges[ei,2]]=fPSN[edges[ei,2],edges[ei,1]]=0
    }
  }
  return(fPSN)
}

#' Determine a cohesive centrality score for each node a patient similarity network
#'
#' Function which computes how much a node is connected to the members of its group, its non-members and then
#' it produces a centrality score. Higher is the score and more a node is strongly connected to the members of its
#' group and lowly connected to the non-members.
#'
#' @param PSN Patient Similarity Network in Adjacent Matrix format
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @return Matrix of values, first column indicates the median of the intraclass similarities,
#' @importFrom Rfast rowMedians
#' @importFrom scales rescale
#' @export
#'
get_cohesive_score = function(PSN, groups){
  #Find indexes of each group
  seq1cl=which(groups==unique(groups)[1])
  seq2cl=which(groups==unique(groups)[2])
  #Compute how much a node is connected inside its group
  INw=c(Rfast::colMedians(PSN[seq1cl,seq1cl]),Rfast::colMedians(PSN[seq2cl,seq2cl]))
  #Compute how much a node is connected outside its group
  ABw=c(1-Rfast::colMedians(PSN[seq2cl,seq1cl]), 1-Rfast::colMedians(PSN[seq1cl,seq2cl]))
  #Finalize the vectors with the names of the nodes and create a matrix
  names(ABw)=names(INw)=colnames(PSN)
  eigs=cbind(INw,ABw)
  #Rank scale
  eigs=apply(eigs,2,rank)

  #Adjusted score of eigencentrality
  eigs=do_standard(eigs)
  eigs_adj <- apply(eigs, 1, function(x){
    get_AUtri(x[1],x[2])
  })
  eigs_adj = scales::rescale(eigs_adj, to=c(0,1))

  #Update
  eigs=cbind(eigs,eigs_adj)
  eigs[,2]=1-eigs[,2]
  colnames(eigs)=c("avg_intra_sim","avg_inter_sim","centrality")
  return(eigs)
}

#' Converts adjacency matrix to edge list
#'
#' Converts a graph represented as adjacency matrix to edge list
#'
#' @param m adjacency matrix of a graph
#' @param symmetric default TRUE, the adjacency matrix is symmetric
#' @param keep_diagonal default FALSE, do not keep the diagonal values in the edge list
#' @param text default FALSE, do not set the column of the edge list as character vector
#' @return Return edge list of the original network as matrix
#' @importFrom reshape2 melt
#' @export
#'
convert_adj2edg<-function(m, symmetric=TRUE, keep_diagonal=FALSE, text=FALSE){
  m<-as.matrix(m)
  if(symmetric){m[lower.tri(m)]=NA} # use to allow missing values
  if(!keep_diagonal){diag(m)=NA}
  el<-reshape2::melt(m)
  colnames(el)<-c("source","target","weight")
  el<-el[!is.na(el$weight),]
  el<-el[el$weight!=0,]
  if(!text){el$weight<-as.numeric(as.character(el$weight))}
  return(el)
}

#' Converts edge list to adjacency matrix
#'
#' Converts a graph represented as edge list to adjacency matrix
#'
#' @param el edge list of a graph
#' @return Return djacency matrix of the original network as matrix
#' @importFrom igraph graph.data.frame
#' @importFrom igraph get.adjacency
#' @export
#'
convert_edg2adj<-function(el){
  mygraph <- igraph::graph.data.frame(el,directed = F)
  adj_mat=igraph::get.adjacency(mygraph, sparse = FALSE, attr=colnames(el)[3])
  return(adj_mat)
}
