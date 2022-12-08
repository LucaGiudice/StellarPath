#' Find which pathway/set of elements are represented by a significant similarity network
#'
#' Remove the pathways/sets which size is bigger than max_size.
#' For each set of elements (pathways_l), keep only the significant elements (sign_df).
#' Keep the (n_top_pathways) having the highest number of significant elements (sign_df).
#' For each set of significant elements (pathways_l), use their values (counts) to compute how much each pair
#' of samples is similar and build the similarity network.
#' Compute the Power of separation for each network.
#'
#' @param counts count matrix of samples at the columns and molecules at the row (e.g. gene expression matrix)
#' @param groups character or factor vector of labels s.t. each label indicates the belonging of a sample to one group
#' @param pathways_l list of pathways or generically set of elements
#' @param sign_df dataframe s.t. rownames are molecules's names belonging to the list of pathways, column logFoldChange
#' with double values indicating log fold change of the molecule's count value between two groups, column Significant
#' which is TRUE or FALSE based on if the row molecule has been tested and resulted signficant, column Barcode has a
#' 0 value for not significant molecules and a 9 for significant values.
#' @param max_size Default 150, integer value filtering out the pathways/sets with more elements than this threshold
#' @param min_elements Default 4, integer value filtering out the pathways/sets containing less than this threshold of significant molecules
#' @param n_top_sets Default 50, integer that set how many best sets/PSNs to retrieve
#' @param n_cores Default 2, integer that set the number of cores to use for running the function in parallel
#' @param keep_PSN Default FALSE, do not keep the PSN generated from a set of element
#' @return list of two elements: top_pathways_l and top_pathways_df.
#' top_pathways_df details the pathways/sets which have a significant PSN.
#' It describes them with the following column features:
#' - pathway name
#' - freq refers to the number of significant molecules/elements belonging to the pathway
#' - regulation refers to the direction of enrichment, in case the molecules belonging to the pathway have count values
#' which are greater in the strongest group than in the opposite group of the network then it is "activated"
#' - groups indicates the group of nodes being significantly more similar than the opposite group in comparison
#' - power indicates from 1 to 10 how much the strongest group is more similar than the non-members
#' - min_intra_STRONG is the minimum similarity between the members of the strongest group which is greater than the opposite group's intrasimilarities and the intergroup similarities
#' - max_intra_WEAK is the maximum similarity between the members of the weakest group which is lower than the opposite group's intrasimilarities and the intergroup similarities
#' - max_inter_similarities is the maximum similarity between members of the two groups
#' While, top_pathways_l contains the data associated to each pathway/set of element that has been tested.
#' - set contains the significant molecules/elements belonging to the pathway/set
#' - set_counts is the count matrix (e.g. genes per samples) with only the count values of the significant molecules/elements belonging to the pathway/set
#' - stats is a dataframe indicating the statistics of the pathway/set
#' - centrality is a dataframe s.t. for each sample indicates how much is strongly similar in its groups vs its not group
#' - orig_set contains all the original elements of the pathway/set
#' @import matrixStats
#' @import data.table
#' @import scales
#' @import stats
#' @import Rfast
#' @import doParallel
#' @import foreach
#' @export
#'
find_PSN = function(counts, groups, pathways_l, sign_df, max_size=200, min_elements=6, n_top_pathways=100, n_cores=2, keep_PSN=FALSE){
  groups_name=unique(groups)

  #Keep pathways/sets having a name that does not include a word of the ban list
  np=names(pathways_l)
  data("ban_list")
  ban_names <- unique(grep(paste(ban_list,collapse="|"), np, value=FALSE, ignore.case = TRUE))
  if(length(ban_names)>0){
    pathways_l=pathways_l[-ban_names]
  }

  #Keep pathways/sets having a number of elements lower than max_size
  keep=sapply(pathways_l,length)<=max_size
  pathways_l=pathways_l[keep]
  rm(max_size)

  #Get which kind of significant molecules are present:
  #Type 9: molecules are used collectively to understand which pathway fit them the most
  #Type 1: molecules which have a higher role and more cohesive in group 1 than group 2 are used
  #Type 2: molecules which have a higher role and more cohesive in group 2 than group 1 are used
  #Type -1: molecules which have a lower role and more cohesive in group 1 than group 2 are used
  #Type -2: molecules which have a lower role and more cohesive in group 1 than group 2 are used
  sign_df=sign_df[sign_df$Barcode!=0,]
  sign_types=unique(sign_df$Barcode)

  #Keep for each type-specific molecules, the pathways/sets which fit them the most
  top_pathways_l=list()
  for(si in 1:length(sign_types)){
    #Get the type of DR to use
    sign_type=sign_types[si]
    query_elements=rownames(sign_df[sign_df$Barcode==sign_type,])

    check_exist_pathways=T
    min_elements1=min_elements
    while(check_exist_pathways & min_elements1>2){
      #Keep pathways with at least n elements
      feat_freq=sapply(pathways_l,function(x){sum(query_elements %in% x)})
      keep=feat_freq>=min_elements1
      feat_freq=feat_freq[keep]
      query_pathways_l=pathways_l[keep]

      #Keep n pathways with the most number of DR elements
      sel_path_l=query_pathways_l
      sel_path_l=lapply(sel_path_l,function(x){x[x %in% query_elements]})
      sel_path_l=sel_path_l[sapply(sel_path_l,length)!=0]

      if(length(sel_path_l)<=10){
        min_elements1=min_elements1-1
      }else{
        check_exist_pathways=FALSE
      }
    }

    if(length(sel_path_l)>n_top_pathways){
      feat_lfc=abs(sapply(sel_path_l,function(x){quantile(sign_df$logFoldChange[rownames(sign_df)%in%x], probs=0.9)}))
      feat_sta=abs(sapply(sel_path_l,function(x){quantile(sign_df$Stability[rownames(sign_df)%in%x], probs=0.9)}))
      new_ord=order(feat_lfc,feat_sta,decreasing = TRUE)
      sel_path_l=sel_path_l[new_ord]
      sel_path_l=sel_path_l[1:n_top_pathways]
    }

    rm(check_exist_pathways,query_pathways_l,min_elements1)

    if(length(sel_path_l)==0){
      top_pathways_df=data.frame(pathways=NA,freq=0,sign_type=0)
      stats=data.frame(regulation=NA,group=NA,power=0,min_intra_STRONG=0,max_intra_WEAK=0,max_inter_similarities=0)
      top_pathways_df=cbind(top_pathways_df,stats)
      PSNs_l=list()
    }else{
      #Split the pathway list into chunks ----
      n_pathways=length(sel_path_l)
      pathways_indxs=seq(1,n_pathways)
      max <- n_pathways/n_cores
      x <- seq_along(pathways_indxs)
      pathway_sets_indxs <- split(pathways_indxs, ceiling(x/max))

      #Create dataframe with the name of the top selected pathways and the number of elements inside
      top_pathways_df=data.frame(pathways=names(sel_path_l),freq=sapply(sel_path_l,length),sign_type=sign_type)
      rownames(top_pathways_df)=seq(1,nrow(top_pathways_df))
      rm(x,max,n_pathways,pathways_indxs)

      #One core one chunk of pathways to process ----
      if(.Platform$OS.type == "unix") {
        #message(">>parallel for linux")
        cl <- makeCluster(n_cores,type="FORK");
      } else {
        #message(">>parallel for windows")
        cl <- makeCluster(n_cores);
      }
      registerDoParallel(cl);

      PSNs_l = list()
      PSNs_l = foreach(
        k_chunk = 1:length(pathway_sets_indxs),
        .inorder = FALSE,
        .noexport = c(),
        .verbose = F,
        #.errorhandling = "pass",
        .packages = c("matrixStats", "data.table", "scales", "stats","Rfast"),
        .export = c("build_PSNs", "get_power","remove_outliers","get_cohesive_score")
      ) %dopar%
      {

        #Get the indexes of the chunck of pathways that the core has to process
        pathway_set_indxs=pathway_sets_indxs[[k_chunk]]
        PSNs_chunk_l=list()

        #Score pathways
        for(qi in pathway_set_indxs){
          gs=sel_path_l[[qi]]
          name_pathway=names(sel_path_l)[qi]
          m_gs=counts[gs,]
          PSNs_l=build_PSNs(m_gs)
          best_stats=data.frame(regulation=NA,group=NA,power=0,min_intra_STRONG=0,max_intra_WEAK=0,max_inter_similarities=0)

          stats_PSN_high=get_power(PSNs_l$PSN_high, groups=groups)
          stats_PSN_low=get_power(PSNs_l$PSN_low, groups=groups)

          if(sign_type==9){
            if(stats_PSN_high$power==0 & stats_PSN_low$power==0){
              PSNs_l$PSN_high=remove_outliers(PSN=PSNs_l$PSN_high, groups=groups, dir_sign = sign_type)
              PSNs_l$PSN_low=remove_outliers(PSN=PSNs_l$PSN_low, groups=groups, dir_sign = abs(sign_type))
            }
          }else{
            if(sign_type>0 & stats_PSN_high$power==0){
              PSNs_l$PSN_high=remove_outliers(PSN=PSNs_l$PSN_high, groups=groups, dir_sign = sign_type)
            }
            if(sign_type<0 & stats_PSN_low$power==0){
              PSNs_l$PSN_low=remove_outliers(PSN=PSNs_l$PSN_low, groups=groups, dir_sign = abs(sign_type))
            }
          }

          stats_PSN_high=get_power(PSNs_l$PSN_high, groups=groups)
          stats_PSN_low=get_power(PSNs_l$PSN_low, groups=groups)

          if(sign_type==9){
            if(stats_PSN_high$power>=(stats_PSN_low$power-1) & stats_PSN_high$power!=0){
              best_stats=cbind(regulation="activated",stats_PSN_high)
              best_PSN=PSNs_l$PSN_high
              centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
            }
            if(stats_PSN_low$power>stats_PSN_high$power & stats_PSN_low$power!=0){
              best_stats=cbind(regulation="inhibited",stats_PSN_low)
              best_PSN=PSNs_l$PSN_low
              centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
            }
          }else{
            if(sign_type>0 & stats_PSN_high$power!=0){
              if(sign_type==1 & stats_PSN_high$group==groups_name[1]){
                best_stats=cbind(regulation="activated",stats_PSN_high)
                best_PSN=PSNs_l$PSN_high
                centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
              }
              if(sign_type==2 & stats_PSN_high$group==groups_name[2]){
                best_stats=cbind(regulation="activated",stats_PSN_high)
                best_PSN=PSNs_l$PSN_high
                centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
              }
            }

            if(sign_type<0 & stats_PSN_low$power!=0){
              if(sign_type==-1 & stats_PSN_low$group==groups_name[1]){
                best_stats=cbind(regulation="inhibited",stats_PSN_low)
                best_PSN=PSNs_l$PSN_low
                centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
              }
              if(sign_type==-2 & stats_PSN_low$group==groups_name[2]){
                best_stats=cbind(regulation="inhibited",stats_PSN_low)
                best_PSN=PSNs_l$PSN_low
                centrality=get_cohesive_score(PSN=best_PSN, groups=groups)[,3]
              }
            }
          }

          if(best_stats$power==0){
            best_PSN=matrix(0,2,length(groups))
            centrality=rep(0,length(groups))
            names(centrality)=names(m_gs)
          }

          centrality=as.data.frame(t(centrality))
          rownames(centrality)=name_pathway
          colnames(centrality)=colnames(best_PSN)

          if(keep_PSN==FALSE){
            PSNs_chunk_l[[name_pathway]]=list(set=gs, set_counts=m_gs, stats=best_stats, centrality=centrality, orig_set=pathways_l[[name_pathway]])
          }else{
            PSNs_chunk_l[[name_pathway]]=list(set=gs, set_counts=m_gs, stats=best_stats, PSN=best_PSN, centrality=centrality, orig_set=pathways_l[[name_pathway]])
          }
        }

        return(PSNs_chunk_l)

      }
      parallel::stopCluster(cl)

      #Format resulting list
      PSNs_l=do.call(c, PSNs_l)
      PSNs_l=PSNs_l[match(top_pathways_df$pathways,names(PSNs_l))]

      #Extract stats about the PSNs created for each corresponding pathway
      stats=lapply(PSNs_l,function(x){x$stats})
      stats=as.data.frame(data.table::rbindlist(stats))
      rownames(stats)=names(PSNs_l)

      #Merge pathway meta information
      top_pathways_df=cbind(top_pathways_df,stats)
    }

    if(sum(is.na(top_pathways_df$pathways))==nrow(top_pathways_df)){
      PSNs_l=list(NAN="NAN")
      names(PSNs_l)[1]=paste("run",sign_type,sep="")
    }

    if(si==1){
      top_pathways_dfs=top_pathways_df
      top_pathways_l=PSNs_l
    }else{
      top_pathways_dfs=rbind(top_pathways_dfs,top_pathways_df)
      top_pathways_l=c(top_pathways_l,PSNs_l)
    }

  }

  #Remove pathways which resulted not significant
  keep=(!is.na(top_pathways_dfs$pathways))
  top_pathways_dfs=top_pathways_dfs[keep,]
  top_pathways_l=top_pathways_l[keep]

  keep=top_pathways_dfs$power!=0
  top_pathways_dfs=top_pathways_dfs[keep,]
  top_pathways_l=top_pathways_l[keep]

  if(sum(duplicated(top_pathways_dfs$pathways))>0){
    #Order by power and minimum similarity for the signature class
    new_ord=order(top_pathways_dfs$power,top_pathways_dfs$min_intra_STRONG,decreasing = T)
    top_pathways_dfs=top_pathways_dfs[new_ord,]
    top_pathways_l=top_pathways_l[new_ord]

    #Find duplicated pathways
    keep=!duplicated(top_pathways_dfs$pathways)
    top_pathways_dfs=top_pathways_dfs[keep,]
    top_pathways_l=top_pathways_l[keep]
  }

  new_order=order(top_pathways_dfs$power,decreasing = T)
  top_pathways_dfs=top_pathways_dfs[new_order,]
  top_pathways_l=top_pathways_l[new_order]
  if(length(new_order)!=0){
    rownames(top_pathways_dfs)=seq(1,nrow(top_pathways_dfs))
  }

  #Return
  res_l=list(top_pathways_l=top_pathways_l,top_pathways_df=top_pathways_dfs)
  return(res_l)
}

#' Combines multiple results of analysis of training omics
#'
#' Combines multiple results of analysis of training omics
#'
#' @param l2merge list to merge s.t. each element of the list includes one top_pathways_l and one top_pathways_df
#' @return list of two elements: top_pathways_l and top_pathways_df
merge_multi_omics_analysis = function(l2merge){
  l2merge=l2merge[sapply(l2merge,length)!=0]
  res=list()
  if(length(l2merge)!=0){
    top_pathways_l=list()
    for(k in 1:length(l2merge)){
      li=l2merge[[k]]
      top_pathways_l=c(top_pathways_l,li$top_pathways_l)
      if(k==1){
        top_pathways_df=li$top_pathways_df
      }else{
        top_pathways_df=rbind(top_pathways_df,li$top_pathways_df)
      }
    }

    new_order=order(top_pathways_df$power,decreasing = T)
    top_pathways_df=top_pathways_df[new_order,]
    top_pathways_l=top_pathways_l[new_order]

    new_names=make.unique(top_pathways_df$pathways,sep=".multi")
    top_pathways_df$pathways=new_names
    names(top_pathways_l)=new_names

    if(nrow(top_pathways_df)>1){
      rownames(top_pathways_df)=seq(1,nrow(top_pathways_df))
    }

    res=list(top_pathways_l=top_pathways_l,top_pathways_df=top_pathways_df)
  }
  return(res)
}

