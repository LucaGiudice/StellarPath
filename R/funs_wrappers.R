#' Check and Prepare the input data
#'
#' This function checks that the input data satisfy all the criteria for being used in the software
#' It creates an object list that includes all the checked and prepared user data
#'
#' @param info Dataframe, two character columns, first has sample's IDs (have to match with column names of the count matrices),
#' second has sample's groups (only two as a pairwise DE analysis, for example c(AD,AD,AD,HEALTHY,HEALTHY,HEALTHY))
#' @param groups_name Character vector of length two, it must contain the same labels of unique(info[,2]), tip: first the case group like AD,
#' while the second name refers to the control group like HEALTHY
#' @param mRNA Numeric matrix of raw count values from mRNA sequencing (e.g. genes x samples)
#' @param miRNA Numeric matrix of raw count values from miRNA sequencing (e.g. miRNA x samples)
#' @param mutations Numeric matrix of binary values from somatic mutation data (e.g. gene X samples, 1 per mutated gene, 0 otherwise)
#' @param tax_id Default 9606, integer value indicating taxonomic id of the samples, either 9606 (human) or 10090 (mouse)
#' @param simplify_names Default TRUE, simplify names of the samples in the info and count matrices
#' @param cv_probs Default c(0.8,0.1,0.1), vector of doubles indicating the proportion of the samples to use as training, validation and testing set
#' @param n_iter_cv Default 2, integer value indicating the number of iterations to perform the cross validation
#' @param n_cores Default 2, integer value greater than 0 indicating the number of cores to use to parallelize and speed up the software's operations
#' @param seed Default 5, integer value indicating the seed to replicate the sampling of the cross validation
#' @param winz Default FALSE, not adjusting outliers by Winsorizing
#' @param min Default 0.1, double value indicating the minimum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and max value
#' @param max Default 0.9, double value indicating the maximum value of the distribution for the Winsorizing,
#' it must be smaller than 1 and greater than min value
#' @param type Default expression, character string which defines the type of transformation applied to the count data, options in the help of the get_rank_01 function
#' @return A list s.t. each element contains the user-input prepared data for a run of cross-validation.
#' The run's data are divided into training and testing sets.
#' Each set contains a specific subset of sample's information (info) and count values (omics).
#' @importFrom utils combn
#' @importFrom splitTools partition
#' @importFrom miRBaseConverter miRNA_PrecursorToMature
#' @importFrom stats aggregate
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph ego
#' @importFrom igraph as_ids
#' @export
#'
prepare_data = function(info, groups_name=NULL, mRNA=NULL, miRNA=NULL, mutations=NULL,
                        tax_id=9606, simplify_names=TRUE,
                        cv_probs=c(0.7,0.1,0.2), n_iter_cv=2, n_cores=2, seed=5,
                        winz=FALSE, min=0.1, max=0.9, type="total_rank"){

  #Initiate data variables
  dataset=list()
  omics=list()

  #Checking user provided omics
  message("Checking input omics/raw count matrices")
  if(is.null(mRNA) & is.null(mutations)){
    stop("The mRNAseq (genes per samples) or somatic mutations matrix (genes per samples) is needed");
  }

  if(!is.null(mRNA)){
    message(">the mRNAseq (mRNA per samples) has been provided");
    omics[["mRNA"]]=mRNA
    rm(mRNA)
  }

  if(!is.null(miRNA)){
    message(">the miRNAseq (miRNA per samples) has been provided");
    omics[["miRNA"]]=miRNA
    if(is.null(omics[["mRNA"]])){
      stop("The mRNAseq (genes per samples) in needed to link miRNAs to targets and pathways");
    }
    rm(miRNA)
  }

  if(!is.null(mutations)){
    omics[["mutations"]]=mutations
    message(">the mutations (genes per samples) have been provided");
    rm(mutations)
  }

  #Checking user provided info about samples
  message("Checking input info")
  if(is.data.frame(info)){
    message(">dataframe format")
    if(ncol(info)==2){
      message(">two columns")
      if(sum(is.na(info))==0){
        message(">does NOT contain NA values")
        groups_freq=as.data.frame(table(info[,2]))
        if(nrow(groups_freq)==2){
          message(">groups in comparison are two")
          if(!is.null(groups_name)){
            if(length(groups_name)!=2){
              stop("hey, the vector of names related to the two groups in comparison is longer than two");
            }
            if(sum(groups_freq$Var1 %in% groups_name) == 2){
              message(">second column of info contains only labels which belong to the user defined group's names")
            }else{
              stop("second column of info dataframe contains labels which don't match with user defined group's names");
            }
          }
        }else{
          stop("the groups are more than two");
        }
      }else{
        stop("info dataframe contains NA values");
      }
    }else{
      stop("info dataframe has more than two columns");
    }
  }else{
    stop("sample's info is NOT a dataframe");
  }
  rm(groups_freq)

  #Checking matching samples info and omics
  message("Checking same sample's IDs between info first column and omics column names")
  pat_ids=list()
  pat_ids[["info"]]=info[,1]
  pat_ids=c(pat_ids,lapply(omics,function(x){colnames(x)}))
  pat_ids_freqs=as.data.frame(table(unlist(pat_ids)))
  check=sum(pat_ids_freqs$Freq!=length(pat_ids))
  if(check==0){
    pairs=t(utils::combn(1:length(pat_ids), 2))
    for(i in 1:nrow(pairs)){
      bin=check_vectors(x=pat_ids[[pairs[i,1]]], y=pat_ids[[pairs[i,2]]])
      if(bin==0){
        stop("samples are missing or are in different order between info dataframe and one omic/raw count matrix");
      }
    }
    message(">same order of existing samples between info dataframe and omics data")
  }else{
    stop("samples are missing between info dataframe and omics data");
  }
  colnames(info)=c("IDs","Groups")
  rm(pat_ids,pat_ids_freqs,pairs,bin)

  #Checking and Prepare format of omics
  message("Checking format of input omics")
  class_omics=as.data.frame(table(sapply(omics,function(x){class(x)})[1,]))
  if(class_omics$Freq[1]==length(omics)){
    message(">omics are in the same format")
    if(class_omics$Var1[1]=="matrix"){
      message(">omics are matrices")
    }else{
      stop("omics are NOT matrices");
    }
  }else{
    stop("omics are NOT all matrices");
  }
  rm(class_omics)

  message("Ordering samples based on groups")
  if(is.null(groups_name)){
    ord=order(info$Groups,decreasing = TRUE)
    info=info[ord,]
    rownames(info)=info$IDs
    omics=lapply(omics,function(x){
      x=x[,ord]
      return(x)
    })
  }else{
    ord=order(info$Groups,decreasing = TRUE)
    tmp_grs=info$Groups[ord]
    if(tmp_grs[1]==groups_name[1]){
      info=info[ord,]
      rownames(info)=info$IDs
      omics=lapply(omics,function(x){
        x=x[,ord]
        return(x)
      })
    }else{
      ord=order(info$Groups,decreasing = FALSE)
      info=info[ord,]
      rownames(info)=info$IDs
      omics=lapply(omics,function(x){
        x=x[,ord]
        return(x)
      })
    }
    groups_freq=as.data.frame(table(info$Groups))
    message(">The comparison and classification is with:")
    rownames(groups_freq)=groups_freq[,1]
    colnames(groups_freq)=c("group","freq")
    message(">>first group : ",groups_freq[groups_name[1],1]," of samples: ",groups_freq[groups_name[1],2])
    message(">>second group : ",groups_freq[groups_name[2],1]," of samples: ",groups_freq[groups_name[2],2])
  }
  rm(ord,groups_freq)

  if(simplify_names){
    message("Simplifying sample's names in info and omics data")
    info$IDs=make.unique(info$Groups)
    omics=lapply(omics,function(x){
      colnames(x)=info$IDs
      return(x)
    })
    message(">>from: ",rownames(info)[2]," to: ",info$IDs[2])
    message(">>original sample's IDs remain as rownames of the info dataframe")
  }

  #Load and set databases ----
  if(tax_id==9606 | tax_id==10090){
    if(tax_id==9606){
      message("Loading human databases")
      data("dbs_hsa_net_l")
      data("dbs_hsa_paths_l")
      dbs_l=list(HUMAN=c(dbs_hsa_paths_l,dbs_hsa_net_l))
      databases=list(mirna_net=dbs_l$HUMAN$MIRWALK,
                     lncrna_net=dbs_l$HUMAN$RNARNA,
                     pathway_db=dbs_l$HUMAN$PATHWAYS,
                     gene_net=dbs_l$HUMAN$NETWORK)
    }else{
      message("Loading mouse databases")
      data("dbs_mmu_net_l")
      data("dbs_mmu_paths_l")
      dbs_l=list(MOUSE=c(dbs_mmu_paths_l,dbs_mmu_net_l))
      databases=list(mirna_net=dbs_l$MOUSE$MIRWALK,
                     lncrna_net=dbs_l$MOUSE$RNARNA,
                     pathway_db=dbs_l$MOUSE$PATHWAYS,
                     gene_net=dbs_l$MOUSE$NETWORK)
    }
  }else{
    stop("the software does not support other organisms out of human and mouse: choose 9606 or 10090");
  }

  #Cleaning and Preparing omics
  message("Processing and Normalizing omics")
  #Process mRNA -----
  if(!is.null(omics[["mRNA"]])){
    #Remove duplicated rows
    message(">preparing mRNA")
    mRNA=omics[["mRNA"]]
    if(sum(duplicated(rownames(mRNA)))!=0){
      mRNA=stats::aggregate(miRNA, by=list(row.names(mRNA)), FUN=median)
      rownames(mRNA)=mRNA[,1];mRNA=mRNA[,-1]
      mRNA=as.matrix(mRNA)
    }

    #Divide the mRNA in the gene expression and in lncRNA expression
    lncRNA_ids=unique(c(names(databases$lncrna_net$lncrna_gene_targ_l),
                        names(databases$lncrna_net$lncrna_mirna_targ_l)))
    omics=omics[-grep("mRNA",names(omics))]
    gex=mRNA[!(rownames(mRNA) %in% lncRNA_ids),]
    lncRNA=mRNA[rownames(mRNA) %in% lncRNA_ids,]

    #Check if you have proper gene and lncRNA names
    gXp=sum(rownames(gex) %in% unique(unlist(databases$pathway_db$pathways_l)))
    if(gXp>0){
      message(">>number of genes falling into pathways: ", gXp)
    }else{
      message(">>gene name example: ",databases$pathway_db$pathways_l[[1]][1])
      stop("no gene as rowname of the mRNA count matrix belongs to any literature pathway");
    }

    lXp=sum(rownames(lncRNA) %in% unique(names(databases$lncrna_net$lncrna_gene_targ_l)))
    if(lXp>0){
      message(">>number of lncRNAs having target sets of genes: ", lXp)
    }else{
      message(">>lncRNA name example: ",unique(names(databases$lncrna_net$lncrna_gene_targ_l))[1])
      stop("no lncRNA as rowname of the mRNA count matrix has a literature target set of genes");
    }

    #Library and Rank normalization
    omics[["gex"]]=prepare_RNAseq(counts=gex, groups=info$Groups, winz=winz, min=min, max=max, type=type)
    omics[["lncRNA"]]=prepare_RNAseq(counts=lncRNA, groups=info$Groups, winz=winz, min=min, max=max, type=type)
    rm(mRNA,gex,lncRNA,gXp,lXp)
  }

  #Process miRNA -----
  if(!is.null(omics[["miRNA"]])){
    message(">preparing miRNA")
    miRNA=omics[["miRNA"]]

    #Finds if there are pre miRNA
    mirnas=rownames(miRNA)
    mirna_map=miRBaseConverter::miRNA_PrecursorToMature(mirnas, version = "v22")
    del=is.na(mirna_map$Mature1) & is.na(mirna_map$Mature2)
    mirna_map=mirna_map[!del,]

    #Add the pre miRNAs as Mature if they are not already
    miRNA_mat=miRNA[1:2,]
    rownames(miRNA_mat)=c("start","start")
    mirna_pre=mirna_map$OriginalName[1]
    for(mirna_pre in mirna_map$OriginalName){
      ex=miRNA[mirna_pre,]
      mirnas_mat=unlist(mirna_map[mirna_map$OriginalName==mirna_pre,c(2,3)])

      miRNA_mat1=t(matrix(ex,length(ex),length(mirnas_mat)))
      rownames(miRNA_mat1)=mirnas_mat
      colnames(miRNA_mat1)=colnames(miRNA)

      miRNA_mat=rbind(miRNA_mat,miRNA_mat1)

    }
    miRNA_mat=miRNA_mat[-c(1,2),]
    miRNA_mat=miRNA_mat[!is.na(rownames(miRNA_mat)),]
    miRNA_mat=miRNA_mat[!(rownames(miRNA_mat) %in% rownames(miRNA)),]
    miRNA=miRNA[!(rownames(miRNA) %in% mirna_map$OriginalName),]
    miRNA=rbind(miRNA,miRNA_mat)
    miRNA=stats::aggregate(miRNA, by=list(row.names(miRNA)), FUN=median)
    rownames(miRNA)=miRNA[,1];miRNA=miRNA[,-1]
    miRNA=as.matrix(miRNA)

    miXp=sum(rownames(miRNA) %in% unique(names(databases$mirna_net)))
    if(miXp>0){
      message(">>number of miRNAs having target sets of genes: ", miXp)
    }else{
      message(">>miRNA name example: ",unique(names(databases$mirna_net))[1])
      stop("no miRNA as rowname of the mRNA count matrix has a literature target set of genes");
    }

    #Library and Rank normalization
    miRNA=prepare_RNAseq(counts=miRNA, groups=info$Groups, winz=winz, min=min, max=max, type=type)
    omics[["miRNA"]]=miRNA
    rm(miRNA,miRNA_mat,miXp)
  }


  #Transform somatic mutation omic by propagation -----
  if(!is.null(omics[["mutations"]])){
    message(">preparing somatic mutations")
    muts=omics[["mutations"]]

    #Get the gene interaction network
    g=igraph::graph_from_adjacency_matrix(
      adjmatrix=databases$gene_net$net_adj,
      mode=c("undirected"),
      weighted=TRUE,
      diag=TRUE
    )

    #Sets to 1 the gene directly connected to the sample's mutations
    for(icol in 1:ncol(muts)){
      mut_genes=rownames(muts)[muts[,icol]==1]
      mut_genes=intersect(colnames(databases$gene_net$net_adj),mut_genes)
      neighs=igraph::ego(graph=g, order=1, nodes=mut_genes, mode="all")
      neighs=unlist(lapply(neighs,function(x){igraph::as_ids(x)}))
      freq_muts=as.data.frame(table(neighs))
      neighs=freq_muts$neighs[freq_muts$Freq>=quantile(freq_muts$Freq,0.75)]
      muts[rownames(muts) %in% neighs,icol]=1
    }

    #Check if you have proper gene and lncRNA names
    gXp=sum(rownames(muts) %in% unique(unlist(databases$pathway_db$pathways_l)))
    if(gXp>0){
      message(">>number of mutated genes falling into pathways: ", gXp)
    }else{
      message(">>mutated gene name example: ",databases$pathway_db$pathways_l[[1]][1])
      stop("no mutated  gene as rowname of the mRNA count matrix belongs to any literature pathway");
    }

    #Propagate mutation information
    prop_muts=get_propagated(net=databases$gene_net$net_adj, counts=muts, n_cores=n_cores, r=0.8, keep_no_nodes=T)
    #Rank transform
    prop_muts=get_rank_01(m=prop_muts)
    omics[["mutations"]]=prop_muts
    rm(muts,g,mut_genes,neighs,freq_muts,gXp,prop_muts)
  }

  #Check cross validation setting
  if(sum(cv_probs)!=1){
    stop("the vector cv_probs has user-defined probabilites which sum is not 1")
  }else{
    message("Preparing partitions of data for cross validation setting")
  }

  #Partitioning of samples for cross validation -----
  partitions_l=list()
  cv_probs[1]=cv_probs[1]+cv_probs[2]
  for(iter_i in 1:n_iter_cv){
    #Get training and testing
    indx <- splitTools::partition(info$Groups, type="stratified", n_bins=2,
                                  p=c(train = cv_probs[1], test=cv_probs[3]),
                                  seed=seed)

    #Set training and testing labels
    set_type=c(rep("TRAINING",length(indx$train)),rep("TESTING",length(indx$test)))
    tmp_info=info[c(indx$train,indx$test),]
    tmp_info$set_type=set_type
    tmp_omics=lapply(omics,function(x){x[,c(indx$train,indx$test)]})

    #Get similarity between training and testing
    if("gex" %in% names(tmp_omics)){
      PSN=build_adjWJ_PSN(tmp_omics$gex)
    }else{
      if("mutations" %in% names(tmp_omics)){
        PSN=build_adjWJ_PSN(tmp_omics$mutations)
      }
    }

    #Find the samples in each group that are similar to the testing samples to compose the validation set
    #This operation is always possible also with real blind testing samples
    PSN1=PSN[tmp_info$set_type=="TRAINING" & tmp_info$Groups==unique(tmp_info$Groups)[1],tmp_info$set_type=="TESTING"]
    n_valid=round(nrow(PSN1)*cv_probs[2])
    if(n_valid<1){n_valid=1}
    rm1=abs(scale(rowMeans(PSN1))[,1])
    rm1=rm1[order(rm1,decreasing = F)]
    val1=names(rm1)[1:n_valid]

    PSN2=PSN[tmp_info$set_type=="TRAINING" & tmp_info$Groups==unique(tmp_info$Groups)[2],tmp_info$set_type=="TESTING"]
    n_valid=round(nrow(PSN2)*cv_probs[2])
    if(n_valid<1){n_valid=1}
    rm2=abs(scale(rowMeans(PSN2))[,1])
    rm2=rm2[order(rm2,decreasing = F)]
    val2=names(rm2)[1:n_valid]

    #Set validation set
    indx$valid=which(info$IDs %in% c(val1,val2))
    indx$train=setdiff(indx$train,indx$valid)
    #colnames(info_df)=c("IDs","Groups","val_indxs","test_indxs","train_indxs")

    #Set binary labels
    info$train_indxs=info$test_indxs=info$val_indxs=0
    info$train_indxs[indx$train]=1
    info$test_indxs[indx$test]=1
    info$val_indxs[indx$valid]=1

    #Generate two lists:
    #one containing only training data for the feature selection and the omics analysis
    #one containing all the samples to build and test the model
    tr_info=info[info$train_indxs==1,]
    tr_omics=lapply(omics,function(x){x[,indx$train]})
    #te_omics=lapply(omics,function(x){x[,c(indx$train,indx$valid,indx$test)]})
    te_omics=lapply(omics,function(x){x})


    TRAINING=list(info=tr_info,
                  omics=tr_omics,
                  databases=databases,
                  tax_id=tax_id,
                  n_cores=n_cores)

    TESTING=list(info=info,
                 omics=te_omics,
                 databases=databases,
                 tax_id=tax_id,
                 n_cores=n_cores)

    partitions_l[[iter_i]]=list(TRAINING=TRAINING,TESTING=TESTING)
  }

  names(partitions_l)=seq(1,length(partitions_l))
  return(partitions_l)
}

#' Analyse the training set to get the best PSNs out of the omics
#'
#' For each omic, it determines the best molecules to build the PSNs.
#' The molecules are fit into the pathways or regulatory target sets.
#' A PSN is built for each pathway, is tested if separate the groups in comparison.
#' If yes, the PSN is kept and used for the training and classification.
#'
#' @param data_l The output list obtained from the processing of the user-defined data with the function: prepare_data
#' @param max_size Default 150, integer value filtering out the pathways/sets with more elements than this threshold
#' @param min_elements Default 4, integer value filtering out the pathways/sets containing less than this threshold of significant molecules
#' @param n_top_sets Default 50, integer that set how many best sets/PSNs to retrieve
#' @param n_cores Default 2, integer that set the number of cores to use for running the function in parallel
#' @param keep_PSN Default FALSE, do not keep the PSN generated from a set of element
#' @return For each unique omic (e.g. gene expression), it provides a dataframe (SDR) which is a result of the function
#' find_SDR, for example $TRAINING$gex$SDR. For each omic type (e.g. transcriptomics which can be combination of genes, lncRNAs and miRNAs), it adds
#' a list of two elements in the $TRAINING$pathway_analysis$transcriptomics section of the data list: top_pathways_l and top_pathways_df
#' top_pathways_df details the pathways/sets which have a significant PSN. The list is the result of the function find_PSN
#' @export
analyse_training = function(data_l, max_size=200, min_elements=6, n_top_pathways=50, n_cores=2, keep_PSN=FALSE){
  #Iterate over training sets for pathway and sets selection based on their patient similarity network
  for(iter_cv in 1:length(data_l)){
    message("Feature selection of training data belonging to: ",iter_cv," run of cross validation")
    #Extract the data of training set of a run of cross validation
    i_tr_data_l=data_l[[iter_cv]]$TRAINING
    groups_freq=as.data.frame(table(i_tr_data_l$info$Groups))
    groups_name=unique(i_tr_data_l$info$Groups)
    rownames(groups_freq)=groups_freq[,1]
    message(">training samples:")
    message(">>first group: ",groups_freq[groups_name[1],1]," of training samples: ",groups_freq[groups_name[1],2])
    message(">>second group: ",groups_freq[groups_name[2],1]," of training samples: ",groups_freq[groups_name[2],2])

    #Instantiate the resulting lists about the pathway and target sets analysis
    i_tr_data_l[["pathway_analysis"]]=list()
    gex_pa_l=mut_pa_l=miRNA2g_pa_l=lncRNA2g_pa_l=lncRNA2mi_pa_l=miRNA2lnc_pa_l=list()

    #Gene expression analysis of training data
    if("gex" %in% names(i_tr_data_l$omics)){
      message(">finding significant group specific genes based on gene expression and update $TRAINING$omics$gex")
      #Detect specific differential role genes and update the gex list of data
      i_tr_data_l$omics$gex = list(norm_counts=i_tr_data_l$omics$gex,
                                   SDR=find_SDR_directional(i_tr_data_l$omics$gex, i_tr_data_l$info$Groups))

      message(">>finding enriched pathways and keep the ones with best PSNs")
      #Detect which pathway fit most of the signifant genes, create and test the PSN
      #Keep the pathways which have the best patient similarity networks
      gex_pa_l = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                          groups = i_tr_data_l$info$Groups,
                          pathways_l = i_tr_data_l$databases$pathway_db$pathways_l,
                          sign_df = i_tr_data_l$omics$gex$SDR,
                          n_cores = i_tr_data_l$n_cores,
                          max_size = max_size,
                          min_elements = min_elements,
                          n_top_pathways = n_top_pathways,
                          keep_PSN = keep_PSN)
      if(length(gex_pa_l$top_pathways_l)!=0){
        gex_pa_l$top_pathways_df$type="gene_canonical_pathways"
        gex_pa_l$top_pathways_df$source="gene_expression"
      }
    }

    #Somatic mutation analysis of training data
    if("mutations" %in% names(i_tr_data_l$omics)){
      message(">finding significant group specific mutated genes based on propagated mutations and update $TRAINING$omics$mutations")
      i_tr_data_l$omics$mutations = list(norm_counts=i_tr_data_l$omics$mutations,
                                         SDR=find_SDR(i_tr_data_l$omics$mutations, i_tr_data_l$info$Groups))

      message(">>finding enriched pathways and keep the ones with best PSNs")
      mut_pa_l = find_PSN(counts = i_tr_data_l$omics$mutations$norm_counts,
                          groups = i_tr_data_l$info$Groups,
                          pathways_l = i_tr_data_l$databases$pathway_db$pathways_l,
                          sign_df = i_tr_data_l$omics$mutations$SDR,
                          n_cores = i_tr_data_l$n_cores,
                          max_size = max_size,
                          min_elements = min_elements,
                          n_top_pathways = n_top_pathways,
                          keep_PSN=keep_PSN)

      if(length(mut_pa_l$top_pathways_l)!=0){
        mut_pa_l$top_pathways_df$type="gene_canonical_pathways"
        mut_pa_l$top_pathways_df$source="mutation_propagation"
      }
    }

    #miRNA analysis
    if("miRNA" %in% names(i_tr_data_l$omics)){
      message(">finding significant group specific miRNAs based on miRNA expression and update $TRAINING$omics$miRNA")
      i_tr_data_l$omics$miRNA = list(norm_counts=i_tr_data_l$omics$miRNA,
                                     SDR=find_SDR_directional(i_tr_data_l$omics$miRNA,
                                                              i_tr_data_l$info$Groups,
                                                              top_molecules_th = 0.15))
      #Filter the list of miRNA's targets to keep only the sets of the significant miRNAs
      SDR_miRNA_genes_l=list()
      #For each significant miRNA keep only the anti-correlated target significant genes and lncRNAs
      miRNA_tT=i_tr_data_l$omics$miRNA$SDR;miRNA_tT=miRNA_tT[miRNA_tT$Significant==TRUE,]
      #Genes
      gex_tT=i_tr_data_l$omics$gex$SDR;gex_tT=gex_tT[gex_tT$Significant==TRUE,]

      #Extract the directions of the miRNAs
      barcodes=unique(miRNA_tT$Barcode)
      for(barcode in barcodes){
        #Extract the anti-correlated
        miRNA_tT_dir=miRNA_tT[miRNA_tT$Barcode==barcode,]
        gex_tT_dir=gex_tT[gex_tT$Barcode==(-barcode),]

        #Filter the target list
        miRNA_genes_l=i_tr_data_l$databases$mirna_net
        SDR_miRNA_genes_l1=miRNA_genes_l[which(names(miRNA_genes_l) %in% rownames(miRNA_tT_dir))]
        SDR_miRNA_genes_l1=lapply(SDR_miRNA_genes_l1,function(x){
          y=x[x %in% rownames(gex_tT_dir)]
          return(y)
        })

        SDR_miRNA_genes_l=c(SDR_miRNA_genes_l,SDR_miRNA_genes_l1)
      }
      #lncRNAs
      miRNA_lncRNAs_l=i_tr_data_l$databases$lncrna_net$mirna_lncrna_targ_l
      SDR_miRNA_lncRNAs_l=miRNA_lncRNAs_l[which(names(miRNA_lncRNAs_l) %in% rownames(miRNA_tT_dir))]

      #Clean
      rm(miRNA_tT,gex_tT,miRNA_lncRNAs_l,miRNA_genes_l)

      message(">>finding enriched miRNA's gene targets and keep the ones with best PSNs")
      miRNA2g_pa_l = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                              groups = i_tr_data_l$info$Groups,
                              pathways_l = SDR_miRNA_genes_l,
                              sign_df = i_tr_data_l$omics$gex$SDR,
                              n_cores = i_tr_data_l$n_cores,
                              max_size = 9999,
                              min_elements=min_elements,
                              n_top_pathways = n_top_pathways,
                              keep_PSN=keep_PSN)
      if(length(miRNA2g_pa_l$top_pathways_l)!=0){
        miRNA2g_pa_l$top_pathways_df$type="miRNA_targets"
        miRNA2g_pa_l$top_pathways_df$source="gene_expression"
      }
    }

    #lncRNA analysis
    if("lncRNA" %in% names(i_tr_data_l$omics)){
      message(">finding significant group specific lncRNAs based on lncRNA expression and update $TRAINING$omics$lncRNA")
      i_tr_data_l$omics$lncRNA = list(norm_counts=i_tr_data_l$omics$lncRNA,
                                      SDR=find_SDR(i_tr_data_l$omics$lncRNA,
                                                   i_tr_data_l$info$Groups,
                                                   top_molecules_th = 0.15))
      message(">finding top enriched gene targets from significant lncRNAs")
      #Extract the significant lncRNAs
      lncRNA_tT=i_tr_data_l$omics$lncRNA$SDR;lncRNA_tT=lncRNA_tT[lncRNA_tT$Significant==TRUE,]
      SDR_lncRNA=rownames(lncRNA_tT)
      #Keep the sets of gene targets and miRNAs belonging to only the significant lncRNAs
      lncRNA_genes_l=i_tr_data_l$databases$lncrna_net$lncrna_gene_targ_l
      lncRNA_miRNA_l=i_tr_data_l$databases$lncrna_net$lncrna_mirna_targ_l
      SDR_lncRNA_genes_l=lncRNA_genes_l[which(names(lncRNA_genes_l) %in% SDR_lncRNA)]
      SDR_lncRNA_miRNAs_l=lncRNA_miRNA_l[which(names(lncRNA_miRNA_l) %in% SDR_lncRNA)]

      message(">>finding enriched lncRNA's gene targets and keep the ones with best PSNs")
      lncRNA2g_pa_l = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                               groups = i_tr_data_l$info$Groups,
                               pathways_l = SDR_lncRNA_genes_l,
                               sign_df = i_tr_data_l$omics$gex$SDR,
                               n_cores = i_tr_data_l$n_cores,
                               max_size=9999,
                               min_elements=min_elements,
                               n_top_pathways = n_top_pathways,
                               keep_PSN = keep_PSN)
      if(length(lncRNA2g_pa_l$top_pathways_l)!=0){
        lncRNA2g_pa_l$top_pathways_df$type="lncRNA_targets"
        lncRNA2g_pa_l$top_pathways_df$source="gene_expression"
      }

      message(">>finding enriched lncRNA's miRNA targets and keep the ones with best PSNs")
      lncRNA2mi_pa_l = find_PSN(counts = i_tr_data_l$omics$miRNA$norm_counts,
                                groups = i_tr_data_l$info$Groups,
                                pathways_l = SDR_lncRNA_miRNAs_l,
                                sign_df = i_tr_data_l$omics$miRNA$SDR,
                                n_cores = i_tr_data_l$n_cores,
                                max_size=9999,
                                min_elements=min_elements,
                                n_top_pathways = n_top_pathways,
                                keep_PSN = keep_PSN)
      if(length(lncRNA2mi_pa_l$top_pathways_l)!=0){
        lncRNA2mi_pa_l$top_pathways_df$type="lncRNA_targets"
        lncRNA2mi_pa_l$top_pathways_df$source="miRNA_expression"
      }

      if("miRNA" %in% names(i_tr_data_l$omics)){
        message(">>finding enriched miRNA's lncRNA targets and keep the ones with best PSNs")
        miRNA2lnc_pa_l = find_PSN(counts = i_tr_data_l$omics$lncRNA$norm_counts,
                                  groups = i_tr_data_l$info$Groups,
                                  pathways_l = SDR_miRNA_lncRNAs_l,
                                  sign_df = i_tr_data_l$omics$lncRNA$SDR,
                                  n_cores = i_tr_data_l$n_cores,
                                  max_size=9999,
                                  min_elements=min_elements,
                                  n_top_pathways = n_top_pathways,
                                  keep_PSN = keep_PSN)
        if(length(miRNA2lnc_pa_l$top_pathways_l)!=0){
          miRNA2lnc_pa_l$top_pathways_df$type="miRNA_targets"
          miRNA2lnc_pa_l$top_pathways_df$source="lncRNA_expression"
        }
      }
    }

    #Merge pathway results
    l2merge=list(gex_pa_l=gex_pa_l,
                 miRNA2g_pa_l=miRNA2g_pa_l,
                 lncRNA2g_pa_l=lncRNA2g_pa_l,
                 lncRNA2mi_pa_l=lncRNA2mi_pa_l,
                 miRNA2lnc_pa_l=miRNA2lnc_pa_l)
    transcriptomics_pa_l=merge_multi_omics_analysis(l2merge)

    #Combine all the pathway results, update training object and return
    message(">adding significant pathways/targets/PSNs into $TRAINING$pathway_analysis")
    pathway_res_l=list(transcriptomics=transcriptomics_pa_l,mutations=mut_pa_l)
    pathway_res_l=pathway_res_l[sapply(pathway_res_l,length)!=0]
    i_tr_data_l[["pathway_analysis"]]=pathway_res_l
    data_l[[iter_cv]]$TRAINING=i_tr_data_l
  }
  return(data_l)
}

#' Classify testing samples based on training and PSNs
#'
#' For each run of cross-validation, extract the pathways/target sets/PSNs that are significant for training data.
#' Generate the PSN, format it as edge list and pass the networks to the graph convolutional network for the training and prediction.
#'
#' @param data_l The output list obtained from the processing of the user-defined data with the function: prepare_data and analyse_training
#' @param n_top_PSNs Default 30, cutoff that defines how many PSNs to use for each omic (e.g. transcriptomics) to classify the testing samples
#' @param py_path Default package python script
#' @return Updates the info dataframe in $TESTING$info based on the predictions made for the testing and validation samples.
#' It adds a new list to $TESTING called classification which contains the dataframe perfsXpathway_df and the matthew correlation coefficent mcc_score.
#' perfsXpathway_df is a dataframe that shows the performances of classification obtained with each pathway/target set/PSN.
#' mcc_score is of double value between -1 and 1 representing the performance of classification.
#' @importFrom reticulate source_python
#' @importFrom scales rescale
#' @export
classify_testing = function(data_l,n_top_PSNs=30,py_path=NULL){
  if(is.null(py_path)){
    py_path=system.file("python", "sage_classification_individual.py", package = "StellarPath")
  }

  #Iterate over the run of cross validation
  for(iter_cv in 1:length(data_l)){
    i_tr_data_l=data_l[[iter_cv]]$TRAINING
    i_te_data_l=data_l[[iter_cv]]$TESTING
    info_df=i_te_data_l$info

    #Initiate the lists that will contain the PSNs and their names
    all_edges <- list()
    all_names <- list()

    #Iterate over the pathways selected with the multi-omics training patients
    for(iter_pth_l in 1:length(i_tr_data_l$pathway_analysis)){
      sets_type=names(i_tr_data_l$pathway_analysis)[iter_pth_l]

      if(sets_type=="transcriptomics"){
        transcripts_types=names(i_tr_data_l$omics)

        m=i_te_data_l$omics$gex
        if("lncRNA" %in% transcripts_types){
          m=rbind(m,i_te_data_l$omics$lncRNA)
        }
        if("miRNA" %in% transcripts_types){
          m=rbind(m,i_te_data_l$omics$miRNA)
        }
      }

      if(sets_type=="mutations"){
        m=rbind(i_te_data_l$omics$mutations)
      }

      pths_l=i_tr_data_l$pathway_analysis[[iter_pth_l]]
      if(nrow(pths_l$top_pathways_df)<5){next;}

      if(n_top_PSNs>nrow(pths_l$top_pathways_df)){
        n_top_PSNs=nrow(pths_l$top_pathways_df)
      }

      use2classify=1:n_top_PSNs
      pths_l$top_pathways_l=pths_l$top_pathways_l[use2classify]
      pths_l$top_pathways_df=pths_l$top_pathways_df[use2classify,]

      #Iterate over the pathways selected from the training set of a run
      for(iter_pth in 1:length(pths_l$top_pathways_l)){
        #Recover stats from the training data
        name_psn=paste(sets_type,": ",names(pths_l$top_pathways_l)[iter_pth],sep="")
        gs=pths_l$top_pathways_l[[iter_pth]]$set
        dir=pths_l$top_pathways_l[[iter_pth]]$stats$regulation

        #Extract the single molecule's values and build the PSN with training, validation and testing samples
        m_gs=m[rownames(m) %in% gs,]
        if(dir=="activated"){
          PSN=build_PSNs(m_gs)$PSN_high
        }else{
          PSN=build_PSNs(m_gs)$PSN_low
        }
        #Convert to edge list
        el_PSN=convert_adj2edg(PSN)
        #Filter low similarities and 01 standardize their values
        el_PSN=el_PSN[el_PSN$weight>=quantile(el_PSN$weight,0.6),]
        el_PSN$weight=scales::rescale(el_PSN$weight,to=c(0,1))

        #Save
        edges <- list(el_PSN)
        all_edges <- c(all_edges,edges)
        all_names <- c(all_names,name_psn)
      }

    }

    #Load and Compile function in Python script
    reticulate::source_python(py_path, convert = TRUE)
    #Function for node classification with Sage that trains on all over the database
    res <- node_classification(all_names,info_df,all_edges)
    #Extract and Format results from the Sage classification
    perfsXpathway_df=res[[1]]
    validation_predictions=as.data.frame(apply(res[[2]],2,unlist))
    consensus_validation_pred_labels=apply(validation_predictions,2,function(x){
      freqs_df=as.data.frame(table(x))
      pred_label=as.character(freqs_df[which.max(freqs_df$Freq),1])
      return(pred_label)
    })

    test_predictions=as.data.frame(apply(res[[3]],2,unlist))
    consensus_test_pred_labels=apply(test_predictions,2,function(x){
      freqs_df=as.data.frame(table(x))
      pred_label=as.character(freqs_df[which.max(freqs_df$Freq),1])
      return(pred_label)
    })

    #Update info dataframe and compute MCC classification performance
    info_df$val_preds=info_df$test_preds="not_used"
    info_df$val_preds[info_df$val_indxs==1]=consensus_validation_pred_labels
    info_df$test_preds[info_df$test_indxs==1]=consensus_test_pred_labels
    mcc_score=get_mcc(info_df$Groups[info_df$test_indxs==1],consensus_test_pred_labels, unique(info_df$Groups)[1])

    #Save
    message("Classification Peformance of the run: ",iter_cv ," is: ", mcc_score)
    message("Updating $TESTING$classification in TESTING data list")
    i_te_data_l[["info"]]=info_df
    i_te_data_l[["classification"]]=list(perfsXpathway_df=perfsXpathway_df,mcc_score=mcc_score)
    data_l[[iter_cv]]$TESTING=i_te_data_l
  }

  #Retrieve MCC scores from the CV runs
  mcc_scores=sapply(data_l[1:iter_cv],function(x){
    x$TESTING$classification$mcc_score
  })
  perfs=data.frame(name_dataset=as.character(Sys.time()),mcc_scores=mcc_scores,runs=seq(1,length(mcc_scores)))

  #Retrieve performances of the model on each training feature selected pathway
  perfsXpathway_dfs=lapply(data_l[1:iter_cv],function(x){
    x$TESTING$classification$perfsXpathway_df
  })

  #Save performances of the multiple runs of cross validations
  data_l[["performances"]]=list(performances=perfs,
                                perfsXpathway_dfs=perfsXpathway_dfs)

  return(data_l)
}

#' Enrichment analysis of the results of classification
#'
#' It takes in input the resulting list of the classification function. It finds the run of cross validation in which the model classified
#' the samples with the best performance. It finds the similarity networks that separate the two sample's groups.
#' It adds annotations to the significant gene pathways/non-coding target sets/PSNs:
#' - ID and database source for canonical gene pathways
#' - Regulation pathway targeted by the non-coding significant elements
#'
#' It creates a similarity network between the pathways/sets of targets that have been found significant with their PSNs.
#' For each pair of pathway/target set, there is a weighted edge which measures how much is strong their overlap and so their similarity.
#' This network is save in an object called EnrichmentMap because it is analog to the idea presented with Cytoscape.
#'
#' @param out_l The output list obtained from the classification of the user-defined data with the function: prepare_data, analyse_training and classification
#' @param extended Default FALSE, the overlap between pathways is computed with only the significant features falling into that pathways
#' @return For each unique omic (e.g. gene expression), it provides a dataframe (SDR) which is a result of the function
#' find_SDR, for example $TRAINING$gex$SDR. For each omic type (e.g. transcriptomics which can be combination of genes, lncRNAs and miRNAs), it adds
#' a list of two elements in the $TRAINING$pathway_analysis$transcriptomics section of the data list: top_pathways_l and top_pathways_df
#' top_pathways_df details the pathways/sets which have a significant PSN. The list is the result of the function find_PSN.
#' - EnrichmentMap similarity network between the pathways/sets of targets that have been found significant with their PSNs
#' @importFrom plyr mapvalues
#' @export
enrichment_analysis = function(out_l,extended=FALSE){
  #Select the best run of cross validation
  best_run=which.max(out_l$performances$performances$mcc_scores)

  #Get the samples info and replace the original groups with the model's predictions
  message("Enrichment with predictions of the run: ", best_run, " of cross validation with best performances")
  info=out_l[[best_run]]$TESTING$info
  info$Groups[info$test_preds!="not_used"]=info$test_preds[info$test_preds!="not_used"]
  #Create a list object to fit in the analyse_training function
  enr_l=list(info=info,
             omics=out_l[[best_run]]$TESTING$omics,
             databases=out_l[[best_run]]$TESTING$databases,
             n_cores=out_l[[best_run]]$TESTING$n_cores)

  data_l=list()
  data_l[["1"]]=list(TRAINING=enr_l)
  #It will run with all the samples and not only the training set
  data_l=analyse_training(data_l, keep_PSN = TRUE)
  #Extract the results and return
  enr_l=data_l$`1`$TRAINING
  #Get literature pathways
  db_info=enr_l$databases$pathway_db$metainfo

  #For each type of omics (e.g. transcriptomics)
  for(si in 1:length(names(enr_l$pathway_analysis))){
    #Recover the stats about the pathways differentiating the two classes in study
    stats=enr_l$pathway_analysis[[si]]$top_pathways_df
    paths_l=enr_l$pathway_analysis[[si]]$top_pathways_l

    if(nrow(stats)!=0){
      #Find the code ID and the database source of the significant pathways of genes based on their PSNs
      code=plyr::mapvalues(stats$pathways, from=db_info$pathway_name, to=db_info$term, warn_missing = F)
      source=plyr::mapvalues(stats$pathways, from=db_info$pathway_name, to=db_info$source, warn_missing = F)
      code[stats$type!="gene_canonical_pathways"]=NA
      source[stats$type!="gene_canonical_pathways"]=NA
      stats$code_db=code
      stats$name_db=source
    }else{
      next;
    }


    #Get significant literature gene pathways of regulation
    bioc_paths_l=enr_l$databases$pathway_db$pathways_l[
      grep("REGULATION",names(enr_l$databases$pathway_db$pathways_l))]

    #Get the significant PSNs created with the gene targets of non coding elements
    indxs=c(which(stats$type=="miRNA_targets" & stats$source=="gene_expression"),
            which(stats$type=="lncRNA_targets" & stats$source=="gene_expression"))

    #Add the pathways that miRNA and lncRNA targets target
    if(length(indxs)!=0){
      for(mi in 1:length(indxs)){
        indx=indxs[mi]
        feature=stats$pathways[indx]
        gs_targs=paths_l[[feature]]$set

        overlap=sapply(bioc_paths_l,function(x){
          tot=sum(gs_targs %in% x)
          if(tot!=0){
            score1=tot/length(gs_targs)
            score2=tot/length(x)
            score=min(score1,score2)
          }else{
            score=0
          }
          return(score)
        })
        overlap=2000-overlap;rnk=rank(overlap);
        paths_overlapping=names(rnk)[rnk<=5]
        if(length(paths_overlapping)!=0){
          codes=plyr::mapvalues(paths_overlapping, from=db_info$pathway_name, to=db_info$term, warn_missing = F)
          paths_overlapping=paste(paths_overlapping,collapse = ",")
          codes=paste(codes,collapse = ",")
          df_extra=data.frame(pathways=feature,pathways_targets=paths_overlapping,code_db_targets=codes)
        }else{
          df_extra=data.frame(pathways=feature,pathways_targets=NA,code_db_targets=NA)
        }

        if(mi==1){
          df_extras=df_extra
        }else{
          df_extras=rbind(df_extras,df_extra)
        }
      }
      stats1=merge(stats,df_extras,by = "pathways",all.x = T)
      stats1=stats1[match(stats$pathways,stats1$pathways),]

      #Update and Return
      stats=enr_l$pathway_analysis[[si]]$top_pathways_df=stats1
    }
  }

  #For each type of omics (e.g. transcriptomics)
  dss=names(enr_l$pathway_analysis)
  #Get the name of the pathways from the different pathway analysis and associate them an UI
  paths="start"
  for(k in 1:length(dss)){
    ds=dss[k]
    stats=enr_l$pathway_analysis[[ds]]$top_pathways_df
    paths_l=enr_l$pathway_analysis[[ds]]$top_pathways_l
    if(nrow(stats)==0){next;}
    paths1=paste(k,stats$pathways,sep="_")
    stats$pathways=paths1
    paths=c(paths,paths1)
    if(k == 1){
      ann_extra=stats[,c("pathways","type","source")]
    }else{
      ann_extra=rbind(ann_extra,stats[,c("pathways","type","source")])
    }
  }
  paths=paths[-1]

  #Create the network
  combs=as.data.frame(t(combn(x=as.character(paths),m=2)))
  colnames(combs)=c("node1","node2")

  #Decompose to get the original source of the pathways and the propert name of the pathways
  analysis_source1=as.numeric(substr(combs$node1,1,1))
  analysis_source2=as.numeric(substr(combs$node2,1,1))
  combs$analysis_source1=analysis_source1
  combs$analysis_source2=analysis_source2
  combs$analysis_source1=dss[combs$analysis_source1]
  combs$analysis_source2=dss[combs$analysis_source2]

  #Add extra annotation about pathways
  type1=plyr::mapvalues(combs$node1,from=ann_extra$pathways, to=ann_extra$type, warn_missing = FALSE)
  source1=plyr::mapvalues(combs$node1,from=ann_extra$pathways, to=ann_extra$source, warn_missing = FALSE)
  type2=plyr::mapvalues(combs$node2,from=ann_extra$pathways, to=ann_extra$type, warn_missing = FALSE)
  source2=plyr::mapvalues(combs$node2,from=ann_extra$pathways, to=ann_extra$source, warn_missing = FALSE)
  combs$type1=type1
  combs$source1=source1
  combs$type2=type2
  combs$source2=source2

  #Return to the original names of the pathways
  combs$node1=substr(combs$node1,3,9000)
  combs$node2=substr(combs$node2,3,9000)

  #Compute jaccard similarity of overlap between pathways
  score_v=c(0)
  for(ri in 1:nrow(combs)){
    f1=combs$node1[ri]
    f2=combs$node2[ri]
    type1=combs$type1[ri]
    type2=combs$type2[ri]
    analysis_source1=combs$analysis_source1[ri]
    analysis_source2=combs$analysis_source2[ri]

    if(length(grep("targets",type1))!=0){
      set1=enr_l$pathway_analysis[[analysis_source1]]$top_pathways_l[[f1]]$set
    }else{
      set1=enr_l$pathway_analysis[[analysis_source1]]$top_pathways_l[[f1]]$orig_set
    }

    if(length(grep("targets",type2))!=0){
      set2=enr_l$pathway_analysis[[analysis_source2]]$top_pathways_l[[f2]]$set
    }else{
      set2=enr_l$pathway_analysis[[analysis_source2]]$top_pathways_l[[f2]]$orig_set
    }

    #Measure the overlap
    score=(length(intersect(set1,set2)))/(min(length(set1),length(set2)))
    score_v=c(score_v,score)
  }
  score_v=score_v[-1]
  combs$weight=score_v

  enr_l[["EnrichmentMap"]]=combs
  return(enr_l)
}
