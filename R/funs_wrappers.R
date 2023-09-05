#' Check and Prepare the input data
#'
#' This function checks that the input data satisfy all the criteria for being used by the method.
#' It creates an object list that includes all the checked, normalized and prepared data.
#'
#' @param info Dataframe, two character columns, first has sample's IDs (have to match with column names of the count matrices),
#' second has sample's groups (only two as for a pairwise DE analysis, for example c(AD,AD,AD,HEALTHY,HEALTHY,HEALTHY))
#' @param groups_name Character vector of length two, it must contain the same labels of unique(info[,2]), tip: first the case group like AD,
#' while the second name refers to the control group like HEALTHY
#' @param gex Numeric matrix of raw count values from bulk RNA sequencing (e.g. genes x samples) (it may include also lncRNAs)
#' @param lncRNA Numeric matrix of raw count values from bulk RNA sequencing (e.g. lncRNA x samples)
#' @param miRNA Numeric matrix of raw count values from miRNA sequencing (e.g. miRNA x samples)
#' @param mutations Numeric matrix of binary values from somatic mutation data (e.g. gene X samples, 1 per mutated gene, 0 otherwise)
#' @param tax_id Default 9606, integer value indicating taxonomic id of the samples, either 9606 (human) or 10090 (mouse)
#' @param simplify_names Default TRUE, simplify names of the samples in the info and count matrices
#' @param cv_probs Default c(0.7,0.1,0.2), vector of doubles indicating the proportion of the samples to use as training, validation and testing set
#' @param user_sets Default NULL, list of elements (one element is used to compose the training and testing set of a cv run) s.t. each element is a list containing two numeric vectors, one called train and one called test, each vector indicate the column index of those patients falling into training or testing set
#' @param n_iter_cv Default 2, integer value indicating the number of iterations to perform the cross validation
#' @param n_cores Default 2, integer value greater than 0 indicating the number of cores to use to parallelize and speed up the operations
#' @param seed Default 5, integer value indicating the seed to replicate the sampling of the cross validation
#' @param winz Default FALSE, boolean to indicate if to adjuste outliers with winsorizing
#' @param type Default total rank, character string which defines the type of transformation applied to the count data, options in the help of the get_rank_01 function
#' @return A list s.t. each element contains the user-input prepared data for a run of cross-validation.
#' The run's data are divided into training and testing sets.
#' Each set contains a specific subset of sample's information (i.e. info) and normalized count matrices (i.e. omics).
#' @importFrom utils combn
#' @importFrom splitTools partition
#' @importFrom miRBaseConverter miRNA_PrecursorToMature
#' @importFrom stats aggregate
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph ego
#' @importFrom igraph as_ids
#' @export
#'
prepare_data = function(info, groups_name=NULL,
                        gex=NULL, lncRNA=NULL, miRNA=NULL, mutations=NULL,
                        tax_id=9606, simplify_names=TRUE,
                        cv_probs=c(0.7,0.1,0.2), n_iter_cv=2, user_sets=NULL,
                        n_cores=2, seed=5,
                        winz=FALSE, type="sample_rank"){

  #Initiate data variables
  omics=list()

  #Checking user provided omics
  #In case an omic exists, put the omic into the list of available omics
  #gene expression matrix is always needed because the method needs to reach the pathway level
  message("Checking input omics/raw count matrices")
  if(is.null(gex)){
    stop("gene expression matrix (genes per samples) is needed to reach the pathways");
  }

  if(!is.null(gex)){
    message(">gene expression data (genes per samples) found");
    omics[["gex"]]=gex
    rm(gex)
  }

  if(!is.null(lncRNA)){
    message(">lncRNA data (lncRNA per samples) found");
    omics[["lncRNA"]]=lncRNA
    rm(lncRNA)
  }

  if(!is.null(miRNA)){
    message(">miRNA data (miRNA per samples) found");
    omics[["miRNA"]]=miRNA
    rm(miRNA)
  }

  if(!is.null(mutations)){
    omics[["mutations"]]=mutations
    message(">mutation data (genes per samples) found");
    rm(mutations)
  }

  #Checking user-provided dataframe info about samples:
  #a dataframe with two columns, ID and Groups' labels in character format
  message("Checking input info")
  if(is.data.frame(info)){
    message(">dataframe format")
    if(ncol(info)==2){
      message(">two columns")
      if(sum(is.na(info))==0){
        message(">does NOT contain NA values")
        if(class(info[,1])=="character" & class(info[,2])=="character"){
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
          stop("info dataframe contains columns in non-character clas");
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

  #Checking the samples' IDs match between info dataframe and column names in the omics
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

  #Checking that omics follow the same matrix format
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
    stop("omics are NOT ALL matrices");
  }
  rm(class_omics)

  #Ordering the samples based on the user's defined groups of samples in comparison
  message("Ordering samples based on groups")
  if(is.null(groups_name)){
    ord=order(info$Groups,decreasing = TRUE)
    info=info[ord,]
    rownames(info)=info$IDs
    omics=lapply(omics,function(x){
      x=x[,ord]
      return(x)
    })
    groups_name=unique(info$Groups)
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
  }

  groups_freq=as.data.frame(table(info$Groups))
  message(">The comparison and classification is with:")
  rownames(groups_freq)=groups_freq[,1]
  colnames(groups_freq)=c("group","freq")
  message(">>first group : ",groups_freq[groups_name[1],1]," of samples: ",groups_freq[groups_name[1],2])
  message(">>second group : ",groups_freq[groups_name[2],1]," of samples: ",groups_freq[groups_name[2],2])
  rm(ord,groups_freq)

  #Convert samples' names which can be articulated and complicated into names that are more easily to recognize
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

  #Load and set networks and lists from databases: human and mouse
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

  #Cleaning, Normalizing and Preparing the omics to be analyzed (DE, Stability, Pathway analysis)
  message("Processing and Normalizing omics")
  #Preparing gene expression matrix
  if(!is.null(omics[["gex"]])){

    #aggregate duplicated rows (usually are non-annotated genomics regions which match in a preliminary gene symbol)
    message(">preparing gex")
    gex_unf=omics[["gex"]]
    gex_unf=gex_unf[setdiff(seq(1,nrow(gex_unf)),grep("\\.[1-999]",rownames(gex_unf))),]
    gex_unf=gex_unf[setdiff(seq(1,nrow(gex_unf)),grep("Gm[1-999]",rownames(gex_unf))),]
    if(sum(duplicated(rownames(gex_unf)))!=0){
      gex_unf=remove_duplicates(gex_unf)
    }

    #find and separate lncRNA appearing in the gene expression matrix (can easily happen)
    lncRNA_ids=unique(c(names(databases$lncrna_net$lncrna_gene_targ_l),
                        names(databases$lncrna_net$lncrna_mirna_targ_l)))
    omics=omics[-grep("gex",names(omics))]

    #in case of existing expressed lncRNA in the gex, then separate
    if(length(lncRNA_ids)>0){
      gex=gex_unf[!(rownames(gex_unf) %in% lncRNA_ids),]
      gex2lncRNA=gex_unf[rownames(gex_unf) %in% lncRNA_ids,]
    }else{
      gex=gex_unf
    }

    #in case the user has also provide a lncRNA count matrix, then bind
    if(("lncRNA" %in% names(omics)) & !is.null(gex2lncRNA)){
      lncRNA=omics[["lncRNA"]]
      lncRNA=rbind(lncRNA,gex2lncRNA)
      rm(gex2lncRNA)
    }
    if(("lncRNA" %in% names(omics)) & is.null(gex2lncRNA)){
      lncRNA=omics[["lncRNA"]]
    }
    if(!("lncRNA" %in% names(omics)) & !is.null(gex2lncRNA)){
      lncRNA=gex2lncRNA
    }

    #Check if genes and lncRNAs belong to the pathways, if yes then normalize and save
    gXp=sum(rownames(gex) %in% unique(unlist(databases$pathway_db$pathways_l)))
    if(gXp>0){
      message(">>number of genes falling into pathways: ", gXp)
      omics[["gex"]]=prepare_RNAseq(counts=gex, groups=info$Groups, winz=winz, type=type)
    }else{
      message(">>gene name example: ",databases$pathway_db$pathways_l[[1]][1])
      stop("no gene as rowname of the gex matrix belongs to any literature pathway");
    }

    if(!is.null(lncRNA) & nrow(lncRNA)>1){

      lncRNA=lncRNA[setdiff(seq(1,nrow(lncRNA)),grep("\\.[1-999]",rownames(lncRNA))),]
      lncRNA=lncRNA[setdiff(seq(1,nrow(lncRNA)),grep("Gm[1-999]",rownames(lncRNA))),]
      if(sum(duplicated(rownames(lncRNA)))!=0){
        lncRNA=remove_duplicates(lncRNA)
      }

      lXp=sum(rownames(lncRNA) %in% unique(names(databases$lncrna_net$lncrna_gene_targ_l)))
      if(lXp>0){
        message(">>number of lncRNAs having target sets of genes: ", lXp)
        omics[["lncRNA"]]=prepare_RNAseq(counts=lncRNA, groups=info$Groups, winz=winz, type=type)
      }else{
        message(">>lncRNA name example: ",unique(names(databases$lncrna_net$lncrna_gene_targ_l))[1])
      }
    }

    #remove temporary variables
    rm(gex,gex_unf,lncRNA,gXp,lXp)
  }

  #Process miRNA
  if(!is.null(omics[["miRNA"]])){
    message(">preparing miRNA")
    miRNA=omics[["miRNA"]]

    #Finds if there are pre miRNA
    mirnas=rownames(miRNA)
    mirna_map=miRBaseConverter::miRNA_PrecursorToMature(mirnas, version = "v22")
    del=is.na(mirna_map$Mature1) & is.na(mirna_map$Mature2)
    mirna_map=mirna_map[!del,]

    #Add the pre miRNAs as mature if they are not there already
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
      stop("no miRNA as rowname of the miRNA count matrix has a literature target set of genes");
    }

    #Library and Rank normalization
    miRNA=remove_duplicates(miRNA)
    miRNA=miRNA[setdiff(seq(1,nrow(miRNA)),grep("\\.[1-999]",rownames(miRNA))),]
    miRNA=prepare_RNAseq(counts=miRNA, groups=info$Groups, winz=winz, type=type)
    omics[["miRNA"]]=miRNA
    rm(miRNA,miRNA_mat,mirna_map,miXp)
  }


  #Transform somatic mutation omic by propagation
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

    #Improve the mutational signal of each sample (protein that are directly connected to mutated protein are likely to be dysfunctional)
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

    #Check if you have proper gene's names for mutations
    gXp=sum(rownames(muts) %in% unique(unlist(databases$pathway_db$pathways_l)))
    if(gXp>0){
      message(">>number of mutated genes falling into pathways: ", gXp)
    }else{
      message(">>mutated gene name example: ",databases$pathway_db$pathways_l[[1]][1])
      stop("no mutated gene as rowname of the mutation matrix belongs to any literature pathway");
    }

    #Propagate mutation information
    prop_muts=get_propagated(net=databases$gene_net$net_adj, counts=muts, n_cores=n_cores, r=0.8, keep_no_nodes=T)
    #Rank transform
    prop_muts=get_rank_01(m=prop_muts)
    omics[["mutations"]]=prop_muts
    rm(muts,g,mut_genes,neighs,freq_muts,gXp,prop_muts)
  }

  #Join omics for future analysis
  available_omics=setdiff(names(omics),"mutations")
  for(ome in available_omics){
    m1=omics[[ome]]

    if(ome == available_omics[1]){
      ms=m1
    }else{
      ms=rbind(ms,m1)
    }

    if(ome == available_omics[length(available_omics)]){
      omics[["all"]]=ms
    }
  }

  #Check cross validation setting
  if(sum(cv_probs)!=1){
    stop("the vector cv_probs has user-defined probabilites which sum is not 1")
  }else{
    message("Preparing partitions of data for cross validation setting")
  }

  #Over-riding the number of default cv based on the number of runs defined by the user
  if(!is.null(user_sets)){
    n_iter_cv=length(user_sets)
  }

  #Partitioning samples for cross validation
  partitions_l=list()
  cv_probs[1]=cv_probs[1]+cv_probs[2]
  for(iter_i in 1:n_iter_cv){
    #Get training and testing
    if(is.null(user_sets)){
      seed=seed+1
      indx <- splitTools::partition(info$Groups, type="stratified", n_bins=2,
                                    p=c(train = cv_probs[1], test=cv_probs[3]),
                                    seed=seed)

    }else{
      indx=user_sets[[iter_i]]
      train_check="train" %in% names(indx)
      test_check="test" %in% names(indx)
      if(train_check & test_check){
        user_sets=lapply(user_sets,function(x){
          res=list(train=(which(rownames(info) %in% unlist(x$train))),
                   test=(which(rownames(info) %in% unlist(x$test))))
          return(res)
        })
        indx=user_sets[[iter_i]]
        train_indxs=indx$train
        test_indxs=indx$test
        indxs=unique(c(train_indxs,test_indxs))
        indxs=indxs[order(indxs,decreasing = FALSE)]
        check1=check_vectors(indxs,seq(1,nrow(info)))==1
        if(check1){
          message(">>User's defined trained and testing/blind patients set correctly")
          indx=user_sets[[iter_i]]
        }else{
          stop("joining the user's train and test sets does not associate all the patient's columns")
        }
      }else{
        stop("a list in user_sets does not contein either the train or test set of indexes")
      }
      rm(train_check,test_check,train_indxs,test_indxs,indxs,check1)
    }

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
    #Pros: the model evaluates how much the training set is able to predict correctly samples that are similar to the testing ones
    #Cons: the model does not train/learn samples that are similar to the unknown ones
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
    te_omics=lapply(omics,function(x){x})


    TRAINING=list(info=tr_info,
                  omics=tr_omics,
                  databases=databases,
                  tax_id=tax_id,
                  n_cores=n_cores,
                  type=type)

    TESTING=list(info=info,
                 omics=te_omics,
                 databases=databases,
                 tax_id=tax_id,
                 n_cores=n_cores,
                 type=type)

    name_iter_i=paste("run",iter_i,sep="")
    partitions_l[[name_iter_i]]=list(training=TRAINING,testing=TESTING)
  }

  return(partitions_l)
}

#' Analyse training patients to find significant PSNs
#'
#' For each omic, it determines the best molecules to build the PSNs.
#' The molecules are fit into the pathways or regulatory target sets.
#' StellarPath finds the most enriched pathways.
#' A PSN is built for each pathway and is tested if separates the groups in comparison.
#' If yes, the PSN is kept and used for the training and classification.
#'
#' @param data_l The output list obtained from processing of the user's data with the function: prepare_data
#' @param max_size Default 150, integer value filtering out the pathways with more elements than this threshold
#' @param min_elements Default 4, integer value filtering out the pathways containing less SIGNIFICANT molecules than this threshold
#' @param n_top_sets Default 50, integer value that defines how many best pathways-specific PSNs to return as result
#' @param n_cores Default 2, integer value that defines the number of cores to use for running the function in parallel
#' @param keep_PSN Default FALSE, boolean that specifies to not keep the PSN generated from significantly enriched pathways
#' @param ban Default TRUE, remove excessive context-specific pathways (e.g. "SOUND", "STABILIZATION", "ODONTOGENESIS", "PREGNANCY")
#' @return
#' For each unique omic (e.g. gene expression), it provides:
#'
#' - the normalized count matrix
#'
#' - the SDR dataframe (obtained with find_SDR()) including the statistics of the molecules tested for separating the patient's groups.
#'
#' It also provide the result of the pathway analysis that combines the information of the different omics.
#' Precisely, a list of two elements accessible with $training$pathway_analysis: top_pathways_l and top_pathways_df.
#'
#' - top_pathways_df contains meta-information about the significantly enriched pathways which have been used to build a significant PSN.
#'
#' - top_pathways_l is the list s.t. each element refers to a significantly enriched pathway and contains:
#'
#' 1) the vector of significant molecules enriching the pathway
#'
#' 2) the count matrix used to build the PSN
#'
#' 3) statistics related to the PSN used to model the relationships between patients based on the specific pathway info
#'
#' @export
analyse_training = function(data_l, max_size=200, min_elements=5, n_top_pathways=100, n_cores=NULL, keep_PSN=FALSE, ban=TRUE){
  if(!is.null(n_cores)){
    data_l[[1]]$training$n_cores=n_cores
  }else{
    n_cores=data_l[[1]]$training$n_cores
  }

  message(">using the following cores: ", n_cores)
  #One core one chunk of pathways to process
  if(.Platform$OS.type == "unix") {
    #message(">>parallel for linux")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    #message(">>parallel for windows")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);

  #Iterate over training sets for finding significantly deregulated molecules, enriched pathways and patient similarity networks
  for(iter_cv in 1:length(data_l)){
    message("Feature selection of training data belonging to: ",iter_cv," run of cross validation")
    #Extract the data of training set of a run of cross validation
    i_tr_data_l=data_l[[iter_cv]]$training
    groups_freq=as.data.frame(table(i_tr_data_l$info$Groups))
    groups_name=unique(i_tr_data_l$info$Groups)
    rownames(groups_freq)=groups_freq[,1]
    path2data=paste(names(data_l[iter_cv]),"training",sep="$")
    message(">training samples:")
    message(">>first group: ",groups_freq[groups_name[1],1]," of training samples: ",groups_freq[groups_name[1],2])
    message(">>second group: ",groups_freq[groups_name[2],1]," of training samples: ",groups_freq[groups_name[2],2])
    message(">analysis:")
    #Gene expression analysis of training data
    if("gex" %in% names(i_tr_data_l$omics)){
      message(">>finding differentially expressed genes")
      #Detect significantly different genes in expression and stability, finally update the omic's data
      i_tr_data_l$omics$gex = list(norm_counts = i_tr_data_l$omics$gex,
                                   molecule_analysis = find_SDR(counts = i_tr_data_l$omics$gex,
                                                                groups = i_tr_data_l$info$Groups))

      message(">>finding enriched pathways represented by a significant PSN")
      #Detect which pathway fit most of the significant genes, then create and test its PSN
      #Keep the pathways which have the best patient similarity networks
      i_tr_data_l$omics$gex$pathway_analysis = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                                                        groups = i_tr_data_l$info$Groups,
                                                        pathways_l = i_tr_data_l$databases$pathway_db$pathways_l,
                                                        sign_df = i_tr_data_l$omics$gex$molecule_analysis,
                                                        pathway_type = "gene_canonical_pathway",
                                                        count_source = "gene_expression",
                                                        n_cores = n_cores,
                                                        max_size = max_size,
                                                        min_elements = min_elements,
                                                        n_top_pathways = n_top_pathways,
                                                        keep_PSN = keep_PSN,
                                                        ban = ban)

    }

    #miRNA expression analysis of training data
    SDR_miRNA_lncRNAs_l=NULL
    if("miRNA" %in% names(i_tr_data_l$omics)){
      message(">finding differentially expressed miRNAs")
      i_tr_data_l$omics$miRNA = list(norm_counts = i_tr_data_l$omics$miRNA,
                                     molecule_analysis=find_SDR(counts = i_tr_data_l$omics$miRNA,
                                                                groups = i_tr_data_l$info$Groups,
                                                                top_molecules_th = 0.15))

      #Find functional target sets of SDR miRNA with SDR genes
      #Filter the list of miRNA's targets to keep only the sets of the significant miRNAs
      #with only the significant targets
      SDR_miRNA_genes_l=list()
      #For each significant miRNA
      miRNA_tT=i_tr_data_l$omics$miRNA$molecule_analysis;miRNA_tT=miRNA_tT[miRNA_tT$significant==TRUE,]
      #keep only the anti-correlated significant gene targets
      gex_tT=i_tr_data_l$omics$gex$molecule_analysis;gex_tT=gex_tT[gex_tT$significant==TRUE,]

      #Extract the directions of the miRNAs
      barcodes=unique(miRNA_tT$sign_type)
      for(barcode in barcodes){
        #Extract the anti-correlated
        miRNA_tT_dir=miRNA_tT[miRNA_tT$sign_type==barcode,]
        gex_tT_dir=gex_tT[gex_tT$sign_type==(-barcode),]

        #Filter the target list
        miRNA_genes_l=i_tr_data_l$databases$mirna_net
        SDR_miRNA_genes_l1=miRNA_genes_l[which(names(miRNA_genes_l) %in% rownames(miRNA_tT_dir))]
        SDR_miRNA_genes_l1=lapply(SDR_miRNA_genes_l1,function(x){
          y=x[x %in% rownames(gex_tT_dir)]
          return(y)
        })

        SDR_miRNA_genes_l=c(SDR_miRNA_genes_l,SDR_miRNA_genes_l1)
      }
      SDR_miRNA_genes_l=SDR_miRNA_genes_l[sapply(SDR_miRNA_genes_l,length)>0]

      message(">>finding enriched miRNA's gene targets passing ecdf test")
      miRNA2gene_ecdf_df = find_cdf_pathways(nc_l = i_tr_data_l$databases$mirna_net,
                                             nc_sign_l = SDR_miRNA_genes_l,
                                             nc_tT_ov = i_tr_data_l$omics$miRNA$molecule_analysis,
                                             target_tT_ov = i_tr_data_l$omics$gex$molecule_analysis,
                                             pathway_type = "miRNA_targets",
                                             count_source = "gene_expression",
                                             anti_corr = TRUE)

      message(">>finding enriched miRNA's gene targets represented by a significant PSN")
      miRNA2gene_pa_l = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                                 groups = i_tr_data_l$info$Groups,
                                 pathways_l = SDR_miRNA_genes_l,
                                 sign_df = i_tr_data_l$omics$gex$molecule_analysis,
                                 pathway_type = "miRNA_targets",
                                 count_source = "gene_expression",
                                 n_cores = n_cores,
                                 max_size = 9999,
                                 min_elements=min_elements,
                                 n_top_pathways = n_top_pathways,
                                 keep_PSN=keep_PSN,
                                 ban=ban)

      miRNA2gene_pa_l$top_pathways_el=combine_test_results(miRNA2gene_pa_l$top_pathways_el,miRNA2gene_ecdf_df)
      i_tr_data_l$omics$miRNA$pathway_analysis$miRNA2gene=miRNA2gene_pa_l

      #Find functional target sets of SDR miRNA with SDR lncRNAs
      miRNA_lncRNAs_l=i_tr_data_l$databases$lncrna_net$mirna_lncrna_targ_l

      if(!is.null(miRNA_tT_dir)){
        #For each significant miRNA
        SDR_miRNA_lncRNAs_l=miRNA_lncRNAs_l[which(names(miRNA_lncRNAs_l) %in% rownames(miRNA_tT_dir))]
        SDR_miRNA_lncRNAs_l=SDR_miRNA_lncRNAs_l[sapply(SDR_miRNA_lncRNAs_l,length)>0]
      }

      #Clean
      rm(miRNA_tT,gex_tT,miRNA_lncRNAs_l,miRNA_genes_l,miRNA2gene_pa_l)
    }

    #lncRNA analysis
    if("lncRNA" %in% names(i_tr_data_l$omics)){
      message(">finding differentially expressed lncRNAs")
      i_tr_data_l$omics$lncRNA = list(norm_counts = i_tr_data_l$omics$lncRNA,
                                      molecule_analysis=find_SDR(counts = i_tr_data_l$omics$lncRNA,
                                                                 groups = i_tr_data_l$info$Groups,
                                                                 top_molecules_th = 0.15))

      #Find functional target sets of SDR lncRNAs with SDR genes/miRNAs
      #Filter the list of lncRNA's targets to keep only the sets of the significant lncRNAs
      #with only the significant targets
      lncRNA_genes_l=i_tr_data_l$databases$lncrna_net$lncrna_gene_targ_l
      lncRNA_miRNA_l=i_tr_data_l$databases$lncrna_net$lncrna_mirna_targ_l

      #For each significant lncRNA
      lncRNA_tT=i_tr_data_l$omics$lncRNA$molecule_analysis;lncRNA_tT=lncRNA_tT[lncRNA_tT$significant==TRUE,]
      #Keep the sets of gene targets and miRNA targets belonging to only the significant lncRNAs
      SDR_lncRNA_genes_l=lncRNA_genes_l[which(names(lncRNA_genes_l) %in% rownames(lncRNA_tT))]
      SDR_lncRNA_miRNAs_l=lncRNA_miRNA_l[which(names(lncRNA_miRNA_l) %in% rownames(lncRNA_tT))]
      #No need to remove non-significant genes and miRNAs because will not be used to create the PSNs
      #and because lncRNA (on the contrary of miRNAs) can target both in corr and anti-corr.
      SDR_lncRNA_miRNAs_l=SDR_lncRNA_miRNAs_l[sapply(SDR_lncRNA_miRNAs_l,length)>0]
      SDR_lncRNA_genes_l=SDR_lncRNA_genes_l[sapply(SDR_lncRNA_genes_l,length)>0]

      message(">>finding enriched lncRNA's gene targets passing ecdf test")
      lncRNA2gene_ecdf_df = find_cdf_pathways(nc_l = i_tr_data_l$databases$lncrna_net$lncrna_gene_targ_l,
                                              nc_sign_l = SDR_lncRNA_genes_l,
                                              nc_tT_ov = i_tr_data_l$omics$lncRNA$molecule_analysis,
                                              target_tT_ov = i_tr_data_l$omics$gex$molecule_analysis,
                                              pathway_type = "lncRNA_targets",
                                              count_source = "gene_expression")

      message(">>finding enriched lncRNA's gene targets represented by a significant PSN")
      lncRNA2gene_pa_l = find_PSN(counts = i_tr_data_l$omics$gex$norm_counts,
                                  groups = i_tr_data_l$info$Groups,
                                  pathways_l = SDR_lncRNA_genes_l,
                                  sign_df = i_tr_data_l$omics$gex$molecule_analysis,
                                  pathway_type = "lncRNA_targets",
                                  count_source = "gene_expression",
                                  n_cores = n_cores,
                                  max_size=9999,
                                  min_elements=min_elements,
                                  n_top_pathways = n_top_pathways,
                                  keep_PSN = keep_PSN,
                                  ban=ban)

      lncRNA2gene_pa_l$top_pathways_el=combine_test_results(lncRNA2gene_pa_l$top_pathways_el,lncRNA2gene_ecdf_df)
      i_tr_data_l$omics$lncRNA$pathway_analysis$lncRNA2gene=lncRNA2gene_pa_l
      rm(lncRNA2gene_pa_l,lncRNA2gene_ecdf_df)

      message(">>finding enriched lncRNA's miRNA targets passing ecdf test")
      lncRNA2miRNA_ecdf_df = find_cdf_pathways(nc_l = i_tr_data_l$databases$lncrna_net$lncrna_mirna_targ_l,
                                               nc_sign_l = SDR_lncRNA_miRNAs_l,
                                               nc_tT_ov = i_tr_data_l$omics$lncRNA$molecule_analysis,
                                               target_tT_ov = i_tr_data_l$omics$miRNA$molecule_analysis,
                                               pathway_type = "lncRNA_targets",
                                               count_source = "miRNA_expression")

      message(">>finding enriched lncRNA's miRNA targets represented by a significant PSN")
      lncRNA2miRNA_pa_l = find_PSN(counts = i_tr_data_l$omics$miRNA$norm_counts,
                                   groups = i_tr_data_l$info$Groups,
                                   pathways_l = SDR_lncRNA_miRNAs_l,
                                   sign_df = i_tr_data_l$omics$miRNA$molecule_analysis,
                                   pathway_type = "lncRNA_targets",
                                   count_source = "miRNA_expression",
                                   n_cores = n_cores,
                                   max_size=9999,
                                   min_elements=1,
                                   n_top_pathways = n_top_pathways,
                                   keep_PSN = keep_PSN,
                                   ban=ban)

      lncRNA2miRNA_pa_l$top_pathways_el=combine_test_results(lncRNA2miRNA_pa_l$top_pathways_el,lncRNA2miRNA_ecdf_df)
      i_tr_data_l$omics$lncRNA$pathway_analysis$lncRNA2miRNA=lncRNA2miRNA_pa_l
      rm(lncRNA2miRNA_pa_l,lncRNA2miRNA_ecdf_df)

      if(!is.null(SDR_miRNA_lncRNAs_l)){
        message(">>finding enriched miRNA's lncRNA targets passing ecdf test")
        miRNA2lncRNA_ecdf_df = find_cdf_pathways(nc_l = i_tr_data_l$databases$lncrna_net$mirna_lncrna_targ_l,
                                                 nc_sign_l = SDR_miRNA_lncRNAs_l,
                                                 nc_tT_ov = i_tr_data_l$omics$miRNA$molecule_analysis,
                                                 target_tT_ov = i_tr_data_l$omics$lncRNA$molecule_analysis,
                                                 pathway_type = "miRNA_targets",
                                                 count_source = "lncRNA_expression")

        message(">>finding enriched miRNA's lncRNA targets represented by a significant PSN")
        miRNA2lncRNA_pa_l = find_PSN(counts = i_tr_data_l$omics$lncRNA$norm_counts,
                                     groups = i_tr_data_l$info$Groups,
                                     pathways_l = SDR_miRNA_lncRNAs_l,
                                     sign_df = i_tr_data_l$omics$lncRNA$molecule_analysis,
                                     pathway_type = "miRNA_targets",
                                     count_source = "lncRNA_expression",
                                     n_cores = n_cores,
                                     max_size=9999,
                                     min_elements=1,
                                     n_top_pathways = n_top_pathways,
                                     keep_PSN = keep_PSN,
                                     ban=ban)

        miRNA2lncRNA_pa_l$top_pathways_el=combine_test_results(miRNA2lncRNA_pa_l$top_pathways_el,miRNA2lncRNA_ecdf_df)
        i_tr_data_l$omics$miRNA$pathway_analysis$miRNA2lncRNA=miRNA2lncRNA_pa_l
        rm(miRNA2lncRNA_pa_l,miRNA2lncRNA_ecdf_df)
      }
    }

    #Combine all the pathway results, update training object and return
    message(">adding significant pathways/targets/PSNs into $training$pathway_analysis")
    i_tr_data_l$pathway_analysis=merge_pathway_analysis(omics = i_tr_data_l$omics, start = path2data)

    #Somatic mutation analysis of training data
    if("mutations" %in% names(i_tr_data_l$omics)){
      message(">>finding differentially mutated genes")
      i_tr_data_l$omics$mutations = list(norm_counts = i_tr_data_l$omics$mutations,
                                         molecule_analysis = find_SDR(counts = i_tr_data_l$omics$mutations,
                                                                      groups = i_tr_data_l$info$Groups))

      message(">>finding altered pathways by mutations")
      #There is no need to test the PSNs because they will be all tested already with the gene expression matrix
      i_tr_data_l$pathway_analysis = find_mutated_pathways(pathways_l = i_tr_data_l$databases$pathway_db$pathways_l,
                                                           sign_df = i_tr_data_l$omics$mutations$molecule_analysis[i_tr_data_l$omics$mutations$molecule_analysis$sign_type>0,],
                                                           gex = i_tr_data_l$pathway_analysis,
                                                           pathway_type = "gene_canonical_pathway",
                                                           count_source = "mutation")

    }

    data_l[[iter_cv]]$training=i_tr_data_l
  }

  parallel::stopCluster(cl)
  return(data_l)
}

#' Classify testing patients based on training ones in top pathway-based PSNs
#'
#' For each run of cross-validation, extract the pathways/PSNs that are significant for training data.
#' Generate the PSN, format it as edge list and pass the network to the graph convolutional network for the training and prediction.
#'
#' @param data_l output list obtained from the processing of the user-defined data with the function: prepare_data and analyse_training
#' @param n_top_PSNs Default 30, cutoff that defines how many PSNs to use to classify
#' @param py_path character string of the path pointing to the python script with the implementation of the GCN
#' @param speed Default fast, character string: fast/medium/slow that changes the parameters of the GCN to run on the PSNs, select according to how much time you want to allow a GCN to train
#' @return Updates the info dataframe in $testing$info based on the predictions made for the testing and validation samples.
#' It adds a new list to $testing called classification which contains the dataframe perfsXpathway_df and the matthew correlation coefficent mcc_score.
#' perfsXpathway_df is a dataframe that shows the performances of classification obtained with each pathway/PSN.
#' mcc_score is a double value between -1 and 1 representing the performance of classification.
#' @importFrom reticulate source_python
#' @importFrom scales rescale
#' @importFrom plyr ddply
#' @importFrom plyr .
#' @export
classify_testing = function(data_l,n_top_PSNs=100,py_path=NULL,speed="fast"){
  if(is.null(py_path)){
    if(nrow(data_l[[1]]$testing$info)<=30){
      speed="slow"
    }
    if(speed=="fast"){
      py_path=system.file("python", "sage_classification_individual_fast.py", package = "StellarPath")
    }
    if(speed=="medium"){
      py_path=system.file("python", "sage_classification_individual_medium.py", package = "StellarPath")
    }
    if(speed=="slow"){
      py_path=system.file("python", "sage_classification_individual_slow.py", package = "StellarPath")
    }
  }

  #Get how many times a significant pathway has been selected during the cv
  pathway_id=c("pathways","sign_type","group")
  for(iter_cv in 1:length(data_l)){
    i_tr_data_l=data_l[[iter_cv]]$training
    i_te_data_l=data_l[[iter_cv]]$testing
    info_df=i_te_data_l$info

    #Initiate the lists that will contain the PSNs and their names
    all_edges <- list()
    all_names <- list()

    #Get the stats of pathways selected with the multi-omics training patients
    top_pathways_df=i_tr_data_l$pathway_analysis$top_pathways_df[,c(pathway_id,"path2pathway")]
    colnames(top_pathways_df)[4]=paste(names(data_l[iter_cv]),"path2pathway",sep="_")
    if(iter_cv==1){
      top_pathways_freq=top_pathways_df[,pathway_id]
      top_pathways_dfs=top_pathways_df
    }else{
      top_pathways_freq=rbind(top_pathways_freq,top_pathways_df[,pathway_id])
      top_pathways_dfs=merge(top_pathways_dfs,top_pathways_df,by=pathway_id,all=TRUE)
    }
  }
  top_pathways_freq <- plyr::ddply(top_pathways_freq,
                                   plyr::.(top_pathways_freq[,pathway_id[1]],
                                           top_pathways_freq[,pathway_id[2]],
                                           top_pathways_freq[,pathway_id[3]]),
                                   nrow)
  colnames(top_pathways_freq)=c(pathway_id,"percent_cv_runs")
  top_pathways_freq$percent_cv_runs=(top_pathways_freq$percent_cv_runs*100)/iter_cv
  top_pathways_tr_stats=merge(top_pathways_dfs,top_pathways_freq,by=pathway_id,all=TRUE)
  rm(top_pathways_freq,top_pathways_dfs,top_pathways_df)

  #For each cv run, build PSN with testing patients in the best pathways, train and predict
  for(iter_cv in 1:length(data_l)){
    i_tr_data_l=data_l[[iter_cv]]$training
    i_te_data_l=data_l[[iter_cv]]$testing
    info_df=i_te_data_l$info

    #Initiate the lists that will contain the PSNs and their names
    all_edges <- list()
    all_names <- list()

    #Get the stats of pathways selected with the multi-omics training patients
    top_pathways_df=i_tr_data_l$pathway_analysis$top_pathways_df
    #Get the matrix with all the omics combined
    m=i_te_data_l$omics$all

    if(nrow(top_pathways_df)<n_top_PSNs){
      n_top_PSNs=nrow(top_pathways_df)
    }

    #Order the pathways based on the best parameters
    top_pathways_df=top_pathways_df[order(top_pathways_df$pathway_type,
                                          top_pathways_df$regulation,
                                          -top_pathways_df$power),]
    #Iterate over the first 100 best pathways with the multi-omics training patients
    for(i_path in seq(1,n_top_PSNs)){

      path2pathway=top_pathways_df$path2pathway[i_path]
      path_l=get_pathway_data(data_l,path2pathway)

      #Recover stats from the training data
      name_psn=top_pathways_df$pathways[i_path]
      gs=path_l$set
      dir=path_l$stats$regulation

      #Extract the single molecule's values and build the PSN with training, validation and testing samples
      m_gs=m[rownames(m) %in% gs,]
      if(dir=="activated"){
        PSN=build_PSN(m_gs, up=TRUE)
      }else{
        PSN=build_PSN(m_gs, up=FALSE)
      }
      #Filter low similarities, convert to edge list and 01 standardize their values
      PSN=set_weak_edges_to_zero(adj_matrix = PSN)
      el_PSN=convert_adj2edg(PSN)
      el_PSN=el_PSN[el_PSN$weight!=0,]
      el_PSN$weight=scales::rescale(el_PSN$weight,to=c(0,1))

      #Save
      edges <- list(el_PSN)
      all_edges <- c(all_edges,edges)
      all_names <- c(all_names,name_psn)

    }

    #Load and Compile function in Python script
    reticulate::source_python(py_path, convert = TRUE)
    #Function for node classification with Sage that trains on all over the database
    res <- node_classification(all_names,info_df,all_edges)

    #Extract and Format results from the Sage classification
    perfsXpathway_df=res[[1]][,c(1,2)]
    colnames(perfsXpathway_df)[1:2]=c("pathways","val_mcc")
    keep=which(perfsXpathway_df$val_mcc>0)
    perfsXpathway_df=perfsXpathway_df[keep,]
    pathway_names=perfsXpathway_df$pathways
    perfsXpathway_df=merge(top_pathways_df[,pathway_id],perfsXpathway_df,by="pathways")
    perfsXpathway_df=perfsXpathway_df[match(pathway_names,perfsXpathway_df$pathways),]

    #For the validation patients
    validation_predictions=as.data.frame(apply(res[[2]],2,unlist))
    validation_predictions=validation_predictions[keep,]
    consensus_validation_pred_labels=apply(validation_predictions,2,function(x){
      freqs_df=as.data.frame(table(x))
      pred_label=as.character(freqs_df[which.max(freqs_df$Freq),1])
      return(pred_label)
    })

    #Compute performance metrics per pathway
    val_true_labels=info_df$Groups[info_df$val_indxs==1]
    for(i in 1:nrow(validation_predictions)){
      predicted_labels=validation_predictions[i,]
      val_perfs_df=get_perf_metrics(val_true_labels, predicted_labels)
      if(i==1){
        val_perfs_dfs=val_perfs_df
      }else{
        val_perfs_dfs=rbind(val_perfs_dfs,val_perfs_df)
      }
    }

    #For the testing patients
    test_predictions=as.data.frame(apply(res[[3]],2,unlist))
    test_predictions=test_predictions[keep,]
    consensus_test_pred_labels=apply(test_predictions,2,function(x){
      freqs_df=as.data.frame(table(x))
      pred_label=as.character(freqs_df[which.max(freqs_df$Freq),1])
      return(pred_label)
    })

    #Compute performance metrics per pathway
    test_true_labels=info_df$Groups[info_df$test_indxs==1]
    for(i in 1:nrow(test_predictions)){
      predicted_labels=test_predictions[i,]
      test_perfs_df=get_perf_metrics(test_true_labels, predicted_labels)
      if(i==1){
        test_perfs_dfs=test_perfs_df
      }else{
        test_perfs_dfs=rbind(test_perfs_dfs,test_perfs_df)
      }
    }
    colnames(test_perfs_dfs)=paste("test_",colnames(test_perfs_dfs),sep="")

    likelihoodsXpat=c(0)
    for(i_lab in 1:length(test_true_labels)){
      true_lab=test_true_labels[i_lab]
      likelihoodsXpat=c(likelihoodsXpat,(sum(test_predictions[,i_lab]==true_lab)*100)/nrow(test_predictions))
    }
    likelihoodsXpat=likelihoodsXpat[-1]

    case_group=unique(info_df$Groups)[1]
    pred_probs=apply(test_predictions,2,function(x){
      freqs_df=as.data.frame(table(x))
      pred_label=as.character(freqs_df[which.max(freqs_df$Freq),1])

      pred_prob=freqs_df$Freq[freqs_df$x==pred_label]/sum(freqs_df$Freq)
      if(pred_label!=case_group){
        pred_prob=1-pred_prob
      }
      return(pred_prob)
    })

    roc_l=get_auc(test_true_labels,pred_probs)

    #Update info dataframe and compute MCC classification performance
    info_df$val_preds="not_used"
    info_df$val_preds[info_df$val_indxs==1]=consensus_validation_pred_labels
    info_df$test_preds="not_used"
    info_df$test_preds[info_df$test_indxs==1]=consensus_test_pred_labels
    info_df$likelihood=NA
    info_df$likelihood[info_df$test_indxs==1]=likelihoodsXpat

    #Get final and overall results
    perfsXpathway_df=cbind(perfsXpathway_df,test_perfs_dfs)
    perfsXpathway_df$name_run_cv=names(data_l[iter_cv])
    perfs_score=get_perf_metrics(info_df$Groups[info_df$test_indxs==1],consensus_test_pred_labels)
    perfs_score$auc_roc=roc_l$auc

    #Save
    message("MCC score of the run: ",iter_cv ," is: ", perfs_score[1])
    message("Updating $testing$classification in testing data list")
    i_te_data_l[["info"]]=info_df
    i_te_data_l[["classification"]]=list(perfsXpathway_df=perfsXpathway_df,
                                         perfs_score=perfs_score,
                                         roc_l=roc_l)
    data_l[[iter_cv]]$testing=i_te_data_l
  }

  #Retrieve MCC scores from the CV runs
  perfs_scores=sapply(data_l[1:iter_cv],function(x){
    x$testing$classification$perfs_score
  })
  perfs_scores=t(perfs_scores)
  perfs_scores=as.data.frame(perfs_scores)
  perfs_scores$name_dataset=as.character(Sys.time())
  perfs_scores$runs=rownames(perfs_scores)
  perfs_scores=perfs_scores[,c(6,7,seq(1,5))]

  roc_ls=lapply(data_l[1:iter_cv],function(x){
    x$testing$classification$roc_l
  })

  #Retrieve performances of the model on each training feature selected pathway
  perfsXpathway_dfs=lapply(data_l[1:iter_cv],function(x){
    x$testing$classification$perfsXpathway_df
  })
  perfsXpathway_dfs=as.data.frame(rbindlist(perfsXpathway_dfs))
  perfsXpathway_dfs=perfsXpathway_dfs[!duplicated(perfsXpathway_dfs[,-ncol(perfsXpathway_dfs)]),]
  perfsXpathway_dfs=merge(perfsXpathway_dfs,top_pathways_tr_stats,by=pathway_id,all=TRUE)
  dups=perfsXpathway_dfs$pathways[duplicated(perfsXpathway_dfs$pathways)]
  for(dup in dups){
    dup_df=perfsXpathway_dfs[perfsXpathway_dfs$pathways==dup,]

    i_dup=1
    for(i_dup in 1:nrow(dup_df)){
      name_run_cv=dup_df$name_run_cv[i_dup]
      dup_df[i_dup,setdiff(grep("run[1-999]",colnames(dup_df)),grep(name_run_cv,colnames(dup_df)))]=NA
    }

    perfsXpathway_dfs[perfsXpathway_dfs$pathways==dup,]=dup_df
  }
  perfsXpathway_dfs=perfsXpathway_dfs[,-grep("name_run_cv",colnames(perfsXpathway_dfs))]

  #Save performances of the multiple runs of cross validations
  data_l[["performances"]]=list(performances=perfs_scores,
                                roc_ls=roc_ls,
                                perfsXpathway_dfs=perfsXpathway_dfs)

  return(data_l)
}

#' Enrichment analysis with predicted data from cross-validation
#'
#' Find the run of cross-validation in which the model predicted at best.
#' Assign the unknown/testing patients to their predicted class.
#' Find the pathways-based PSNs that separate the patient's classes with all the patients.
#' Combine the info about the significant PSNs coming from this analysis with the info coming from the cross-validation phase.
#'
#' @param data_l data list updated from the function "classify_testing"
#' @return data list with the new element "enrichment" that contains resulting dataframe of pathways
#' @export
enrichment_analysis = function(data_l){
  #Select the best run of cross validation
  best_run=which.max(data_l$performances$performances$mcc)

  #Get the samples info and replace the original groups with the model's predictions
  message("Enrichment with predictions of the run: ", best_run, " of cross validation with best performances")
  info=data_l[[best_run]]$testing$info
  info$Groups[info$test_preds!="not_used"]=info$test_preds[info$test_preds!="not_used"]

  #Create a list object to fit in the analyse_training function
  enr_l=list(info=info,
             omics=data_l[[best_run]]$testing$omics,
             databases=data_l[[best_run]]$testing$databases,
             tax_id=data_l[[best_run]]$testing$tax_id,
             n_cores=data_l[[best_run]]$testing$n_cores,
             type=data_l[[best_run]]$testing$type)
  data_l1=list()
  data_l1[["run1"]]=list(training=enr_l)

  #Get the significant pathway-based PSNs with the predicted patients
  data_l1=analyse_training(data_l1, keep_PSN = TRUE)
  #Extract the results
  enr_l=data_l1[[1]]$training

  #Change the pointer of the pathways selected in the enrichment phase
  #to not overlap with the pointers refering to the pathways' data used in the cv runs
  top_pathways_df_in_enr=enr_l$pathway_analysis$top_pathways_df
  top_pathways_df_in_enr$path2pathway=gsub("run1\\$training","enrichment",top_pathways_df_in_enr$path2pathway)

  #Update the enrichment results with the info of the pathways-based PSN
  #from the real training and testing phase
  top_pathways_df_in_cv=data_l$performances$perfsXpathway_dfs[,-c(2,3)]
  top_pathways_df_in_enr=merge(top_pathways_df_in_enr,top_pathways_df_in_cv,by=c("pathways"),all.x=TRUE)
  top_pathways_df_in_enr=top_pathways_df_in_enr[order(top_pathways_df_in_enr$power,top_pathways_df_in_enr$freq,decreasing = TRUE),]
  enr_l$pathway_analysis$top_pathways_df=top_pathways_df_in_enr
  data_l[["enrichment"]]=enr_l

  return(data_l)
}


#' Infer the class of new patients and find significant PSNs
#'
#' Infer the class of new patients and find significant PSNs
#'
#' @param info Dataframe, two character columns, first has sample's IDs (have to match with column names of the count matrices),
#' second has sample's groups (only two as for a pairwise DE analysis, for example c(AD,AD,AD,HEALTHY,HEALTHY,HEALTHY))
#' @param groups_name Character vector of length two, it must contain the same labels of unique(info[,2]), tip: first the case group like AD,
#' while the second name refers to the control group like HEALTHY
#' @param gex Numeric matrix of raw count values from bulk RNA sequencing (e.g. genes x samples) (it may include also lncRNAs)
#' @param lncRNA Numeric matrix of raw count values from bulk RNA sequencing (e.g. lncRNA x samples)
#' @param miRNA Numeric matrix of raw count values from miRNA sequencing (e.g. miRNA x samples)
#' @param mutations Numeric matrix of binary values from somatic mutation data (e.g. gene X samples, 1 per mutated gene, 0 otherwise)
#' @param tax_id Default 9606, integer value indicating taxonomic id of the samples, either 9606 (human) or 10090 (mouse)
#' @param simplify_names Default TRUE, simplify names of the samples in the info and count matrices
#' @param cv_probs Default c(0.7,0.1,0.2), vector of doubles indicating the proportion of the samples to use as training, validation and testing set
#' @param user_sets Default NULL, list of elements (one element is used to compose the training and testing set of a cv run) s.t. each element is a list containing two numeric vectors, one called train and one called test, each vector indicate the column index of those patients falling into training or testing set
#' @param data_l Default NULL, results list of data produced with StellarPath workflow on a previous patient dataset
#' @param n_cores Default 2, integer value greater than 0 indicating the number of cores to use to parallelize and speed up the operations
#' @param seed Default 5, integer value indicating the seed to replicate the sampling of the cross validation
#' @param winz Default FALSE, boolean to indicate if to adjuste outliers with winsorizing
#' @param type Default total rank, character string which defines the type of transformation applied to the count data, options in the help of the get_rank_01 function
#' @return A list s.t. each element contains the user-input prepared data for a run of cross-validation.
#' The run's data are divided into training and testing sets.
#' Each set contains a specific subset of sample's information (i.e. info) and normalized count matrices (i.e. omics).
#' @export
#'
infer_new_data = function(info, groups_name=NULL,
                          gex=NULL, lncRNA=NULL, miRNA=NULL, mutations=NULL,
                          tax_id=9606, simplify_names=TRUE,
                          cv_probs=c(0.7,0.1,0.2),
                          user_sets=NULL, data_l=NULL,
                          n_cores=2, seed=5,
                          winz=FALSE, type="sample_rank"){

  #Prepare data with the new patients
  new_data_l=prepare_data(info = info, groups_name = groups_name,
                          gex = gex, miRNA = miRNA, lncRNA=lncRNA, mutations = mutations,
                          tax_id = tax_id,
                          user_sets = user_sets,
                          n_cores = n_cores, seed=seed,
                          winz=winz, type=type);tmp=gc();

                          #Analyse the data
                          new_data_l=analyse_training(new_data_l)

                          #Filter the significant pathway-based PSN to keep only those used during the training
                          new_top_pathways_df=new_data_l$run1$training$pathway_analysis$top_pathways_df
                          cv_top_pathways_df=data_l$enrichment$pathway_analysis$top_pathways_df
                          new_top_pathways_df=new_top_pathways_df[new_top_pathways_df$pathways %in% cv_top_pathways_df$pathways,]
                          new_data_l$run1$training$pathway_analysis$top_pathways_df=new_top_pathways_df

                          #Classify the new patients
                          new_data_l=classify_testing(new_data_l)

                          #Enrichment
                          new_data_l=enrichment_analysis(new_data_l)
                          data_l[["inference"]]=new_data_l

                          return(data_l)
}

