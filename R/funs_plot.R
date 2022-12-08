#' Plot a Patient Similarity Network
#'
#' Plot a Patient Similarity Network
#'
#' @param info Dataframe, two character columns, first has sample's IDs (have to match with column names of the count matrices),
#' second has sample's groups (only two as a pairwise DE analysis, for example c(AD,AD,AD,HEALTHY,HEALTHY,HEALTHY))
#' @param PSN Patient Similarity Network in Adjacent matrix format
#' @param k_clusters Default 10, number of nodes to visualize for each group
#' @param edge_threshold Default 0.3, similarities lower than 0.3 will be hide in order to simplify the plot
#' @param edge_size Default 5, size of the edges
#' @param label.cex Default 3, size of the node's labels
#' @param legend.cex Default 0.5, size of the legend text
#' @param edge_cols Defulat c("#03c800","#96b5a3","#9ab8ff"), vector of three exadecimal colours, first colour is used for first group's
#' intrasimilarities, second colour for intersimilarities and third colour for second group's intrasimilarities
#' @param node_cols Default c("#9affb3","#9aebff"), vector of two exadecimal colours, first colour is used for first group's nodes
#' @param file_path Default "PSN", indicates the path and name of the plot file
#' @param file_type Default "pdf", indicates the format of the plot file
#' @param file_width Default 14, integer that indicates the width of the plot file
#' @param file_height Default 16, integer that indicates the height of the plot file
#' @param GLratio Relative size of the graph compared to the layout. Defaults to 2.5
#' @return Save the plot on the disk and return a list with the qgraph object and the list of the patient's clusters
#' @importFrom stats as.dist
#' @importFrom stats kmeans
#' @importFrom FNN get.knnx
#' @importFrom scales rescale
#' @importFrom plyr mapvalues
#' @import qgraph qgraph
#' @export
#'
plot_PSN = function(info, PSN, k_clusters=10,
                    edge_threshold = 0.3, edge_size=5, label.cex = 3, legend.cex=0.5,
                    edge_cols=c("#03c800","#96b5a3","#9ab8ff"),
                    node_cols=c("#9affb3","#9aebff"),
                    file_path="PSN", file_type="pdf",
                    file_width=14,file_height=16,GLratio=2
){
  colnames(info)[1:2] = c("IDs", "Groups")
  group_freq_df = as.data.frame(table(info$Groups))
  group_freq_df = group_freq_df[unique(match(info$Groups,
                                             group_freq_df$Var1)), ]

  IDGroups = group_freq_df$Var1
  Freqs = group_freq_df$Freq
  seq1cl = seq(1, Freqs[1])
  seq2cl = seq(Freqs[1] + 1, length(info$Groups))
  PSN2 = build_adjWJ_PSN(PSN)
  nodes_v = c("start")
  for (i in 1:2) {
    if (i == 1) {
      n_nodes = length(seq1cl)
      PSN2G1 = PSN2[seq1cl, seq1cl]
    }
    else {
      n_nodes = length(seq2cl)
      PSN2G1 = PSN2[seq2cl, seq2cl]
    }
    if (n_nodes > k_clusters) {
      k = k_clusters
    }
    else {
      k = n_nodes - 1
    }
    PSN2G1 = stats::as.dist(PSN2G1)
    kres = stats::kmeans(PSN2G1, k)
    kcentroids <- FNN::get.knnx(PSN2G1, kres$centers, 1)$nn.index
    colnames(kcentroids) = "centroids"
    rownames(kcentroids) = paste("cl", seq(1, length(unique(kres$cluster))),
                                 sep = "")
    cls_v = roots = "start"
    for (k_cl in 1:length(unique(kres$cluster))) {
      cluster = names(kres$cluster)[kres$cluster == k_cl]
      root = names(kres$cluster)[kcentroids[k_cl, 1]]
      for (l1 in seq(5, 100, 5)) {
        if (length(cluster) > l1) {
          cluster = append(cluster, "\n", after = l1)
        }
      }
      cluster = paste(cluster, collapse = ", ")
      cls_v = c(cls_v, cluster)
      roots = c(roots, root)
    }
    cls_v = cls_v[-1]
    roots = roots[-1]
    names(cls_v) = roots
    nodes_v = c(nodes_v, cls_v)
  }
  nodes_v = nodes_v[-1]
  keep = which(colnames(PSN) %in% names(nodes_v))
  info = info[keep, ]
  PSN = PSN[keep, keep]

  names_in_l=strsplit(nodes_v,split=",",fixed = TRUE)
  names_in_l=sapply(names_in_l,function(x){gsub(" ","",x)})
  names_in_k_l=sapply(names_in_l,function(x){y=x[nchar(x)>2];return(y)})

  #Recover stats of the new compact PSN
  group_freq_df = as.data.frame(table(info$Groups))
  group_freq_df = group_freq_df[unique(match(info$Groups, group_freq_df$Var1)), ]
  IDGroups = group_freq_df$Var1
  Freqs = group_freq_df$Freq
  seq1cl = seq(1, Freqs[1])
  seq2cl = seq(Freqs[1] + 1, length(info$Groups))
  eigs = get_cohesive_score(PSN, info$Groups)
  centrality = scales::rescale(eigs[, "centrality"], to = c(2,5))

  #Prepare adjacent matrix for plot
  grs_l=list(x=seq1cl,y=seq2cl)
  names(grs_l)=IDGroups

  col_m=PSN
  col_m[seq1cl,seq1cl]=edge_cols[1]
  col_m[seq2cl,seq2cl]=edge_cols[3]
  col_m[seq1cl,seq2cl]=edge_cols[2]
  col_m[seq2cl,seq1cl]=edge_cols[2]

  v_size=colnames(col_m)
  v_size=centrality[v_size]

  qgraph::qgraph(PSN, threshold = edge_threshold, edgelist = FALSE, weighted = TRUE, directed = FALSE,
                 layout = "spring", repulsion = 0.7,
                 legend = TRUE, legend.cex = legend.cex, GLratio = GLratio,
                 groups = grs_l,
                 edge.color=col_m, esize=edge_size,
                 color=node_cols, vsize=v_size, label.cex = label.cex,
                 filetype = file_type, filename = file_path, width = file_width,
                 height = file_height, node.resolution = 300)

  qp=qgraph::qgraph(PSN, threshold = edge_threshold, edgelist = FALSE, weighted = TRUE, directed = FALSE,
                    layout = "spring", repulsion = 0.7,
                    legend = TRUE, legend.cex = legend.cex, GLratio = GLratio,
                    groups = grs_l,
                    edge.color=col_m, esize=edge_size,
                    color=node_cols, vsize=v_size, label.cex = label.cex)

  res_l=list(qp=qp,names_in_k_l=names_in_k_l)
  return(res_l)
}
