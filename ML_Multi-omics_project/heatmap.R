library(pheatmap)

rownames(GE_new) <- paste0("p_", seq(nrow(GE_new)))

# create heatmap using pheatmap
pheatmap(GE_new)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

GE_norm <- t(apply(GE_new, 1, cal_z_score))
pheatmap(GE_norm, cluster_cols = FALSE)

# load package
library(dendextend)
my_heatmap <- pheatmap(GE_norm, silent = TRUE, cluster_cols = FALSE)

my_heatmap$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)

GE_col <- cutree(tree = as.dendrogram(my_heatmap$tree_row), k = 3)
sw_fun <- function(x) switch(x, '1'='cluster_1', '2'='cluster_2', '3'='cluster_3')
GE_col <- data.frame(cluster = sapply(GE_col, sw_fun))

pheatmap(GE_norm, annotation_row = GE_col,
         cutree_rows = 3, cluster_cols = FALSE)
## from icluster
clusters  = output2[[2]]$fit[[306]]$clusters

my_clusters <- as.factor(clusters)
GE_col$iclusters <- my_clusters 

pheatmap(GE_norm, annotation_row = GE_col,
         cutree_rows = 3, cluster_cols = FALSE)

### analysis why there is a huge gap


### CNA analysis
rownames(CNA_new) <- paste0("p_", seq(nrow(CNA_new)))
pheatmap(CNA_new, cluster_cols = FALSE, cluster_rows = FALSE)

CNA_heatmap <- pheatmap(CNA_new, clustering_distance_rows="manhattan", silent = TRUE, cluster_cols = FALSE)
CNA_col <- cutree(tree = as.dendrogram(CNA_heatmap$tree_row), k = 3)
CNA_col <- data.frame(cluster = sapply(CNA_col, sw_fun))

CNA_col$iclusters <- my_clusters 

pheatmap(CNA_new, annotation_row = CNA_col,
         cutree_rows = 3, cluster_cols=FALSE)


### Somatic mutation
rownames(SM_new) <- paste0("p_", seq(nrow(SM_new)))
#pheatmap(SM_new, cluster_cols = FALSE, cluster_rows = FALSE)

SM_heatmap <- pheatmap(SM_new, clustering_distance_rows="manhattan", silent = TRUE, cluster_cols = FALSE)
SM_col <- cutree(tree = as.dendrogram(SM_heatmap$tree_row), k = 2)
SM_col <- data.frame(cluster = sapply(SM_col, sw_fun))

SM_col$iclusters <- my_clusters 

pheatmap(SM_new, annotation_row = SM_col,
         cutree_rows = 2, cluster_cols=FALSE)
















