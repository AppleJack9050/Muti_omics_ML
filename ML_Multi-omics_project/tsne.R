library(dplyr)
library(cluster)
library(ggplot2)
library(Rtsne)

GE_dist <- daisy(GE_norm)

tsne_obj <- Rtsne(GE_dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(GE_col$cluster),
         name = rownames(GE_col))

ggplot(aes(x = X, y = Y), data = tsne_data) + geom_point(aes(color = cluster))


