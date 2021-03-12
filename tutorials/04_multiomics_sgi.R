#' ---
#' title: "SGI tutorial 4: Multi-omics SGI"
#' author: 
#' date: 
#' output: github_document
#' ---

#+include = F
rm(list = ls())

knitr::opts_chunk$set(
  # out.width = 103*7,
  # out.height = 65*7,
  fig.width = 7.21,
  fig.height = 4.25,
  warning = F, 
  message = F, 
  include = T
)

#' This tutorial showcases the multi-omics capabilities of SGI. We will use plasma, urine, and saliva metabolomics data as three different levels of omics.

library(sgi)
library(magrittr)
library(ggplot2)

#' 
#' ## Multi-omics SGI

# find samples measured on all 3 platforms
ids = intersect( intersect(rownames(sgi::qmdiab_plasma), 
                           rownames(sgi::qmdiab_saliva)),
                 rownames(sgi::qmdiab_urine))

# extract data frames with overlapping samples
plasma = qmdiab_plasma[ids,]
urine = qmdiab_urine[ids,]
saliva = qmdiab_saliva[ids,]

# clinical variables
clins = qmdiab_clin[ids, ]

# generate list of distance matrices using dist() function
ldist = lapply( list(plasma = plasma, urine = urine, saliva =  saliva), dist)

# integrate distance matrices into one using SGI's ndist function
d = ndist(ldist) 

# create hierarchical clustering with multi-omics distance matrix 
hc = hclust(d, method = "ward.D2")
plot(hc, labels = F)


# run SGI on this multi-omics-derived tree
# initialize SGI structure; minsize is set to 5% of sample size
sg = sgi_init(hc, minsize = ceiling( length(hc$height)/20 ), outcomes = clins)
# run SGI
as = sgi_run(sg)

(gg_multiomics = plot(as, padj_th =0.01) + 
    theme(legend.position = "none")+ggtitle(label = "multi-omics"))


#' ## Compare with regular, single-omics SGI
#' 
# also run SGI on each individual omics layer
nods <- lapply( structure(names(ldist), names = names(ldist)), function(i){
  hc  = hclust(ldist[[i]], method = "ward.D2")
  # initialize SGI structure; minsize is set to 5% of sample size
  sg = sgi_init(hc, minsize = ceiling( length(hc$height)/20 ), outcomes = clins)
  # run SGI
  sgi_run(sg)
})

# plot each result
ggs = lapply(names(nods), function(i) 
  plot(nods[[i]],  padj_th = 0.01,  
       draw_cluster_ids = F, 
       plot_only_valid_cluster_pairs = T, 
       overlay = F, 
       tree_pruning = 2) + 
    theme(legend.position = "none") + 
    ggtitle(label = i))

#+out.width = NULL, out.height = NULL, fig.width = 15, fig.height = 10
# plot trees together
library(patchwork)
gg_multiomics/(ggs[[1]] + ggs[[2]] + ggs[[3]]) 

# results are improved by integration, as shown by diabetes enrichment at first split
options(digits = 3)
c(list(multi_omics = as), nods) %>% 
  sapply(function(as) as$results$Diabetes[1,1:4] ) %>% t

# HbA1c enrichments at first split is also improved
c(list(multi_omics = as), nods) %>% 
  sapply(function(as) as$results$HbA1c[1,1:4] ) %>% t

# same for sodium 
c(list(multi_omics = as), nods) %>% 
  sapply(function(as) as$results$Sodium[1,1:4] ) %>% t




















