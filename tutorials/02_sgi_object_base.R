#' ---
#' title: "SGI tutorial 2: Objects and functions"
#' author:
#' date:
#' output: github_document
#' ---

#+include = F
rm(list = ls())

knitr::opts_chunk$set(
  out.width = 95*7,
  out.height = 65*7,
  # fig.width = 9.5,
  # fig.height = 6.5,
  warning = F, 
  message = F, 
  include = T
)

#'

#' This tutorial explains the functionality of SGI objects for result extraction and plotting.

library(sgi)
library(magrittr)
library(ggplot2)

# QMdiab clinical variables and plasma metabolomics data 
plasma = sgi::qmdiab_plasma
clins = sgi::qmdiab_clin

#' ## Generate hierarchical clustering tree 
hc = hclust(  dist(plasma), method = "ward.D2")
plot(hc, labels = F)

#' ## Generate and inspect SGI object
# initialize SGI structure, minsize is set to 0.05 of total sample size
# no analysis is performed yet
sg = sgi_init(hc, minsize = round(max(hc$order)/20), outcomes = clins)
class(sg)

# print the sgi object
sg
# summarize the valid cluster pairs
summary(sg)
# plot SGI object to visualize valid cluster pairs to be tested
plot(sg) 


#' ## Run association analysis
as = sgi_run(sg)
class(as)

# print significant associations, default adjusted p-value threshold is 0.05
as
# print results with different p-value threshold
print(as, padj_th = 0.001)
# print with respect to outcomes 
print(as, padj_th = 0.1, by_clusters = F)

# summary of statistical results 
summary(as, padj_th = 0.01)
# print summary with respect to outcomes
(sm = summary(as, padj_th = 0.01, by_clusters = F))
# list-based objects for further processing
sm$Diabetes
# access adjusted p-values
sm$Diabetes$padj

#' ## Visualization

# SGI produces ggplot objects
gg_tree = plot(as, padj_th = 0.05)
gg_tree

# the user can access all functionality of the ggplot package
gg_tree + coord_flip()

gg_tree + 
  scale_y_reverse() +  
  scale_color_viridis_d(na.value = "black") + 
  scale_fill_viridis_d()

gg_tree + coord_polar()

gg_tree + 
  scale_y_reverse() + 
  coord_polar() 


