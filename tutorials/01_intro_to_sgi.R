#' ---
#' title: "SGI Tutorial 1: Quickstart guide"
#' author:
#' date:
#' output: github_document
#' ---

#+include = F
rm(list = ls())

knitr::opts_chunk$set(
  out.width = 95*7,
  out.height = 65*7,
  fig.width = 11,
  fig.height = 7,
  warning = F, 
  message = F, 
  include = T
)

#'

#' This tutorial provides a quick introduction to the functionality of the SGI package.
library(sgi)
library(magrittr)
library(ggplot2)

#' 
#' ## The QMdiab data

# QMdiab clinical variables and plasma metabolomics data 
# see '?sgi::qmdiab_plasma' for details
plasma = sgi::qmdiab_plasma
# see '?sgi::qmdiab_clin' for details
clins = sgi::qmdiab_clin
  
# data type and size
plasma %>% class(); plasma %>% dim()
# first 4 rows and columns of the data matrix
plasma[1:4,1:4]

# clinical parameters/outcomes
clins %>% str

#' ## Generate hierarchical clustering tree 
hc = hclust(  dist(plasma), method = "ward.D2")
plot(hc, labels = F)

#' ## Run SGI & Visualize results
# initialize SGI structure; minsize is set to 5% of sample size
sg = sgi_init(hc, minsize = 18, outcomes = clins)
# run SGI
as = sgi_run(sg)

#+fig.width = 9.5, fig.height = 6.5
# plot, show results for adjusted p-values <0.05
(gg_tree = plot(as, padj_th = 0.05))

# show overview, including clinical data and metabolomics data matrix
plot_overview( gg_tree = gg_tree, as = as, 
               outcomes = clins, xdata = plasma, 
               ncut_tree = 4, data_title = "plasma")

# generate distribution plots of outcomes
ggs = plot_outcomes(sg, as, padj_th = 0.05)
# show that ggs is a list of plots 
str(ggs, max.level = 1)

library(patchwork) # needed to combine plots, overwrites + operator between ggplots

# look at 4vs5 subgroups
ggs$Diabetes$`4vs5`+ ggs$HbA1c$`4vs5`+ ggs$BMI$`4vs5`



