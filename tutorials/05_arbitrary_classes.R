#' ---
#' title: "SGI tutorial 5: SGI with user-defined statistical tests"
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

#'
#' This tutorial explains how user-defined tests can be use in SGI.
#' The SGI algorithm accesses statistical test functions by the classes of the input variables.
#' The user can control SGI's behavior by overwriting the test function for a given variable class, and
#' by setting custom variable classes.
#'
#' Setup SGI: 

# Initialize
library(sgi)
library(magrittr)
library(ggplot2)

# Numeric output formatting
options(digits = 3)

# load QMDiab plasma metabolomics
plasma = sgi::qmdiab_plasma
# hierarchical based on plasma data
hc = hclust(dist(plasma), method = "ward.D2")

# load QMDiab clinical variables
clins = sgi::qmdiab_clin

# verify sample order
stopifnot(identical(rownames(plasma), rownames(clins)))


#' ## Change default test for 'numeric' from t.test to wilcox.test

# SGI's default test for numeric is t.test.
# --> we here change it to wilcox.test.

# first, inspect the template function body
sgi::create_stattest_function_skeleton()

# implement wilcox.test based on this template
my_wilcox_test <- function(v, y){

  # try to run wilcox.test  
  res <- tryCatch({
    wfit <- wilcox.test(v~y)
    # return results if test
    list(pval = wfit$p.value, stat = wfit$statistic)

  }, error = function(e){
    # something didn't work with the wilcox.test (probably no samples in one of the groups); return NA vector
    list(pval = NA, stat = NA)
  })

  # return results
  return(res)

}

# inspect classes of each outcome
sapply(clins, class)

# set user-defined test as the test to be used for numeric 
sg = sgi_init(hc, minsize = 18, outcomes = clins, 
              user_defined_tests = list(numeric = my_wilcox_test))

# run SGI
as = sgi_run(sg)

# check if wilcox.test has been applied properly and p-values are correct
# to test this, manually extract each cluster pair and perform the tests

# get valid cluster pairs, each column is a pair of valid clusters
mm = get_vcps(sg)
# run wilcox.test manually and record p-values for, here for HbA1c as an example
wpvals = apply(mm, 2, function(vcp) wilcox.test(clins$HbA1c~vcp)$p.value)

# confirm that all p-values are equal
cbind(wpvals, sgi = as$results$HbA1c$pval)


#' ## User-defined test for an arbitrary class

#' Change default test for some of the numeric variables to wilcox.test.
#' We need to assign a new class to those variables.

# show all variable classes
sapply(clins, class)

# define a new class called "ordinal" for sodium and hemobglobin
# (could be any arbitrary string, R allows any strings as class names)
class(clins$Sodium) = c("ordinal", class(clins$Sodium))
class(clins$Hemoglobin) = c("ordinal", class(clins$Hemoglobin))

# set user-defined test as the test to be performed for 'ordinal'
sg = sgi_init(hc, minsize = 18, outcomes = clins, 
              user_defined_tests = list(ordinal = my_wilcox_test))

# run SGI
as = sgi_run(sg)

# verify that SGI produced the correct p-values

# wilcox.test for hemoglobin ('ordinal')
cbind( wilcox = apply(mm, 2, function(vcp)
  wilcox.test(clins$Hemoglobin~vcp)$p.value),
  sgi = as$results$Hemoglobin$pval)

# t.test for age ('numeric')
cbind(
  wilcox = apply(mm, 2, function(vcp)
    wilcox.test(clins$AGE~vcp)$p.value),
  t.test = apply(mm, 2, function(vcp)
    t.test(clins$AGE~vcp)$p.value),
  sgi = as$results$AGE$pval)



#' ## User-defined test for an outcome that consists of more than one variable

#' This example showcases a scenario where a column in a data frame contains another data frame with two columns.
#' This can be used to define multivariable outcomes.
#' The application here is to correct a diabetes association for sex effects inside the association function.

# create data frame (that will become a single column later)
my_df = data.frame(diab = clins$Diabetes, sex = clins$SEX)

# define class
class(my_df) = c("diabetes_and_sex", class(my_df))
class(my_df)

# add this newly created variable to outcomes data frame
clins$diab_corrected_by_sex <- my_df
class(clins$diab_corrected_by_sex)

# Implement user-defined test function.
# Note that in this case, the argument v will be a data.frame.
# In order to asses the significance of the clusers, we compare the following two models:
# diab~sex and diab~sex+clusters. Significance is determined by an ANOVA likelihood ratio test
my_diab_test <- function(v, y){
  # ----- input -----------------------
  # v: outcome
  # y: valid cluster pairs as a vector of factors
  # -----------------------------------

  v$cluster_pair = y
  # null model with only sex
  mnull = glm(diab~as.numeric(sex), data = v, family = binomial())
  # model testing the cluster pair
  m1 = glm(diab~cluster_pair+as.numeric(sex), data = v, family = binomial())
  anv = anova(mnull, m1, test = "LRT")
  pval = anv$`Pr(>Chi)`[2]
  stat = anv

  # return
  return(list(pval = pval, stat = stat))
}

# set user-defined test as the test to be performed for 'diabetes_and_sex'
sg = sgi_init(hc, minsize = 18, outcomes = clins, 
              user_defined_tests = list(numeric = my_wilcox_test, 
                                        diabetes_and_sex = my_diab_test ))
# run SGI
as = sgi_run(sg)
# inspect results
as$results$diab_corrected_by_sex


#' ## Check clustering aggrement w.r.t a reference data

#' This is an advanced example of using user-defined functions, which checks the clustering agreement between
#' plasma-based clusters in comparison with urine and saliva clusters. In this case, each "outcome" is an entire
#' data matrix. The testing function will receive the respective matching subsets of urine and saliva data matrices
#' from the SGI engine. We cluster these with hierarchical clustering inside the test function and derive two
#' clusters. The function then tests whether the two obtained clusters in urine/saliva overlap with the given clusters
#' in plasma using a Fisher's exact test. In other words, this approach checks for each plasma-based branching point whether urine and saliva
#' cluster similarly for the same samples. The lower the resulting p-values, the higher the agreement.

# identify samples measure on all three platforms
ids = intersect( intersect(rownames(sgi::qmdiab_plasma), 
                           rownames(sgi::qmdiab_saliva)),
                 rownames(sgi::qmdiab_urine))

# subset to those common samples
plasma = qmdiab_plasma[ids,]
urine = qmdiab_urine[ids,]
saliva = qmdiab_saliva[ids,]

# reference hc is plasma
hc = hclust(dist(plasma), method = "ward.D2")

# generate outcome data.frame that contains the entire data matrices
df_platforms = data.frame( row.names = rownames(urine))
df_platforms$urine = urine
df_platforms$saliva = saliva

# verify classes
sapply(df_platforms, class)


# implement user-defined test function,
# in this case v is a matrix
my_cool_test <- function(v, y){
  # ----- inputs -----------------------
  # v: outcome
  # y: valid cluster pairs as a vector of factors
  # -----------------------------------

  # y is plasma clusters
  y_plasma = y
  # generate 2 clusters based on v (which is the matching subset of urine/saliva samples)
  yh = cutree( hclust(dist(v), method = "ward.D2"), k = 2)
  # compare with plasma clusters
  tb = table(y_plasma,yh)
  pval = fisher.test(tb)$p.value
  stat = tb

  # return results
  return(list(pval = pval, stat = stat))
}

# set user-defined test as the test to be performed for 'matrix'
sg = sgi_init(hc, minsize = 18, outcomes = df_platforms, 
              user_defined_tests = list(matrix = my_cool_test))
# run SGI
as = sgi_run(sg)

# inspect results, urine shows some agreement with plasma clusters, saliva doesn't agree at all
summary(as, padj_th = 1, by_clusters = F)

# plot plasma tree annotated with other platforms based on clustering aggrement
plot(as, plot_only_valid_cluster_pairs = T,
     tree_pruning = 2,
     overlay_only_sgs = F,
     overlay.color = "black",
     overlay.alpha = 0.2,
     overlay.fill = "wheat")

# inspect the confusion tables for cluster pair 22vs23

# urine
as$results$urine$stat[[8]]

# saliva
as$results$saliva$stat[[8]]



















