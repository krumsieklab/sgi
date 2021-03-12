#' Template for user-defined statistical test
#'
#' Use this template to define customized statistical tests for a specific variable class.
#'
#' @return returns function template
#' @export
#'
#' @examples
#' \dontrun{
#' create_stattest_function_skeleton()
#' }
create_stattest_function_skeleton <- function(){

  function(v, y){
    # ----- input -----------------------
    # v: outcome
    # y: valid cluster pairs as a vector of factors
    # -----------------------------------

    # your code comes here,

    # store p-value in 'pval' and
    #        test statistic in 'stat'

    # return results
    return(list(pval = pval, stat = stat))
  }

}


#' Combine multiple distance matrices into one for integrated hierarchical clustering
#'
#' @param dlist List of distance matrices of class \code{dist}
#' @param w Weights for distance matrices in \code{dlist}, if weighted combination is desired; all matrices have equal weight by default
#'
#' @return Combined distance matrix
#' @export
#' @author mubu
#'
#' @examples
#' \dontrun{
#' d1 = dist(rnorm(10))
#' d2 = dist(rnorm(10))
#' d3 = dist(rnorm(10))
#'
#' # list of distance matrices
#' dlist = list(d1,d2,d3)
#' # function call
#' ndist(dlist)
#' }
#'
ndist <- function(dlist, w = rep(1, length(dlist))){

  if(!is.list(dlist)) stop("dlist has to be a list of distance matrices")
  if(any(sapply(dlist, function(x) class(x)[1])!="dist"))
    stop("every object in dlist has to be of class 'dist'!")
  if(any(w<0)) stop("weights must be positive")

  if(length(dlist)<2) return(dlist[[1]])

  # map all to same scale
  dlist = lapply(dlist, function(d) d/max(d))

  re = dlist[[1]]*w[1]
  for(i in 2:length(dlist))
    re = re + dlist[[i]]*w[2]

  re/sum(w)
}


#' Get valid cluster pairs
#'
#' @param sg \code{sgi.object} created by \code{sgi::sgi()}
#'
#' @return Returns a matrix where each column is a valid cluster pair and each row corresponds to a sample
#' @export
#'
#' @examples
#' \dontrun{
#' # for a given hclust object 'hc' and outcome data frame 'outcomes'
#' sg = sgi_init(hc, minsize = 18, outcomes)
#' get_vcps(sg)
#' }
get_vcps <- function(sg){

  # valid cluster pairs
  vcps = na.omit(sg$clinds)
  structure(
    apply(vcps, 1, sgi:::ginds_vcp, n = length(sg$hc$order))[order(sg$hc$order),],
    dimnames = list( rownames(sg$outcomes), paste0("l", vcps$level))
  )
}


# get valid cluster pair indices as factors
ginds_vcp  <- function(x, n=NULL){
  a = if(is.data.frame(x)) unlist(x[,c("c1_1", "c1_2", "c2_1", "c2_2", "cid1", "cid2")]) else x[c("c1_1", "c1_2", "c2_1", "c2_2", "cid1", "cid2")]
  if(is.null(n)) return(
    list( y=c(rep( unname(a[5]), a[2]-a[1]+1), rep(unname(a[6]), a[4]-a[3]+1)), loc = a[1]:a[4] )
  )
  y = rep(NA, n)
  y[a[1]:a[4]] = c(rep(a[5], a[2]-a[1]+1), rep(a[6], a[4]-a[3]+1))
  y
}











