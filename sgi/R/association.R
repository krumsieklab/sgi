
#' Subgroup Identification Analysis
#'
#' This function runs the subgroup identification analysis initialized by \code{sgi::sgi_init()} and returns statistics of subgroup enrichment for all valid cluster pairs.
#' Statistical tests are decided based on the type of outcomes. T-test, Fisher's exact tests/Chi-square tests, and log-rank test are employed for continuous, dichotomous/factor, and survival outcomes, respectively.
#'
#' @param sg an \code{sgi.object} created by \code{sgi::sgi_init()}
#'
#' @return association results as a list which contains
#' \item{results}{statistical results correponding to each outcome}
#' \item{cluster_pairs}{summary of cluster pairs, e.g. cluster ids, sample sizes, tree levels etc.}
#' \item{hc}{ \code{hclust} object given to \code{sgi} }
#' \item{minsize}{minimum cluster size threshold}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' as = sgi::sgi_run(sg)
#' as
#' }
#'
#' @author mubu
#'
sgi_run <- function(sg){

  # hiertest = T # WY is deprecated

  if(class(sg)!="sgi.object") stop("sg has to be an sgi.object")

  # tree parameters
  hc = sg$hc
  minsize = sg$minsize
  # outcomes
  outcomes = sg$outcomes
  # outc_classes = sapply(lapply(outcomes,class),`[`,1)

  # valid cluster pairs
  clinds = sg$clinds
  # tests defined by user
  tests_user = sg$user_defined_tests

  if(length(tests_user)>0 & any(!sapply(tests_user, is.function)))
    stop("user_defined_tests should be function!")


  # fix here after package is created
  # default test options
  tests_default =  list( numeric = ff_numeric, integer = ff_numeric, factor = ff_factor, Surv = ff_Surv )

  # final list of tests to be applied
  # check if a new test supercedes the older one
  f_tests = c(tests_default[setdiff(names(tests_default), names(tests_user))], tests_user)

  # decide which classes should be used
  outc_classes = sapply(outcomes, function(x){
    xc = class(x)
    i = which(xc %in% names(f_tests))[1]
    if(is.na(i)) return(NA)
    xc[i]
  })

  # hierarchical test
  cmpset <- na.omit(clinds)
  rownames(cmpset) = paste(cmpset[,"cid1"], cmpset[,"cid2"], sep = "vs")

  # with new update taking functions to higher level no need for this trick
  # # keep classes as it is
  # classes = lapply(outcomes, class)
  outcomes = outcomes[hc$order, ,drop = F]

  re <- list( minsize = minsize,
              hc = hc,
              cluster_pairs = cmpset,
              results = lapply( structure(seq(outcomes), names = names(outcomes)) , function(i)
                ftc( .ftest = f_tests[[outc_classes[i]]], x = outcomes[[i]], cmpset = cmpset)
              ))

  class(re)<-"sgi.assoc"

  return(re)
}



#' Summary of subgroup identification analysis
#'
#' Summary function that returns a result summary for a \code{sgi.assoc} object created by \code{sgi::sgi_run()}
#'
#' @param as \code{sgi.assoc} object created by \code{sgi::sgi_run()}
#' @param padj_th Adjusted p-value threshold for filtering results
#' @param by_clusters \code{TRUE}: Summarize results at by clusters; \code{FALSE}: Summarize results by outcomes
#'
#' @return Summary of SGI results
#'
#' @method summary sgi.assoc
#' @export
#'
#' @examples
#' \dontrun{
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' as = sgi::sgi_run(sg)
#' summary(as, padj_th = 0.01)
#' }
summary.sgi.assoc <- function( as, padj_th = 0.05, by_clusters = NULL){

  if(is.null(by_clusters))
    by_clusters = ifelse(length(as$results)<2,F,T)

  if(by_clusters){
    cp = rownames(as$cluster_pairs)
    names(cp) = cp
    ll = lapply(cp, function(x){
      # browser()
      # r = t(sapply(as$results, function(y) y[x,]))
      # r = r[, !(colnames(r) %in% "level"), drop = F]
      # # filter significant clusters
      # r = r[r[,"padj"] <= padj_th,, drop = F]

      r = do.call(rbind, lapply(as$results, function(y) subset(y,rownames(y) %in% x)) )
      r = cbind(outcome = names(as$results), r) #[,c("cid1","cid2","level"):=NULL]
      # filter significant clusters
      r = r[padj<= padj_th, ]

      # add cluster info as attr
      attr(r, "vcp") = as.matrix(as$cluster_pairs[x,])
      class(r) = c("assoc_table",class(r))
      r
    })
  }else{
    ll = lapply(as$results, function(r){
      #browser()
      r = cbind(cluster_pair = rownames(as$cluster_pairs), r)
      r = r[padj <= padj_th, ]
      if(nrow(r)>0){
        vcp = as$cluster_pairs[r$cluster_pair,]
        attr(r, "hvcp") = as.matrix(vcp)
      }
      class(r) = c("assoc_table",class(r))
      r
    })
    # ll = Filter(function(x) nrow(x)>0, ll)
  }
  class(ll) = c("listofx", "listof","summary.sgi.assoc",class(ll))
  ll
}

# modified listof to print assoc results properly
#' @export
print.listofx<- function( x, drop_list = T, ...){
  cat("\nSummary of SGI associations...\n\n")

  if(drop_list) x = x[sapply(x, nrow)>0]
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll)
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)) {
    cat(nn[i], ": ")
    print(x[[i]], ...)
    cat("\n")
  }
  invisible(x)
}

# result of association for each outcome, summarized table
#' @export
print.assoc_table <- function( x, cluster_info = T){
  vcp = unlist(attr(x,"vcp")[1,])

  if(is.null(vcp)) cluster_info = F

  if(nrow(x)>0){
    if(cluster_info){
      # hiertest
      cat(vcp["cid1"],"(n=", vcp["n_c1"] ,") vs ", vcp["cid2"],"(n=", vcp["n_c2"] ,")",sep = "")

      cat(" at L=", vcp["level"], ", h=", round(vcp["h"],2), "\n", sep = "")

    }else{cat("\n")}
  }


  if(nrow(x)>0) print(data.table::as.data.table(x[,.SD, .SDcols = !c("cid1","cid2","test")]),
                      row.names = F, digits = 3)
}

#' Print subgroup identification analysis results
#'
#' Print results contained in an \code{sgi.assoc} object created by \code{sgi::sgi_run()}
#'
#' @param as \code{sgi.assoc} object created by \code{sgi::sgi_run()}
#' @param padj_th Adjusted p-value threshold for filtering results
#' @param by_clusters \code{TRUE}: Summarize results at by clusters; \code{FALSE}: Summarize results by outcomes
#'
#' @method print sgi.assoc
#' @export
#'
#' @examples
#' \dontrun{
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' as = sgi::sgi_run(sg)
#' print(as, padj_th = 0.01)
#' }
print.sgi.assoc <- function( as, padj_th = 0.05, by_clusters = T, ...){
  
  # ellipsis(...) is to support for silent option which can be helpful to return 
  # printed text as a list
  silent = list(...)$silent
  if(is.null(silent)) silent = FALSE

  if(!silent){
    cat("\nSGI associations...")
    cat("\n--- Significant results ---\n\n")
  }

  if(is.null(by_clusters))
    by_clusters = ifelse(length(as$results)<2,F,T)

  # filtered results with padj_th
  stable = sapply(as$results, function(x) x[,"padj"] < padj_th)

  # l = apply( stable,2-by_clusters, function(x)
  #   if(length(names(which(x)))>0) return(names(which(x))))

  #-----
  # if only one valid cluster
  if(is.vector(stable)) stable = t(stable)

  rownames(stable) = rownames(as$cluster_pairs)
  # browser()

  l = apply( stable,2-by_clusters, function(x)
    if(length(names(which(x)))>0) return(list(names(which(x)))))
  l = lapply(l, unlist)
  #-----

  l = Filter( Negate(is.null), l)
  if(is.null(l)) return(l)

  if(silent) return(invisible(l))

  nn <- names(l)
  ll <- length(l)
  for (i in seq_len(ll)) {
    cat(nn[i], ": ")
    cat(paste(l[[i]],collapse = ", "),"\n")
    cat("\n")
  }
  invisible(l)
}

