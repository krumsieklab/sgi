#' Initialize subgroup identification analysis
#'
#' This is a constructor function which initializes the subgroup identification process
#'
#' @param hc Hierarchical clustering object of class \code{hclust}
#' @param minsize Minimum size for a cluster to be considered
#' @param outcomes \code{data.frame} of outcomes to be tested; if this data frame has different rownames than \code{hc$labels}, a warning will be thrown
#' @param user_defined_tests [optional] Named list of functions for user-defined tests. The keys of the list are the classes for which the functions are defined, the values are the functions. See SGI tutorial for more details.
#' @return Object of class \code{sgi.object} which contains all required information to run subgroup identification.
#' @export
#' @examples
#' \dontrun{
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' sg
#' }
#'
#' @author mubu
#'
sgi_init <- function(hc, minsize, outcomes, user_defined_tests=c()){
  
  # default is 5% of sample size
  if(missing(minsize)) minsize = length(hc$height)*0.05

  if(class(hc)!="hclust"){stop("hc has to be 'hclust' object!")}
  if(minsize>length(hc$height) | minsize<2){stop("invalid minsize!")}
  #if(levels>(length(hc$height)-1)|| levels<2){stop("invalid levels.")}

  if(is.data.frame(outcomes)){
    defined_test_types = unique( c("integer","numeric","factor", "Surv",
                                   names(user_defined_tests)) )

    # decide which classes should be used
    outc_classes = sapply(outcomes, function(x){
      xc = class(x)
      i = which(xc %in% defined_test_types)[1]
      if(is.na(i)) return(NA)
      xc[i]
    })

    if(any(is.na(outc_classes)))
      stop("Test for class(es)  of ", paste(names(outcomes)[is.na(outc_classes)], collapse = ", "), " not defined!")
  }else{
    stop("outcomes has to be data.frame!")
  }

  if(!identical(hc$labels, rownames(outcomes)))
    warning("be careful, 'hc$labels' and 'rownames(outcomes)' not identical.")

  ruc <- get_uc_vc(hc, minsize)
  clinds = ruc$clinds

  if(nrow( na.omit(clinds) ) < 1 ) stop("no valid cluster pairs for given minsize!")

  # if type of the test is not defined via attr, set it "user-defined"
  if(length(user_defined_tests)>0)
    user_defined_tests = lapply(user_defined_tests, function(f){
      test_type = attr(f, "type")
      if(is.null(test_type)) return(structure(f, type = "user-defined"))
      else f
    })

  # user defined stat tests for arbitrary outcomes
  new_stat_tests = intersect(names(user_defined_tests), outc_classes)

  # for some classes if no available test, check and warn or throw an error

  # if some tests are defined but never used
  if(length(user_defined_tests) != length(new_stat_tests))
    warning(paste0(setdiff(names(user_defined_tests), new_stat_tests),collapse = ", "),
            " tests are defined but no outcomes of these types applied for!")

  re=list(hc=hc, minsize=minsize, outcomes=outcomes, clinds = clinds,
          user_defined_tests = user_defined_tests[new_stat_tests])
  class(re)<-"sgi.object"
  return(re)
}

#' Summary of \code{sgi.object}
#'
#' Summary function that returns information of valid clusters of \code{sgi.object} created by \code{sgi::sgi_init()}
#'
#' @param sg \code{sgi.object} created by \code{sgi::sgi_init()}
#' @param summary Whether unpaired valid clusters should be summarized or not. \code{TRUE} (default)
#' @param silent \code{TRUE} returns the object but does not print, default is \code{FALSE}
#'
#' @return summary of valid clusters
#'
#' @method summary sgi.object
#'
#' @export
#' @examples
#' \dontrun{
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' summary( sg )
#' }
#'
#' @author mubu
#'
summary.sgi.object<- function(sg, summary = T, silent = F){

  # work around for older version with summary =1 which returns
  # unpaired valid clusters
  if(summary) summary = 2

  cat("\nSummary of SGI clusters...\n")

  # remove cluster locations
  sm = sg$clinds[,-(1:4)]

  # summary comparable pairs
  smcp =na.omit(sm)
  # summary of valid clusters only
  smvc = sm[attr(smcp,"na.action"),]

  # cat("minsize: ", sg$minsize, "\n")
  if(!silent){
    cat(" - # of comparable valid cluster pairs: ", nrow(smcp), "\n")
    print(smcp,row.names = F)
    if(summary < 2){
      cat("\n - - # of unpaired valid clusters: ", nrow(smvc), "\n")
      print(smvc,row.names = F)
    }
  }

  invisible(list(valid_cluster_pairs = smcp,  valid_clusters_unpaired = if(summary < 2) smvc else NULL))
}

#' Print an \code{sgi.object}
#'
#' Print function for \code{sgi.object} created by \code{sgi::sgi_init()}
#'
#' @param sg \code{sgi.object} created by \code{sgi::sgi_init()}
#' @param summary Print summary of valid clusters? \code{FALSE} by default
#' @param outcomes Print outcomes? \code{TRUE} by default
#'
#' @return NULL
#' @method print sgi.object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' sg
#' print(sg, summary = T, outcomes = F)
#' }
#'
#' @author mubu
#'
print.sgi.object<- function(sg, summary = F, outcomes = T){

  cat("\nSGI initialization...\n")
  cat("\n-------- hierarchical clustering --------------")
  print(sg$hc)
  cat("-------- min valid cluster size ---------------\n")
  cat("minsize: ", sg$minsize, "\n")

  if(summary){
    cat("\n-------- summary ------------------------------\n")
    summary.sgi.object(sg, summary)
  }

  if(outcomes){
    cat("\n-------- outcomes -----------------------------\n")
    str( sg$outcomes, list.len = 10, vec.len = 2 )
    cat("\n")
  }


  # user defined tests
  if(length(sg$user_defined_tests)>0){
    cat("\n----- user defined tests for outcome types ----\n")
    cat( names(sg$user_defined_tests), sep = ", " )
    cat("\n")
  }
}
