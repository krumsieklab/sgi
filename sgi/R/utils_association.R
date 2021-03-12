# TO-DOs:
# 1. add covariates z
# 2. returning rich set of stats, i.e. dev and df
# - - list may be
# 3. collect and print/summary rich stats principly
#

# function to get unique, valid clusters etc
get_uc_vc <-  function(hc, th){
  # # of samples = # of possible clustering including 1 cluster
  k = length(hc$order)

  # cutting the tree every possible level possible given th
  omat <- cutree(hc, k = seq(length(hc$order)-th-1))[hc$order,]
  omat = omat[,!apply(omat,2, function(x) all(table(x) < th) ),drop = F]

  # giving a hashed id to every split at every level as hashed.id = i1_i2
  chah <-
    apply(omat,2, function(x)
      sapply(unique(x), function(i){
        inds = which(x==i)
        paste(min(inds), max(inds), sep = "_")
      })
    )

  if(length(chah)==1){
    if(length(unlist(chah))==1)
      stop("no valid cluster pairs!\ntry decreasing threshold.")
  }

  # get unique splits on the tree with level
  csm = sapply(2:length(chah), function(i) c(setdiff(chah[[i]],chah[[i-1]]),i) )
  rm(chah)
  # set cluster ids
  mmat = t(matrix(as.numeric( factor(c(csm[1:2,]), levels = c(csm[1:2,]))), ncol = ncol(csm)))+1
  # locations of clusters
  sz = sapply(apply(csm[1:2,],2,strsplit,split="_"), function(x) as.numeric(unlist(x)))
  l = as.numeric(csm[3, ])
  mmat = cbind(t(sz), l = l, h = rev(hc$height)[l-1], mmat)
  # sizes of clusters
  sz = apply(sz,2,diff)[-2, ,drop = F]+1
  rm(csm,l)
  # to keep, at least 1 valid clusters at level h
  mmat = mmat[!(colSums(sz<th)==2), , drop = F]
  sz = sz[,!(colSums(sz<th)==2), drop = F]
  colnames(mmat) = c("c1_1", "c1_2", "c2_1", "c2_2", "level", "h", "cid1", "cid2")
  # unvalid cluster of the pair
  mmat[sz[1,]<th,"cid1"] = NA
  mmat[sz[2,]<th,"cid2"] = NA
  rownames(sz) = c("n_c1","n_c2")
  # set sized for non-hiertest
  sz[1,is.na(mmat[,"cid1"])] = k - sz[2,is.na(mmat[,"cid1"])]
  sz[2,is.na(mmat[,"cid2"])] = k - sz[1,is.na(mmat[,"cid2"])]

  list( clinds = data.frame(cbind(mmat,t(sz))), minsize = th )
}

# multiple testing correction
# z: covariates
# .ftest : function that does the statistical test
ftc<-function(.ftest , x, cmpset){

  test0 = tst(.ftest, x, cmpset)
  # browser()
  # pvals = test0[,"p"]
  # ttype = attr(test0, "type")
  pvals = sapply( test0$test_results, `[[`, "pval" )
  ttype = test0$type
  lvels = test0$level

  if(is.null(ttype)) ttype = NA

  padj = p.adjust(pvals, method = "bonferroni")

  re <- data.table::data.table(
    cid1 = cmpset[,"cid1"], cid2 = cmpset[,"cid2"],
    padj = padj, pval = pvals, level = lvels,
    stat = lapply(test0$test_results, `[[`, "stat" ),
    test = test0$type)
  rownames(re) = rownames(cmpset)
  attr(re,"test") <- ttype
  return(re)
}


# workhorse function for association testing
# z: covariates
# here y and x are reversed
# y: clusters so ~y
# x: outcomes so x~
#
# convention to return a test result object
# list(pval = pvalue, # numeric value
#      stat = stats: list(each stat element),
#      ...: will not be summarized)
#
tst <- function(fft, x, cmpset){
  # if(class(x)[1]=="integer") x=as.numeric(x);
  # n = ifelse( class(x) != "Surv", length(x), nrow(x) )

  # n = nrow(data.frame(x))
  # stat function, old way of doing this
  # fft <- eval(as.name(paste0("ff_",class(x)[1])))

  # collect all of the results
  test_results <- apply(cmpset[,c("c1_1", "c1_2", "c2_1", "c2_2", "cid1", "cid2")], 1, function(a){
    i1 = a[1]:a[2]
    i2 = a[3]:a[4]
    inds = c(i1,i2)
    y = factor( c(i1*0+1,i2*0+2), labels = a[5:6] )

    # if the outcome consists of more than 1 column
    return( fft( if( length(dim(x)) >1 ) x[inds,] else x[inds], y) )
  })

  # structure( cbind(t(test_results), level = cmpset[,"level"]), type = attr(fft, "type") )
  list( test_results = test_results, level = cmpset[,"level"], type = attr(fft, "type") )
}

# tests for each outcomes
ff_numeric <- function(x, y){

  tt=try(t.test(x~y))
  if(inherits(tt, "try-error")){
    pval=NA
    stat=NA
  }else{
    pval=tt$p.value
    stat=tt$statistic
  }
  return(list(pval = unname(pval), stat = stat))

}
attr(ff_numeric, "type")= "t.test"

ff_factor <- function(x, y){

  # no level
  if( length( levels(x) )<2 ) return( list(pval = NA, stat = NA) )

  # too many level
  if(length(levels(x)) == length(x)){
    warning("# of factor levels equals number of samples: length(levels(x)) = length(x)")
    return(list(pval = NA, stat = NA))
  }

  # 2 level: fisher test
  if(length(levels(x)) == 2){
    fu=try(fisher.test(x, y=y))

    if(inherits(fu, "try-error")){
      fu=try(fisher.test(x, y=y, simulate.p.value = T))
      if(inherits(fu, "try-error")) return(list(pval = NA, stat = NA))
    }
    pval=fu$p.value
    stat=fu$estimate
  }else{ # more than 2 levels
    cht = try( tryCatch( chisq.test(x, y=y), warning=function(w) w) )
    if(inherits(cht, "try-error")) return(list(pval = NA, stat = NA))

    if(inherits(cht, "warning")){ # if warning for simulate p val
      cht = try( tryCatch( chisq.test(x, y=y, simulate.p.value = T), warning=function(w) w) )
      if(inherits(cht, "try-error") | inherits(cht, "warning")) return(list(pval = NA, stat = NA))
    }
    pval=cht$p.value
    stat=cht$statistic
  }

  if(is.null(pval)) {pval=NA}
  if(is.null(stat)) {stat=NA}

  return( list(pval = unname(pval), stat = stat) )
}
attr(ff_factor, "type")= "fisher|chisq"

ff_Surv<-function(x, y){
  # does not make any difference in this case
  y=as.numeric(y)
  re <-
    try({
      cp <- survival::coxph(x~., data.frame(x=x, y=y) )
      coef = cp$coefficients["y"]
      # logrank p-value if no covariates, otherwise wald test with coef
      pval = -expm1(pchisq(cp$score, 1, log.p = T))
      list(pval = unname(pval), stat = c(HR = unname(exp(coef))))
    })

  if(inherits(re, "try-error")) return(list(pval = NA, stat = NA))
  return(re)
}
attr(ff_Surv, "type")= "log-rank"
