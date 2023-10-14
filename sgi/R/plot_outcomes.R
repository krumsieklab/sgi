
#' Outcome plots for subgroups
#'
#' Generates plots that visualize the differences of a list of outcomes with their associating cluster splits.
#'
#' @param sg SGI object created by \code{sgi::sgi()}
#' @param as SGI association object created by \code{sgi::sgi_run()}
#' @param outcome_names Names of outcomes to plot
#' @param cluster_pairs Names of cluster pairs to plot, e.g. \code{c("2vs3","4vs5")}
#' @param padj_th Adjusted p-value threshold
#' @param drop_nulls Whether clinical outcomes with no associations will be excluded from the result list (\code{TRUE}) or included as a \code{NULL} (\code{FALSE})
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # initialize SGI structure, minsize is set to 0.05 of total sample size
#' sg = sgi_init(hc, minsize = 18, outcomes = df_clinicalvars)
#' # run SGI
#' as = sgi_run(sg)
#' # inspect distributions of outcomes in subgroups
#' ggs = plot_outcomes(sg, as, padj_th = 0.05)
#' }
#'
plot_outcomes <- function(sg, as, outcome_names = NULL, cluster_pairs = NULL, padj_th = 0.05, drop_nulls =T,... ){

  params = list(...)

  # outcomes has to be chracter vector
  if(is.null(outcome_names)){
    df = sg$outcomes[sg$hc$order,names(as$results), drop = F]
  }else{
    df = sg$outcomes[sg$hc$order, outcome_names, drop = F]
  }

  ous = colnames(df)

  if( is.null(cluster_pairs) ){
    sgs = rownames( as$cluster_pairs )
  } else{
    sgs = cluster_pairs
  }
  #browser()
   # fix the names
  names(sgs) = sgs
  ous = structure(make.names(ous), names = ous)
  colnames(df) = make.names(colnames(df))
  as.res = as$results
  names(as.res) = make.names(names(as.res))

  N = length( sg$hc$order )

  gg_factor <- function(df_, nonsg = F){

    if(nonsg) df_$sub_group[is.na(df_$sub_group)]<-"NA_"
    df <-  reshape2::melt(table(df_))
    df[,2]<- factor(df[,2])
    df$vcount <- ave(df$value,df$sub_group, FUN = sum)
    df$prop = df$value/ df$vcount

    xbreaks = c(0,df$vcount) %>% sort %>% unique
    ybreaks = c(0,unique(df$sub_group) %>% lapply(function(x) df[df$sub_group==x, "prop"] %>% rev %>% cumsum) %>% unlist) %>%
      sort %>% unique %>% round(2)

    ggplot(df, aes(x = vcount/2, y = prop, width =vcount)) +
      geom_bar(stat = "identity", position = "fill",aes(colour = factor(sub_group)), lwd = 8) +
      scale_x_continuous(breaks = xbreaks, minor_breaks = xbreaks )+
      scale_y_continuous(breaks = ybreaks, minor_breaks = ybreaks )+
      facet_grid(as.formula("~sub_group"), scales = "free_x", space = "free_x") +
      geom_bar(stat = "identity", position = "fill", colour = NA, aes_string(fill = colnames(df)[2]), width = 0.5) +
      scale_fill_brewer(palette = 8,type = "seq") +scale_color_viridis_d() +
      theme(legend.position = "bottom",
            panel.background = element_rect(fill = NA),
            panel.grid = element_line(colour = "gray", linetype = 2),
            axis.text.x = element_text(vjust = 0.5),
            axis.text.y = element_text(angle = 30, hjust = 1) ) +
      labs(y = "distribution within sub-group",x = "sub-group size", color = "sub-groups") +
      scale_color_viridis_d(end = 0.6) +guides(color = F)
  }

  gg_numeric <- gg_integer <- function(df_ , nonsg = F, jitter = T){

    if(!nonsg) df_ =  na.omit(df_)
    df_$sub_group = factor(df_$sub_group)

    gg = ggplot(df_,aes_string(x = "sub_group", y = colnames(df_)[2]))
    if(jitter) gg = gg + geom_jitter(aes_string(color = "sub_group"), width = 0.25)
    gg +
      geom_boxplot(notch = T, alpha = 0.5, aes_string(fill = "sub_group"),
                   outlier.size = if(!jitter) 1 else -1, width = 0.8)+
      labs( x = "sub-group", fill = "sub-group", color  = "sub-group") +
      scale_color_viridis_d(end = 0.6) + scale_fill_viridis_d(end = 0.6)
  }

  gg_Surv <- function(df_ , nonsg = F){
    lvls = c( sort( unique( na.omit(df_$sub_group) ) ), "NA" )

    if(!nonsg) df_ =  na.omit(df_)
    else df_$sub_group[is.na(df_$sub_group)]= "NA"


    df_$sub_group = factor(df_$sub_group, lvls)
    Sname = colnames(df_)[2]
    sfit = survfit(df_[[Sname]]~df_$sub_group)

    gg = ggfortify:::autoplot.survfit(sfit, conf.int.alpha = 0.15) +
      labs(fill = "sub-group", color  = "sub-group", y = names(df_[2]) )

    gg + theme_minimal()
  }


  re<-
  lapply(ous, function(ou){
    res = as.res[[ou]]
    lapply(sgs, function(sg){
      # filtering based on pth
      cpadj = res[paste0(cid1,"vs",cid2)==sg,]$padj
      if (is.na(cpadj) | cpadj > padj_th) return(NULL)

      cpa = t(as$cluster_pairs[sg,])[,1]
      x = rep(NA, N)
      x[ cpa["c1_1"]:cpa["c1_2"] ] =  cpa["cid1"]
      x[ cpa["c2_1"]:cpa["c2_2"] ] =  cpa["cid2"]

      y = df[[ou]]
      df_ = data.frame(x,y)
      colnames(df_) = c("sub_group", ou)
      # if(is.Surv(y)) browser()
      try( do.call(paste0("gg_", class(y)), c(list(df_),params)) + theme_minimal() )

  })})

  # browser()
  if(drop_nulls){
    f <- function(x){
      x = x[!sapply(x, is.null)]
      if(length(x)<1) return(NULL)
      x
    }
    f(lapply(re, f))
  }else re

}



#' to_do:
#' 1. add results to subtitle/title i.e. pvals
#' 2. keep names with space as they are
#'
#'

