
#' Overview plot
#'
#' Plots all parts of an SGI analysis in one combined overview plot. Shows the annotated SGI association tree on top, and outcomes and data matrix heatmaps underneath.
#'
#' @param gg_tree ggplot object created by \code{plot(association_object)} (see \code{plot.sgi.assoc()})
#' @param as SGI association results, an object of class \code{sgi.assoc} created by \code{sgi::sgi_run()}
#' @param outcomes Names of the outcomes to be plotted in heatmap underneath the tree
#' @param xdata Data matrix of hierarchical clustering to be plotted underneath the tree (usually the original matrix used for clustering).
#' @param h_tree Height of the tree in the plot, relative to 1 plotted outcome strip
#' @param h_xdata Height of the data matrix, relative to 1 plotted outcome strip
#' @param padj_th Adjusted p-value threshold to decide which splits to mark; if \code{NULL}, set to \code{padj_th} with which \code{gg_tree} was created
#' @param outcomes_nmax Maximum number of outcomes, i.e. strip tiles, to be plotted, default is first 10 variables in outcome data.frame.
#' @param ncut_tree Cutting the tree to visually emphasize some clusters i.e. number of clusters to highlight \code{ncut_tree>1}
#' @param split_markers_for_subgroups Should split markers for subgroups be added into the strip tiles? \code{TRUE} by default
#' @param split_markers_only_for_cuts Mark splits only at the level of \code{ncut_tree} if given. FALSE by default.
#' @param split_marker.size Size of split marker
#' @param split_marker.color Color of split marker
#' @param split_marker.fill Filling color of split marker
#' @param split_marker.alpha Alpha of split marker
#' @param split_marker.shape Shape of the split marker point, one of \code{pch} options
#' @param data.color Color palette of data matrix heatmap, see color parameter of \code{pheatmap::pheatmap()}
#' @param data_cluster_cols Should data matrix columns be clustered? \code{0}: data not clustered, \code{1} data clustered, but cluster tree not shown, \code{2} data clustered and tree shown
#' @param data_title Title to be added to the right of data matrix
#' @param draw_legends Draw legends? \code{FALSE} by default; Note: Legends are returned as a list of \code{grid} objects within a ggplot object
#' @param ...
#'
#' @return Overview plot of SGI analysis
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # hierarchical clustering
#' hc = hclust(data_matrix)
#' # SGI initialization
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' # association analysis
#' as = sgi::sgi_run(sg)
#' # plot the association analysis tree
#' gg_tree = plot(as, padj_th = 0.001)
#' # plot overview
#' plot_overview( gg_tree, as, outcomes, data_matrix, draw_legends = T, ncut_tree = 3)
#'
#' }
#'
#' @author mubu
#'
#'
plot_overview <- function(gg_tree,
                          as,
                          outcomes = NULL,
                          xdata = NULL,
                          h_tree = 25,
                          h_xdata = 10,
                          padj_th = NULL,
                          outcomes_nmax =10,

                          ncut_tree = 0,
                          split_markers_for_subgroups = TRUE,
                          split_markers_only_for_cuts = FALSE,
                          split_marker.size = 2.5,
                          split_marker.color = "red",
                          split_marker.fill = "red",
                          split_marker.alpha = 1,
                          split_marker.shape = 19,

                          data.color = colorRampPalette(c("blue", "floralwhite", "red"))(100),
                          data_cluster_cols = 2,#if(h_xdata<15) T else 2, # tree not shown yet
                          data_title = "data",
                          draw_legends = F,
                          ...
                          ){

  # deprecated color palatte
  #colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdYlBu")))(30),
  if(is.null(padj_th)) padj_th = gg_tree$padj_th
  if(is.null(padj_th) & split_markers_for_subgroups) stop("padj_th should be given if split_markers_for_subgroups=T")

  add_ys = !is.null(outcomes)
  add_xs = !is.null(xdata)

  if( (add_ys+add_xs) < 1 ) stop("nothing to add underneath sgi tree, try plot(as)!")

  check_consistency_with_as(as = as, outcomes = outcomes, xdata = xdata)

  # # scale max height to 24, add some space above for ggrepel
  # # max will be 25
  # as$hc = fix_height(as$hc, hscale = h_tree*0.95)
  hc= as$hc

  # sgi tree
  # gg_tree = plot(as, relavance.tree.pruning = relavance.tree.pruning, padj_th = padj_th, overlay. = overlay.) + ylim(c(0, h_tree))

  if(gg_tree$padj_th != padj_th) warning("padj_th with which tree is plotted does not match with the padj_th given!")

  # change max y value of gg_tree so that some room for repels

  # phetotype data
  # max outcomes to be ploted
  if(add_ys){
    gg_phes = p_yheal(outcomes[hc$order, , drop=FALSE], nmax =  outcomes_nmax)
    if (outcomes_nmax < ncol(outcomes)) warning("outcomes displayed in phenotype tiles is less than total #outcomes")
    n_tiles = length(gg_phes)
  }else{
    gg_phes = NULL
    n_tiles = 0
  }

  if(add_xs){
    gg_data = p_xhea(xdata[hc$order, ], color = data.color, cluster_rows = data_cluster_cols>0)
  }else{
    h_xdata = 0
  }


  # p_add_ytiles should handle if no tiles to be plotted
  obj = p_add_ytiles(gg_tree, gg_phes, ty = h_tree)
  legends = obj$legends
  gg = obj$gg


  # ggplot(data.frame(x=rnorm(10), y = rnorm(10)), aes(x,y))+
  #   geom_point() + annotation_custom(aa, 0,1,0,1) +
  #   scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  #   scale_y_continuous(expand = c(0,0), limits = c(0,1))

  if(add_xs){
    # extract the legend object from gtable
    legends$xdata = gtable::gtable_filter(gg_data$gtable, "legend")$grob[[1]]
    gg = gg +
      annotation_custom(grob = gg_data$gtable$grobs[[2]],
                        xmax = 100,
                        xmin = 0,
                        ymax = -n_tiles -h_tree/50,
                        ymin = -n_tiles -h_tree/50 - h_xdata)
  }
  gg = gg + xlim(c(-10,101)) + ylim(c(-n_tiles -2*h_tree/50 - h_xdata, h_tree*1.04))

  #browser()

  if(add_xs){
    if(!is.null(data_title)){
      gg= gg + annotate('text', x = 101, y = -n_tiles -h_tree/50 - h_xdata/2,
                        label = data_title, angle = -90, size = min(12,h_xdata*0.8), vjust = 0)
    }

    if(data_cluster_cols >1){
      gg= gg +
        annotation_custom(
          grob = gg_data$gtable$grobs[[1]],
          xmin = -10, xmax = -1,
          ymax = -n_tiles -h_tree/50,
          ymin = -n_tiles -h_tree/50 - h_xdata)
    }
  }

  if(ncut_tree>1){
    # cluster points
    v =  get_cluster_points(hc, ncut_tree)
    gg = gg + geom_segment(data = data.frame(x = 100*v, xend = 100*v, y = 0, yend = -Inf),
                           aes(x = x , y = y, xend = xend, yend = yend),
                           color = "white",  lwd=1.5, lty = 1)
  }

  if(split_markers_for_subgroups & add_ys){
    # browser()
    # significant hits
    ih = summary(as, padj_th = padj_th, by_clusters = F)
    ys = intersect(names(ih), names(gg_phes))
    #browser()
    if(length(ys)>0){
      ih = ih[ys]
      # get where to put red dots
      max_cut = if(!split_markers_only_for_cuts)  Inf else ncut_tree
      if(ncut_tree<2){
        #warning("split_markers_only_for_cuts = TRUE but ncut_tree < 2, thus no tree cut to place the hits on!")
        max_cut = Inf
      }

      which_clusters = lapply(lapply(ih,`[`,,level), function(x) x[x<(max_cut+2)]-1)
      which_clusters = which_clusters[sapply(which_clusters, length)>0]

      # if ever there is to be annotated
      if(length(which_clusters)>0){
        which_clusters = do.call(rbind, lapply(seq(which_clusters), function(i)
          data.frame(yname = names(which_clusters)[i],
                     cluster = which_clusters[[i]],
                     stringsAsFactors = F)))

        which_clusters$yi = -structure(seq(gg_phes),names = names(gg_phes))[ which_clusters$yname] + h_tree/50
        if(!(split_markers_only_for_cuts & ncut_tree>1)){
          # we need more clusters
          v =  get_cluster_points(hc, max(which_clusters$cluster)+1)
        }
        which_clusters$xj = v[which_clusters$cluster]*100

        gg = gg +geom_point(data = which_clusters, aes(x=xj,y=yi),
                            color = split_marker.color,
                            fill = split_marker.fill,
                            alpha = split_marker.alpha,
                            shape = split_marker.shape,
                            size = split_marker.size)
      }
    }

  }


  if(!draw_legends) return(structure(gg, legends_all = legends))

  # combine legends and add into the plot
  # browser()

  gg <- gg +
    annotation_custom(grob = legends$tree_legend,
                      xmin = 105, xmax = 115,
                      ymin = 0, ymax = h_tree)+
    xlim(c(-12, 115))
  if(add_xs){
    gg = gg +
      annotation_custom(grob = legends$xdata,
                        xmin = 107, xmax = 115,
                        ymax = -n_tiles -h_tree/50,
                        ymin = -n_tiles -h_tree/50 - h_xdata)
  }

  if(add_ys){
    xmax = 115+ceiling(length(legends$ytiles)/4)*10
    gg = gg+
      annotation_custom(
        grob = gridExtra::arrangeGrob(grobs= legends$ytiles, nrow = 4),
        xmin = 115, xmax = xmax,
        ymax = h_tree,
        ymin = -n_tiles -h_tree/50 - h_xdata) +
      xlim(c(-12, xmax))
  }
  gg
}


# heatmap of data X
# later add other parameters, like type of clustering etc
p_xhea <- function(dt, color=NULL, cluster_rows = T ){
  # for aes_string, could be more cleaver solution with aes_q
  cnames = colnames(dt)
  colnames(dt) = make.names(cnames)

  pheatmap::pheatmap(t(dt),
                     silent = T,
                     cluster_rows = cluster_rows,
                     clustering_distance_rows = "correlation",
                     cluster_cols = F,
                     #labels_row = "", labels_col = "",
                     color = if(is.null(color)) colorRampPalette(c("blue", "white", "red"))(100)
                     else color,
                     scale = "none")
}

# dont support survival yet
# tiles of outcomes as list
p_yheal <- function(dp, nmax =  10){

  # for aes_string, could be more cleaver solution with aes_q
  inds = sapply(dp, function(x) class(x) %in% c("integer", "numeric", "factor") )
  if(any(!inds)) warning("only outcomes of types: factor, numeric are supported!")
  dp = dp[,inds, drop = F]
  if(ncol(dp)<1) return(NULL)

  dp = dp[, seq(min(nmax, ncol(dp))), drop = F]
  cnames = colnames(dp)
  colnames(dp) = paste0( make.names(colnames(dp)), seq(cnames) )
  df = data.frame(x = seq(nrow(dp)), dp)

  structure( lapply(structure(seq(cnames),names = cnames), function(i)try({
    ggplot( df[,c(1,i+1)], aes(x=x)) +
      geom_tile(aes_string(y=1,  fill= colnames(df)[i+1])) +
      {if(is.factor(df[[i+1]])) scale_fill_viridis_d() else scale_fill_viridis_c()} +
      theme_void() +
      labs(fill = cnames[i])
  })), class = c("pgg_list", "list"))
}

# fix tree height for controlling height ratio, you may change this later to reflect anno custom trick
fix_height <- function(hcx, hscale = 25){
  h = hcx$height
  h = h/max(h) * hscale
  hcx$height = h
  hcx
}

get_legend_from_grobs <- function(grobs){
  l = grobs[sapply(grobs, function(x) x$name == "guide-box")]
  if(length(l)<1) return(NULL)
  l
}

# add tiles underneath to tree.
p_add_ytiles <- function(ggt, tlist, ty = 25){

  # ratios
  # tree heigh: 25,
  # one tile h: 1,
  # convert tree to grob
  obj = ggplotGrob(ggt +
                       scale_x_continuous(expand = c(0, 0))+
                       scale_y_continuous(expand = expansion(mult = c(1/50, 1/10)))
                   )

  tree_grob =  obj$grobs[[5]]
  legends = list()
  legends$tree_legend = get_legend_from_grobs(obj$grobs)[[1]]

  # place to 0-100, 0-ty
  gg_ = ggplot() +
    xlim(c(-1,101)) + ylim(-length(tlist)-1, ty*1.03)+
    theme_void()+
    annotation_custom(tree_grob, xmin = 0, xmax = 100, ymin = 0, ymax = ty)
  if(is.null(tlist))  return(list(gg = gg_, legends = legends))

  k = length(tlist)
  for(i in seq(k)){
    bobj =  ggplotGrob(tlist[[i]]+
                       scale_x_continuous(expand = c(0, 0))+ #expansion(mult = c(0.02, 0.02)) )
                       scale_y_continuous(expand = c(0,0 )) )
    bb = bobj$grobs[[5]]
    legends$ytiles[names(tlist)[i]] = get_legend_from_grobs(bobj$grobs)
    gg_ = gg_+ annotation_custom(bb, ymax =1-i, ymin = -i, xmin = 0, xmax = 100) +
      annotate("text", x=-1, y=0.5-i, label= names(tlist)[i],  hjust = 1)
  }
  list(gg = gg_, legends = legends)
}

# given hc and k gives the cluster points where gaps should be
# mid_d added to mid of two clusters i.e.  (c1_2+c2_1)/2+mid_d
get_cluster_points <- function(hc, k, mid_d = 0){
  uc = get_uc_vc(hc, -1)$clinds
  uc = uc[uc$level <= k,]
  v = 0.5*(uc$c1_2 + uc$c2_1) + mid_d
  v/max(uc$c2_2)
}


# check if consistend objects are supplied,
# order of patients labels etc.
check_consistency_with_as <- function(as, sg=NULL, outcomes=NULL, xdata=NULL){
  TRUE
}





# #  change the implementation to get precalculated gg object for tree
# # plot overview of analysis
# #  tree + outcomes+ data
# #
# #  add aes params like ratios
# #
# #  refrain from annotation of results
# #  to-do:
# #  1. option for adding metaboliets names
# #  2. rich set of cutting functionality
# #  3. sgi outcomes annotations like p values, size with log10(p)
# #  4. better tree cutting.
# #  5. overlapping legend and data title
# #
# #  eger number of cut verilirse, as ile padj_th cek et ve ona gore kirmizi nokta ekle
# #
