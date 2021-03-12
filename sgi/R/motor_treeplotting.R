#' Plotting and annotating hierarchical clustering tree with cluster information
#'
#' This functions offer various options to visualize, manipulate and annotate the hierarchical tree with the clustering information such as with the cluster pairs.
#' Some of the capabilities of the function are pruning the tree according to cluster info, annotation of the tree with cluster ids, coloring tree parts accordingly and
#' adding overlaying rectangles representing the cluster information. This is designed to be an utility function for the visualization of association results in the SGI package.
#'
#' @param hc hierarchical clustering tree, \code{hclust} object
#' @param uvc unique valid clusters , \code{data.frame} returned by \code{get_uc_vc()}, or \code{sg$clinds}
#' @param sgs names of subgroup pairs, i.e \code{c("2vs3","6vs7")} that tree will be annotated with. If missing, will be set to all of the cluster pairs in \code{uvc}.
#'
#' @param min_height height threshold for plotting the tree
#' @param tree_opt shape of the tree: \cr
#' \code{"full_rect"} rectangular tree (default) \cr
#' \code{"up_tri_down_rect"} valid cluster pairs as triangles and rest as rectangles \cr
#' \code{"up_rect_down_tri"} valid cluster pairs as rectangles and rest as triangles \cr
#' \code{"full_tri"} triangular tree
#' @param plot_only_vcs plotting only vcs \code{FALSE} (default), overly simplified
#' @param sgs.color color of subgroups, if \code{NULL} (default), each subgroups is colored discreate
#' @param ivcs.color color of part of the tree which is not valid cluster pairs
#' @param vcs.color color of valid cluster pairs which are not subgroups
#' @param color_only_sgs_branches  if only subgroups are colored discreate, \code{TRUE} (default)
#' @param lwds  line widths, vector of length two as  \code{c(lwd.ivcs, lwd.vcs)} widths of invalid clusters and valid cluster pairs
#' @param tree_pruning option for pruning the tree: \cr
#' \code{0}: no pruning \cr
#' \code{1}: prune with \code{min_height} \cr
#' \code{2}: adaptive pruning based on height of the highest invalid cluster segment below the lowest valid cluster segment \cr
#' \code{3}: adaptive pruning based on height of the highest segment of the tree leaves but pruned part left blank \cr
#' \code{4}: same as option \code{3}, but tree pulled down to the bottom, so pruned part is removed without left blank
#' @param lty.ivcs linetype of the segments which are not of valid clusters
#' @param lty.pruned linetype of pruned part of the tree
#'
#' @param overlay should valid cluster pairs or subgroups be overlayed with rectangles, \code{TRUE} (default)
#' @param overlay.only_sgs overlaying rectangles for only subgroups, \code{TRUE} (default), \code{FALSE} for all valid cluster pairs to be overlayed i.e. subgroups or not
#' @param overlay.fill filling color of the overlaying rectangles, \code{NULL} (default) to fill discreate with the cluster pairs ids
#' @param overlay.alpha alpha of overlaying rectangles
#' @param overlay.color color of the overlaying rectangles
#' @param overlay.lwd line width of the overlaying rectangles
#' @param overlay.lty linetype of the rectangles
#'
#' @param draw_cluster_ids  should tree be annotated with cluster ids, \code{TRUE} {default}
#' @param cluster_ids_for_only_sgs cluster ids printed only for subgroups, \code{TRUE} {default}, \code{FALSE} ids printed for all valid cluster pairs
#' @param cluster_ids.color color of cluster ids
#' @param cluster_ids_text.size size of text of the cluster ids
#' @param cluster_ids_as_label \code{FALSE} (default), should ids be printed with \code{geom_label} or  \code{geom_text} (default)
#'
#' @param vcs.lineend lineend type of segments of the valid cluster pairs, default is \code{"round"}
#'
#' @author mubu
#' @return ggplot object of the manipulated tree
#'
#'
#' @examples
#' \dontrun{
#'
#' # hierarchical clustering
#' hc = hclust( dist( sgi::kirc.rppa ), "ward.D2")
#' # cluster information
#' clinds = get_uc_vc(hc, length(hc$height)/20)$clinds
#'
#' gg <- plot_tree_with_nested_clusters(
#'    hc = hc, uvc =  clinds, sgs = c("2vs3", "4vs5", "22vs23"),
#'    tree_opt = "up_rect_down_tri",
#'    cluster_ids_as_label = T,
#'    cluster_ids.color = "blue",
#'    cluster_ids_for_only_sgs = F,
#'    lty.ivcs = 1,
#'    tree_pruning = 4,
#'    vcs.color = "black",
#'    plot_only_vcs = F,
#'    overlay = T,
#'    overlay.color = NA,
#'    overlay.lwd = 2,
#'    overlay.alpha = 0.1,
#'    overlay.only_sgs = T
#' )
#'
#' # plot the tree with different aesthetics
#' gg
#' gg+coord_polar()
#' gg+scale_y_reverse()+coord_polar()
#' }
#'
#' @noRd
#'
plot_tree_with_nested_clusters <- function(
  #---------------------- sgi parameters
  hc,          # hierarchical tree from hclust
  uvc,         # clusters info supplied by get_uc_vc: unique valid clusters
  sgs,         # names of sub group pairs

  #---------------------- tree visualisation parameters
  min_height = NULL,               # min length to visualize hclust
  tree_opt = c("full_rect",        # full rectangle vcs & ivcs
               "up_tri_down_rect", # vcs triangle ivcs rectangle
               "up_rect_down_tri", # vcs rectangle ivcs triangle
               "full_tri"),        # full triangle vcs & ivcs
  plot_only_vcs = FALSE,           # plotting only vcs, overly simplified
  sgs.color = NULL,                # if NULL, each (non) sgs colored discreate,
  ivcs.color = "gray",             # if color, all (non) sgs will be colored same with the color, ivcs: invalid clusters
  vcs.color = ivcs.color,          #
  color_only_sgs_branches = TRUE,  # should only sg branches be colored
  lwds = c(0.5, 3.0),              # line weights, lwds = (lwd.ivcs, lwd.vcs)
  tree_pruning = 4,                # no pruning(0) or prune with min_height(1) or adaptive(2), or drop to min leaf height(3),
                                   # squize to bottom after (3): (4)
  lty.ivcs = 2,                    # lty of ivcs segments ie. dashed lines
  lty.pruned = lty.ivcs,           # lty of pruned leafs i.e dotted lines

  #---------------------- overlaying parameters
  overlay = TRUE,          # should vcs be overlayed with rectangle
  overlay.only_sgs = TRUE, # overlay only sgs
  overlay.fill = NULL,     # filling of overlay, if NA, filled discreate
  overlay.alpha = 0.085,   # alpha of overlay rectangles
  overlay.color = "white", # color of rectangle lines
  overlay.lwd =  0.1,
  overlay.lty = 1,

  #---------------------- cluster_id parameters
  draw_cluster_ids  = TRUE,        # should cluster ids be added on top of the tree parts
  cluster_ids_for_only_sgs = TRUE,  # do that for only subgroups
  cluster_ids.color = "black",      # color of cluster ids
  cluster_ids_text.size = 4,        # size of cluster ids
  cluster_ids_as_label = FALSE,         # should they be plotted as label or text

  #---------------------- specs
  vcs.lineend = c("round", "butt", "square"), # valid clusters line end type

  #---------------------- if more
  ...
){

  tree_opt <- match.arg(tree_opt)
  vcs.lineend <- match.arg(vcs.lineend)

  if(!is.null(min_height)){
    if(tree_pruning!=1){
      warning("tree_pruning option set to (1) with given min_height.")
      tree_pruning = 1
    }
  }else{
    if(tree_pruning==1) stop("for tree_pruning option (1), min_height should be given!")
  }

  if(tree_pruning==0) min_height =0

  if(missing(sgs)){ # set all of the vcp in uvc
    mm = na.omit(uvc)
    sgs = paste( mm[,"cid1"], mm[,"cid2"], sep = "vs")
  }

  if(is.null(sgs)) sgs = c()

  # # of samples = # of possible clustering including 1 cluster
  k = length(hc$order)

  # height(y axis) intervals of every split
  h = c(Inf, rev(hc$height),0)
  h = cbind(h[-length(h)], h[-1] )
  colnames(h) = c("ymax", "ymin")

  # clusters information for valid pairs
  cui = na.omit(uvc)

  # get coordinates of tree segments, for each level
  dfaptree <- function(df, h, cui, add_sm = T){

    df$g = NA # id of vcp, NA means not vc
    if(add_sm) sm = c()  # valid cluster summary vector

    for(i in seq(nrow(cui)) ){
      y = h[i, ,drop = F]
      lux = unlist(cui[i,])
      x1 = lux[c("c1_1", "c1_2")]
      x2 = lux[c("c2_1", "c2_2")]

      # if valid
      # if( !is.na(lux["cid1"]) & !is.na(lux["cid2"]) ){
      inds = (df$y <= y[1]) & (df$y > y[2]) &
        (df$x <= x2[2]) & (df$x >= x1[1]) &
        (df$xend <= x2[2]) & (df$xend >= x1[1])

      # give the vcp id
      tag =  paste(lux[c("cid1","cid2")], collapse = "vs")
      df$g[inds] <- tag

      if(add_sm){
        yends = sort(unique(c(df$yend[inds], df$y[inds])))
        # tree ....
        #   |-----|     ymax
        #   |    y_smax
        # y_smin

        # idk if it ever happen but, this is for if left and right same heights
        if(length(yends)<3) yends = c(rep(yends[1], 3-length(yends)),yends)
        names(yends) = c("y_smin", "y_smax", "ymax")
        xends = df$xend[inds]
        xends = c(xmin = min(xends), xmax = max(xends))
        # rowname the summary vector with vcp id
        sm = rbind(c(x1[1],x2[2], cmid = median(c(x1,x2)),  yends, xsplit =  mean(df$x[inds]), xends), sm)
        rownames(sm)[1] <- tag
      }
      # }
    }
    if(add_sm) attributes(df)$sm <- sm
    df
  }

  # get the triangle and rectangele segment data of the tree
  # if(tree_opt != "full_tri"){
  #   df.rect = ggdendro::dendro_data(hc, type = "rectangle")$segments
  #   df.rect = dfaptree(df = df.rect, h = h[cui$level,], cui = cui)
  # }
  # do it all the time, otherwise I couldn't solve the summary data of triangle problem

  # summary of segments of valid cluster pairs, and all other segments
  sm_tree <- df.rect <- df.tria <- NULL

  if(tree_opt != "full_tri"){
    df.rect = ggdendro::dendro_data(hc, type = "rectangle")$segments
    df.rect = dfaptree(df = df.rect, h = h[cui$level,,drop = F], cui = cui)
    sm_tree = attr(df.rect, "sm")
  }

  if(tree_opt != "full_rect"){
    df.tria = ggdendro::dendro_data(hc, type = "triangle")$segments
    df.tria = dfaptree(df = df.tria, h = h[cui$level,,drop=F], cui = cui)
    if(is.null(sm_tree)) sm_tree = attr(df.tria, "sm")
  }

  # different tree visualisation options
  if(tree_opt == "up_tri_down_rect"){
    df = rbind(df.tria[!is.na(df.tria$g),], df.rect[is.na(df.rect$g),])
  }else if(tree_opt == "up_rect_down_tri"){
    df = rbind(df.tria[is.na(df.tria$g),], df.rect[!is.na(df.rect$g),])
  }else if(tree_opt == "full_rect"){
    df = df.rect
  }else if(tree_opt == "full_tri"){
    df = df.tria
  }else{
    stop("invalid tree_opt!")
  }

  # h of min leaf
  min_leaf_h = min(df$y)

  # keep only valid clusters, overly simplified
  if(plot_only_vcs){
    df = df[!is.na(df$g),]
    # no invalid cluster segments so adjust lwd accordingly
    lwds = rev(lwds)[1]
    lwds = c(lwds, lwds)
  }

  # summary of line segments of valid pairs
  sdf = data.frame(sm_tree)

  # browser()
  # colnames(sdf) <- c("xmin","xmax","mid","ymax","ymin","yend1","yend2","xbranch","xend1","xend2")

  # cluster pairs ids
  # sgs: subgroups, vcp: valid clusters
  sdf$vcp <- sdf$sgs <- rownames(sdf)

  # if sgs NULL, then all NA
  sdf$sgs[sdf$sgs %in% setdiff(sdf$vcp, sgs)] = NA

  # overlaying ids
  sdf$ovl <- if(!overlay.only_sgs) sdf$vcp else  sdf$sgs

  df$lty = 1 # lty of line segments
  # lty of ivcs
  df$lty[is.na(df$g)] = lty.ivcs

  # should tree be pruned as keeping only relavent part
  h.prune = NA
  if(tree_pruning > 0){

    h.prune = min_height # default is prune with min_height
    # adaptive pruning based on height of the highest ivcs below the lowest vcp
    if(tree_pruning == 2){
      tmp = min( df[!is.na(df$g), c("y","yend")] )
      h.prune = mean(h[max(which(h[,2] >= tmp))+1,])
      # print(h.prune)
    }
    if(tree_pruning >2){
      # min h of leafs
      h.prune = min_leaf_h
      # # min_height is the pruned point
      # if(tree_pruning != 3) min_height = h.prune
      # else min_height = 0
    }
    df = df[ !(df$y <=  h.prune & df$yend <= h.prune), ]
    linds = df$yend < h.prune # pruned leaves
    df$yend[linds] = h.prune
    df$lty[linds] = lty.pruned

  }

  # is valid cluster pairs, technically we don't need this because of NAs below
  # if NA, no vcs, T: for subgroups, F: for vcs
  df$vc = !is.na(df$g) # flag for valid clusters for size
  # is subgroup
  df$sg = FALSE # flag for subgroups
  df$sg[df$g %in% sgs] = TRUE
  df$sg[is.na(df$g)] = NA

  #--- re order factors for coloring
  df$g = factor(df$g, levels = sgs)
  sdf$vcp = factor(sdf$vcp, levels = rev(sdf$vcp))
  sdf$ovl = factor(sdf$ovl, levels = levels(sdf$vcp))
  # sdf = sdf[order(sdf$ovl),]
  #---

  # color only sgs not all vcs
  df$coloring = df$g
  if(color_only_sgs_branches){ df$coloring[-which(df$sg)] = NA }

  # fix the issue of nothing to color, etc
  # browser()
  # anything_to_color = length(levels(df$coloring)) > 0

  #---gg plotting
  # overlaying with rectangle
  gg <-
  if(overlay){

    # if I am not mistaken min of hs is on left and max is on the right by desing for hc
    #---- data for overlay rectangles

    d1 = sdf[,c("c1_1", "cmid", "y_smin", "ovl")]
    d2 = sdf[,c("cmid", "c2_2", "y_smax", "ovl")]
    colnames(d1) <- colnames(d2) <- c("xmin", "xmax", "hmax", "ovl")
    df_overlay = na.omit(rbind( d1, d2))
    #----

    # update rectangles according to pruning
    if(tree_pruning %in% c(2,4)) min_height = h.prune
    if(tree_pruning == 3) min_height = 0
    #

    ggplot(data = df_overlay) +
    if( is.null(overlay.fill) ){ # fill with overlay levels sgs or vcs
      geom_rect( aes(xmin = xmin, xmax = xmax, ymin = min_height, ymax = hmax, fill = ovl),
                 alpha = overlay.alpha, size = overlay.lwd, lty = overlay.lty, color = overlay.color)
    }else{ # fill with given color
      geom_rect( aes(xmin = xmin, xmax = xmax, ymin = min_height, ymax = hmax),
                 alpha = overlay.alpha, size = overlay.lwd, lty = overlay.lty, color = overlay.color, fill = overlay.fill)
    }

  }else{
    ggplot()
  }

  # add valid clusters which are not subgroups
  inds_only_vc = df$vc & !df$sg
  if(any(inds_only_vc)){
    gg = gg +
      geom_segment( data = df[inds_only_vc,],
                    aes(x = x, y = y, xend = xend, yend =yend, size = vc, lty = factor(lty)),
                    lineend = vcs.lineend, color = vcs.color) +
      scale_linetype_manual(values = sort(unique(df[inds_only_vc, "lty"])))
  }

  # add invalid clusters and subgroups

  gg<-
  if( is.null(sgs.color) ){ # color sgs with factor levels
    gg + geom_segment( data = df[!inds_only_vc, ],
                       aes(x = x, y = y, xend = xend, yend =yend, color = coloring, size = vc, lty = factor(lty)),
                       lineend = vcs.lineend)+
      scale_colour_discrete(na.value = ivcs.color) +
      scale_linetype_manual(values = sort(unique(df[!inds_only_vc, "lty"])))
  }else{ # color sgs and valid clusters only
    gg + geom_segment( data = df,
                       aes(x = x, y = y, xend = xend, yend =yend, color = factor(sg, levels = c(FALSE,TRUE)), size = vc, lty = factor(lty)),
                       lineend = vcs.lineend) +
      scale_colour_manual(values = c(vcs.color, sgs.color), na.value = ivcs.color, drop = F)+
      scale_linetype_manual(values = sort(unique(df[, "lty"])))
  }
  gg <- gg + theme_void() + scale_size_discrete(range = lwds)

  # leave the pruned part blank
  if(tree_pruning==3) gg = gg + ylim(c(0, NA))

  attributes(gg)$sdf <- sdf

  # should cluster ids be added
  if(draw_cluster_ids ){

    df.names <-
      data.frame( cluster_id = as.numeric(c(do.call(rbind,strsplit(rownames(sdf),"vs")))),
                  x = c(sdf$xmin,sdf$xmax),y = c(sdf$ymax, sdf$ymax), pairs = rep(sdf$vcp,2),
                  is_sg = rep(sdf$ovl,2))

    if(cluster_ids_for_only_sgs) df.names = na.omit(df.names)

    gg <- gg +
    if(cluster_ids_as_label){
      geom_label(data = df.names, aes(x,y,label=cluster_id),fontface = "bold",
                color = cluster_ids.color, size = cluster_ids_text.size)#, vjust = -0.25)
    }else{
      geom_text(data = df.names, aes(x,y,label=cluster_id),fontface = "bold",
                color = cluster_ids.color, size = cluster_ids_text.size, vjust = -0.5)
    }
  }

  gg$h.prune = h.prune
  return(gg)
}

# to.dos
# 1. beta testing every combinations of input args
# 2. cluster size depiction over leafs
#    a. branch segment size
#    b. scaled point
# 3. optimizing adaptive pruning
# 4. implementing annotative overlayed outcome colored trees
#
# if only vcs be plotted then give as sgs

