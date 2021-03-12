#' Plot subgroup identification results on a tree
#'
#' Plots the results of a subgroup identification run from \code{sgi::sgi_run()}.
#' Offers various options for association highlighting and tailored visualizations.
#'
#' @param as Association object returned by \code{sgi::sgi_run()}
#' @param padj_th Adjusted p-value threshold to filtering out significant enrichments
#'
#' @param outcomes_text.size Font size of outcomes on tree
#' @param outcomes_at_middle Position of outcomes on the tree, TRUE(default) for middle of samples or FALSE for at branching point on the tree
#' @param branching_point.size Size of branching point marking for significant outcomes
#' @param branching_point.color Color of branching point marking for significant outcomes
#'
#' @param tree_opt Shape of the tree branches: \cr
#' \code{"full_rect"} Rectangles (default) \cr
#' \code{"up_tri_down_rect"} Valid cluster pairs as triangles, rest as rectangles \cr
#' \code{"up_rect_down_tri"} Valid cluster pairs as rectangles, rest as triangles \cr
#' \code{"full_tri"} Triangles
#' @param plot_only_valid_cluster_pairs Plot only valid cluster pairs \code{FALSE} (default)
#' @param color.subgroups Single color for all subgroups, if \code{NULL} (default), each subgroup is colored differently
#' @param color.valid_cluster_pairs Color of valid cluster pairs, default is black
#' @param color.invalid_clusters Color of invalid clusters, default is gray
#' @param color_only_sgs_branches  Whether only significant subgroups are colored differently, \code{TRUE} (default), otherwise all valid clusters will get different colors.
#' @param lwds Line widths, vector of length two containing widths of invalid cluster lines and valid cluster lines
#' @param tree_pruning Options for tree pruning (i.e. removing bottom branches in the visualization): \cr
#' \code{0}: no pruning \cr
#' \code{1}: prune with \code{prune_height} (separate function argument) \cr
#' \code{2}: adaptive pruning, cut underneath lowest valid cluster pair  \cr
#' \code{3}: adaptive pruning, cut at highest leaf branching point (cuts away long leaf branches), pruned part is left blank \cr
#' \code{4}: same as option \code{3}, but pruned part is cut out.
#' @param prune_height Pruning height if \code{tree_pruning==1}
#'
#' @param lty.invalid_clusters Line type for invalid clusters, see linetypes options in R
#' @param lty.pruned Line type of cut branches after pruning (indicating that pruning has taken place)
#'
#' @param overlay Whether valid cluster pairs or subgroups are overlayed with rectangles, \code{TRUE} (default)
#' @param overlay_only_sgs \code{TRUE} (default): Overlay subgroups with rectangles,\code{FALSE}: overall all valid cluster pairs with rectangles
#' @param overlay.fill Single fill color of overlay rectangles for subgroups. \code{NULL} (default) to fill each subgroups differently.
#' @param overlay.alpha Alpha of overlay rectangles
#' @param overlay.color Border color of overlay rectangles
#' @param overlay.lwd Line width of overlay rectangles
#' @param overlay.lty Line type of overlay rectangles
#'
#' @param draw_cluster_ids Whether tree is annotated with cluster IDs, \code{TRUE} {default}
#' @param cluster_ids_for_only_sgs  \code{TRUE} {default}: Cluster IDs are only printed for subgroups, \code{FALSE}: cluster IDs are printed for all valid cluster pairs
#' @param cluster_ids.color Color of cluster ID text
#' @param cluster_ids_text.size Size of cluster ID text
#' @param cluster_ids_as_label Should IDs be printed with \code{geom_label} (\code{TRUE}), or \code{geom_text} (\code{FALSE}, default)
#'
#' @param vcs.lineend Line end type of segments of the valid cluster pairs, default is \code{"round"}, options: \code{c("round", "butt", "square")}
#'
#' @author mubu
#' @return ggplot object of annotated tree
#'
#' @method plot sgi.assoc
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # initialize sgi
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' # run sgi
#' as = sgi::sgi_run(sg)
#'
#' # plot sgi results with default options
#' gg = plot(as)
#' gg
#' }
plot.sgi.assoc <- function( as,          # association object from sgi::sgi_run
                            padj_th = 0.05, # p-value threshold for outcomes

                            #---------------------- outcome annotation on the tree
                            outcomes_text.size = 3.7,
                            outcomes_at_middle = TRUE,
                            branching_point.size = 3,
                            branching_point.color = "white",
                            # repel_geom = ggrepel::geom_label_repel,

                            #---------------------- tree visualisation parameters
                            tree_opt = "full_rect",          # full rectangle vcs & ivcs
                            plot_only_valid_cluster_pairs = FALSE, # plotting only vcs, overly simplified
                            color.subgroups = NULL,                # if NULL, each (non) sgs colored discreate,
                            color.valid_cluster_pairs = "gray",    # color of valid cluster pairs but not sgs, should be null when discreate
                            color.invalid_clusters = "gray",       # color of ivcs: invalid clusters
                            color_only_sgs_branches = TRUE,  # should only sg branches be colored
                            lwds = c(0.5, 3.0),              # line weights, lwds = (lwd.ivcs, lwd.vcs)
                            prune_height = NULL,             # min length to visualize hclust
                            tree_pruning = 4,                # no pruning(0) or prune with prune_height(1) or adaptive(2),
                                                             # or drop to min leaf height(3), squize to bottom after (3): (4)
                            lty.invalid_clusters = 2,        # lty of ivcs segments ie. dashed lines
                            lty.pruned = 2,                  # lty of pruned leafs i.e dotted lines

                            #---------------------- overlaying parameters
                            overlay = TRUE,          # should vcs be overlayed with rectangle
                            overlay_only_sgs = TRUE, # overlay only sgs
                            overlay.fill = NULL,     # filling of overlay, if NA, filled discreate
                            overlay.alpha = 0.085,   # alpha of overlay rectangles
                            overlay.color = "white", # color of rectangle lines
                            overlay.lwd =  0.1,
                            overlay.lty = 1,

                            #---------------------- cluster_id parameters
                            draw_cluster_ids = TRUE,          # should cluster ids be added on top of the tree parts
                            cluster_ids_for_only_sgs = TRUE,  # do that for only subgroups
                            cluster_ids.color = "black",      # color of cluster ids
                            cluster_ids_text.size = 4,        # size of cluster ids
                            cluster_ids_as_label = FALSE,     # should they be plotted as label or text

                            #---------------------- specs
                            vcs.lineend = "round",            # valid clusters line end type

                            #---------------------- if more
                            ...
                            ){

  # pr = print(as, padj_th = padj_th, by_clusters = T, silent = T)
  # if(is.null(pr) | length(pr) < 1){
  #   warning("no significant clusters for given p-value threshold!")
  #   pr = c()
  # }

  pr = do.call( rbind, lapply( summary(as, padj_th, by_clusters = T),
                               `[`,,c("outcome", "cid1", "cid2", "padj")) )

  if(is.null(pr) | nrow(pr) < 1){
    stop("no significant clusters for given p-value threshold!
         still want to plot the tree, try plotting sgi.object created by sgi::sgi()!")
  }
  sgs = paste( pr$cid1, pr$cid2, sep = "vs")

  # gg0 = plot(sg, sgs = names(pr), col.sgs = NA, col.non.sg = "gray90")
  gg0 <- suppressMessages(
    plot_tree_with_nested_clusters(
      #---------------------- sgi parameters
      hc = as$hc, uvc = as$cluster_pairs, sgs = unique(sgs) ,
      #---------------------- sgi parameters
      min_height = prune_height,
      tree_opt = tree_opt,
      plot_only_vcs = plot_only_valid_cluster_pairs,
      sgs.color = color.subgroups,
      ivcs.color = color.invalid_clusters,
      vcs.color = color.valid_cluster_pairs,
      color_only_sgs_branches = color_only_sgs_branches,
      lwds = lwds,
      tree_pruning = tree_pruning,
      lty.ivcs = lty.invalid_clusters,
      lty.pruned = lty.pruned,
      overlay = overlay,
      overlay.only_sgs = overlay_only_sgs,
      overlay.fill = overlay.fill,
      overlay.alpha = overlay.alpha,
      overlay.color = overlay.color,
      overlay.lwd =  overlay.lwd,
      overlay.lty = overlay.lty,
      draw_cluster_ids = draw_cluster_ids,
      cluster_ids_for_only_sgs = cluster_ids_for_only_sgs,
      cluster_ids.color = cluster_ids.color,
      cluster_ids_text.size = cluster_ids_text.size,
      cluster_ids_as_label = cluster_ids_as_label,
      vcs.lineend = vcs.lineend,
      ...
    )
  )

  # get cluster splits points and add outcomes to the enrichment points
  sdf = attr(gg0,"sdf")
  df.clusters = na.omit(sdf)
  df.clusters = df.clusters[sgs, ,drop = F]
  df.clusters$outcome = pr$outcome

  if(outcomes_text.size>0){
    # decide where to put outcomes labels 
    df.clusters$xpos <- if(outcomes_at_middle) df.clusters$cmid else df.clusters$xsplit
    
    gg0<-gg0 +
      geom_point(data = df.clusters,
                 aes(x = xpos, y = ymax),
                 size = branching_point.size,
                 color = branching_point.color ) +
      ggrepel::geom_label_repel( inherit.aes = F,
                                 data = df.clusters,
                                 box.padding = 0.7,
                                 size = outcomes_text.size,
                                 aes(x = xpos, y = ymax, label = outcome),
                                 segment.size = 0.3,
                                 fill = NA,
                                 segment.alpha = 0.3,
                                 max.overlaps = Inf)
  }
 gg0 = gg0 +
  labs(color ="subgroups",
       fill ="subgroups",
       lty = "valid cluster pairs") +
  guides(lwd = F, fill = F, lty = F)

  gg0$padj_th = padj_th
  gg0

}

