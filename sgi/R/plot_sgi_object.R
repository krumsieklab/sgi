#'Plot \code{SGI} object
#'
#'This function provides detailed visualization options for \code{SGI} objects before association analysis.
#'Plotting an \code{SGI} initialization object before performing the association analysis can be used to guide the statistical design.
#'
#' @param sg \code{SGI} initialization object created by \code{sgi::sgi()}
#' @param tree_opt Shape of the tree branches: \cr
#' \code{"full_rect"} Rectangles (default) \cr
#' \code{"up_tri_down_rect"} Valid cluster pairs as triangles, rest as rectangles \cr
#' \code{"up_rect_down_tri"} Valid cluster pairs as rectangles, rest as triangles \cr
#' \code{"full_tri"} Triangles
#' @param color.valid_cluster_pairs Color of valid cluster pairs, default is black
#' @param color.invalid_clusters Color of invalid clusters, default is gray
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
#' 
#' @param overlay Whether valid cluster pairs are overlayed with rectangles, \code{FALSE} (default)
#' @param overlay.fill Single fill color of overlay rectangles. \code{NULL} (default) to fill each cluster pair differently.
#' @param overlay.alpha Alpha of overlay rectangles
#' @param overlay.color Border color of overlay rectangles
#' @param overlay.lwd Line width of overlay rectangles
#' 
#' @param draw_cluster_ids Whether tree is annotated with cluster IDs, \code{TRUE} {default}
#' @param cluster_ids.color Color of cluster ID text
#' @param cluster_ids_text.size Size of cluster ID text
#' @param cluster_ids_as_label Should IDs be printed with \code{geom_label} (\code{TRUE}), or \code{geom_text} (\code{FALSE}, default)
#' 
#' @param vcs.lineend Line end type of segments of the valid cluster pairs, default is \code{"round"}, options: \code{c("round", "butt", "square")}
#'
#' @author mubu
#' @return ggplot object of the manipulated tree where valid cluster pairs are depicted
#' @method plot sgi.object
#' @export
#'
#' @examples
#' \dontrun{
#'
#' sg = sgi::sgi_init(hc, minsize, outcomes)
#' plot(sg)
#' }
plot.sgi.object <- function( sg,                                 # sg object
                             tree_opt = "full_rect",            # full triangle vcs & ivcs
                             color.valid_cluster_pairs = "black",# if null colored discreate
                             color.invalid_clusters = "gray",
                             lwds = c(0.5, 3.0),              # line weights, lwds = (lwd.ivcs, lwd.vcs)
                             tree_pruning = 0,                # no pruning(0) or prune with prune_height(1) or adaptive(2), or drop to min leaf height(3),
                             prune_height = 0,               # min length to visualize hclust
                             lty.invalid_clusters = 2,        # lty of ivcs segments ie. dashed lines

                             #---------------------- overlaying parameters
                             overlay = FALSE,          # should vcs be overlayed with rectangle
                             overlay.fill = NULL,     # filling of overlay, if NA, filled discreate
                             overlay.alpha = 0.085,   # alpha of overlay rectangles
                             overlay.color = "white", # color of rectangle lines
                             overlay.lwd =  0.1,

                             #---------------------- cluster_id parameters
                             put_cluster_ids_on = TRUE,        # should cluster ids be added on top of the tree parts
                             cluster_ids.color = "black",      # color of cluster ids
                             cluster_ids_text.size = 4,        # size of cluster ids
                             cluster_ids_as_label = FALSE,     # should they be plotted as label or text

                             #---------------------- specs
                             vcs.lineend = "round", # valid clusters line end type
                             ...
                             ){
  type = "tree"

  mm = na.omit(sg$clinds)
  vcp = paste0(mm[,"cid1"],"vs", mm[,"cid2"])

  if(type == "tree"){
    gg <- suppressWarnings( plot_tree_with_nested_clusters(
      hc = sg$hc, uvc = sg$clinds, sgs = vcp,
      tree_opt = tree_opt,
      ivcs.color = color.invalid_clusters,
      sgs.color =  color.valid_cluster_pairs,
      color_only_sgs_branches = TRUE,
      lwds = lwds,
      min_height = prune_height,
      tree_pruning = tree_pruning,
      lty.ivcs =  lty.invalid_clusters,
      overlay = overlay,
      overlay.only_sgs = TRUE,
      overlay.fill = overlay.fill,
      overlay.alpha = overlay.alpha,
      overlay.color = overlay.color,
      overlay.lwd =  overlay.lwd,
      put_cluster_ids_on = put_cluster_ids_on,
      cluster_ids_for_only_sgs = TRUE,
      cluster_ids.color = cluster_ids.color,
      cluster_ids_text.size = cluster_ids_text.size,
      cluster_ids_as_label = cluster_ids_as_label,
      vcs.lineend = vcs.lineend,
      ...
    ) + guides(lty = F) + labs(lwd = "valid cluster pairs") ) + theme(legend.position = "top")

    if(!is.null(color.valid_cluster_pairs)) gg = gg+ guides(color = F)

    # gg <- gg + theme(
    #   legend.position = c(.95, .95),
    #   legend.justification = c(1, 1)
    # )

  }else{
    # to do: add arc and rect plot
    stop("type != 'tree' is not supported yet!")
  }

  return(gg)
}

