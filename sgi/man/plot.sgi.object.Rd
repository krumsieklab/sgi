% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sgi_object.R
\name{plot.sgi.object}
\alias{plot.sgi.object}
\title{Plot \code{SGI} object}
\usage{
\method{plot}{sgi.object}(
  sg,
  tree_opt = "full_rect",
  color.valid_cluster_pairs = "black",
  color.invalid_clusters = "gray",
  lwds = c(0.5, 3),
  tree_pruning = 0,
  prune_height = 0,
  lty.invalid_clusters = 2,
  overlay = FALSE,
  overlay.fill = NULL,
  overlay.alpha = 0.085,
  overlay.color = "white",
  overlay.lwd = 0.1,
  put_cluster_ids_on = TRUE,
  cluster_ids.color = "black",
  cluster_ids_text.size = 4,
  cluster_ids_as_label = FALSE,
  vcs.lineend = "round",
  ...
)
}
\arguments{
\item{sg}{\code{SGI} initialization object created by \code{sgi::sgi()}}

\item{tree_opt}{Shape of the tree branches: \cr
\code{"full_rect"} Rectangles (default) \cr
\code{"up_tri_down_rect"} Valid cluster pairs as triangles, rest as rectangles \cr
\code{"up_rect_down_tri"} Valid cluster pairs as rectangles, rest as triangles \cr
\code{"full_tri"} Triangles}

\item{color.valid_cluster_pairs}{Color of valid cluster pairs, default is black}

\item{color.invalid_clusters}{Color of invalid clusters, default is gray}

\item{lwds}{Line widths, vector of length two containing widths of invalid cluster lines and valid cluster lines}

\item{tree_pruning}{Options for tree pruning (i.e. removing bottom branches in the visualization): \cr
\code{0}: no pruning \cr
\code{1}: prune with \code{prune_height} (separate function argument) \cr
\code{2}: adaptive pruning, cut underneath lowest valid cluster pair  \cr
\code{3}: adaptive pruning, cut at highest leaf branching point (cuts away long leaf branches), pruned part is left blank \cr
\code{4}: same as option \code{3}, but pruned part is cut out.}

\item{prune_height}{Pruning height if \code{tree_pruning==1}}

\item{lty.invalid_clusters}{Line type for invalid clusters, see linetypes options in R}

\item{overlay}{Whether valid cluster pairs are overlayed with rectangles, \code{FALSE} (default)}

\item{overlay.fill}{Single fill color of overlay rectangles. \code{NULL} (default) to fill each cluster pair differently.}

\item{overlay.alpha}{Alpha of overlay rectangles}

\item{overlay.color}{Border color of overlay rectangles}

\item{overlay.lwd}{Line width of overlay rectangles}

\item{cluster_ids.color}{Color of cluster ID text}

\item{cluster_ids_text.size}{Size of cluster ID text}

\item{cluster_ids_as_label}{Should IDs be printed with \code{geom_label} (\code{TRUE}), or \code{geom_text} (\code{FALSE}, default)}

\item{vcs.lineend}{Line end type of segments of the valid cluster pairs, default is \code{"round"}, options: \code{c("round", "butt", "square")}}

\item{draw_cluster_ids}{Whether tree is annotated with cluster IDs, \code{TRUE} {default}}
}
\value{
ggplot object of the manipulated tree where valid cluster pairs are depicted
}
\description{
This function provides detailed visualization options for \code{SGI} objects before association analysis.
Plotting an \code{SGI} initialization object before performing the association analysis can be used to guide the statistical design.
}
\examples{
\dontrun{

sg = sgi::sgi_init(hc, minsize, outcomes)
plot(sg)
}
}
\author{
mubu
}
