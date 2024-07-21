#' @importFrom ggnewscale new_scale
#' @importFrom ggtext geom_textbox
#' @include fun.R gg.R
#' @noRd
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Plot Cosine Theta and Like Value
#' 
#' @description Plot cosine theta and like value of a candidate sub-cluster to determine the 
#' similarity of a candidate sub-cluster.
#' 
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param plot.scell character.\cr Candidate sub-cluster name. Default is "C0:0002".
#' @param point.size numeric.\cr Size of points. Default is 0.5.
#' @param label.subc.number integer.\cr Number of sub-cluster label. Default is 20.
#' @param label.subc.size numeric.\cr Size of sub-cluster label. Default is 1.
#' @param label.like.size numeric.\cr Size of like value. Default is 2.5.
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' 
#' @return Return a plot of cosine theta and like value of a candidate sub-cluster.
#' 
#' @export
#' 
plotThetaLike <- function(hemato.subc, plot.scell = "C0:0002", 
                          point.size = 0.5,
                          label.subc.number = 20, label.subc.size = 1, label.like.size = 2.5,
                          line.size = 0.35, text.size = 7){
  alist <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = rownames(hemato.subc@meta.data))
  theta <- alist$theta; like.mat <- alist$like.mat
  subc.ref <- hemato.subc@dist$cosine@subc.ref
  plot.info <- data.frame(
    Subc = subc.ref,
    CellType = remove.subc.suffix(subc.ref)
  )
  plot.info$Like <- as.numeric(like.mat[plot.info$Subc, plot.scell])
  plot.info$Theta <- as.numeric(theta[plot.info$Subc, plot.scell])
  label.num <- round(max(plot.info$Like), 3)
  coord.x <- plot.info$Theta[which(plot.info$Like == max(plot.info$Like))[1]]
  
  sub <- head(plot.info$Subc[order(plot.info$Like, decreasing = T)], label.subc.number)
  plot.info$Label = ""
  plot.info$Label[plot.info$Subc %in% sub] <- plot.info$Subc[plot.info$Subc %in% sub]
  
  ggplot() +
    geom_point(data = plot.info, aes_(x = ~Theta, y = ~Like, color = ~CellType), size = point.size) +
    geom_text(data = plot.info, aes_(x = ~Theta, y = ~Like, color = ~CellType, label = ~Label), size = label.subc.size) +
    geom_hline(yintercept = label.num, linetype="dashed", color = "#222222", size = 0.5 * 0.47) +
    geom_vline(xintercept = coord.x, linetype="dashed", color = "#222222", size = 0.5 * 0.47) +
    annotate("text", x = coord.x, y = label.num + diff(range(plot.info$Like))*0.01, 
             label = paste0("Max like = ", label.num), 
             size = label.like.size, vjust = 0) +
    scale_color_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% plot.info$CellType]) +
    labs(subtitle = plot.scell) +
    theme_standard(text.size = text.size, line.size = line.size) +
    theme(legend.position = "none") 
}







#' Plot Different Mean of Like Value between Normal BMMC and Candidate Group.
#' 
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param top.subc integer.\cr Number of top like sub-clusters for each candidate 
#'   sub-cluster. Default is 5.
#' @param group.name character.\cr Name of candidate group. Default is "Group".
#' @param color character.\cr Color for two groups. Default is c("blue","red").
#' @param label.size numeric.\cr Size of label. Default is 2.45.
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' 
#' @return Return a plot of different mean of like value between normal BMMC and 
#'   candidate group.
#' 
#' @importFrom stats wilcox.test
#' @export
#' 
plotMeanLike <- function(hemato.subc, top.subc = 5, group.name = "Group", 
                         color = c("blue","red"), label.size = 2.45,
                         line.size = 0.35, text.size = 7){
  dat <- data.frame(
    x = factor(rep(c("Normal", group.name), each = length(hemato.subc@dist[["cosine"]]@subc.ref)), levels = c("Normal", group.name)),
    y = c(computeMeanLike(hemato.subc = hemato.subc, cell.subc = hemato.subc@dist[["cosine"]]@subc.ref, top.subc = top.subc),
          computeMeanLike(hemato.subc = hemato.subc, cell.subc = rownames(hemato.subc@meta.data), top.subc = top.subc))
    )
  y_range <- range(dat$y, na.rm = T)
  y <- max(y_range) + diff(y_range)*0.03
  dat_label <- data.frame(x = "Normal", y = y, xend = group.name, yend = y)
  fit <- wilcox.test(y ~ x, data = dat)
  
  ggplot(data = dat, aes_(x = ~x, y = ~y, color = ~x)) +
    # geom_point(position = position_jitter(width = 0.1), alpha = 0.2, size = 0.5, show.legend = F) +
    geom_violin(size = 0.23, show.legend = F) +
    geom_boxplot(size = 0.23, outlier.colour = NA, width = 0.2, fill = "#FFFFFF", show.legend = F) +
    scale_color_manual(values = color, name = NULL) +
    geom_bracket(data = dat_label, 
                 mapping = aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend, annotation = format.pval(fit$p.value)), show.legend = F, 
                 linewidth = 0.35, label.size = label.size, color = "#000000") +
    labs(x = NULL, y = "Mean Like") +
    theme_standard(text.size = text.size, line.size = line.size)
}





#' @title Plot Circle Tree with Lineage Mapping in scRNA-seq.
#' @description Circle tree is designed to show lineage info of acute leukemia cells.
#' 
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param group.subc charater.\cr c("external", "internal", "reference"). 
#'   "external" means plot the whole sub-clusters and cell clusters of the 
#'   candidate sample. "internal" means plot the whole sub-clusters and cell 
#'   clusters of the reference in the same cosine similarity matrix. "reference"
#'   means directly plot the reference sub-clusters and cell clusters info.
#' @param plot.scell character.\cr Candidate sub-cluster name. If provided, it will
#'   only plot the like value of a sub-cluster to reference sub-cluster and cell cluster.
#' @param rm.nkt logical.\cr Weather to remove NK and T cell when plotting the like
#'   value of a sub-cluster to reference sub-cluster and cell cluster.
#' @param top.subc integer.\cr Number of top like sub-clusters for each candidate 
#'   sub-cluster.
#' @param color.mapping character.\cr c("cell.type", "cell.percentage"). "cell.type"
#'   means the color mapping of points is based on different cell types. "cell.percentage"
#'   means the color mapping of points is based on the different percentages of cell types.
#' @param size.mode character.\cr c("absolute.pct", "relative.pct"). "absolute.pct"
#'   means the size of points is based on actual percentages. "relative.pct" means
#'   the size of points is based on percentages relative to normal BMMC.
#' @param relative.compute character.\cr c("ratio", "difference"). Method for calculate 
#'   relative to normal BMMC. It will only take effect when parameter size.mode selects 
#'   "relative.pct". "ratio" means the method is cell percentages of acute leukemia 
#'   divided by that of normal BMMC. "difference" means the method is cell percentages of 
#'   acute leukemia minus that of normal BMMC.
#' @param label.population logical.\cr Whether to label population names. Default is TRUE.
#' @param label.size numeric.\cr Text size of label. Default is 3.
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param ribbon.size numeric.\cr Size of lineage ribbon. Default is 5.
#' @param point.size numeric.\cr Size of points. Default is 10.
#' @param curve.size numeric.\cr Size of curve. Default is 0.05.
#' @param title character.\cr Title of the plot. Default is NULL.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' 
#' @return Return a plot of circle tree with lineage info of acute leukemia cells.
#' 
#' @export
#' 
plotCircleTree <- function(hemato.subc, 
                           group.subc = c("external", "internal", "reference"), 
                           plot.scell = NULL, rm.nkt = T,
                           top.subc = 5,
                           color.mapping = c("cell.type", "cell.percentage"), 
                           size.mode = c("absolute.pct", "relative.pct"),
                           relative.compute = c("ratio", "difference"),
                           label.population = T, label.size = 3,
                           line.size = 0.35, ribbon.size = 5,
                           point.size = 10,
                           curve.size = 0.05,
                           title = NULL, text.size = 7){
  circle.label <- data.circle$label
  circle.node <- data.circle$node
  circle.edge <- data.circle$edge
  colnames(circle.node)[colnames(circle.node) %in% "CellName"] <- "Cell"
  circle.edge$size <- line.size
  circle.edge$size[circle.edge$Type %in% c(0, -2)] <- curve.size # outside - edge
  circle.edge$alpha <- 0.8
  
  group.subc <- match.arg(group.subc)
  color.mapping <- match.arg(color.mapping)
  size.mode <- match.arg(size.mode)
  relative.compute <- match.arg(relative.compute)
  
  if(!is.null(plot.scell)){
    if(!plot.scell %in% rownames(hemato.subc@meta.data)) stop(sprintf("cell.subc(%s) is not in hemato.subc", cell.subc[1]))
    like.mat <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = rownames(hemato.subc@meta.data))$like.mat
    like <- setNames(as.numeric(like.mat[hemato.subc@dist[["cosine"]]@subc.ref, plot.scell]), hemato.subc@dist[["cosine"]]@subc.ref)
    like <- c(like, tapply(like, remove.subc.suffix(names(like)), max))
    like <- like[match(circle.node$Cell, names(like))]
    if(rm.nkt) like[remove.subc.suffix(names(like)) %in% c("CD8 Tnaive", "CD8 Teff", "CD8 Tex", "CD8 Tdpe", "CD8 Tmpe", "CD4 Tnaive", "CD4 Tem", "CD4 Treg", "NK", "NK-XCL1", "cycling NK", "cycling NK/T")] <- 0.5
    like[is.na(like)] <- 0
    circle.node$Like <- like
    if(is.null(title)) title <- plot.scell
  }else{
    
    if(!group.subc %in% "reference"){
      cell.subc <- switch(group.subc, "external" = rownames(hemato.subc@meta.data), "internal" = hemato.subc@dist[["cosine"]]@subc.ref)
      # Summary Percentage
      like.mat <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = cell.subc)$like.mat
      ref.subc <- hemato.subc@dist[["cosine"]]@subc.ref
      subc.percentage <- computePredictSubClustePercentage(like.mat = like.mat, top.subc = top.subc, ref.subc = ref.subc)
      subc.percentage <- c(subc.percentage, tapply(subc.percentage, remove.subc.suffix(names(subc.percentage), correct.cycling = T), sum))
      subc.percentage <- subc.percentage[match(circle.node$Cell, names(subc.percentage))]
      subc.percentage[is.na(subc.percentage)] <- 0
      circle.node$Percentage <- subc.percentage
      circle.node$size <- circle.node$Percentage
      if(size.mode == "relative.pct"){
        if(relative.compute == "ratio"){
          circle.node$size <- circle.node$size / circle.node$PercentageRef
          circle.node$size[which(circle.node$size > 20)] <- 20 # too much
          circle.node$size[circle.node$Value == 1] <- circle.node$size[circle.node$Value == 1] / 10
        }else{
          circle.node$size <- circle.node$size - circle.node$PercentageRef
          if(color.mapping == "cell.type") color.mapping <- "cell.percentage"
        }
      }
      ## circle.edge - fixed line-size
      # circle.edge <- circle.edge %>%
      #   left_join(., circle.node %>% select(Root = Cell, "line.size1" = size), by = "Root")  %>%
      #   left_join(., circle.node %>% select(Cell = Cell, "line.size2" = size), by = "Cell") %>%
      #   mutate(line.size = (line.size1 + line.size2)/2, size = size + line.size2/15)
    }else{
      if(size.mode == "relative.pct"){
        if(relative.compute == "ratio"){
          circle.node$size <- 1
          circle.node$size[circle.node$Value == 1] <- circle.node$size[circle.node$Value == 1] / 10
        }else{
          cluster.node$size <- 0.1
        }
      }else{
        circle.node$size <- circle.node$PercentageRef
      }
    }
  }
  
 
  p <- ggplot() +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, color = ~CellLink), 
                 linewidth = ribbon.size, lineend = "round",
                 data = circle.edge[rev(which(!is.na(circle.edge$CellLink))),], na.rm = T) +
    geom_rect(aes(xmin = -95, xmax = 95, ymin = -95, ymax = 95), fill = "#FFFFFF", alpha = 0.85) +
    scale_colour_manual(values = color.panel.cell.population, guide = "none") +
    new_scale("color") +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size), alpha = ~I(alpha)), 
                 color = "#212121", show.legend = F, lineend = "round",
                 data = circle.edge[which(circle.edge$Type == 1), ]) +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size), alpha = ~I(alpha)), 
                 color = "#212121", show.legend = F, lineend = "round", linewidth = line.size, # fixed size
                 linetype = 2, data = circle.edge[which(circle.edge$Type == -1), ]) +
    geom_curve(mapping = aes_(xend = ~from.x, yend = ~from.y, x = ~to.x, y = ~to.y),
               linewidth = curve.size, curvature = 0.2, lineend = "round",
               inherit.aes = FALSE, show.legend = F,
               color = "#D0D0D0", data = circle.edge[which(circle.edge$Type == 0), ] ) +
    geom_curve(mapping = aes_(xend = ~from.x, yend = ~from.y, x = ~to.x, y = ~to.y),
               linewidth = curve.size, curvature = -0.1, lineend = "round",
               inherit.aes = FALSE, show.legend = F,
               color = "#D0D0D0", data = circle.edge[which(circle.edge$Type == -2), ] ) +
    new_scale("size") +
    labs(x = NULL, y = NULL, title = title) +
    theme_tree(text.size = text.size) +
    coord_fixed(clip = "off") 
  
  
  if(!is.null(plot.scell)){
    if(color.mapping %in% "cell.type"){
      p <- p + 
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~CellType), shape = 16, size = 0.5,
                   data = circle.node[circle.node$Value == 1,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~CellType), shape = 16, size = 0.5,
                   stroke = 0.23,
                   data = circle.node[circle.node$Like < 0.8 & circle.node$Value == 10,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~CellType, size = ~Like), shape = 16, 
                   stroke = 0.23,
                   data = circle.node[circle.node$Like >= 0.8 & circle.node$Value == 10,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~Like), shape = 21, 
                   stroke = 0.23, color = "#333333",
                   data = circle.node[which.max(circle.node$Like),]) +
        geom_text(mapping = aes_(x = ~pos.x, y = ~pos.y, 
                                 label = sprintf("Max like: (%s: %.3f)", circle.node[which.max(circle.node$Like),"CellType"], circle.node[which.max(circle.node$Like),"Like"])),
                  size = label.size,
                  data = circle.node[which.max(circle.node$Like),]) +
        scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType], name = "CellType",
                          guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top")) +
        scale_color_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType], name = "CellType",
                           guide = "none") +
        new_scale("color")
      
    }else{
      p <- p + 
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~Like), color = "#E5E5E5",shape = 16, size = 1,
                   data = circle.node[circle.node$Value == 1,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~Like), 
                   fill = "#E5E5E5",shape = 21, size = 0.5,
                   stroke = 0.23, color = "#333333",
                   data = circle.node[circle.node$Like < 0.8 & circle.node$Value == 10,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~Like, size = ~Like), 
                   shape = 21,
                   stroke = 0.23, color = "#333333",
                   data = circle.node[circle.node$Like >= 0.8 & circle.node$Value == 10,]) +
        scale_fill_gradient(low = "#E5E5E5", high = "#003366", limits = c(0.8, 1), breaks = c(0.8, 0.85, 0.9, 0.95, 1),
                             guide = guide_colorbar( direction = "vertical")) 
    }
    p <- p +
      scale_size(range = c(point.size/20, point.size), limits = c(0.8, 1),
                 breaks = c(0.8, 0.85, 0.9, 0.95),
                 guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), direction = "vertical")) 
  }else{
    if(color.mapping %in% "cell.type"){
      p <- p + 
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~CellType, size = ~abs(size)), shape = 16,
                   data = circle.node[which(circle.node$Value == 1),]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~abs(size)), shape = 21, 
                   stroke = 0.23, color = "#333333",
                   data = circle.node[which(circle.node$Value == 10),]) +
        scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType], name = "CellType",
                          guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top")) +
        scale_color_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType], guide = "none") +
        new_scale("color")
    }else if(color.mapping %in% "cell.percentage"){
      if(size.mode == "absolute.pct" || (size.mode == "relative.pct" & relative.compute == "ratio")){
        p <- p + 
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~size, size = ~size), 
                     shape = 21, 
                     stroke = 0.23, color = "#333333",
                     data = circle.node[circle.node$size < 20 & circle.node$Value == 10,]) +
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~size), 
                     shape = 21, fill = "#003366",
                     stroke = 0.23, color = "#333333",
                     data = circle.node[circle.node$size >= 20 & circle.node$Value == 10,]) +
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~size, size = ~size),
                     shape = 16, 
                     data = circle.node[circle.node$size < 20 & circle.node$Value == 1,]) +
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~size), 
                     shape = 16, color = "#003366",
                     data = circle.node[circle.node$size >= 20 & circle.node$Value == 1,]) +
          scale_fill_gradient(low = "#E5E5E5", high = "#003366", limits = c(0, 20), 
                              name = if(size.mode == "relative.pct") "Relative Ratio" else "Percentage (%)",
                              guide = guide_colorbar( direction = "vertical")) +
          scale_color_gradient(low = "#E5E5E5", high = "#003366", limits = c(0, 20), guide = "none") +
          new_scale("color") 
      }else{
        p <- p + 
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~size, size = ~abs(size)), 
                     shape = 21, 
                     stroke = 0.23, color = "#333333",
                     data = circle.node[circle.node$Value == 10,]) +
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~size, size = ~abs(size)),
                     shape = 16, 
                     data = circle.node[circle.node$Value == 1,]) +
          scale_fill_gradientn(colours = c("#00398e","#00398e","#00599F","#54a2dd","#a8dcff","#FFFFFF", "#ffca79","#ef6e51","#db1b18","#b50600","#b50600"),
                               values = seq.default(0, 1, by = 0.1), limits = c(-100, 100), 
                               name = "Relative (%)", guide = guide_colorbar( direction = "vertical")) +
          scale_color_gradientn(colours = c("#00398e","#00398e","#00599F","#54a2dd","#a8dcff","#FFFFFF", "#ffca79","#ef6e51","#db1b18","#b50600","#b50600"),
                                values = seq.default(0, 1, by = 0.1), limits = c(-100, 100), 
                                name = "Relative (%)", guide = "none") +
          new_scale("color")
      }
    }
    if(size.mode == "absolute.pct"){
      size.max <- max(circle.node$size, na.rm = T)
      size.min <- min(circle.node$size, na.rm = T)
      size.breaks <- c(0, 2, 5, 10, 20, 40, 80, 100)
      p <- p +  
        scale_size(range = c(point.size/20, point.size), limits = c(0, 100), name = "Percentage (%)",
                   breaks = size.breaks[1:(min(which(size.breaks >= size.max)))],
                   guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), 
                                        direction = "vertical", ncol = 2))
    }else{
      if(relative.compute == "ratio"){
        size.max <- max(circle.node$size, na.rm = T)
        size.min <- min(circle.node$size, na.rm = T)
        size.breaks <- c(0, 0.5, 1, 1.5, 2, 5, 10, 20)
        p <- p +
          scale_size(range = c(point.size/20, point.size/1.5), limits = c(0, 20), name = "Relative Ratio",
                     breaks = size.breaks[1:(min(which(size.breaks >= size.max)))],
                     guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"),
                                          direction = "vertical", ncol = 2)) 
      }else{
        value.max <- max(abs(circle.node$size), na.rm = T)
        p <- p +  
          scale_size(range = c(point.size/20, point.size/1.2), limits = c(0, 100), name = "|Relative| (%)",
                     breaks = base::pretty(c(-value.max, value.max)),
                     guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"),
                                          direction = "vertical", ncol = 2)) 
      }
    }
  }
  if(label.population){
    p <- p + 
      geom_text(aes_(x = ~x, y = ~y, label = ~Population, color = ~Population, angle = ~angle, vjust = ~vjust), 
                data = circle.label, size = label.size) +
      scale_colour_manual(values = color.panel.cell.population, guide = "none") +
      new_scale("color")
  }
  
  return(p)
}
  




#' Plot Circle Tree with Feature Expression Mapping in Sub-cluster.
#' 
#' @param feature character.\cr A feature in sub-cluster expression matrix.
#' @param scell.data matrix.\cr Sub-cluster expression matrix.
#' @param color.mapping character.\cr c("cell.type", "cell.value"). "cell.type"
#'   means the color mapping of points is based on different cell types. "cell.value"
#'   means the color mapping of points is based on the different feature expression 
#'   values of cell types.
#' @param label.population logical.\cr Whether to label population names. Default is TRUE.
#' @param title character.\cr Title of the plot. Default is NULL.
#' 
#' @export
#' 
plotCircleTreeRefFeature <- function(feature, scell.data, 
                                     color.mapping = c("cell.type", "cell.value"), 
                                     label.population = T,
                                     title = NULL){
  circle.label <- data.circle$label
  circle.node <- data.circle$node
  colnames(circle.node)[colnames(circle.node) %in% "CellName"] <- "Cell"
  circle.edge <- data.circle$edge
  circle.edge$size <- 0.4
  circle.edge$size[circle.edge$Type %in% c(0, -2)] <- 0.05 # outside - edge
  circle.edge$alpha <- 0.8
  
  color.mapping <- match.arg(color.mapping)
  # Extract feature expression
  if(length(feature) == 0) stop("Should provide 1 feature")
  if(length(feature) > 1) feature <- feature[1]
  if(!feature %in% rownames(scell.data)) stop(sprintf("Feature(%s) is not in scell.data", feature))
  value <- scell.data[feature,]
  value <- c(value, tapply(value, remove.subc.suffix(names(value)), mean))
  value <- value[match(circle.node$Cell, names(value))]
  value[is.na(value)] <- 0
  circle.node$value <- value
  
  p <- ggplot() +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, color = ~CellLink), 
                 linewidth = 6, lineend = "round",
                 data = circle.edge[rev(which(!is.na(circle.edge$CellLink))),], na.rm = T) +
    geom_rect(aes(xmin = -95, xmax = 95, ymin = -95, ymax = 95), fill = "#FFFFFF", alpha = 0.85) +
    scale_colour_manual(values = color.panel.cell.population, guide = "none") +
    new_scale("color") +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size), alpha = ~I(alpha)), 
                 color = "#212121", show.legend = F,
                 # linewidth = circle.edge$size[which(circle.edge$Type == 1)],
                 data = circle.edge[which(circle.edge$Type == 1), ]) +
    geom_segment(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size), alpha = ~I(alpha)), 
                 color = "#212121", show.legend = F,
                 linetype = 2, data = circle.edge[which(circle.edge$Type == -1), ]) +
    geom_curve(mapping = aes_(xend = ~from.x, yend = ~from.y, x = ~to.x, y = ~to.y),
               linewidth = 0.05, curvature = 0.2,
               inherit.aes = FALSE, show.legend = F,
               color = "#D0D0D0", data = circle.edge[which(circle.edge$Type == 0), ] ) +
    geom_curve(mapping = aes_(xend = ~from.x, yend = ~from.y, x = ~to.x, y = ~to.y),
               linewidth = 0.05, curvature = -0.1, 
               inherit.aes = FALSE, show.legend = F,
               color = "#D0D0D0", data = circle.edge[which(circle.edge$Type == -2), ] ) +
    new_scale("size") +
    labs(x = NULL, y = NULL, title = title) +
    theme_tree() +
    coord_fixed(clip = "off") 
  
  if(color.mapping %in% "cell.type"){
    p <- p + 
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, color = ~CellType, size = ~value/5), shape = 16,
                 data = circle.node[circle.node$Value == 1,]) +
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~value), shape = 21, 
                 stroke = 0.23, color = "#333333",
                 data = circle.node[circle.node$Value == 10,]) +
      scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType],
                        guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top")) +
      scale_color_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% circle.node$CellType], guide = "none") +
      new_scale("color") 
  }else{
    p <- p + 
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~value, size = ~value), shape = 21, 
                 stroke = 0.23, color = "#333333",
                 data = circle.node[circle.node$Value == 10,]) +
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~value/5), shape = 16, color = "#E5E5E5",
                 data = circle.node[circle.node$Value == 1,]) +
      scale_fill_gradient(low = "#E5E5E5", high = "#003366", name = feature,
                           guide = guide_colorbar( direction = "vertical")) 
  }
  p <- p +
    scale_size(range = c(0.3, 5), name = feature,
               guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), direction = "vertical")) 
  if(label.population){
    p <- p + 
      geom_text(aes_(x = ~x, y = ~y, label = ~Population, color = ~Population, angle = ~angle, vjust = ~vjust), 
                data = circle.label, size = 3) +
      scale_colour_manual(values = color.panel.cell.population, guide = "none") +
      new_scale("color")
  }
  return(p)
}




#' @title Plot Cluster Dendrogram Tree with Lineage Mapping in scRNA-seq.
#' @description Circle dendrogram tree is designed to show lineage info of acute 
#'   leukemia cells.
#' 
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param group.subc charater.\cr c("external", "internal", "reference"). 
#'   "external" means plot the whole sub-clusters and cell clusters of the 
#'   candidate sample. "internal" means plot the whole sub-clusters and cell 
#'   clusters of the reference in the same cosine similarity matrix. "reference"
#'   means directly plot the reference sub-clusters and cell clusters info.
#' @param plot.scell character.\cr Candidate sub-cluster name. If provided, it will
#'   only plot the like value of a sub-cluster to reference sub-cluster and cell cluster.
#' @param rm.nkt logical.\cr Weather to remove NK and T cell when plotting the like
#'   value of a sub-cluster to reference sub-cluster and cell cluster.
#' @param top.subc integer.\cr Number of top like sub-clusters for each candidate 
#'   sub-cluster.
#' @param color.mapping character.\cr c("cell.type", "cell.percentage"). "cell.type"
#'   means the color mapping of points is based on different cell types. "cell.percentage"
#'   means the color mapping of points is based on the different percentages of cell types.
#' @param size.mode character.\cr c("absolute.pct", "relative.pct"). "absolute.pct"
#'   means the size of points is based on actual percentages. "relative.pct" means
#'   the size of points is based on percentages relative to normal BMMC.
#' @param relative.compute character.\cr c("ratio", "difference"). Method for calculate 
#'   relative to normal BMMC. It will only take effect when parameter size.mode selects 
#'   "relative.pct". "ratio" means the method is cell percentages of acute leukemia 
#'   divided by that of normal BMMC. "difference" means the method is cell percentages of 
#'   acute leukemia minus that of normal BMMC.
#' @param label.cell logical.\cr Whether to label cell names. Default is FALSE.
#' @param label.population logical.\cr Whether to label population names. Default is TRUE.
#' @param label.pop.size logical.\cr Text size of population names. Default is 3.
#' @param line.alpha numeric.\cr Alpha of line. Default is 0.7
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param line.size.ratio numeric.\cr Ratio of line size when line size is mapping. 
#'   Default is 0.05
#' @param ribbon.size numeric.\cr Size of lineage ribbon. Default is 5.
#' @param point.size numeric.\cr Size of points. Default is 10.
#' @param title character.\cr Title of the plot. Default is NULL.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' 
#' @return Return a plot of dendrogram tree with lineage info of acute leukemia cells.
#' 
#' @export
#' 
plotClusterTree <- function(hemato.subc, 
                            group.subc = c("external", "internal", "reference"), 
                            plot.scell = NULL, rm.nkt = T, 
                            top.subc = 5,
                            color.mapping = c("cell.type", "cell.percentage"), 
                            size.mode = c("absolute.pct", "relative.pct"),
                            relative.compute = c("ratio", "difference"),
                            label.cell = F, label.population = T, label.pop.size = 3,
                            line.alpha = 0.7,
                            line.size = 0.35, line.size.ratio = 0.05,
                            ribbon.size = 5,
                            point.size = 10,
                            title = NULL, text.size = 7){
  cluster.node <- data.cluster$node
  cluster.node$size <- 1
  cluster.edge <- data.cluster$edge
  cluster.edge$size <- line.size
  cluster.label <- data.cluster$label
  
  group.subc <- match.arg(group.subc)
  color.mapping <- match.arg(color.mapping)
  size.mode <- match.arg(size.mode)
  relative.compute <- match.arg(relative.compute)
  
  if(!is.null(plot.scell)){
    if(!plot.scell %in% rownames(hemato.subc@meta.data)) stop(sprintf("cell.subc(%s) is not in hemato.subc", cell.subc[1]))
    like.mat <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = rownames(hemato.subc@meta.data))$like.mat
    like <- setNames(as.numeric(like.mat[hemato.subc@dist[["cosine"]]@subc.ref, plot.scell]), hemato.subc@dist[["cosine"]]@subc.ref)
    like <- c(like, tapply(like, remove.subc.suffix(names(like)), max))
    like <- like[match(cluster.node$Cell, names(like))]
    if(rm.nkt) like[remove.subc.suffix(names(like)) %in% c("CD8 Tnaive", "CD8 Teff", "CD8 Tex", "CD8 Tdpe", "CD8 Tmpe", "CD4 Tnaive", "CD4 Tem", "CD4 Treg", "NK", "NK-XCL1", "cycling NK", "cycling NK/T")] <- 0.5
    like[is.na(like)] <- 0
    cluster.node$Like <- like
    if(is.null(title)) title <- plot.scell
  }else if(color.mapping %in% c("cell.type", "cell.percentage")){
    if(!group.subc %in% "reference"){
      cell.subc <- switch(group.subc, "external" = rownames(hemato.subc@meta.data), "internal" = hemato.subc@dist[["cosine"]]@subc.ref)
      # Summary Percentage
      like.mat <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = cell.subc)$like.mat
      ref.subc <- hemato.subc@dist[["cosine"]]@subc.ref
      subc.percentage <- computePredictSubClustePercentage(like.mat = like.mat, top.subc = top.subc, ref.subc = ref.subc)
      subc.percentage <- c(subc.percentage, tapply(subc.percentage, remove.subc.suffix(names(subc.percentage), correct.cycling = T), sum))
      subc.percentage <- subc.percentage[match(cluster.node$Cell, names(subc.percentage))]
      subc.percentage[is.na(subc.percentage)] <- 0
      cluster.node$Percentage <- subc.percentage
      cluster.node$size <- cluster.node$Percentage
      if(size.mode == "relative.pct"){
        if(relative.compute == "ratio"){
          cluster.node$size <- cluster.node$size / cluster.node$PercentageRef
          cluster.node$size[which(cluster.node$size > 20)] <- 20 # too much
          line.size.ratio <- line.size.ratio * 3
        }else{
          cluster.node$size <- cluster.node$size - cluster.node$PercentageRef
          if(color.mapping == "cell.type") color.mapping <- "cell.percentage"
        }
      }
      cluster.edge <- cluster.edge %>%
        left_join(cluster.node %>% select(from = "CellType", line.size1 = "size"), by = "from")  %>%
        left_join(cluster.node %>% select(to = "CellType", line.size2 = "size"), by = "to") %>%
        mutate(line.size = (.data$line.size1 + .data$line.size2)/2, size = .data$size + .data$line.size2*line.size.ratio)
      cluster.edge$size[which(cluster.edge$size < line.size)] <- line.size
    }else{
      if(size.mode == "relative.pct"){
        if(relative.compute == "ratio"){
          cluster.node$size <- 1
        }else{
          cluster.node$size <- 0.1
        }
      }else{
        cluster.node$size <- cluster.node$PercentageRef
      }
    }
  }
  
  p <- ggplot() +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, color = ~CellLink), 
                   size = ribbon.size, lineend = "round",
                   data = cluster.edge[which(!is.na(cluster.edge$CellLink)),], na.rm = T) +
    geom_rect(aes(xmin = 0.5, xmax = 12.5, ymin = 0.5, ymax = 9.5), fill = "#FFFFFF", alpha = 0.9) +
    scale_colour_manual(values = color.panel.cell.population, guide = "none") +
    new_scale("color") +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size)), 
                   color = "#212121", lineend = "round", alpha = line.alpha,
                   data = cluster.edge[cluster.edge$value == 1,], linetype = 1) +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size)), 
                   color = "#212121", lineend = "round", alpha = line.alpha, size = 0.4,
                   data = cluster.edge[cluster.edge$value == -1,], linetype = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    new_scale("size") +
    labs(x = NULL, y = NULL, title = title) +
    theme_tree(text.size = text.size)
  
  if(!is.null(plot.scell)){
    if(color.mapping %in% "cell.type"){
      p <- p +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType), 
                   shape = 21, size = 0.5, stroke = 0.23, color = "#333333",
                   data = cluster.node[cluster.node$Like <= 0.8,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~Like), 
                   shape = 21, stroke = 0.23, color = "#333333",
                   data = cluster.node[cluster.node$Like > 0.8,]) +
        scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% cluster.node$CellType], name = "CellType",
                          guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top"))
    }else{
      p <- p + 
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~Like, size = ~Like),
                   shape = 21, stroke = 0.23, color = "#333333", fill = "E5E5E5", size = 1,
                   data = cluster.node[cluster.node$Like < 80,]) +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~Like), shape = 21,
                   stroke = 0.23, color = "#333333",
                   data = cluster.node[cluster.node$Like >= 80,]) +
        scale_fill_gradient(low = "#E5E5E5", high = "#003366", limits = c(0.8, 1), breaks = c(0.8, 0.85, 0.9, 0.95, 1),
                            guide = guide_colorbar( direction = "vertical"), name = "Like") 
    }
    p <- p +
      scale_size(range = c(0.5, 7), limits = c(0.8, 1),
                 breaks = c(0.8, 0.85, 0.9, 0.95),
                 guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), direction = "vertical")) 
  }else if(color.mapping %in% c("cell.type", "cell.percentage")){
    if(color.mapping %in% "cell.type"){
      p <- p +
        geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~abs(size)), 
                   shape = 21, stroke = 0.23, color = "#333333",
                   data = cluster.node) +
        scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% cluster.node$CellType],
                          guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top"))
    }else{
      if(size.mode == "absolute.pct" || (size.mode == "relative.pct" & relative.compute == "ratio")){
        p <- p + 
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~size, size = ~size),
                     shape = 21, stroke = 0.23, color = "#333333",
                     data = cluster.node[cluster.node$size < 20,]) +
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~size), 
                     shape = 16, color = "#003366",
                     data = cluster.node[cluster.node$size >= 20,]) +
          scale_fill_gradient(low = "#E5E5E5", high = "#003366", limits = c(0, 20), 
                              name = if(size.mode == "relative.pct") "Relative Ratio" else "Percentage (%)", 
                              guide = guide_colorbar( direction = "vertical")) 
      }else{
        p <- p + 
          geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~size, size = ~abs(size)), 
                     shape = 21, stroke = 0.23, color = "#333333",
                     data = cluster.node) +
          scale_fill_gradientn(colours = c("#00398e","#00398e","#00599F","#54a2dd","#a8dcff","#FFFFFF", "#ffca79","#ef6e51","#db1b18","#b50600","#b50600"),
                               values = seq.default(0, 1, by = 0.1), limits = c(-100, 100), 
                               name = "Relative (%)", guide = guide_colorbar( direction = "vertical"))
      }
    }
    if(size.mode == "absolute.pct"){
      size.max <- max(cluster.node$size, na.rm = T)
      size.min <- min(cluster.node$size, na.rm = T)
      size.breaks <- c(0, 2, 5, 10, 20, 40, 80, 100)
      p <- p +
        scale_size(range = c(point.size/10, point.size), limits = c(0, 100), name = "Percentage (%)",
                   breaks = size.breaks[1:(min(which(size.breaks >= size.max)))],
                   guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), 
                                        direction = "vertical", ncol = 2)) 
    }else if(size.mode == "relative.pct"){
      if(relative.compute == "ratio"){
        size.max <- max(cluster.node$size, na.rm = T)
        size.min <- min(cluster.node$size, na.rm = T)
        size.breaks <- c(0, 0.5, 1, 1.5, 2, 5, 10, 20)
        p <- p +
          scale_size(range = c(point.size/20, point.size), limits = c(0, 20), name = "Relative Ratio",
                     breaks = size.breaks[1:(min(which(size.breaks >= size.max)))],
                     guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"),
                                          direction = "vertical", ncol = 2))
      }else{
        value.max <- max(abs(cluster.node$size), na.rm = T)
        p <- p +  
          scale_size(range = c(point.size/10, point.size), limits = c(0, 100), name = "|Relative| (%)",
                     breaks = base::pretty(c(-value.max, value.max)),
                     guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"),
                                          direction = "vertical", ncol = 2)) 
      }
    }
  }
  
  if(label.cell){
    p <- p + 
      geom_textbox(data = cluster.node, aes_(x = ~pos.x, y = ~pos.y, label = ~gsub("-", "- ", CellType)),
                   maxwidth = unit(0.7, "cm"), size = 2.1,
                   halign = 0.5, fill = NA,
                   box.colour = NA,
                   box.padding = unit(c(0,0,0,0), "pt"), box.margin = unit(c(0,0,0,0), "pt"))
  }
  
  if(label.population){
    p <- p + 
      geom_text(aes_(x = ~x, y = ~y, label = ~Population, color = ~Population, angle = ~angle, vjust = ~vjust, hjust = ~hjust), 
                data = cluster.label, size = label.pop.size) +
      scale_colour_manual(values = color.panel.cell.population, guide = "none") +
      new_scale("color")
  }
  
  return(p)
}








#' Plot Lasso Dendrogram Tree with Lineage Mapping of Acute Leukemia Cells in RNA-seq.
#' 
#' @param mat.score matrix.\cr Lasso score matrix.
#' @param color.mapping character.\cr c("value", "cell.type"). "value" means 
#'   the color mapping of points is based on lasso score value. "cell.type" means 
#'   the color mapping of points is based on different cell types.
#' @param label.cell logical.\cr Whether to label cell names. Default is FALSE.
#' @param label.population logical.\cr Whether to label population names. Default is TRUE.
#' @param label.pop.size logical.\cr Text size of population names. Default is 3.
#' @param line.alpha numeric.\cr Alpha of line. Default is 0.7
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param line.size.ratio numeric.\cr Ratio of line size when line size is mapping. 
#'   Default is 1.5.
#' @param ribbon.size numeric.\cr Size of lineage ribbon. Default is 5.
#' @param point.size numeric.\cr Size of points. Default is 1.
#' @param point.color.limits numeric.\cr Limits of point color. Default is c(-2.2, 2.2)
#' @param title character.\cr Title of the plot. Default is NULL.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' @param legend.barheight unit.\cr Bar height of legend. Default is unit(2.5, "mm").
#' 
#' @return Return a plot of dendrogram tree with lasso score of acute leukemia cells.
#' 
#' @export
#' 
plotLassoTree <- function(mat.score,
                          color.mapping = c("value", "cell.type"), 
                          label.cell = F, label.population = T, label.pop.size = 3,
                          line.alpha = 0.7,
                          line.size = 0.35, line.size.ratio = 1.5,
                          ribbon.size = 5,
                          point.size = 1, point.color.limits = c(-2.2, 2.2),
                          title = NULL, text.size = 7,
                          legend.barheight = unit(2.5, "mm")){
  cluster.node <- data.cluster$node
  cluster.node$size <- 1
  cluster.edge <- data.cluster$edge
  cluster.edge$size <- line.size
  cluster.label <- data.cluster$label
  
  color.mapping <- match.arg(color.mapping)
  
  if(!is.matrix(mat.score)) stop("mat.score should be matrix")
  if(nrow(mat.score) > 1){
    mat.score <- t(as.matrix(colMeans(mat.score)))
  }
  score <- setNames(as.vector(mat.score), colnames(mat.score))
  # align
  score <- score[match(cluster.node$Cell, names(score))]
  score[is.na(score)] <- 0
  cluster.node$score <- score
  cluster.edge <- cluster.edge %>%
    left_join(cluster.node %>% select(from = "CellType", line.size1 = "score"), by = "from")  %>%
    left_join(cluster.node %>% select(to = "CellType", line.size2 = "score"), by = "to") %>%
    mutate(line.size = (.data$line.size1) + (.data$line.size2)/2, 
           size = .data$size + .data$line.size2*line.size.ratio)
  cluster.edge$size[cluster.edge$size < line.size] <- line.size
  cluster.node$size <- case_when(cluster.node$score <= -0.5 ~ "< -0.50",
                                 cluster.node$score > -0.5 & cluster.node$score <= 0~ "-0.50 ~ 0.00",
                                 cluster.node$score > 0 & cluster.node$score <= 0.32~ "0.00 ~ 0.32",
                                 cluster.node$score > 0.32 & cluster.node$score <= 0.5~ "0.32 ~ 0.50",
                                 cluster.node$score > 0.5 & cluster.node$score <= 0.75~ "0.50 ~ 0.75",
                                 cluster.node$score > 0.75 & cluster.node$score <= 1~ "0.75 ~ 1.00",
                                 cluster.node$score > 1 & cluster.node$score <= 1.25~ "1.00 ~ 1.25",
                                 cluster.node$score > 1.25 & cluster.node$score <= 1.5~ "1.25 ~ 1.50",
                                 cluster.node$score > 1.50 ~ "> 1.50")
  size_value <- c("< -0.50" = 1, "-0.50 ~ 0.00" = 1.5, "0.00 ~ 0.32" = 2, "0.32 ~ 0.50" = 2.25, "0.50 ~ 0.75" = 2.5, "0.75 ~ 1.00" = 2.75, 
                  "1.00 ~ 1.25" = 3, "1.25 ~ 1.50" = 3.25,  "> 1.50" = 3.5)
  cluster.node$size <- factor(cluster.node$size, levels = names(size_value))
  
  
  p <- ggplot() +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, color = ~CellLink), 
                   linewidth = ribbon.size, lineend = "round",
                   data = cluster.edge[which(!is.na(cluster.edge$CellLink)),], na.rm = T) +
    geom_rect(aes(xmin = 0.5, xmax = 12.5, ymin = 0.5, ymax = 9.5), fill = "#FFFFFF", alpha = 0.9) +
    scale_colour_manual(values = color.panel.cell.population, guide = "none") +
    new_scale("color") +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size)), 
                   color = "#212121", lineend = "round", alpha = line.alpha,
                   data = cluster.edge[cluster.edge$value == 1,], linetype = 1) +
    geom_stairstep(mapping = aes_(x = ~from.x, y = ~from.y, xend = ~to.x, yend = ~to.y, size = ~I(size)), 
                   color = "#212121", lineend = "round", alpha = line.alpha, linewidth = 0.4,
                   data = cluster.edge[cluster.edge$value == -1,], linetype = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    new_scale("size") +
    labs(x = NULL, y = NULL, title = title) +
    theme_tree(text.size = text.size)
  
  if(color.mapping %in% "cell.type"){
    p <- p +
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, fill = ~CellType, size = ~size), 
                 shape = 21, stroke = 0.23, color = "#333333",
                 data = cluster.node, na.rm = T) +
      scale_fill_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% cluster.node$CellType], name = "CellType",
                        guide = guide_legend(order = 1, override.aes = list(size = 2.5), title.position = "top"))
  }else{
    p <- p + 
      geom_point(mapping = aes_(x = ~pos.x, y = ~pos.y, size = ~size, fill = ~score), shape = 21,
                 stroke = 0.23, color = "#333333",
                 data = cluster.node, na.rm = T) +
      scale_fill_gradient2(low = "#00599F", high = "#db1b18", limits = point.color.limits,
                           guide = guide_colorbar( direction = "horizontal", title.position = "top",
                                                   barheight = legend.barheight), name = "Normalized Lasso Score")
  }
  
  p <- p +
    scale_size_manual(values = size_value*point.size, name = "Normalized Lasso Score",
               guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#CCCCCC", color = "#000000"), ncol = 2, direction = "vertical")) 
  
  if(label.cell){
    p <- p + 
      geom_textbox(data = cluster.node, aes_(x = ~pos.x, y = ~pos.y, label = ~gsub("-", "- ", CellType)),
                   maxwidth = unit(0.7, "cm"), size = 2.1,
                   halign = 0.5, fill = NA,
                   box.colour = NA,
                   box.padding = unit(c(0,0,0,0), "pt"), box.margin = unit(c(0,0,0,0), "pt"))
  }
  
  if(label.population){
    p <- p + 
      geom_text(aes_(x = ~x, y = ~y, label = ~Population, color = ~Population, angle = ~angle, vjust = ~vjust, hjust = ~hjust), 
                data = cluster.label, size = label.pop.size) +
      scale_colour_manual(values = color.panel.cell.population, guide = "none") +
      new_scale("color")
  }
  return(p)
}


#' Plot Lollipop with Lasso Score of Acute Leukemia Cells in RNA-seq.
#' 
#' @param mat.score matrix.\cr Lasso score matrix.
#' @param point.size numeric.\cr Size of points. Default is 3.
#' @param line.size numeric.\cr Size of line. Default is 0.35.
#' @param text.size numeric.\cr Size of text. Default is 7.
#' 
#' @return Return a plot of lollipop with lasso score of acute leukemia cells.
#' 
#' @export
#' 
plotLassoLollipop <- function(mat.score, point.size = 3, line.size = 0.35, text.size = 7){
  
  if(!is.matrix(mat.score)) stop("mat.score should be matrix")
  mat.score <- colMeans(mat.score)
  plot.data <- data.frame(Cell = factor(names(mat.score), levels = rev(names(mat.score))), value = mat.score)
  
  
  ggplot() +
    geom_vline(xintercept = 0, size = line.size, color = "#000000") +
    geom_linerange(data = plot.data, aes_(y = ~Cell, xmin = 0, xmax = ~value), linewidth = line.size) +
    geom_point(data = plot.data, aes_(x = ~value, y = ~Cell, color = ~Cell, size = ~abs(value)), 
               stat = "identity", show.legend = F) +
    geom_vline(xintercept = c(-0.32, 0.32), linetype = "dashed", linewidth = line.size, color = "#555555") +
    scale_color_manual(values = color.panel.cell.type[names(color.panel.cell.type) %in% plot.data$Cell], name = NULL) +
    scale_size(range = c(0.01, point.size)) +
    labs(y = NULL, x = "Normalized Lasso Score") +
    theme_standard(text.size = text.size, line.size = line.size)
  
}






