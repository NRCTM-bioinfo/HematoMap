#' @import dplyr
#' @import ggplot2
#' @importFrom graphics text
#' @noRd
NULL



GeomStairstep <- ggproto(
  'GeomStairstep', ggplot2::GeomSegment,
  # draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
  #                       lineend = "butt", linejoin = "round", na.rm = FALSE){
  #   data <- ggplot2::remove_missing(data, na.rm = na.rm, c(ggplot2:::GeomSegment$required_aes),
  #                                   name = "geom_segment")
  #   data <- self$stairstep(data)
  #   ggplot2:::GeomSegment$draw_panel(data, panel_params, coord, arrow, arrow.fill, lineend, linejoin, na.rm = na.rm)
  # },
  # stairstep = function(data){
  #   data <- data[order(data$y, data$yend, data$group, decreasing = T), ]
  #   n <- nrow(data)
  #   xs <- as.vector(aperm(matrix(c(1:(2*n)), ncol = 2)))
  #   ys <- rep(1:n, each = 2)
  #   data.frame(x = c(data$x, data$xend)[xs], 
  #              y = c(data$y, data$y)[xs],
  #              xend = c(data$xend, data$xend)[xs],
  #              yend = c(data$y, data$yend)[xs],
  #              data[ys, setdiff(names(data), c("x", "xend", "y", "yend"))])
  # }
  draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                        lineend = "butt", linejoin = "round", na.rm = FALSE){
    data <- remove_missing(data, na.rm = na.rm, c(ggplot2:::GeomSegment$required_aes),
                           name = "geom_segment")
    data <- self$stairstep(data)
    ggplot2:::GeomPath$draw_panel(data, panel_params, coord, arrow, lineend, linejoin, na.rm = na.rm)
  },
  stairstep = function(data){
    data <- data[order(data$y, data$yend, decreasing = T), ]
    data$group <- 1:nrow(data)
    rbind(
      data[, setdiff(colnames(data), c("xend", "yend"))],
      data[, setdiff(colnames(data), c("x", "yend"))] %>% rename(x = xend),
      data[, setdiff(colnames(data), c("x", "y"))] %>% rename(x = xend, y = yend)
    ) %>% arrange(group) %>% distinct()
  }
)


geom_stairstep <- function(mapping = NULL, data = NULL, stat = "identity",
                           position = "identity", show.legend = NA, inherit.aes = T, na.rm = FALSE, ...) {
  layer(mapping = mapping, data = data, stat = stat,
        geom = GeomStairstep, position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}


#' Add Brackets with Labels to a GGPlot
#' 
#' @description add brackets with label annotation to a ggplot. Helpers for
#'   adding p-value or significance levels to a plot. **It is based on `ggpubr::geom_bracket`**
#'  
#' @param mapping Set of aesthetic mappings created by [aes()]. If specified and
#'   `inherit.aes = TRUE` (the default), it is combined with the default mapping
#'   at the top level of the plot. You must supply `mapping` if there is no plot
#'   mapping.
#' @param data The data to be displayed in this layer. There are three
#'    options:
#'
#'    If `NULL`, the default, the data is inherited from the plot
#'    data as specified in the call to [ggplot()].
#'
#'    A `data.frame`, or other object, will override the plot
#'    data. All objects will be fortified to produce a data frame. See
#'    [fortify()] for which variables will be created.
#'
#'    A `function` will be called with a single argument,
#'    the plot data. The return value must be a `data.frame`, and
#'    will be used as the layer data. A `function` can be created
#'    from a `formula` (e.g. `~ head(.x, 10)`).
#' @param stat The statistical transformation to use on the data for this
#'    layer, either as a `ggproto` `Geom` subclass or as a string naming the
#'    stat stripped of the `stat_` prefix
#' @param position Position adjustment, either as a string naming the adjustment, 
#'   or the result of a call to a position adjustment function. Use the latter 
#'   if you need to change the settings of the adjustment.
#' @param show.legend logical. Should this layer be included in the legends?
#'   `NA`, the default, includes if any aesthetics are mapped.
#'   `FALSE` never includes, and `TRUE` always includes.
#'   It can also be a named logical vector to finely select the aesthetics to
#'   display.
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics,
#'   rather than combining with them. This is most useful for helper functions
#'   that define both data and aesthetics and shouldn't inherit behaviour from
#'   the default plot specification.
#' @param label Character vector with alternative label, if not null test is
#'   ignored.
#' @param xmin Numeric vector with the positions of the left sides of the
#'   brackets.
#' @param xmax Numeric vector with the positions of the right sides of the
#'   brackets.
#' @param size The width of the lines of the bracket
#' @param label.size The size of the label text
#' @param family The font used for the text
#' @param vjust Move the text up or down relative to the bracket
#' @param coord.flip If \code{TRUE}, flip x and y coordinates so that
#'   horizontal becomes vertical, and vertical, horizontal.
#' @param ... other arguments passed on to \code{\link{layer}}.
#' 
#' @export
geom_bracket <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE,
                         label = NULL, xmin = NULL, xmax = NULL,
                         size = 0.3, label.size = 2.5, family="", vjust = 0,
                         coord.flip = FALSE,
                         ...) {

  ggplot2::layer(
    stat = stat, geom = GeomBracket, mapping = mapping,  data = data,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      size = size, label.size = label.size,
      family = family, na.rm = na.rm, coord.flip = coord.flip,
      ...
    )
  )
}


CalculateGroupPath <- function(df){
  path <- droplevels(df[, 1])
  theGroupName <- colnames(df)[1]
  nPathPoints <- ncol(df) - 1
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi)/nPathPoints)
  nDataPoints <- ncol(df) * length(levels(path))
  graphData <- data.frame(seg = rep("", nDataPoints), x = rep(0, nDataPoints), y = rep(0, nDataPoints))
  colnames(graphData)[1] <- theGroupName
  rowNum <- 1
  for (i in 1:length(levels(path))) {
    pathData <- subset(df, df[, 1] == levels(path)[i])
    for (j in c(2:ncol(df))) {
      graphData[rowNum, theGroupName] <- levels(path)[i]
      graphData$x[rowNum] <- pathData[, j] * sin(angles[j - 1])
      graphData$y[rowNum] <- pathData[, j] * cos(angles[j - 1])
      rowNum <- rowNum + 1
    }
    graphData[rowNum, theGroupName] <- levels(path)[i]
    graphData$x[rowNum] <- pathData[, 2] * sin(angles[1])
    graphData$y[rowNum] <- pathData[, 2] * cos(angles[1])
    rowNum <- rowNum + 1
  }
  graphData[, 1] <- factor(graphData[, 1], levels = levels(path))
  graphData
}

CalculateAxisPath <- function (var.names, min, max){
  n.vars <- length(var.names)
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi)/n.vars)
  min.x <- min * sin(angles)
  min.y <- min * cos(angles)
  max.x <- max * sin(angles)
  max.y <- max * cos(angles)
  axisData <- list()
  for (i in 1:n.vars) {
    a <- c(i, min.x[i], min.y[i])
    b <- c(i, max.x[i], max.y[i])
    axisData[[i]] <- rbind(a, b)
  }
  axisData <- do.call(rbind, axisData)
  colnames(axisData) <- c("axis.no", "x", "y")
  rownames(axisData) <- seq(1:nrow(axisData))
  as.data.frame(axisData)
}


funcCircleCoords <- function(center = c(0, 0), r = 1, npoints = 100){
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

ggradar <- function(dat, axis.x.expand = 0.05,
                    axis.label.offset = 1.1,
                    axis.label.size = 2.5,
                    point.size = 1,
                    colours =c("#AAAAAA","#F6A7DF","#4076EB")){
  
  dat.offset <- dat <- as.data.frame(dat)
  y <- as.numeric(as.matrix(dat[, 2:ncol(dat)]))
  y.range <- range(y, na.rm = T)
  y.breaks <- pretty(y.range)
  if(y.range[1] > min(y.breaks)) y.range[1] <- min(y.breaks)
  if(y.range[2] < max(y.breaks)) y.range[2] <- max(y.breaks)
  y.center = min(y.range) - (1/9) * diff(y.range)
  x.center = 0.02 * (max(y.range) - y.center)
  
  var.names <- colnames(dat.offset)[-1]
  dat.offset[, 2:ncol(dat.offset)] <- dat.offset[, 2:ncol(dat.offset)] - y.center
  group <- list()
  group$path <- CalculateGroupPath(dat.offset)
  
  # outer - label
  axis <- NULL
  axis$path <- CalculateAxisPath(var.names, min(y.range) - y.center, max(y.range) - y.center)
  axis$label <- data.frame(text = var.names, x = NA, y = NA)
  n.vars <- length(var.names)
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi)/n.vars)
  axis$label$x <- sapply(1:n.vars, function(i, x)
    ((max(y.range) - y.center) * axis.label.offset) * sin(angles[i])
  )
  axis$label$y <- sapply(1:n.vars, function(i, x)
    ((max(y.range) - y.center) * axis.label.offset) * cos(angles[i])
  )
  axis$label$hjust <- 0.5
  axis$label$vjust <- 0.5
  axis$label$hjust[axis$label$x < (-x.center)] <- 1
  axis$label$vjust[abs(axis$label$x) <= x.center & axis$label$y > 0] <- 0
  axis$label$hjust[axis$label$x > x.center] <- 0
  axis$label$color <- color.panel.cell.type[match(axis$label$text, names(color.panel.cell.type))]
  
  
  # data - circle - grid
  axis_circ <- function(axis.breaks){
    pathList <- labelList <- list()
    for(i in 1:length(axis.breaks)){
      pathList[[i]] <- funcCircleCoords(c(0, 0), axis.breaks[i] - y.center, npoints = 360)
      pathList[[i]]$group = as.character(axis.breaks[i])
      labelList[[i]] <- data.frame(x = -0.05 * (max(y.range) - y.center), 
                                   y = axis.breaks[i] - y.center, 
                                   text = as.character(axis.breaks[i]))
    }
    list(path = do.call(rbind, pathList), label = do.call(rbind, labelList))
  }
  gridline <- axis_circ(axis.breaks = y.breaks)
  
  
  p <- ggplot() + 
    geom_text(data = axis$label, 
              aes_(x = ~x, y = ~y, label = ~text, hjust = ~hjust, vjust = ~vjust),
              size = axis.label.size, color = axis$label$color) +
    geom_path(data = gridline$path, aes_(x = ~x, y = ~y, group = ~group),
              colour = "grey90", linewidth = 0.5*0.47, linetype = 1) +
    geom_path(data = axis$path, aes_(x = ~x, y = ~y, group = ~axis.no), 
              colour = "grey90", linewidth = 0.5*0.47, linetype = 1)
  
  p <- p + geom_path(data = group$path, aes_(x = ~x, y = ~y, group = group$path[,1], colour = group$path[,1]), 
                     linewidth =0.5*0.47) +
    geom_point(data = group$path, aes_(x = ~x, y = ~y, fill = group$path[,1], colour = group$path[,1]), 
               shape = 21, size = point.size) +
    geom_polygon(data = group$path, aes_(x = ~x, y = ~y, fill = group$path[,1]), alpha = 0.1, show.legend = F) +
    scale_colour_manual(values = colours, name = colnames(group$path)[1]) +
    scale_fill_manual(values = colours, guide = "none")
  
  
  # axis-tick-label
  p <- p + geom_text(aes_(x = ~x, y = ~y, label = ~text), data = gridline$label, 
                     size = axis.label.size-0.4, hjust = 1)
  p <- p + 
    coord_equal(clip = "off") + 
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
    theme(legend.position = "right",
          legend.key.height = unit(0.25,'cm'), 
          legend.key.width = unit(0.35,'cm'),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.spacing.x = unit(0.1,'cm'),
          legend.spacing.y =  unit(0.05,'cm'),
          legend.margin=margin(r = 0.1, t = 0, b = 0, l = 0.1, unit='cm'),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
  p <- p + 
    scale_x_continuous(expand = expansion(c(axis.x.expand, axis.x.expand)))
  return(p)
}


GeomBracket <- ggplot2::ggproto("GeomBracket", ggplot2::Geom,
                                required_aes = c("x", "xend", "y", "yend", "annotation"),
                                default_aes = ggplot2::aes(
                                  shape = 19, colour = "black", label.size = 3.88, angle = NULL, hjust = 0.5,
                                  vjust = 0, alpha = NA, family = "", fontface = 1, lineheight = 1.2, linetype=1, size = 0.3,
                                  xmin = NULL, xmax = NULL, label = NULL,
                                ),
                                draw_key = ggplot2::draw_key_path,
                                draw_group = function(data, panel_params, coord, type = "text",
                                                      coord.flip = FALSE) {
                                  lab <- as.character(data$annotation)
                                  coords <- coord$transform(data, panel_params)
                                  label.x <- mean(c(coords$x[1], tail(coords$xend, n=1)))
                                  label.y <- max(c(coords$y, coords$yend))+0.01
                                  label.angle <- coords$angle
                                  if(coord.flip){
                                    label.y <- mean(c(coords$y[1], tail(coords$yend, n=1)))
                                    label.x <- max(c(coords$x, coords$xend))+0.01
                                    if(is.null(label.angle)) label.angle <- -90
                                  }
                                  if(is.null(label.angle)) label.angle <- 0
                                  grid::gList(
                                    grid::textGrob(
                                      label = lab,
                                      x = label.x,
                                      y = label.y,
                                      default.units = "native",
                                      hjust = coords$hjust, vjust = coords$vjust,
                                      rot = label.angle,
                                      gp = grid::gpar(
                                        col = scales::alpha(coords$colour, coords$alpha),
                                        fontsize = coords$label.size * ggplot2::.pt,
                                        fontfamily = coords$family,
                                        fontface = coords$fontface,
                                        lineheight = coords$lineheight
                                      )
                                    ),
                                    grid::segmentsGrob(
                                      coords$x, coords$y,
                                      default.units = "native",
                                      coords$xend, coords$yend,
                                      gp = grid::gpar(
                                        col = scales::alpha(coords$colour, coords$alpha),
                                        lty = coords$linetype,
                                        lwd = coords$size * ggplot2::.pt
                                      )
                                    )
                                  )
                                }
)



theme_tree <- function(text.size = 7){
  theme(plot.title = element_text(size = text.size, hjust = 0.5), 
        panel.background = element_blank(), panel.grid = element_blank()) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(), 
          axis.ticks.length = unit(0, "pt")) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.25,'cm'), 
          legend.key.width = unit(0.35,'cm'),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.spacing.x = unit(0.1,'cm'),
          legend.spacing.y =  unit(0.05,'cm'),
          legend.margin=margin(r = 0.1, t = 0, b = 0, l = 0.1, unit='cm'),
          legend.background = element_blank(),
          legend.text = element_text(size = text.size-1),
          legend.title = element_text(size = text.size-1)) 
}


theme_standard <- function(text.size = 7, line.size = 0.35){
  theme(plot.subtitle = element_text(size = text.size, hjust = 0.5, vjust = 0,
                                     margin = margin(b = 0.1, t = -0.04, l = 0, r = 0, unit = "cm"))) +
  theme(axis.text = element_text(size = text.size-1), title = element_text(size = text.size)) + 
  theme(axis.line = element_line(linewidth = line.size), axis.ticks = element_line(linewidth = line.size)) +
  theme(legend.text = element_text(size = text.size-1), legend.key.size = unit(3, "mm")) +
  theme(strip.text = element_text(size = text.size),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
        panel.spacing=unit(0.2, "cm")) +
  theme(axis.line = element_line(linewidth = line.size), axis.ticks = element_line(linewidth = line.size)) +
  theme(panel.background = element_blank(), panel.grid = element_blank())
}



color.panel.cell.type <- c(
  "HSC/MPP" = "#e31a00", "LMPP" = "#ffd35f", "CLP" = "#8BCE9D",
  "CMP" = "#FDDBB8",  "MDP" = "#F9A62C","GMP" = "#F6A7DF",
  "CDP" = "#DDA860",
  "MEP" = "#D67B79",
  "pro-Mono" = "#9340C4",
  "CD14 Mono" = "#ACBBE4",
  "CD16 Mono" = "#7A6FD0",
  "pre-DC" = "#E08563",
  "pDC" = "#F8D371", "cDC1" = "#CF5A00", "mo-DC" = "#F0DEC7",
  "MKP" = "#A281C1", "MK" = "#D6C8D1",
  "pro-Ery1" = "#F6999A", "pro-Ery2" = "#E69920", "Ery" = "#F8C5C5",
  "pre-pro-B" = "#BBE7F9",
  "Early pro-B" = "#ADC1F2",
  "Early cycling pro-B" = "#8790E7",
  "Late pro-B" = "#78D3F8",
  "Late cycling pro-B" = "#5F91F3",
  "pre-B" = "#5C7AEA",
  "Immature B" = "#7BA6D6",
  "Naive B" = "#3174F3",
  "Memory B 1" = "#272EB3",
  "Memory B 2" = "#802385",
  "CD8 Tnaive" = "#53B871", "CD8 Teff" = "#D2DFC7", "CD8 Tex" = "#B6D6E3", "CD8 Tdpe" = "#CFF992", "CD8 Tmpe" = "#586B36",
  "CD4 Tnaive" = "#AFCE69", "CD4 Tem" = "#59903F", "CD4 Treg" = "#71CD6D", 
  "NK" = "#80AAB1",
  "NK-XCL1" = "#3B8777",
  "cycling NK" = "#1E4B2B", "cycling NK/T" = "#679163"
  # "Stromal" = "#C1B39B"
  )

color.panel.cell.population = c(
  "HSPCs" = "#D4211A", "Monocytes" = "#C3AEE0", 
  "Dendritic cells"= "#DC7835", "Erythrocytes" = "#F0C7C6", 
  "B cells" = "#1F7EE0","CD8+ T cells" = "#57B170","CD4+ T cells" = "#ADC969",
  "NK cells" = "#7FA8B0"
  # "Stromal" = "#C1B39B"
  )



