# globals shared across markdown docs

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(matrixStats)
library(purrr)
library(R.utils)
library(viridis)
library(ComplexHeatmap)
library(here)
#### Paths ####

project_dir <- here()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")


##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file.
#' Useful for avoid gene names being clobbered to dates and
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}

#' Extract out reduced dimensions and cell metadata to tibble
#'
#' @param obj Seurat Object
#' @param embedding dr slot to extract (defaults to tsne)
#'
get_metadata <- function(obj, embedding = "tsne") {

  mdata <- as_tibble(obj@meta.data, rownames = "cell")

  if(!embedding %in% names(obj@reductions)){
    stop(paste0(embedding, " not found in seurat object"), call. = FALSE)
  }

  embed_dat <- obj@reductions[[embedding]]@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  embed_dat <- left_join(mdata,
                         embed_dat,
                         by = "cell")
  embed_dat
}



#' Download a file if it doesn't exist
#' @param url download url
#' @param filename output file name
#' @param gunzip uncompress file using gzip Default = FALSE, if suffix contains .gz this will
#' be removed.
#' @export
get_file <- function(url, filename, gunzip = FALSE){
  if(!file.exists(filename)){
    download.file(url, filename)
  }

  if(gunzip){
    if(endsWith(filename, ".gz")) {
      R.utils::gunzip(filename)
    } else {
      stop(paste0(filename, " not a gzipped file"))
    }
  }
}


tcols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])[1:6]

# from colorblindr
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

palette_OkabeIto_black <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot())

discrete_palette_default <- c(brewer.pal(12, "Paired"),
                              brewer.pal(9, "Set1"),
                              brewer.pal(8, "Set2"),
                              brewer.pal(8, "Dark2"))

#' Plot cells in reduced dimensionality 2D space
#'
#' @description Cells can be colored by gene or feature in meta.data dataframe
#'
#' @param seurat_obj object of class Seurat
#' @param feature feature to plot, either gene name or column in seurat_obj@meta.data
#' @param plot_dat supplemental data.frame containing feature to plot.
#' Must have a column named cell that contains matching colnames in seurat_obj@data
#' @param pt_size size of points produced by geom_point
#' @param pt_alpha alpha value for points plotted by geom_point
#' @param label_text if TRUE display feature labels on plot
#' @param label_size size of label text
#' @param label_color color of label text
#' @param .cols vector of colors to use for plot.
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' dr method present in seurat_obj@dr (e.g. umap, pca, tsne). defaults to tsne
#'
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = "grey",
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne"){

  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")

  if(!embedding %in% names(seurat_obj@reductions)){
    stop(paste0(embedding, " not found in seurat object"))
  }

  embed_dat <- seurat_obj@reductions[[embedding]]@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  embed_cols <- colnames(embed_dat)
  xcol <- embed_cols[2]
  ycol <- embed_cols[3]

  embed_dat <- left_join(mdata, embed_dat, by = "cell")

  if (!is.null(cell_filter)){
    embed_dat <- dplyr::filter(embed_dat,
                               cell %in% cell_filter)
  }

  meta_data_col <- feature %in% colnames(embed_dat)

  if (!is.null(feature) & !meta_data_col) {
    feature_dat <- FetchData(seurat_obj, feature) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")
    embed_dat <- left_join(embed_dat, feature_dat, by = "cell")
  }

  if (!is.null(plot_dat)){
    embed_dat <- left_join(embed_dat, plot_dat, by = "cell")
  }

  color_aes_str <- feature

  color_aes_str_q <- quo(color_aes_str)
  embed_dat <- embed_dat %>% arrange_at(.vars = color_aes_str)

  p <- ggplot(embed_dat,
              aes_string(xcol, ycol)) +
    geom_point(aes_string(color = paste0("`", color_aes_str, "`")),
               size = pt_size,
               alpha = pt_alpha)

  ## discrete or continuous data?
  if (typeof(embed_dat[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(embed_dat[[feature]])) {
    discrete <- T
  } else {
    discrete <- F
  }

  ## increase legend size
  if (discrete) {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_blank())
  }

  if (label_text) {
    if(discrete) {
      tsne_mean_dat <- embed_dat %>%
        group_by_at(vars(one_of(feature))) %>%
        summarize(med_dim_1 = median(tSNE_1),
                  med_dim_2 = median(tSNE_2))

      p <- p +
        geom_text(data = tsne_mean_dat,
                  aes_string(x = "med_dim_1",
                             y = "med_dim_2",
                             label = feature),
                  size = label_size,
                  color = label_color)
    } else {
      warning("label_text not compatible with continuous features")
    }
  }

  ## handle legend limit
  if (is.null(max_y) & !discrete) {
    max_y <- c(0, max(embed_dat[[color_aes_str]]))
  } else if (discrete & is.null(max_y)){
    max_y <- c(NA, NA)
  }

  # loupe-like colors
  cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])

  #handle legend name
  if(is.null(legend_title)) legend_title <- color_aes_str

  ## handle zero expression
  if (!all(is.na(max_y)) && all(max_y == c(0, 0))){
    p <- p + scale_color_gradient(low = cols[1], high = cols[1], name = legend_title)
    return(p)
  }

  ## handle colors
  if (is.null(.cols) && !discrete){
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y, name = legend_title)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, name = legend_title)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols, name = legend_title)
    }
  } else if (!is.null(.cols) && !discrete){
    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = .cols, name = legend_title)
  } else {

    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title
      )
    }
  }
  p
}


get_high_low <- function(vec, n){
  half <- n / 2
  split <- floor(half)
  len_v <- length(vec)
  vec[c(1:split, (len_v - split):len_v)]
}

importCDS_v3 <- function (otherCDS, seurat_scale=F, import_all = FALSE)
{
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if (!seurat_scale) {
      data <- otherCDS@assays$RNA@counts
    } else {
      data <- otherCDS@assays$RNA@scale.data
    }
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data),
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    #lowerDetectionLimit <- otherCDS@is.expr
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
      expr <- "negbinomial.size"
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
      expr <- "unimormal"
    }
    else {
      expressionFamily <- tobit()
      expr <- "tobit"
    }
    print(paste0("expressionFamily ",expr))
    # valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                                  #lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if ("var.genes" %in% slotNames(otherCDS)) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }
  else if (class(otherCDS)[1] == "SCESet") {
    requireNamespace("scater")
    message("Converting the exprs data in log scale back to original scale ...")
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if ("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else lowerDetectionLimit <- 1
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    if (import_all) {
      mist_list <- otherCDS
    }
    else {
      mist_list <- list()
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    monocle_cds@auxOrderingData$scran <- mist_list
  }
  else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
}
