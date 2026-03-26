`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

pick_first_existing <- function(values, preferred) {
  matches <- preferred[preferred %in% values]
  if (length(matches) == 0L) {
    return(NULL)
  }
  matches[[1]]
}

parse_cli_args <- function(args) {
  if (length(args) == 0L) {
    stop("Expected a command: inspect or convert.")
  }

  command <- args[[1]]
  rest <- args[-1]
  parsed <- list()
  i <- 1L
  while (i <= length(rest)) {
    token <- rest[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token)
    }

    key <- sub("^--", "", token)
    next_token <- if (i < length(rest)) rest[[i + 1L]] else NULL
    if (is.null(next_token) || startsWith(next_token, "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1L
      next
    }

    parsed[[key]] <- next_token
    i <- i + 2L
  }

  list(command = command, options = parsed)
}

read_input_object <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("Input file does not exist: ", input_path)
  }

  rds_result <- tryCatch(readRDS(input_path), error = function(err) err)
  if (!inherits(rds_result, "error")) {
    return(rds_result)
  }

  load_env <- new.env(parent = emptyenv())
  load_result <- tryCatch(load(input_path, envir = load_env), error = function(err) err)
  if (!inherits(load_result, "error")) {
    if (length(load_result) != 1L) {
      stop(
        "Input is an R workspace with ", length(load_result),
        " objects. Expected exactly one object."
      )
    }
    return(get(load_result[[1]], envir = load_env, inherits = FALSE))
  }

  stop(
    "Could not read input as RDS or R workspace. ",
    "readRDS error: ", conditionMessage(rds_result), ". ",
    "load error: ", conditionMessage(load_result), "."
  )
}

split_csv_arg <- function(value) {
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(NULL)
  }
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

sanitize_atomic_column <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x) || is.integer(x) || is.numeric(x)) {
    return(x)
  }
  if (inherits(x, "Date")) {
    return(as.character(x))
  }
  if (inherits(x, "POSIXt")) {
    return(format(x, tz = "UTC", usetz = TRUE))
  }
  as.character(x)
}

sanitize_data_frame <- function(df, rownames_value = NULL) {
  out <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  for (column_name in names(out)) {
    column <- out[[column_name]]
    if (is.list(column) && !is.data.frame(column)) {
      out[[column_name]] <- vapply(column, function(value) {
        if (length(value) == 0L || is.null(value)) {
          return(NA_character_)
        }
        paste(as.character(value), collapse = "|")
      }, character(1))
      next
    }
    out[[column_name]] <- sanitize_atomic_column(column)
  }

  if (!is.null(rownames_value)) {
    rownames(out) <- rownames_value
  } else if (!is.null(rownames(df))) {
    rownames(out) <- rownames(df)
  }
  out
}

is_non_empty_matrix <- function(x) {
  !is.null(x) && length(dim(x)) == 2L && all(dim(x) > 0L)
}

extract_assay_feature_names <- function(assay_obj, fallback_matrix = NULL) {
  if (!is.null(fallback_matrix) && !is.null(rownames(fallback_matrix))) {
    return(as.character(rownames(fallback_matrix)))
  }
  if ("meta.features" %in% slotNames(assay_obj)) {
    meta_features <- assay_obj@meta.features
    if (is.data.frame(meta_features) && nrow(meta_features) > 0L && !is.null(rownames(meta_features))) {
      return(as.character(rownames(meta_features)))
    }
  }
  character()
}

extract_seurat_layer_names <- function(assay_obj) {
  layer_names <- character()
  if ("layers" %in% slotNames(assay_obj)) {
    layer_names <- c(layer_names, names(assay_obj@layers))
  }
  for (slot_name in c("counts", "data", "scale.data")) {
    if (!(slot_name %in% slotNames(assay_obj))) {
      next
    }
    value <- methods::slot(assay_obj, slot_name)
    if (is_non_empty_matrix(value)) {
      layer_names <- c(layer_names, slot_name)
    }
  }
  unique(layer_names)
}

extract_seurat_layer_matrix <- function(assay_obj, layer_name) {
  if ("layers" %in% slotNames(assay_obj) && layer_name %in% names(assay_obj@layers)) {
    return(assay_obj@layers[[layer_name]])
  }
  if (layer_name %in% slotNames(assay_obj)) {
    value <- methods::slot(assay_obj, layer_name)
    if (is_non_empty_matrix(value)) {
      return(value)
    }
  }
  NULL
}

choose_seurat_assay_name <- function(x, requested_assay = NULL) {
  assay_names <- names(x@assays)
  if (length(assay_names) == 0L) {
    stop("Seurat object contains no assays.")
  }
  if (!is.null(requested_assay) && nzchar(requested_assay)) {
    if (!(requested_assay %in% assay_names)) {
      stop(
        "Requested assay not found in Seurat object: ", requested_assay,
        ". Available assays: ", paste(assay_names, collapse = ", ")
      )
    }
    return(requested_assay)
  }
  preferred <- unique(c("SCT", "Spatial", x@active.assay, "RNA", "integrated"))
  pick_first_existing(assay_names, preferred) %||% assay_names[[1]]
}

choose_sce_assay_name <- function(x, requested_assay = NULL) {
  assay_names <- SummarizedExperiment::assayNames(x)
  if (length(assay_names) == 0L) {
    stop("SingleCellExperiment contains no assays.")
  }
  if (!is.null(requested_assay) && nzchar(requested_assay)) {
    if (!(requested_assay %in% assay_names)) {
      stop(
        "Requested assay not found: ", requested_assay,
        ". Available assays: ", paste(assay_names, collapse = ", ")
      )
    }
    return(requested_assay)
  }
  pick_first_existing(assay_names, c("logcounts", "data", "normalized", "counts")) %||% assay_names[[1]]
}

choose_x_name <- function(names_value, preferred) {
  pick_first_existing(names_value, preferred) %||% names_value[[1]]
}

extract_reduced_dim_matrix <- function(x) {
  if (!("cell.embeddings" %in% slotNames(x))) {
    return(NULL)
  }
  matrix_value <- x@cell.embeddings
  if (!is_non_empty_matrix(matrix_value)) {
    return(NULL)
  }
  matrix_value
}

extract_seurat_image_coordinate_df <- function(x) {
  if (!("images" %in% slotNames(x)) || length(x@images) == 0L) {
    return(NULL)
  }

  coord_frames <- list()
  for (image_name in names(x@images)) {
    image_obj <- x@images[[image_name]]
    if ("boundaries" %in% slotNames(image_obj) && length(image_obj@boundaries) > 0L) {
      centroid_name <- pick_first_existing(names(image_obj@boundaries), c("centroids", "centroid"))
      if (!is.null(centroid_name)) {
        centroid_obj <- image_obj@boundaries[[centroid_name]]
        if (all(c("coords", "cells") %in% slotNames(centroid_obj))) {
          centroid_coords <- centroid_obj@coords
          centroid_cells <- centroid_obj@cells
          if (is_non_empty_matrix(centroid_coords) &&
              length(centroid_cells) == nrow(centroid_coords)) {
            centroid_coords <- as.matrix(centroid_coords)
            if (ncol(centroid_coords) >= 2L) {
              x_col <- pick_first_existing(colnames(centroid_coords), c("x", "X", "imagecol", "col"))
              y_col <- pick_first_existing(colnames(centroid_coords), c("y", "Y", "imagerow", "row"))
              if (is.null(x_col) || is.null(y_col)) {
                x_col <- colnames(centroid_coords)[[1]]
                y_col <- colnames(centroid_coords)[[2]]
              }

              coord_frames[[image_name]] <- data.frame(
                x = as.numeric(centroid_coords[, x_col]),
                y = as.numeric(centroid_coords[, y_col]),
                row.names = as.character(centroid_cells),
                check.names = FALSE,
                stringsAsFactors = FALSE
              )
              next
            }
          }
        }
      }
    }

    if (!("coordinates" %in% slotNames(image_obj))) {
      next
    }

    image_coords <- image_obj@coordinates
    if (!is.data.frame(image_coords) || nrow(image_coords) == 0L) {
      next
    }

    x_col <- pick_first_existing(colnames(image_coords), c("imagecol", "col", "x"))
    y_col <- pick_first_existing(colnames(image_coords), c("imagerow", "row", "y"))
    if (is.null(x_col) || is.null(y_col)) {
      numeric_cols <- colnames(image_coords)[vapply(image_coords, is.numeric, logical(1))]
      if (length(numeric_cols) < 2L) {
        next
      }
      x_col <- x_col %||% numeric_cols[[1]]
      y_col <- y_col %||% numeric_cols[[2]]
    }

    coord_frames[[image_name]] <- data.frame(
      x = as.numeric(image_coords[[x_col]]),
      y = as.numeric(image_coords[[y_col]]),
      row.names = rownames(image_coords),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (length(coord_frames) == 0L) {
    return(NULL)
  }

  coord_df <- do.call(rbind, unname(coord_frames))
  coord_df[!duplicated(rownames(coord_df)), , drop = FALSE]
}

extract_spatial_matrix <- function(x) {
  if (inherits(x, "SpatialExperiment")) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
      return(NULL)
    }
    spatial <- SpatialExperiment::spatialCoords(x)
    if (is_non_empty_matrix(spatial)) {
      return(as.matrix(spatial))
    }
  }

  if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      return(NULL)
    }
    rd_names <- SingleCellExperiment::reducedDimNames(x)
    for (candidate in c("spatial", "Spatial")) {
      if (candidate %in% rd_names) {
        spatial <- SingleCellExperiment::reducedDim(x, candidate)
        if (is_non_empty_matrix(spatial)) {
          return(as.matrix(spatial))
        }
      }
    }
    return(NULL)
  }

  if (inherits(x, "Seurat")) {
    coord_df <- extract_seurat_image_coordinate_df(x)
    if (!is.null(coord_df) && nrow(coord_df) > 0L) {
      return(as.matrix(coord_df[, c("x", "y"), drop = FALSE]))
    }
  }

  if (is.list(x) && !is.null(x$coordinates)) {
    coords <- as.matrix(x$coordinates)
    if (is_non_empty_matrix(coords)) {
      return(coords)
    }
  }

  NULL
}

normalize_reduced_dim_names <- function(matrix_value, prefix) {
  if (is.null(colnames(matrix_value))) {
    colnames(matrix_value) <- sprintf("%s_%d", prefix, seq_len(ncol(matrix_value)))
  }
  matrix_value
}

extract_sce_reduced_dims <- function(x, requested = NULL, include_spatial = TRUE) {
  out <- list()
  rd_names <- SingleCellExperiment::reducedDimNames(x)
  if (!is.null(requested)) {
    rd_names <- intersect(requested, rd_names)
  }
  for (rd_name in rd_names) {
    matrix_value <- SingleCellExperiment::reducedDim(x, rd_name)
    if (!is_non_empty_matrix(matrix_value)) {
      next
    }
    out[[rd_name]] <- normalize_reduced_dim_names(as.matrix(matrix_value), rd_name)
  }

  if (include_spatial && !("spatial" %in% names(out))) {
    spatial <- extract_spatial_matrix(x)
    if (!is.null(spatial)) {
      out[["spatial"]] <- normalize_reduced_dim_names(spatial, "spatial")
    }
  }
  out
}

extract_seurat_reduced_dims <- function(x, requested = NULL, include_spatial = TRUE) {
  out <- list()
  reduction_names <- names(x@reductions)
  if (!is.null(requested)) {
    reduction_names <- intersect(requested, reduction_names)
  }

  for (reduction_name in reduction_names) {
    matrix_value <- extract_reduced_dim_matrix(x@reductions[[reduction_name]])
    if (is.null(matrix_value)) {
      next
    }
    out[[reduction_name]] <- normalize_reduced_dim_names(as.matrix(matrix_value), reduction_name)
  }

  if (include_spatial && !("spatial" %in% names(out))) {
    spatial <- extract_spatial_matrix(x)
    if (!is.null(spatial)) {
      out[["spatial"]] <- normalize_reduced_dim_names(spatial, "spatial")
    }
  }
  out
}

build_sce_from_list <- function(x, x_layer_name = NULL, requested_layers = NULL, requested_reduced_dims = NULL, include_spatial = TRUE) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("List conversion requires SingleCellExperiment, SummarizedExperiment, and S4Vectors.")
  }

  obs <- x$obs %||% NULL
  expr <- x$expression %||% x$expr %||% NULL
  if (is.null(obs) || !is.data.frame(obs)) {
    stop("List inputs must include an obs data.frame.")
  }
  if (is.null(expr)) {
    stop("List inputs must include an expression matrix in $expression or $expr.")
  }

  expr <- as.matrix(expr)
  obs <- sanitize_data_frame(obs, rownames(obs) %||% colnames(expr))
  if (nrow(obs) != ncol(expr)) {
    stop("List input obs rows must match expression columns.")
  }

  gene_names <- rownames(expr) %||% x$gene_names %||% sprintf("gene_%d", seq_len(nrow(expr)))
  row_data <- S4Vectors::DataFrame(row.names = as.character(gene_names))
  assays <- list()
  assays[[x_layer_name %||% "X"]] <- expr

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs),
    rowData = row_data
  )

  if (!is.null(x$umap) && (is.null(requested_reduced_dims) || "UMAP" %in% requested_reduced_dims)) {
    SingleCellExperiment::reducedDim(sce, "UMAP") <- as.matrix(x$umap)
  }
  if (include_spatial && !is.null(x$coordinates)) {
    SingleCellExperiment::reducedDim(sce, "spatial") <- as.matrix(x$coordinates)
  }

  list(
    sce = sce,
    x_name = names(assays)[[1]],
    summary = list(
      object_type = "list",
      cells = ncol(expr),
      genes = nrow(expr),
      selected_assay = NULL,
      available_assays = character(),
      available_layers = names(assays),
      reduced_dims = SingleCellExperiment::reducedDimNames(sce)
    )
  )
}

build_sce_from_sce <- function(x, requested_assay = NULL, x_layer_name = NULL, requested_layers = NULL, requested_reduced_dims = NULL, include_spatial = TRUE) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("SingleCellExperiment conversion requires SingleCellExperiment, SummarizedExperiment, and S4Vectors.")
  }

  assay_name <- choose_sce_assay_name(x, requested_assay)
  assay_names <- SummarizedExperiment::assayNames(x)
  selected_layers <- requested_layers %||% assay_names
  selected_layers <- intersect(selected_layers, assay_names)
  if (length(selected_layers) == 0L) {
    stop("No compatible assays selected for export.")
  }

  assays <- list()
  for (layer_name in selected_layers) {
    assays[[layer_name]] <- SummarizedExperiment::assay(x, layer_name)
  }

  x_name <- x_layer_name %||% choose_x_name(selected_layers, c("logcounts", "data", "normalized", assay_name, "counts"))
  if (!(x_name %in% names(assays))) {
    stop("Requested x-layer not found among exported assays: ", x_name)
  }

  obs <- sanitize_data_frame(as.data.frame(SummarizedExperiment::colData(x)))
  row_data_source <- tryCatch(as.data.frame(SummarizedExperiment::rowData(x)), error = function(err) NULL)
  if (is.null(row_data_source) || nrow(row_data_source) != nrow(assays[[x_name]])) {
    row_data <- S4Vectors::DataFrame(row.names = rownames(assays[[x_name]]))
  } else {
    row_data <- S4Vectors::DataFrame(sanitize_data_frame(row_data_source, rownames(assays[[x_name]])))
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs),
    rowData = row_data
  )

  reduced_dims <- extract_sce_reduced_dims(
    x = x,
    requested = requested_reduced_dims,
    include_spatial = include_spatial
  )
  for (rd_name in names(reduced_dims)) {
    SingleCellExperiment::reducedDim(sce, rd_name) <- reduced_dims[[rd_name]]
  }

  list(
    sce = sce,
    x_name = x_name,
    summary = list(
      object_type = paste(class(x), collapse = ","),
      cells = ncol(sce),
      genes = nrow(sce),
      selected_assay = assay_name,
      available_assays = assay_names,
      available_layers = assay_names,
      reduced_dims = SingleCellExperiment::reducedDimNames(sce)
    )
  )
}

build_sce_from_seurat <- function(x, requested_assay = NULL, x_layer_name = NULL, requested_layers = NULL, requested_reduced_dims = NULL, include_spatial = TRUE) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Seurat conversion requires SingleCellExperiment and S4Vectors.")
  }

  assay_name <- choose_seurat_assay_name(x, requested_assay)
  assay_obj <- x@assays[[assay_name]]
  available_layers <- extract_seurat_layer_names(assay_obj)
  if (length(available_layers) == 0L) {
    stop("Selected Seurat assay does not contain any exportable layers.")
  }

  selected_layers <- requested_layers %||% available_layers
  selected_layers <- intersect(selected_layers, available_layers)
  if (length(selected_layers) == 0L) {
    stop("No compatible Seurat assay layers selected for export.")
  }

  x_name <- x_layer_name %||% choose_x_name(selected_layers, c("data", "counts", "scale.data"))
  if (!(x_name %in% selected_layers)) {
    stop("Requested x-layer not found among exported layers: ", x_name)
  }

  x_matrix <- extract_seurat_layer_matrix(assay_obj, x_name)
  if (is.null(x_matrix)) {
    stop("Could not resolve Seurat layer for X: ", x_name)
  }

  cell_names <- colnames(x_matrix)
  gene_names <- extract_assay_feature_names(assay_obj, fallback_matrix = x_matrix)
  if (length(gene_names) == 0L) {
    gene_names <- sprintf("gene_%d", seq_len(nrow(x_matrix)))
  }

  obs <- sanitize_data_frame(x@meta.data, rownames_value = rownames(x@meta.data))
  if (!all(cell_names %in% rownames(obs))) {
    stop("Seurat meta.data row names do not cover all assay columns.")
  }
  obs <- obs[cell_names, , drop = FALSE]

  row_data_source <- NULL
  if ("meta.features" %in% slotNames(assay_obj)) {
    meta_features <- assay_obj@meta.features
    if (is.data.frame(meta_features) && nrow(meta_features) == length(gene_names)) {
      row_data_source <- sanitize_data_frame(meta_features, rownames_value = gene_names)
    }
  }
  row_data <- if (is.null(row_data_source)) {
    S4Vectors::DataFrame(row.names = gene_names)
  } else {
    S4Vectors::DataFrame(row_data_source)
  }

  assays <- list()
  for (layer_name in selected_layers) {
    matrix_value <- extract_seurat_layer_matrix(assay_obj, layer_name)
    if (is.null(matrix_value)) {
      next
    }
    if (!identical(dim(matrix_value), dim(x_matrix))) {
      warning("Skipping layer with incompatible dimensions: ", layer_name, call. = FALSE)
      next
    }
    if (!identical(colnames(matrix_value), cell_names)) {
      matrix_value <- matrix_value[, cell_names, drop = FALSE]
    }
    assays[[layer_name]] <- matrix_value
  }
  if (!(x_name %in% names(assays))) {
    stop("The chosen x-layer was dropped due to incompatible dimensions: ", x_name)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs),
    rowData = row_data
  )

  reduced_dims <- extract_seurat_reduced_dims(
    x = x,
    requested = requested_reduced_dims,
    include_spatial = include_spatial
  )
  for (rd_name in names(reduced_dims)) {
    matrix_value <- reduced_dims[[rd_name]]
    if (!identical(rownames(matrix_value), cell_names) && !is.null(rownames(matrix_value))) {
      if (all(cell_names %in% rownames(matrix_value))) {
        matrix_value <- matrix_value[cell_names, , drop = FALSE]
      }
    }
    if (nrow(matrix_value) == ncol(sce)) {
      SingleCellExperiment::reducedDim(sce, rd_name) <- matrix_value
    }
  }

  list(
    sce = sce,
    x_name = x_name,
    summary = list(
      object_type = paste(class(x), collapse = ","),
      cells = ncol(sce),
      genes = nrow(sce),
      selected_assay = assay_name,
      available_assays = names(x@assays),
      available_layers = available_layers,
      reduced_dims = SingleCellExperiment::reducedDimNames(sce)
    )
  )
}

object_summary <- function(x, requested_assay = NULL) {
  if (inherits(x, "Seurat")) {
    assay_name <- choose_seurat_assay_name(x, requested_assay)
    assay_obj <- x@assays[[assay_name]]
    layer_names <- extract_seurat_layer_names(assay_obj)
    x_name <- choose_x_name(layer_names, c("data", "counts", "scale.data"))
    x_matrix <- extract_seurat_layer_matrix(assay_obj, x_name)
    spatial <- extract_spatial_matrix(x)
    return(list(
      object_type = paste(class(x), collapse = ","),
      selected_assay = assay_name,
      available_assays = names(x@assays),
      available_layers = layer_names,
      reduced_dims = names(x@reductions),
      cells = ncol(x_matrix),
      genes = nrow(x_matrix),
      has_spatial = !is.null(spatial)
    ))
  }

  if (inherits(x, "SingleCellExperiment")) {
    assay_name <- choose_sce_assay_name(x, requested_assay)
    assay_names <- SummarizedExperiment::assayNames(x)
    matrix_value <- SummarizedExperiment::assay(x, assay_name)
    spatial <- extract_spatial_matrix(x)
    return(list(
      object_type = paste(class(x), collapse = ","),
      selected_assay = assay_name,
      available_assays = assay_names,
      available_layers = assay_names,
      reduced_dims = SingleCellExperiment::reducedDimNames(x),
      cells = ncol(matrix_value),
      genes = nrow(matrix_value),
      has_spatial = !is.null(spatial)
    ))
  }

  if (is.list(x)) {
    expr <- x$expression %||% x$expr %||% NULL
    coords <- extract_spatial_matrix(x)
    return(list(
      object_type = "list",
      selected_assay = NULL,
      available_assays = character(),
      available_layers = if (is.null(expr)) character() else "X",
      reduced_dims = Filter(Negate(is.null), c(if (!is.null(x$umap)) "UMAP", if (!is.null(coords)) "spatial")),
      cells = if (is.null(expr)) NA_integer_ else ncol(expr),
      genes = if (is.null(expr)) NA_integer_ else nrow(expr),
      has_spatial = !is.null(coords)
    ))
  }

  stop(
    "Unsupported input class: ", paste(class(x), collapse = ", "),
    ". Supported inputs: Seurat, SingleCellExperiment, SpatialExperiment, or a plain list with obs/expression."
  )
}

build_export_sce <- function(x, requested_assay = NULL, x_layer_name = NULL, requested_layers = NULL, requested_reduced_dims = NULL, include_spatial = TRUE) {
  if (inherits(x, "Seurat")) {
    return(build_sce_from_seurat(x, requested_assay, x_layer_name, requested_layers, requested_reduced_dims, include_spatial))
  }
  if (inherits(x, "SingleCellExperiment")) {
    return(build_sce_from_sce(x, requested_assay, x_layer_name, requested_layers, requested_reduced_dims, include_spatial))
  }
  if (is.list(x)) {
    return(build_sce_from_list(x, x_layer_name, requested_layers, requested_reduced_dims, include_spatial))
  }
  stop("Unsupported input class for conversion: ", paste(class(x), collapse = ", "))
}

write_h5ad <- function(sce, output_path, x_name) {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop(
      "zellkonverter is required to write .h5ad files. ",
      "Install it in R with BiocManager::install(\"zellkonverter\")."
    )
  }
  zellkonverter::writeH5AD(sce, file = output_path, X_name = x_name)
}

json_write <- function(value) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for inspect output.")
  }
  if (!is.null(value$available_assays)) {
    value$available_assays <- unname(as.list(value$available_assays))
  }
  if (!is.null(value$available_layers)) {
    value$available_layers <- unname(as.list(value$available_layers))
  }
  if (!is.null(value$reduced_dims)) {
    value$reduced_dims <- unname(as.list(value$reduced_dims))
  }
  cat(jsonlite::toJSON(value, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  cat("\n")
}

run_inspect <- function(options) {
  input_path <- options$input %||% stop("Missing required --input.")
  requested_assay <- options$assay %||% NULL
  obj <- read_input_object(input_path)
  summary <- object_summary(obj, requested_assay)
  summary$input <- normalizePath(input_path, mustWork = TRUE)
  json_write(summary)
}

run_convert <- function(options) {
  input_path <- options$input %||% stop("Missing required --input.")
  output_path <- options$output %||% stop("Missing required --output.")
  requested_assay <- options$assay %||% NULL
  x_layer_name <- options[["x-layer"]] %||% NULL
  requested_layers <- split_csv_arg(options$layers %||% NULL)
  requested_reduced_dims <- split_csv_arg(options[["reduced-dims"]] %||% NULL)
  include_spatial <- !isTRUE(options[["no-spatial"]])

  obj <- read_input_object(input_path)
  export_result <- build_export_sce(
    x = obj,
    requested_assay = requested_assay,
    x_layer_name = x_layer_name,
    requested_layers = requested_layers,
    requested_reduced_dims = requested_reduced_dims,
    include_spatial = include_spatial
  )

  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  write_h5ad(export_result$sce, output_path = output_path, x_name = export_result$x_name)

  summary <- export_result$summary
  summary$output <- normalizePath(output_path, mustWork = FALSE)
  summary$x_name <- export_result$x_name
  json_write(summary)
}

run_validate <- function(options) {
  input_path <- options$input %||% stop("Missing required --input.")
  if (!file.exists(input_path)) {
    stop("Input file does not exist: ", input_path)
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop(
      "zellkonverter is required to validate .h5ad files. ",
      "Install it in R with BiocManager::install(\"zellkonverter\")."
    )
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment is required to validate .h5ad files.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required to validate .h5ad files.")
  }

  sce <- zellkonverter::readH5AD(input_path)
  assay_names <- SummarizedExperiment::assayNames(sce)
  reduced_dims <- SingleCellExperiment::reducedDimNames(sce)
  obs_cols <- colnames(SummarizedExperiment::colData(sce))
  var_cols <- colnames(SummarizedExperiment::rowData(sce))

  summary <- list(
    input = normalizePath(input_path, mustWork = TRUE),
    cells = ncol(sce),
    genes = nrow(sce),
    assays = assay_names,
    reduced_dims = reduced_dims,
    obs_columns = obs_cols,
    var_columns = var_cols
  )

  if (length(assay_names) > 0L) {
    assay_dims <- lapply(assay_names, function(name) {
      matrix_value <- SummarizedExperiment::assay(sce, name)
      list(name = name, dims = dim(matrix_value))
    })
    summary$assay_dims <- assay_dims
  }

  if ("spatial" %in% reduced_dims) {
    spatial <- SingleCellExperiment::reducedDim(sce, "spatial")
    summary$spatial_dims <- dim(spatial)
    summary$spatial_columns <- colnames(spatial)
  }

  json_write(summary)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- parse_cli_args(args)
  switch(
    parsed$command,
    inspect = run_inspect(parsed$options),
    convert = run_convert(parsed$options),
    validate = run_validate(parsed$options),
    stop("Unknown command: ", parsed$command)
  )
}

tryCatch(
  main(),
  error = function(err) {
    message(conditionMessage(err))
    quit(save = "no", status = 1L)
  }
)
