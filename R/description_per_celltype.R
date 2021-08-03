#' Description of various statistics for each celltype in the Seurat object
#'
#' @param seurat_object Seurat single-cell RNA sequencing object
#'
#' @importFrom SeuratObject FetchData
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%

description_per_celltype <- function(seurat_object,
                                     condition_to_compare,
                                     condition_ident,
                                     celltype_ident,
                                     sample_ident,
                                     celltypes) {

  results <- list()

  # Subset the Seurat object for the right condition
  conditions <- FetchData(object = seurat_object, vars = condition_ident)
  subset_seurat <- seurat_object[, which(x = conditions == condition_to_compare)]

  # Amount of cells and percentage of each celltype for this condition
  results$cells_per_condition <- table(subset_seurat$celltype)
  results$cell_percentage_per_condition <- results$cells_per_condition / nrow(subset_seurat@meta.data)

  # Calculate amount of cells and percentage of each celltype for each sample
  samples <- FetchData(object = subset_seurat, vars = sample_ident) %>% pull(sample_ident) %>% unique()

  cells_per_patient <- c()
  cell_percentages_per_patient <- c()

  # One sample at a time
  for (sample in samples) {
    # Dummy results that contain 0 for all the celltypes
    # 0 for missing values automatically
    result.dummy <- rep(0, length(celltypes))
    names(result.dummy) <- celltypes

    # Subset the seurat object for this sample
    data.samples <- FetchData(object = subset_seurat, vars = sample_ident)
    data.sample <- subset_seurat[, which(x = data.samples == sample)]

    # Amount of cells per patient
    # Order in the same way as result.dummy and add there
    results.preliminary <- table(data.sample@meta.data[[celltype_ident]])
    result.dummy <- results.preliminary[celltypes]

    # Aggregate results and change the rowname at the right spot
    cells_per_patient <- rbind(cells_per_patient, result.dummy)
    rownames(cells_per_patient)[nrow(cells_per_patient)] <- sample

    # Percentages per celltype for the samples individually
    result.dummy <- result.dummy / sum(result.dummy)
    cell_percentages_per_patient <- rbind(cell_percentages_per_patient, result.dummy)
    rownames(cell_percentages_per_patient)[nrow(cell_percentages_per_patient)] <- sample
  }

  results$cells_per_patient <- cells_per_patient
  results$cell_percentages_per_patient <- cell_percentages_per_patient

  # Mean, variance and coefficient of variance per celltype
  results$mean_celltype <- colMeans(cell_percentages_per_patient)
  results$var_celltype <- apply(X = cell_percentages_per_patient, MARGIN = 2, FUN = var)
  results$coefficient_variance_celltype <- results$var_celltype / results$mean_celltype
  results$condition <- condition_to_compare
  return(results)
}
