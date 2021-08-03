#' @title
#' sc_prop
#'
#' @description
#' Dirichlet regression to test for differences in cell proportions between conditions in single-cell RNA sequencing data
#'
#' @return
#' Output consists of:
#' 1. dirichlet_plot.pdf: Ggplot2 visualization of the cell proportions per condition
#' 2. dirichlet_regression_results.csv: Results for all comparisons made
#' 3. dirichlet_regression_results_sig.csv: Significant results if any
#'
#' @param seurat_object Seurat single-cell RNA sequencing object
#' @param conditions Vector of conditions to compare. All conditions must be present in the Seurat object. Order in which conditions are given determines order in the final figure
#' @param celltypes Vector of celltypes to compare in the Seurat object
#' @param condition_ident Name of the metadata column containing the condition. Default: condition.
#' @param sample_ident Name of the metadata column containing the sample information. Default: sample
#' @param celltype_ident Name of the metadata column (Seurat Ident) containing the celltypes. Default: celltype
#' @param outdir Output directory. Will raise error if does not exist. Default: current directory
#'
#' @import Seurat
#' @import DirichletReg
#' @import ggsci
#' @import ggpubr
#' @import ggplot2
#'
#' @importFrom reshape2 melt
#' @importFrom dplyr pull mutate filter group_by summarise arrange
#' @importFrom SeuratObject FetchData
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @export "sc_prop"
sc_prop <- function(seurat_object,
                    conditions,
                    celltypes,
                    outdir,
                    condition_ident = 'condition',
                    sample_ident = 'sample',
                    celltype_ident = 'celltype', 
                    cols = NULL) {
  
  #1. Check input arguments
  if( !dir.exists(outdir) ) {
    stop('Output directory does not exist')
  }
  
  if(is.null(cols)) {
    cols <- pal_nejm()(length(unique(conditions)))
  } else {
    stopifnot(length(cols) == length(unique(conditions)))
  }
  
  stopifnot("Idents must be in the Seurat object" = 
              condition_ident %in% colnames(seurat_object@meta.data) & 
              sample_ident %in% colnames(seurat_object@meta.data) &
              celltype_ident %in% colnames(seurat_object@meta.data))
  
  condition_subset <- FetchData(object = seurat_object, vars = condition_ident)
  seurat_object <- seurat_object[ , which(x = condition_subset[, 1] %in% conditions)]
  
  celltypes_subset <- FetchData(object = seurat_object, vars = celltype_ident)
  seurat_object <- seurat_object[ , which(x = celltypes_subset[, 1] %in% celltypes)]
  
  #2. Descriptive statistics per celltype
  res <- lapply(X = conditions,
                FUN = description_per_celltype,
                seurat_object = seurat_object,
                condition_ident = condition_ident,
                sample_ident = sample_ident,
                celltype_ident = celltype_ident,
                celltypes = celltypes)
  
  names(res) <- conditions

    #3. Combine data for regression and final plot
  data.combined <- c()
  
  for (condition in conditions) {
    data.boxplot <- data.frame(res[[condition]][['cell_percentages_per_patient']],
                               check.names = F) %>%
      mutate(patient = row.names(.)) %>%
      melt(., id.vars = 'patient', value.name = 'percentage') %>%
      mutate(condition = condition)
    data.combined <- rbind(data.combined, data.boxplot)
  }
  
  colnames(data.combined) <- c('patient', 'celltype', 'percentage', 'condition')
  
  #4. Dirichlet Regression to find which comparisons are significant
  regress.results <- data.frame()
  
  for (i in 1:length(conditions)) {
    
    # We prepend each of the conditions with a LETTER
    # This is necessary because we need to iterate over the conditions and supply in a different order each time
    # Ordering the letters is consistent instead of worrying about factors
    # prepend <- LETTERS[unique(c(i:length(as.character(unique(data[[identity.to.compare]])[,1])), 1:i))]
    
    prepend <- LETTERS[unique(c(i:length(conditions), 1:i))]
    
    data.test <- data.combined %>%
      mutate(rem = 1 - percentage)
    
    for (prepending in 1:length(prepend)) {
      data.test[which(data.test$condition == conditions[prepending]), 'condition'] <- paste0(prepend[prepending], conditions[prepending])
    }
    
    conditions.vec <- unique(data.test$condition)
    conditions.vec <- conditions.vec[order(conditions.vec)]
    
    for (cellt in celltypes) {
      
      try( {
        data.celltype <- data.test %>% filter(celltype == cellt)
        data.celltype %<>% arrange(condition)
        
        data.celltype$dr_data <- DR_data(data.celltype[, c('percentage', 'rem')])

        dirich_res <- DirichReg(formula = dr_data ~ condition,
                                data = data.celltype,
                                model = "alternative",
                                base = 1,
                                verbosity = 0,
                                control = list(iterlim = 5000))

        # Weird construction to keep the summary(dirich_res) to spit output to console
        invisible(capture.output(x <- summary(dirich_res)))
        coefs <- x$coef.mat
        coefs <- coefs[which(! rownames(coefs) %in% '(Intercept)'),  , drop=F]
        rownames(coefs) <- substr(rownames(coefs), start = 11, stop = nchar(rownames(coefs)))
        conditions.stripped <- rownames(coefs)

        # Add results to the main dataframe
        for (k in 1:length(rownames(coefs))) {
          tmp <- data.frame('condition1' = substr(conditions.vec[1], start = 2, stop = nchar(conditions.vec[1])),
                            'condition2' = rownames(coefs)[k],
                            'celltype' = cellt,
                            'Z' = coefs[rownames(coefs)[k], 'z value'],
                            'pval' = coefs[rownames(coefs)[k], 'Pr(>|z|)'])
          
          regress.results <- rbind(regress.results, tmp)
          
        }
      })
    }
  }
  
  # 5. Generate the output
  # Filter out the duplicates
  write.csv(x = regress.results, file = paste0(outdir, '/dirichlet_regression_results_all.csv'))
  regress.results <- regress.results[which(regress.results$condition1 < regress.results$condition2), ]
  regress.results$padj <- p.adjust(regress.results$pval, method = 'BH')
  stopifnot(all(celltypes %in% regress.results$celltype))
  
  # Export only sig
  write.csv(x = regress.results %>% filter(padj < 0.05), file = paste0(outdir, '/dirichlet_regression_results_sig.csv'))
  
  # Construct the plot
  # order the celltypes so that higher proportions on average come to the left.
  data.combined %>%
    group_by(celltype) %>%
    summarise(mean = mean(percentage), max = max(percentage)) %>%
    arrange(desc(mean)) -> coordinates
  
  data.combined$condition <- factor(data.combined$condition, levels = conditions)

  # Get the order
  celltypes_ord <- coordinates %>% pull(celltype)
  
  # Construct the figure
  write.csv(paste0(outdir, 'plot_data.csv'), x = data.combined)
  
  pl <- ggplot(data = data.combined,
               aes(x = condition, y = percentage, fill = condition)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2)) +
    labs(x = 'Condition', y = 'Cell proportion', fill = 'Condition') +
    facet_wrap(~factor(celltype, levels = celltypes_ord), nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = cols)

  if(any(regress.results$padj < 0.05)) {
    df.sig <- regress.results %>% filter(padj < 0.05)
    df.sig$celltype <- factor(df.sig$celltype, levels = celltypes_ord)
    df.sig$group1 <- df.sig$condition1
    df.sig$group2 <- df.sig$condition2
    
    df.sig$label <- if (df.sig$padj < 0.001) {
      '***'
    } else if (df.sig$padj < 0.01) {
      '**'
    } else if (df.sig$padj < 0.05) {
      '*'
    }
    
    # Merge ypos and change to 10% above max for this celltype
    df.sig <- merge(x = df.sig, y = coordinates, by.x = 'celltype', by.y = 'celltype')
    df.sig$ypos <- df.sig$max + ((df.sig$max / 100) * 10)
    
    pl <- pl + stat_pvalue_manual(data = df.sig, label = "label",
                                  y.position = 'ypos', tip.length = 0.01)
    
  }
  
  pdf(paste0(outdir, 'dirichlet_plot.pdf'), width = 13, height = 6)
  print(pl)
  dev.off()
}

