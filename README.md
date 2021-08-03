# sc_proportions

Dirichlet regression to test the abundance of cell populations between conditions in single-cell RNA sequencing data.

## Usage

- Pre-process single-cell RNA sequencing data using Seurat as normal. Refer to the Seurat vignettes for these steps. 
```
rm(list = ls())
data <- readRDS('inst/extdata/test_data_HDV.RDS')
```

- Specify the celltypes and conditions to test

```
celltypes <- unique(data$celltype)
conditions <- c('healthy control', 'HBV control', 'HDV neg', 'HDV high')
```
- Optionally, specify the Idents in which the required information is stored (e.g. the colnames in the Seurat object)
```
condition_ident = 'condition'
sample_ident = 'sample'
celltype_ident = 'celltype'
```

- Run the main function `sc_prop`
```
sc_prop(seurat_object   = data,
        conditions      = conditions,
        celltypes       = celltypes,
        condition_ident = condition_ident,
        sample_ident    = sample_ident,
        celltype_ident  = celltype_ident,
        outdir          = './output_dir/')
```

## Output

Visualization is made using ggplot2. Significance stars are added when the comparison is significant after multiple testing correction

![output](https://github.com/MZoodsma/sc_proportions/blob/61c706fb35d36a82d33d62df443a679f049b5f35/inst/images/dirichlet.png)

Furthermore, several output files are written:
1. The Dirichlet regression results, comparing each celltype within each condition
2. The data used to generate the figure
