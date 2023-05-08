
# HematoMap <img src="https://nrctm-bioinfo.github.io/HematoMap/figure/logo.png" align="right" height=170 width=150/>

`HematoMap` is a comprehensive R package that can provide a visual representation of the cellular hematopoietic hierarchy of human bone marrow. It is a mapping tool that integrates information from single-cell based gene expression data and generate a unified model of hematopoietic cell differentiation. `HematoMap` allows researchers to explore the relationships between different hematopoietic cell types (such as HSC/MPP, CML and CLP) and interpret normal hematopoiesis cellular compositions. 

See the tutorial of HematoMap, please visit [Tutorial of HematoMap](https://nrctm-bioinfo.github.io/HematoMap/index.html).

## Installation

### Installation from GitHub

This requires the `devtools` package to be pre-installed first. The `HematoMap` also required for pre-installation of the [`Seurat`](https://CRAN.R-project.org/package=Seurat) package (https://satijalab.org/seurat/). 

``` {r eval = FALSE}

# If Seurat is not already installed, you need to install Seurat first
install.packages("Seurat") 

# If devtools is not already installed, you need to install devtools first
install.packages("devtools") 

```

To install HematoMap

``` {r eval = FALSE}
# Installation of HematoMap
devtools::install_github("NRCTM-bioinfo/HematoMap")

# library HematoMap
library(HematoMap)

```

The source code of `HematoMap` on GitHub can be accessed at https://github.com/NRCTM-bioinfo/HematoMap.


## Quick start

### Load package


``` {r eval = FALSE, out.width='60%', fig.asp=1, fig.align='center'}

# library the Seurat and HematoMap packages
library(Seurat)
library(HematoMap)

```


### Show normal BMMCs

``` {r eval = FALSE, fig.width = 6, fig.height = 6, fig.asp=1, fig.align='center'}

# Show circle plot of normal BMMC
plotCircleTree(group.subc = "reference", 
               color.mapping = "cell.type", 
               title = "Normal BMMCs")


```

``` {r eval = FALSE, fig.width = 5, fig.height = 6, fig.asp=1, fig.align='center'}
# Show cluster-tree plot of normal BMMC
plotClusterTree(group.subc = "reference",
                color.mapping = "cell.type",
                title = "Normal BMMCs",
                point.size = 20, 
                label.cell = TRUE)


```

### Quick-start code


``` {r eval = FALSE, fig.width = 6, fig.height = 6, fig.asp=1, fig.align='center'}

# For scRNA-seq analysis
# the `scrna` is a seurat object
hmap <- HematoMap::CreateSubclusterObject(scrna)

# For bulk RNA-seq analysis
# the `bulk` is a expression matrix of bulk RNA-seq
hmap.bulk <- computeLassoScore(bulk)

```

## Run shiny app

``` {r eval = FALSE, out.width='50%', fig.asp=1, fig.align='center'}

# Run shiny app
HematoMap::runApp()

```


