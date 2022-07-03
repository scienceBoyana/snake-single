args <- commandArgs(trailingOnly = TRUE)

## hardcoded paramters:
min_cells_per_gene <- 25
max_doublet_score <- .3
leiden_resolution <- .6

## for UMAP
n_neighbors <- 25
min_dist <- .15

## testing commandArgs
## args <- list(
##     '../results/mouse/cerebellum/SN104/full_sce.rds',
##     '../results/mouse/cerebellum/SN044/'
## )

sce_input_path <- args[[1]]
sce_out <- args[[2]]
umap_plot_out <- args[[3]]
gene_score_out <- args[[4]]

stopifnot(file.exists(sce_input_path))


## required packages:
suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(Matrix)
    library(sparseMatrixStats)
    library(leidenAlg)
    library(igraph)
    source('scripts/R/helper_functions.R')
})


sce <- readRDS(sce_input_path)
rownames(sce) <- rowData(sce)[,1]
sce <- sce[ rowSums(counts(sce)>0) > min_cells_per_gene, ]

## caluclating doublet score
## based on 100 nearest neighbors and simulated doublets
## adapted from doubletFinder
scr <- calc_doublet_score( counts(sce) )

sce$doublet_score <- scr

sce <- sce[,sce$doublet_score <= max_doublet_score]

# counts per million
#sce <- sce[ , sce$frac_doublet_nn < .3]
normcounts(sce) <- t(t(counts(sce)) / colSums(counts(sce))) * 1e6


## finding overdispersed genes
hvg <- get.info.genes( counts(sce), vmr_factor = 1.2 )
rowData(sce)$is_hvg <- rownames(sce) %in% hvg

## running PCA
set.seed(1234)
pca <- irlba::prcomp_irlba( t(log1p( normcounts(sce)[hvg,] )),
                            center = TRUE,
                            scale. = TRUE,
                            n = 50)
rownames(pca$x) <- colnames(sce)
reducedDim(sce, 'pca') <- pca$x

## since we are also interested in some PCs after the elbow, I hardcoded the 10 following
## PCs after the elbow to use as well.
use_max_pca <- find_inflection(1:length(pca$sdev), pca$sdev, type = 'elbow', do.plot = TRUE) + 10

## running UMAP and returning the approx. nearest neighbor graph
set.seed(1234)
umap <- uwot::umap( reducedDim(sce,'pca')[,1:use_max_pca],
                   metric = 'cosine',
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   ret_nn = TRUE )

## saving the UMAP to the sce object
reducedDim(sce, 'umap') <- umap$embedding

## formatting the nn matrix to use with igraph
ij <- umap$nn$cosine$idx %>%
    apply(1, function(x) {
        cbind(rep(x[1], length(x)), x)
    }, simplify = FALSE) %>%
    do.call(rbind,.)

## generating the graph sparse matrix
nn <- sparseMatrix(i=ij[,1], j=ij[,2])

## running the leiden algorithm
leiden <- leidenAlg::leiden.community( graph_from_adjacency_matrix(nn),
                                      resolution = leiden_resolution)

sce$clusters <- str_pad(leiden$membership, width = 2, side = 'left', pad = '0')


plt <- ggplot(NULL,
       aes( x = umap$embedding[,1],
            y = umap$embedding[,2],
            color = leiden$membership )) +
    geom_point(size=.25) +
        coord_fixed() +
  labs(
    title = paste0('UMAP - Leiden Clusters: ', leiden_resolution),
    x = 'UMAP1',
    y = 'UMAP2'
  )

ggsave(plot = plt, filename = umap_plot_out)


## running gene enrichment per cluster
smry <- score_genes( counts(sce), sce$clusters )

smry %>%
  write_csv(gene_score_out)

saveRDS(sce, sce_out)
