suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(Matrix)
})

## column variance of sparse matrix
colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

## row variances of sparse matrix
rowVars_spm <- function( spm ) {
  colVars_spm( t(spm) )
}

## finding informative genes based on overdispersion
get.info.genes <- function( counts, vmr_factor = 1.25 ) {
  library(Matrix)
  stopifnot(!is.null(dim(counts)))
  counts <- counts[ rowSums(counts > 0) > 0, ]
  norm_counts <- t(t(counts) / Matrix::colSums(counts))

  gene_means <- rowMeans( norm_counts )
  gene_vars <- rowVars_spm( norm_counts )

  counts <- counts[!is.na(gene_vars), ]
  poisson_vmr <- mean( 1 / Matrix::colSums( counts ) )

  informative_genes <- names(which(
    gene_vars / gene_means  >  vmr_factor * poisson_vmr ))

  return(
    informative_genes
  )

}


## knee/elbow detection
## uses a simple principle to find elbow / kneepoints:
## connect the first and last point of the curve using a linear function.
## Then the distance of each observation to this line is calculated.
## The inflection point has the highest distance to the line.
find_inflection <- function(x, y, type=c('kneepoint', 'elbow'), do.plot=FALSE){

  type <- match.arg(type)

  if( type == 'kneepoint'){
    x1 <- min(x)
    y1 <- min(y)
    x2 <- max(x)
    y2 <- max(y)
  } else {
    x1 <- min(x)
    y1 <- max(y)
    x2 <- max(x)
    y2 <- min(y)
  }

  dst <- apply(cbind(x,y), 1, function(z) {
    abs( (x2-x1) * (y1 - z[2]) - (x1 - z[1]) * (y2 - y1) ) / sqrt( (x2-x1)^2 + (y2-y1)^2 )
  })

  if( do.plot ){
    y_scaled <- y - min(y)
    y_scaled <- y_scaled / max(y_scaled)
    plot(x, y_scaled)
    dst_scaled <- dst - min(dst)
    dst_scaled <- dst_scaled / max(dst_scaled)
    points(x, dst_scaled, col='red')
  }

  return( x[dst == max(dst)] )

}

## Probabalistic rounding of floats
##
## scalar function
## principle: the decimal part of a float is the probability of
## rounding up.
## e.g. 1.25 has a .25 chance of being rounded to 2
## and a .75 chance of being rounded to 1
prob_round_scalar <- function( x ){
    x_prob <- x - floor(x)
    x <- floor(x) + sample(c(0,1), 1, prob = c(1-x_prob, x_prob))
    return(x)
}

## Vectorized version of the prob_round_scalar function
prob_round <- Vectorize(prob_round_scalar)


## Calculates doublet score per cell
##
## inspired and partly taken from: https://github.com/chris-mcginnis-ucsf/DoubletFinder
## whom I'd like to thank for the elegant approach
##
## I changed the nearest neighbor identification by using the ANNOY
## NN approximation implemented in uwot::umap function.
##
## Parameters:
##   umi = [sparse matrix] of UMI counts, genes x cells
##   perc_sim = [float] fraction of original number of cells to simulate doublets
##   n_neighbors = [int] number of neighbors to consider
calc_doublet_score <- function( umi, perc_sim = .5, n_neighbors = 100 ){

    ## doublet generation
    cells <- colnames(umi)
    n_doublets <- floor( length(cells) * perc_sim )

    dplt1 <- sample(cells, n_doublets, replace = TRUE)
    dplt2 <- sample(cells, n_doublets, replace = TRUE)

    ## merging two random cells and dividing each UMI by 2
    doublets <- (umi[,dplt1] + umi[,dplt2]) / 2

    ## probabalistic rounding to recover integer
    doublets@x <- prob_round( doublets@x )
    colnames(doublets) <- paste0('D_', colnames(doublets))

    ## merging the matrices
    orig_sim <- cbind( umi, doublets )


    ## and keeping track of the simulated doublets
    ## using a prefix
    is_sim_doublet <- startsWith(colnames(orig_sim), 'D_')

    ## highly variable gene identification
    sim_i_genes <- get.info.genes(orig_sim, vmr_factor = 1.25)
    ## normalization
    sim_nrm <- log1p(t(orig_sim[sim_i_genes,]) / colSums(orig_sim) * 1e6)
    ## running PCA with 15 components
    sim_pca <- irlba::prcomp_irlba(sim_nrm, n = 15, center = TRUE, scale. = TRUE)

    ## and UMAP
    ## this is mainly done to access the ANNOY KNN approximation
    ## (quicker than implementing on my own)
    sim_umap <- uwot::umap( sim_pca$x,
                        metric = 'cosine',
                        n_neighbors = 100,
                        min_dist = .15,
                        ret_nn = TRUE )

    ## nearest neighbor indices to sparse matrix
    sim_nn <- lapply(1:nrow(sim_umap$nn$cosine$idx), function(i){
        idx <- sim_umap$nn$cosine$idx[i,]
        cbind( rep(i, length(idx)), idx)
    }) %>%
        do.call(rbind,.)
    sim_nn <- sparseMatrix( i = sim_nn[,1], j = sim_nn[,2] )
    rownames(sim_nn) <- rownames(sim_nrm)

    ## using the dot product to calculate the fraction
    ## of simulated doublets of NN per cell
    frac_doublets <- as.numeric( sim_nn %*% is_sim_doublet ) / ncol(sim_umap$nn$cosine$idx)
    frac_doublets_cells <- frac_doublets[!is_sim_doublet]

    ## returning a vector per cell
    return( frac_doublets_cells )
}

## score genes for enrichment analysis
## inspired by the batchelor package, quickMarkers function
## https://bioconductor.org/packages/release/bioc/html/batchelor.html
##
## PARAMETERS:
##   x = [sparse matrix] UMI matrix, genes x cells
##   clusters = [vector] vector of cluster labels, order of x columns
##   max_fdr = [float, (0, 1] ] maximum false discovery rate to report per cluster
##   min_enrichment = [float] minimal enrichment of a gene within a cluster to report
score_genes <- function(x, clusters, max_fdr = .01, min_enrichment = 1.25){
    group_counts <- table(clusters)
    n_feature <- tapply( seq(ncol(x)), clusters, function(i) {
        rowSums( x[,i,drop=FALSE] > 0 )
    }) %>%
        do.call(cbind, .)

    n_tot <- rowSums(x > 0)
    tf <- t(t(n_feature) / colSums(n_feature))
    out_tf <- t(t(rowSums(n_feature) - n_feature) / colSums(rowSums(n_feature)-n_feature))
    idf <- log( ncol(x) / n_tot )
    tfidf <- tf*idf

    enrich <- tf / out_tf

    pval <- lapply( colnames(tfidf), function(i) {
        phyper(
            n_feature[,i] - 1,
            n_tot,
            ncol(x) - n_tot,
            group_counts[ i ],
            lower.tail = FALSE
        )
    }) %>%
        do.call(cbind,.)
    colnames(pval) <- colnames(tfidf)

    smry <- tf %>%
        as.data.frame() %>%
        rownames_to_column('gene') %>%
        gather('cluster', 'tf', -gene) %>%
        left_join(
            out_tf %>%
            as.data.frame() %>%
            rownames_to_column('gene') %>%
            gather('cluster', 'out_tf', -gene)
        ) %>%
        left_join(
            tfidf %>%
            as.data.frame() %>%
            rownames_to_column('gene') %>%
            gather('cluster', 'tfidf', -gene)
        ) %>%
        left_join(
            pval %>%
            as.data.frame() %>%
            rownames_to_column('gene') %>%
            gather('cluster', 'pval', -gene)
        ) %>%
        left_join(
            enrich %>%
            as.data.frame() %>%
            rownames_to_column('gene') %>%
            gather('cluster', 'enrichment', -gene)
        ) %>%
        mutate(FDR = p.adjust(pval, method = 'BH'))  %>%
        filter(FDR < max_fdr) %>%
        filter(enrichment > min_enrichment)

    return(smry)
}
