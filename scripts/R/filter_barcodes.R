## how to filter
use_kneepoint = TRUE
use_frac_intronic = TRUE

## parsing arguments
args <- commandArgs(trailingOnly = TRUE)

## for testing
## args <- list(
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/Gene/raw/matrix.mtx',
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/Gene/raw/features.tsv',
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/Gene/raw/barcodes.tsv',
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/GeneFull/raw/matrix.mtx',
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/GeneFull/raw/features.tsv',
##     '../Solo.out/mouse_cerebellum_SN044/Solo.out/GeneFull/raw/barcodes.tsv',
##     '50',  # minimum UMI per cell
##     '25',   # minimum cells per gene
##     'SN044',  # batch name
##     'tmp/results'  # results folder path
## )

## checking for presence of input
for( i in args[1:6] ) stopifnot(file.exists(i))

## sorting the input paths
mm_paths <- list(
    exonic = list(
        mm_path = args[[1]],
        features = args[[2]],
        barcodes = args[[3]]
    ),
    full = list(
        mm_path = args[[4]],
        features = args[[5]],
        barcodes = args[[6]]
    )
)

## parsing the user specified parameters
min_umi <- as.numeric(args[[7]])
min_cells <- as.numeric(args[[8]])
batch_name <- args[[9]]

## and defining the output path
output_path <- args[[10]]
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## needed packages
suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(Matrix)
    library(mclust)
})

## reading the data into SCE objects
mm <- lapply(mm_paths, function(pths){
    umi <- readMM(pths$mm_path)

    features <- read_tsv(pths$features, col_names = FALSE)

    rownames(umi) <- features %>% pull(1)
    colnames(umi) <- read_tsv(pths$barcodes, col_names = FALSE) %>% pull(1)

    sce <- SingleCellExperiment(
        list(counts = umi),
        rowData = features %>% as.data.frame() %>% column_to_rownames("X1")
    )

    sce <- sce[ , colSums(umi) > min_umi ]
    sce <- sce[ rowSums(umi>0) > min_cells , ]
    sce$batch <- batch_name

    colnames(sce) <- paste(batch_name, colnames(sce), sep = '_')

    return(sce)
})

# pairing both counting modes:
shared_cells <- intersect( colnames(mm$exonic), colnames(mm$full) )

stopifnot( length(shared_cells) > 0 )

mm <- lapply(mm, function(m) m[,shared_cells])

## preparation for proportion of intronic UMI plot
n_full = colSums( counts(mm$full) )
n_exonic = colSums( counts(mm$exonic) )

## and for the kneepoint
mdata <- tibble(
    cell_id = shared_cells
) %>%
    add_column( prop_intronic = (n_full - n_exonic) / n_full ) %>%
    mutate( prop_intronic = ifelse(prop_intronic < 0, 0, prop_intronic) ) %>%
    mutate( exonic = n_exonic ) %>%
    mutate( full = n_full ) %>%
    mutate( rank = rank(-full, ties.method = 'first') ) %>%
    mutate( rank = rank - 1 )

## scaling the cummulative UMI sum between 0 and 1
mdata <- mdata %>%
    arrange(rank) %>%
    mutate( umi_cs = cumsum(full) ) %>%
    mutate( umi_cs = umi_cs - min(umi_cs) ) %>%
    mutate( umi_cs = umi_cs / max(umi_cs) )

## kneepoint detection
## this uses the idea that the kneepoint is furthest away
## from a linear function which connects the first and last cell
md <- mdata %>%
    select('rank', 'umi_cs') %>%
    as.data.frame() %>%
    as.matrix()

colnames(md) <- NULL

x1 <- min(md[,1])
y1 <- min(md[,2])
x2 <- max(md[,1])
y2 <- max(md[,2])

## distance calculation
dst <- apply(md, 1, function(x) {
   abs( (x2-x1) * (y1 - x[2]) - (x1 - x[1]) * (y2 - y1) ) / sqrt( (x2-x1)^2 + (y2-y1)^2 )
})

## and scaling of the distance
## for better plotting
dst <- dst - min(dst)
dst <- dst / max(dst)

## and this gives us the kneepoint
kneepoint <- which( dst == max(dst) )

pdf(file.path(output_path, 'kneepoint.pdf'))
plot(md[,2], type ='l',
    xlab = 'cell UMI rank',
    ylab = 'scaled values',
    main = 'Barcode detection')
points(dst, pch = '.', col = 'red')
abline( v = kneepoint, lty = 'dashed')
abline( a = 0, b = y2/x2, lty = 'dotted')
legend(x = .3*x2,
       y = .25,
       col = c('black', 'red', 'black'),
       lty = c('solid', 'solid','dashed'),
       legend = c('UMI cummulative sum',
                  'distance to linear function',
                  'kneepoint'))
dev.off()

## when finding nuclei containing droplets,
## we can use the proportion of intronic UMI to
## cluster the barcodes into two clusters
## just look at the plot and you can clearly see the difference.
mdata <- mdata %>%
    mutate(intronic_cluster = as.factor(Mclust(prop_intronic, G = 2)$classification))

## using the fraction of intronic UMIa
theme_set(theme_classic())
plt <- mdata %>%
    ggplot(aes(x = rank,
               y = prop_intronic,
              color = intronic_cluster,
              shape = rank <= kneepoint)) +
    geom_point() +
    guides( color = guide_legend( override.aes = list(size=2) ) ) +
    labs(
        title = 'nuclei identification',
        subtitle = 'using intronic UMI fraction',
        x = 'cell UMI rank',
        y = 'fraction of intronic UMI'
    )

ggsave(plot = plt, file.path(output_path, 'kneepoint_frac_intronic_umi.pdf'))


## as set at the top of this script
## we can filter based on the kneepoint and/or
## on the fraction of intronic UMI classes.
if( use_kneepoint ) mdata <- mdata %>% filter( rank <= kneepoint )
if( use_frac_intronic ) mdata <- mdata %>% filter( intronic_cluster == 2 )


## data sanatization and saving
## saves:
##   SCE objects
##   matrix market matrices
use_cells <- mdata %>%
    pull('cell_id')

mm <- lapply(mm, function(m) m[ , use_cells] )

for( n in names(mm) ){
    m <- mm[[n]]
    saveRDS( m , file.path(output_path, paste0(n, '_sce.rds')) )

    u <- counts(m)

    tmp <- file.path(output_path, paste0(n, '_counts'))
    dir.create( tmp, showWarnings = FALSE, recursive = TRUE)

    writeMM( u, file.path(tmp, 'matrix.mtx') )

    write_lines( colnames(u), file.path(tmp, 'barcodes.tsv') )
    write_lines( rownames(u), file.path(tmp, 'features.tsv') )
}
