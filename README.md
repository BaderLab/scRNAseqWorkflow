# scRNAseqWorkflow
Brendan's skeleton scRNAseq workflow using scran, Seurat, and scClustViz

This is an RStudio notebook that reflects my opinion of best practices in single-sample processing of scRNAseq data from the 10X Genomics platform.  It is heavily based on the [SimpleSingleCell](https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/tenx.html) and [Seurat](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) tutorials.  Normalization is performed using *scran*, with *Seurat* for clustering.  Clustering is performed iteratively at higher resolutions and stopped when differential expression between clusters is lost, as assessed by [scClustViz](https://baderlab.github.io/scClustViz/) using the wilcoxon rank-sum test.
