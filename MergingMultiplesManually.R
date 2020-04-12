# This is a script for merging any number of 10X Genomics datasets
# from raw counts with the option of subseting each dataset to a 
# predefined set of barcodes.  It feeds into scRNAseqWorkflow.Rmd


library(SingleCellExperiment)
library(DropletUtils)
library(Matrix)

setwd("D:/Dropbox/GDB/MKlab_v1to2/For Brendan/")

# These are the folder names of each raw dataset.  
# They can be relative or full filepaths.
# Name the folders something meaningful though
# (i.e. rename "mm10" to something useful),
# because these names (last folder in the path)
# will be used as dataset labels.
temp_data_file_paths <- c("../For Brendan/10Day/",
                          "D:/Dropbox/GDB/MKlab_v1to2/For Brendan/14Day/",
                          "L11Day/")
# The folder / dataset names are cleaned to 
# remove anything that isn't a letter or number.
names(temp_data_file_paths) <- gsub("[^A-Za-z0-9]","",
                                    sapply(strsplit(temp_data_file_paths,"/"),
                                           function(X) X[length(X)]))

# These MUST be in the same order as the datasets!
# If you want to keep all cells from all datasets, 
# skip this and the subsetting step below.
temp_barcodes_to_keep <- list(cluster_10Day_barcodes,
                              cluster_14Day_barcodes,
                              cluster_L11Day_barcodes)


# Reading raw data, creating a list of SCE objects.
temp_all_sce <- sapply(temp_data_file_paths,
                       read10xCounts,
                       col.names=T,
                       simplify=F)

# Subsetting each SCE object
temp_all_sce <- mapply(function(SCE,BC) SCE[,BC],
                       SCE=temp_all_sce,
                       BC=temp_barcodes_to_keep)

# Collecting all gene names across datasets
temp_all_genes <- sort(unique(unlist(lapply(temp_all_sce,rownames))))

# Making unified rowData
temp_rowdata <- do.call("rbind",
                        lapply(temp_all_sce,function(Y) 
                          rowData(Y)[,Reduce("intersect",
                                             lapply(temp_all_sce,
                                                    function(X) colnames(rowData(X))))]))
temp_rowdata <- temp_rowdata[!duplicated(rownames(temp_rowdata)),]
temp_rowdata <- temp_rowdata[temp_all_genes,]

# Adding missing genes to each count matrix and 
# generating replacement SCEs with missing genes included
for (i in names(temp_all_sce)) {
  temp_missing_genes <- temp_all_genes[!temp_all_genes %in% rownames(temp_all_sce[[i]])]
  temp_new_counts <- rbind(
    counts(temp_all_sce[[i]]),
    Matrix(0,nrow=length(temp_missing_genes),ncol=ncol(temp_all_sce[[i]]),
           dimnames=list(temp_missing_genes,colnames(temp_all_sce[[i]])))
  )
  temp_new_counts <- temp_new_counts[temp_all_genes,]
  colnames(temp_new_counts) <- paste(i,colnames(temp_new_counts),sep="_")
  temp_sce <- SingleCellExperiment(assays=list(counts=temp_new_counts))
  rowData(temp_sce) <- temp_rowdata
  colData(temp_sce) <- colData(temp_all_sce[[i]])[,Reduce("intersect",
                                                          lapply(temp_all_sce,function(X) 
                                                            colnames(colData(X))))]
  colData(temp_sce)$orig.ident <- i
  temp_all_sce[[i]] <- temp_sce
}
sce <- do.call("cbind",temp_all_sce)

# cleaning up
rm(list=c("i",grep("^temp",ls(),value=T)))

# Now (AND NOT BEFORE) can you launch into the scRNAseqWorkflow.Rmd 
# starting at:
# sce <- sce[rowSums(counts(sce)) > 0,]
