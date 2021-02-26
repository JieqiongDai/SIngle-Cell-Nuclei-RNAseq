# merge tumor noraml seurat objects for the same patient

# Check required packages
# function for installing needed packages
installpkg <- function(x){
    r = getOption("repos")
    r["CRAN"] = "http://cran.us.r-project.org"
    options(repos = r)
    if(x %in% rownames(installed.packages())==FALSE) {
        if(x %in% rownames(available.packages())==FALSE) {
            if (!requireNamespace("BiocManager", quietly = TRUE))
               install.packages("BiocManager")
               BiocManager::install(x,ask = FALSE)
        } else {
            install.packages(x)
        }
        paste(x, "package is installed...")
    } else {
        paste(x,"package already installed...")
    }
}

# install necessary packages
required_packages  <- c('Seurat')
lapply(required_packages,installpkg)

# Read in tumor seurat data
library(Seurat)
load(snakemake@input[[1]])
seurat_t <- seurat2
seurat_t
sample_t <- data.frame(rep("T",ncol(seurat_t)))
dim(sample_t)
rownames(sample_t)<- colnames(seurat_t)
colnames(sample_t) <- c("sample")

# Read in normal seurat data
load(snakemake@input[[2]])
seurat_n <- seurat2
seurat_n
sample_n <- data.frame(rep("N",ncol(seurat_n)))
dim(sample_n)
rownames(sample_n)<- colnames(seurat_n)
colnames(sample_n) <- c("sample")

# Merge two seurat objects 
merge <- merge(x=seurat_t,y=seurat_n,add.cell.ids=c("T","N"),merge.data=T)
merge
sample <- rbind(sample_t,sample_n)
dim(sample)
merge@meta.data$sample.type <- sample$sample
write.csv(merge[["SCT"]]@counts,snakemake@output[[1]],quote=F)
write.csv(merge@meta.data,snakemake@output[[2]],quote=F)
rownames(sample)<- colnames(merge)
write.csv(sample,snakemake@output[[3]],quote=F)
save.image(snakemake@output[[4]])

