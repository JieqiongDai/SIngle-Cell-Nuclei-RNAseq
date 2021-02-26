#get filtered counts matrix from cellranger
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
required_packages  <- c('Matrix')
lapply(required_packages,installpkg)

library(Matrix)
mat <- readMM("filtered_feature_bc_matrix/matrix.mtx.gz")
feature.names = read.delim("filtered_feature_bc_matrix/features.tsv.gz", header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim("filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE,stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
mat2 <- as.data.frame(as.matrix(mat))
mat2$gene_id <- feature.names$V1
write.csv(mat2,"filtered_feature_bc_matrix.csv")
