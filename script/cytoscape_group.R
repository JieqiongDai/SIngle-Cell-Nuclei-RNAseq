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
required_packages  <- c('RCy3','argparse')
lapply(required_packages,installpkg)

#load libary for cytoscape and dependent packages
library(RCy3)
cytoscapePing ()

library(argparse)
parser <- ArgumentParser()
parser$add_argument('-s',
                    '--style',
                    nargs='+',
                    required=TRUE,
                    help="cytoscape style file")
args <-parser$parse_args()

style <- c(args$style)

#create cytoscape network
nodes <- data.frame(read.table("../inferCNV/node_cyto.txt",sep="\t",header=T,row.names=NULL,check.names = F))
nodes <- nodes[,1:3]
colnames(nodes)<- c("id","cell.number","infor")
nodes[nrow(nodes)+1,] <- c("root",1,"root")
nodes $ cell.number <- as.numeric(nodes $ cell.number)
nodes $ id <- as.character(nodes $ id)
edges <- read.table("CNV.tree2.txt",header=T, stringsAsFactors = F,check.names = F)
edges $ target <- as.character(edges $ target)
edges $ source <- as.character(edges $ source)
edges $ dist <- as.numeric(edges $ dist)
createNetworkFromDataFrames(nodes,edges)
lockNodeDimensions(TRUE)
#switch styles
if ("scCNV_group" %in% getVisualStyleNames() == "FALSE")
importVisualStyles(style)
setVisualStyle('scCNV_group')
#set node size mapping
setNodeSizeMapping('cell.number', c(min(nodes$cell.number),max(nodes$cell.number)), c(10,60),style.name = "scCNV_group")
##select root node
selectNodes ('root','name')
# export pdf file
saveSession('medalt.group')
full.path=paste(getwd(),'medalt.group',sep='/')
exportImage(full.path, 'PDF')
#change to force directed layout
layoutNetwork('force-directed edgeAttribute=dist')
# export pdf file
saveSession('medalt.group.force.directed')
full.path=paste(getwd(),'medalt.group.force.directed',sep='/')
exportImage(full.path, 'PDF')
