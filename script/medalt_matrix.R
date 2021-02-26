# read in tumor cnv matrix
tu_mt <- read.table("infercnv.observations.txt",header=T,row.names=1,check.names = F)

# read in tumor group table
tu_gp <- read.table("infercnv.observation_groupings.txt",header=T,row.names=1,check.names = F)

# get average cnv data of each tumor group
gp <- unique(tu_gp[,1])
for (i in gp)
{
  cell <- rownames(tu_gp[tu_gp==i,])
  a <- tu_mt[,colnames(tu_mt)%in%cell,drop=F]
  n <- ncol(a)
  t <- data.frame("name"=i,"cellnumber" =n, "infor"=paste0(i,"-",n))
  b <- as.data.frame(rowMeans(a))
  colnames(b) <- i
  if(!exists("tu_mt_gp"))
  {tu_mt_gp = b} else if (exists("tu_mt_gp"))
  {tu_mt_gp <- cbind(tu_mt_gp,b)}
  if(!exists("node"))
  {node = t} else if (exists("node"))
  {node <- rbind(node,t)}
}

# Read in normal cnv data and get average 
nor_mt <- read.table("infercnv.references.txt",header=T,row.names=1,check.names = F)
nor_mt_gp <- as.data.frame(rowMeans(nor_mt))
colnames(nor_mt_gp) <-c("0")
t2 <- data.frame("name"=c("0"),"cellnumber" =ncol(nor_mt), "infor"=paste0("0","-",ncol(nor_mt)))

#merge tumor and normal groupped cnv matrix
cnv <- cbind(nor_mt_gp,tu_mt_gp)
cnv <- as.matrix(cnv)
write.table(cnv,"scRNA.CNV.txt",quote=F,sep = "\t")
node2 <- rbind(t2,node)
write.table(node2,"node_cyto.txt",quote=F,sep = "\t",row.names=F)
