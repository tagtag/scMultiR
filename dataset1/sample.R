#Suppose there are GSMXXXXXXX_DNA_XXX_scXX.GCH.bw files and chr.lst in ./GCH/ directory
#in ./GCH/
require(rtracklayer)
require(Matrix)
bw <- import("GSM4678346_DNA_GO1_sc1.GCH.bw")
chr.lst<- read.csv("chr.lst",sep="\t",header=F)
chr<-seqnames(bw)@values
ID_NDR <- NULL
for (i in c(1:length(chr)))
{
    cat(i, " ")
    breaks <- seq(1,chr.lst[chr.lst[,1]==chr[i],2],by=200)
    ID_NDR <- rbind(ID_NDR,data.frame(chr[i],breaks[1:(length(breaks)-1)]))
}
save(file="ID_NDR",ID_NDR)
#------------------- 
#run batch_1.R to batch_9.R in ./GCH/ directory
#run batch_1p.R to batch_9p.R in ./GCH/ directory
#run batch_1pp.R to batch_9pp.R in ./GCH/ directory
#run batch_1.R to batch_9.R in ./WCG/ directory
#run batch_1p.R to batch_9p.R in ./WCG/ directory
#-----
ID_all_tot <-NULL
files <- list.files("./WCG",pattern="ID_all")
ID_all_LST<-rep(list(NA),9)
for (i in c(1:length(files)))
{
    cat(i," ")
    load(paste("./WCG/",files[i],sep=""))
    ID_all_tot <- union(ID_all_tot,ID_all)
    ID_all_LST[[i]] <- ID_all
}
index<-rep(list(NA),9)
for (i in c(1:length(files)))
{
    cat(i," ")
    index[[i]] <- match(ID_all_LST[[i]],ID_all_tot)
}
files <- list.files("./WCG",pattern="sm_all")
i1 <- NULL
j1<-NULL
x1<- NULL
j0<-0
for(i in c(1:length(files)))
{
    cat(i, " ")
    load(paste("./WCG/",files[i],sep=""))
    i1 <- c(i1,index[[i]][summary(sm_all)$i])
    j1 <- c(j1,summary(sm_all)$j+j0)
    x1 <- c(x1,summary(sm_all)$x)
    j0 <- j0+dim(sm_all)[2] 
}
sm_all <- sparseMatrix(i = i1, j = j1, x = x1, 
                       dims = c(length(ID_all_tot),j0))
save(file="./WCG/sm_all_tot",sm_all)
#---
files <- list.files("./GCH/",pattern="SPLT_all")
SPLT_all_tot<- NULL
for (i in c(1:length(files)))
{
    cat(i," ")
    load(paste("./GCH/",files[i],sep=""))
    SPLT_all_tot <- cbind(SPLT_all_tot,SPLT_all)
}
save(file="./GCH/SPLT_all_tot",SPLT_all_tot)
#--------------------------
#Suppose there is GSE154762_hO_scChaRM_count_matix.txt.gz in current directory
load("./WCG/sm_all_tot")
load("./GCH/SPLT_all_tot")
require(Matrix)
require(rTensor)
require(irlba)
sm_all@x <- 2*sm_all@x-1
A <- t(t(SPLT_all_tot)/colSums(abs(SPLT_all_tot)))
B <- t(t(sm_all)/colSums(abs(sm_all)))
SVD_WCG <- irlba(B,10)
SVD_GCH <- irlba(A,10)
x <- read.csv("GSE154762_hO_scChaRM_count_matix.txt.gz",sep="\t")
SVD <- svd(scale(x[,-1]))
sm_all_10 <- t(SVD_WCG$u) %*% B
SPLT_all_10 <-  t(SVD_GCH$u) %*% A
x_10 <- t(SVD$u[,1:10]) %*% scale(x[,-1])
Z <- array(NA,c(10,dim(SPLT_all_10)[2],3))
Z[,,1] <- data.matrix(sm_all_10) #methylation
Z[,,2] <- data.matrix(SPLT_all_10) #accessibility
Z[,,3] <- data.matrix(x_10) #gene expression
require(rTensor)
HOSVD <- hosvd(as.tensor(Z))
#---------------
files <- list.files("./GCH/",pattern="GSM")
class <- unlist(lapply(strsplit(files,"_"),"[",3))
LM <- lm(HOSVD$U[[2]][,1:30]~class)
SLM <- summary(LM)
fs <- t(data.frame(lapply(SLM,"[",10)))
P <- pf(fs[,1],fs[,2],fs[,3],lower.tail=F)
table(p.adjust(P,"BH")<0.01)
#--------------------- 
apply(HOSVD$Z@data[,1:30,][,p.adjust(P,"BH")<0.01,]^2,c(1,3),sum)
#-----------------------
P <- pchisq(scale(SVD$u[,1:10] %*% HOSVD$U[[1]][,1])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
data.frame(x[p.adjust(P,"BH")<0.01,1])
#------------------- 
require(umap)
custom.config <- umap.defaults
custom.config$n_neighbors<-100
UMAP <- umap(HOSVD$U[[2]][,1:30],config=custom.config)
plot(UMAP$layout,col=as.numeric(factor(class)),pch=as.numeric(factor(class)))
legend(-4,-2,c("FGO","GO1","GO2","Granulosa","Immune","MI","MII","StromaC1","StromaC2"), col=1:9,pch=1:9)
