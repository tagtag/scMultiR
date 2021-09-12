require(rtracklayer)
require(Matrix)
files <- list.files("./",pattern="GSM")
TABLE<-table(unlist(lapply(strsplit(files,"_"),"[",3)))
CL <- names(TABLE)
j=7
cat("\nj=",j,"\n")
files1 <- files[grep(CL[j],files)]
load(paste("ID_all_",CL[j],sep=""))
sm_all <- NULL
for (i in c(1:length(files1)))
{
cat(i," ")
bw <- import(files1[i])
ID <- paste(seqnames(bw),start(bw),sep="_")
ix <- match(ID,ID_all)
iy <- rep(1,length(ix))
value <- score(bw)
sm <- sparseMatrix(ix,iy,x=value,dims=c(length(ID_all),1))
sm_all <- cbind(sm_all,sm)
}
save(file=paste("sm_all_",CL[j],sep=""),sm_all)
