require(rtracklayer)
files <- list.files("./",pattern="GSM")
TABLE<-table(unlist(lapply(strsplit(files,"_"),"[",3)))
CL <- names(TABLE)
j=4
cat("\nj=",j,"\n")
files1 <- files[grep(CL[j],files)]
ID_all <-NULL
for (i in c(1:length(files1)))
{
cat(i," ")
bw <- import(files1[i])
ID <- paste(seqnames(bw),start(bw),sep="_")
ID_all <- c(ID_all,setdiff(ID,ID_all))
}
save(file=paste("ID_all_",CL[j],sep=""),ID_all)
