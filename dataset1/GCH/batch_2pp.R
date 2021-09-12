require(rtracklayer)
require(Matrix)
bw <- import("GSM4678346_DNA_GO1_sc1.GCH.bw")
chr<-seqnames(bw)@values
load("ID_NDR")
files <- list.files("./",pattern="GSM")
TABLE<-table(unlist(lapply(strsplit(files,"_"),"[",3)))
CL <- names(TABLE)
j<-2
cat("j=",j, "\n")
load(paste("ID_all",CL[j],sep="_"))
load(paste("sm_all",CL[j],sep="_"))
SPLT_all <- NULL
for (i in c(1:length(chr)))
{
    cat("i=",i,"\n")
    index <- grep(paste(chr[i],"_",sep=""),ID_all)
    num <- gsub(paste(chr[i],"_",sep=""),"",ID_all[index])
    breaks <- ID_NDR[ID_NDR[,1]==chr[i],2]
    breaks <- c(breaks,breaks[length(breaks)]+200)
    num <- as.numeric(num)
    CUT <-cut(num,breaks)
    SPLT_temp<-NULL
    for (k in c(1:dim(sm_all)[2]))
    {
    cat(k, " ")
    SPLT <- split(sm_all[index,k],CUT)
    SPLT_sum <- lapply(SPLT,sum)
    SPLT_sum <- unlist(SPLT_sum)
    ix <- (1:length(SPLT_sum))[SPLT_sum!=0]
    value <- SPLT_sum[ix]
    iy <- rep(1,length(ix))
    sm <- sparseMatrix(ix,iy,x=value,dims=c(length(breaks)-1,1))
    SPLT_temp <- cbind(SPLT_temp,sm)
    }
    SPLT_all <- rbind(SPLT_all,SPLT_temp)
}
save(file=paste("SPLT_all",CL[j],sep="_"),SPLT_all)
