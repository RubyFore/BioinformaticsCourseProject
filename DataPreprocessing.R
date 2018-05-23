


op <- par()

source("~/School/PHP2620_Bioinformatics/FinalProject/loadData.R")

#### EXPLORATORY ANALYSIS ####

# barplot of lib size by tissue response type
LibSize = colSums(counts) 
plot(LibSize~colData(rse_gene)$responseType)



# outlier in Library size is the very first sample
which(LibSize==min(LibSize))

# check for batch effect in samples
plot(LibSize~colDat$index)
longCounts <- gather(as.data.frame(counts))
longCounts$sampleID <- colDat$index[match(longCounts$key,colDat$run)]
# Individual sample barplot
# individual samples look pretty good
plot(log2(longCounts$value+1)~as.factor(longCounts$sampleID), main="Batch effect looking at individual samples", 
     xlab="Sample Index", ylab="log2 range of counts")

# summing the number of zero counts for each sample
zeros <- c()
for (i in 1:20){
  zeros[[i]] <-sum(counts[,i]==0)
}
#hist(zeros)
percentzeros <- zeros/58037
nonzeros <- 58037-zeros
colors<- c(rep('red', 8), rep('green', 8), rep('blue', 4))
barplot(zeros, names.arg=paste0(round(percentzeros, digits=2)*100, "%"), legend.text=c('normal', 'resistant', 'sensitive'),col=colors, 
        main="Percent 0 count by sample", ylab="Number of samples with zero count", xlab = "% of total",
        las=2, ylim=c(0,50000), args.legend = list(fill=c('red', 'green', 'blue')))


# Average expression ratios comparing the three groups 
mNormal<- rowMeans(counts[,which(colData(rse_gene)$responseType=="normal")])
mResistant<-rowMeans(counts[,which(colData(rse_gene)$responseType=="resistance")])
mSensitive <-rowMeans(counts[,which(colData(rse_gene)$responseType=='sensitive')]) 

hist(log2(mResistant/mNormal), breaks= 30, xlim = c(-15,15))
hist(log2(mSensitive/mNormal), breaks=30,xlim=c(-15,15))
hist(log2(mSensitive/mResistant), breaks=30,xlim=c(-15,15)) 
# They look mostly the same

# MA plots for averages across groups? Can take from hw2

layout(matrix(c(1,2,3), nrow=1, ncol=3))
a1 <- (log2(mNormal+1)+log2(mResistant+1))/2
m1 <- log2(mNormal+1) - log2(mResistant+1)
a2 <- (log2(mResistant+1)+log2(mSensitive+1))/2
m2 <- log2(mResistant+1)-log2(mSensitive+1)
a3 <- (log2(mNormal+1) +log2(mSensitive+1))/2
m3 <- log2(mNormal+1)-log2(mSensitive+1)
plot(a1,m1,pch=16,cex=.3, main='MA plot of normal vs resistant')
abline(h=0,col='red',lwd=1) 
plot(a2,m2, pch=16, cex=.3, main="MA plot of resistant vs. sensitive")
abline(h=0,col='red',lwd=1) 
plot(a3,m3, pch=16, cex=.3, main="MA plot of normal vs. sensitive")
abline(h=0,col='red',lwd=1)


#### GETTING GENE SYMBOLS, NAMES ####
Ids <-str_split_fixed(rowData(rse_gene)$gene_id, "\\.", n=2)[,1]
geneNames <- mapIds(EnsDb.Hsapiens.v86, keys=Ids, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
lengths <- lapply(geneNames, length)
multipleMask<- lengths>1
namedGenes <- geneNames[multipleMask]
singleIds <- Ids[multipleMask]

#### NORMALIZATION ####
# want to use cqn. But requires read length and gc content. 
#gcLengthInfo <- getGeneLengthAndGCContent(Ids,'hsa', mode='biomart')
#df <- as.data.frame(gcLengthInfo)
# https://support.bioconductor.org/p/58846/ page with hack for how to do this. 
#saveRDS(df, file="/Users/rubyfore/School/PHP2620_Bioinformatics/FinalProject/gcLengthInfo.Rdata")
gcLengthInfo <- readRDS("/Users/rubyfore/School/PHP2620_Bioinformatics/FinalProject/gcLengthInfo.Rdata")

zerolengthM <- !is.na(gcLengthInfo$length)
cqnObj <- cqn(counts[zerolengthM,], x=as.data.frame(gcLengthInfo)$gc[zerolengthM], lengths=as.data.frame(gcLengthInfo)$length[zerolengthM], sizeFactors=NULL, verbose=F)


#matching lengt and gc info to cqn obj
# need to trim version number off cqnObj rownames
cqnIds <-str_split_fixed(rownames(cqnObj$counts), "\\.", n=2)[,1]
cqnGCLengthInfo<- subset(gcLengthInfo, rownames(gcLengthInfo)%in%cqnIds)
colors<- c(rep('red', 8), rep('green', 8), rep('blue', 4))

cqnplot(cqnObj, n = 1, xlab = "GC content", lty = 1, col=colors)
legend("topright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))
cqnplot(cqnObj, n = 2, xlab = "length", lty = 1, col=colors)
legend("bottomright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))


#### WHAT AM I DOING HERE??????
RPKM.cqn <- cqnObj$y + cqnObj$offset
zerolengthM <- !is.na(gcLengthInfo$length)
RPM <- sweep(log2(counts[zerolengthM,] + 1), 2, log2(LibSize/10^6))

RPKM.std <- sweep(RPM, 1, log2(as.data.frame(gcLengthInfo[zerolengthM,])$length / 10^3))


normalMask <- colDat$responseType=="normal"
resisMask <- colDat$responseType=="resistance"
sensMask <- colDat$responseType=="sensitive"


# MA plots for normal versus sensitive
 whGenes <- which(rowMeans(RPKML.std) >= 2 & as.data.frame(gcLengthInfo)$length >= 100)
 cqnWhGenes<-which(rowMeans(RPKM.cqn) >= 2 & as.data.frame(cqnGCLengthInfo)$length >= 100) 
 M.std <- rowMeans(RPKM.std[whGenes, normalMask]) - rowMeans(RPKM.std[whGenes, sensMask])
 A.std <- rowMeans(RPKM.std[whGenes,])
 M.cqn <- rowMeans(RPKM.cqn[cqnWhGenes, normalMask]) - rowMeans(RPKM.cqn[cqnWhGenes, sensMask])
 A.cqn <- rowMeans(RPKM.cqn[cqnWhGenes,])


 # mA plots 
  par(mfrow = c(1,2))
 plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
        main = "Standard RPKM (N vs S)", ylim = c(-4,4), xlim = c(0,12),
        col = alpha("black", 0.25)
      )
 points(A.std[whHigh], M.std[whHigh], cex = 0.2, pch = 16, col = "red")
 points(A.std[whLow], M.std[whLow], cex = 0.2, pch = 16, col = "blue")
 plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
        main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12),
        col = alpha("black", 0.25)
        )
 points(A.cqn[cqnWhLow], M.cqn[cqnWhLow], cex = 0.2, pch = 16, col = "blue")
 points(A.cqn[cqnWhHigh], M.cqn[cqnWhHigh], cex = 0.2, pch = 16, col = "red")
 

 
 # MA plots for resistant versus sensitive
 whGenes <- which(rowMeans(RPKM.std) >= 2 & as.data.frame(gcLengthInfo)$length >= 100)
 cqnWhGenes<-which(rowMeans(RPKM.cqn) >= 2 & as.data.frame(cqnGCLengthInfo)$length >= 100) 
 M.std <- rowMeans(RPKM.std[whGenes, resisMask]) - rowMeans(RPKM.std[whGenes, sensMask])
 A.std <- rowMeans(RPKM.std[whGenes,])
 M.cqn <- rowMeans(RPKM.cqn[cqnWhGenes, resisMask]) - rowMeans(RPKM.cqn[cqnWhGenes, sensMask])
 A.cqn <- rowMeans(RPKM.cqn[cqnWhGenes,])
 
 # mA plots 
 par(mfrow = c(1,2))
 plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
      main = "Standard RPKM (R vs S)", ylim = c(-4,4), xlim = c(0,12),
      col = alpha("black", 0.25)
 )
 points(A.std[whHigh], M.std[whHigh], cex = 0.2, pch = 16, col = "red")
 points(A.std[whLow], M.std[whLow], cex = 0.2, pch = 16, col = "blue")
 plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
      main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12),
      col = alpha("black", 0.25)
 )
 points(A.cqn[cqnWhLow], M.cqn[cqnWhLow], cex = 0.2, pch = 16, col = "blue")
 points(A.cqn[cqnWhHigh], M.cqn[cqnWhHigh], cex = 0.2, pch = 16, col = "red")
 
 
 #intermediate GC content. We define high/low GC content as the 10% most extreme genes:
   par(mfrow = c(1,2))
 gccontent <- as.data.frame(gcLengthInfo)$gc[whGenes]
 cqnGcContent <- as.data.frame(cqnGCLengthInfo)$gc[cqnWhGenes]
  whHigh <- which(gccontent > quantile(gccontent, 0.9))
  whLow <- which(gccontent < quantile(gccontent, 0.1))
  cqnWhHigh <- which(cqnGcContent > quantile(cqnGcContent, 0.9))
  cqnWhLow <-which(cqnGcContent < quantile(cqnGcContent, 0.1))
  plot(A.std[whHigh], M.std[whHigh], cex = 0.2, pch = 16, xlab = "A",
         ylab = "M", main = "Standard RPKM",
         ylim = c(-4,4), xlim = c(0,12), col = "red")
  points(A.std[whLow], M.std[whLow], cex = 0.2, pch = 16, col = "blue")
  plot(A.cqn[cqnWhHigh], M.cqn[cqnWhHigh], cex = 0.2, pch = 16, xlab = "A",
         ylab = "M", main = "CQN normalized RPKM",
         ylim = c(-4,4), xlim = c(0,12), col = "red")
  points(A.cqn[cqnWhLow], M.cqn[cqnWhLow], cex = 0.2, pch = 16, col = "blue")


