library("GEOsearch")
library("SummarizedExperiment")
library("tidyr")
library("cqn")
library("EDASeq")
library("stringr")
library('scales')
library("EnsDb.Hsapiens.v86")
library("formattable")

layout(matrix(c(1)))

load(file.path("SRP043085", "rse_gene.Rdata"))


#Putting the counts into their own datafram
counts <- assays(rse_gene)$counts


# making shorter column of response type
colData(rse_gene)$responseType <- c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4))
colData(rse_gene)$responseType <- as.factor(colData(rse_gene)$responseType)
colData(rse_gene)$responseType <- as.factor(colData(rse_gene)$responseType)
colData(rse_gene)$index <- 1:20


# Getting the column data
colDat<- colData(rse_gene)

##### Getting GC Content and length info ###########

# first we need the gene IDs 
Ids <-str_split_fixed(rowData(rse_gene)$gene_id, "\\.", n=2)[,1]
geneNames <- mapIds(EnsDb.Hsapiens.v86, keys=Ids, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
lengths <- lapply(geneNames, length)
multipleMask<- lengths>1
namedGenes <- geneNames[multipleMask]
singleIds <- Ids[multipleMask]


# Now we get the information. 
# Not run, because downloading the info takes a long time. 
# gcLengthInfo <- getGeneLengthAndGCContent(Ids,'hsa', mode='biomart')
# df <- as.data.frame(gcLengthInfo)

# Because I saved the data after downloading it, we can simply load it into R
gcLengthInfo <- readRDS("/Users/rubyfore/School/PHP2620_Bioinformatics/FinalProject/gcLengthInfo.Rdata")

# If there is no gene info, cqn will not be happy
zerolengthM <- !is.na(gcLengthInfo$length)

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


# Examining too-large library sizes
colSums(counts)/1e6
which.max(counts[,1])
which.max(counts[,20])
plot(counts[21484,])
range(counts[,1])
hist(log2(1+counts[,1]))
id=21484;gcLengthInfo$length[id]
counts[id,1]/1e8
plot(counts[,2],pch=3)
id=which(counts[,1]/1e7>1)
plot(counts[21484,]~LibSize)
abline(0,.1)
abline(0,.2)
L=colSums(counts[id,])
plot(L/LibSize, col=colors, pch=10, main = "% of library due to 11 highly expressed genes")
legend("topright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))


##### MA plots with raw RPKM #######

# Average expression ratios comparing the three groups 
mNormal<- rowMeans(counts[,which(colData(rse_gene)$responseType=="normal")])
mResistant<-rowMeans(counts[,which(colData(rse_gene)$responseType=="resistance")])
mSensitive <-rowMeans(counts[,which(colData(rse_gene)$responseType=='sensitive')]) 

# Checking distribution of log2 normalized counts
hist(log2(mResistant/mNormal), breaks= 30, xlim = c(-15,15))
hist(log2(mSensitive/mNormal), breaks=30,xlim=c(-15,15))
hist(log2(mSensitive/mResistant), breaks=30,xlim=c(-15,15)) 

# MA plots
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



### Running CQN and obtaining plots #########
cqnObjno4 <- cqn(counts[zerolengthM,-4], x=as.data.frame(gcLengthInfo)$gc[zerolengthM], lengths=as.data.frame(gcLengthInfo)$length[zerolengthM], sizeFactors=NULL, verbose=F)
cqnObj <- cqn(counts[zerolengthM,], x=as.data.frame(gcLengthInfo)$gc[zerolengthM], lengths=as.data.frame(gcLengthInfo)$length[zerolengthM], sizeFactors=NULL, verbose=F)

#matching length and gc info to cqn obj
# need to trim version number off cqnObj rownames
cqnIds <-str_split_fixed(rownames(cqnObj$counts), "\\.", n=2)[,1]
cqnGCLengthInfo<- subset(gcLengthInfo, rownames(gcLengthInfo)%in%cqnIds)

layout(matrix(c(1)))
# Making cqn plots
colors<- c(rep('red', 8), rep('green', 8), rep('blue', 4))
cqnplot(cqnObj, n = 1, xlab = "GC content", lty = 1, col=colors, main='GC content effect')
legend("topright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))
cqnplot(cqnObj, n = 2, xlab = "length", lty = 1, col=colors, main="Length effect")
legend("bottomright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))

# And making them again without sample 4
colors<- c(rep('red', 7), rep('green', 8), rep('blue', 4))
cqnplot(cqnObjno4, n = 1, xlab = "GC content", lty = 1, col=colors, main='GC content effect')
legend("topright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))
cqnplot(cqnObjno4, n = 2, xlab = "length", lty = 1, col=colors, main='Length effect')
legend("bottomright", legend=c('normal', 'resistant', 'sensitive'),  fill=c('red', 'green', 'blue'))


### Making MA plots demonstrating the impact of cqn normalization ####3
# getting cqn standardized rpkm
RPKM.cqn <- cqnObj$y + cqnObj$offset
# getting un-standardized RPKM
zerolengthM <- !is.na(gcLengthInfo$length)
RPM <- sweep(log2(counts[zerolengthM,] + 1), 2, log2(LibSize/10^6))
RPKM.std <- sweep(RPM, 1, log2(as.data.frame(gcLengthInfo[zerolengthM,])$length / 10^3))

# Creating masks
normalMask <- colDat$responseType=="normal"
resisMask <- colDat$responseType=="resistance"
sensMask <- colDat$responseType=="sensitive"


# MA plots for normal versus sensitive
whGenes <- which(rowMeans(RPKM.std) >= 2 & as.data.frame(cqnGCLengthInfo)$length >= 100)
cqnWhGenes<-which(rowMeans(RPKM.cqn) >= 2 & as.data.frame(cqnGCLengthInfo)$length >= 100) 
M.std <- rowMeans(RPKM.std[whGenes, normalMask]) - rowMeans(RPKM.std[whGenes, sensMask])
A.std <- rowMeans(RPKM.std[whGenes,])
M.cqn <- rowMeans(RPKM.cqn[cqnWhGenes, normalMask]) - rowMeans(RPKM.cqn[cqnWhGenes, sensMask])
A.cqn <- rowMeans(RPKM.cqn[cqnWhGenes,])


whHigh <- which(cqnGCLengthInfo$gc > quantile(cqnGCLengthInfo$gc, 0.9))
whLow <- which(cqnGCLengthInfo$gc < quantile(cqnGCLengthInfo$gc, 0.1))
cqnWhHigh <- which(cqnGCLengthInfo$gc > quantile(cqnGCLengthInfo$gc, 0.9))
cqnWhLow <-which(cqnGCLengthInfo$gc < quantile(cqnGCLengthInfo$gc, 0.1))

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
whGenes <- which(rowMeans(RPKM.std) >= 2 & as.data.frame(cqnGCLengthInfo)$length >= 100)
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

####### DIFFERENTIAL GENE EXPRESSION ################


# Make design matrix 
design <- model.matrix(~0+responseType, data=colDat)
design2 <- model.matrix(~0 + responseType, data=colDat[-4,])
colnames(design) <- c('normal', 'resistance', 'sensitive')
colnames(design2) <- c('normal', 'resistance', 'sensitive')

contrasts <- makeContrasts(trtA = resistance-normal, trtB=resistance-sensitive, trtC = sensitive-normal ,
                           levels=c("normal", 'resistance', 'sensitive'))

# make a DGElist object from the counts matrix

#choosing rows that cqn did not remove
indices <- which(rownames(cqnObj$counts)%in%names(id) )

# making a dge object removing the samples of low quality
dgeObj <- DGEList(counts=cqnObjno4$counts[-indices,], lib.size=NULL, group=c(rep('normal', 7), rep('resistance', 8), rep('sensitive', 4)))
dgeObj <- calcNormFactors(dgeObj, method="TMM")
dgeObj <- estimateDisp(dgeObj, design2)

# using data with sample 4, and genes that have counts that are too high
dgeObj.Normal <- DGEList(counts=counts,group=c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4)))
dgeObj.Normal <- calcNormFactors(dgeObj.Normal)

# estimate dispersion
dgeObj.Normal <- estimateDisp(dgeObj.Normal, design)


#fitting model with complete data 
efit.norm <- glmFit(dgeObj.Normal, design = design, prior.count=1)

# fitting model with selected data 
## NOTE: to run for complete data, comment out this line 
efit.norm <- glmFit(dgeObj, design = design2, prior.count=1)

elrt.norm.A <- glmLRT(efit.norm, contrast = contrasts[,'trtA'])
elrt.norm.B <- glmLRT(efit.norm, contrast = contrasts[,'trtB'])
elrt.norm.C <- glmLRT(efit.norm, contrast = contrasts[,'trtC'])


tta <- topTags(elrt.norm.A, n=Inf, sort.by="logFC", p=0.01)
ttb <- topTags(elrt.norm.B, n=Inf,sort.by="logFC",p=0.1)
ttc <- topTags(elrt.norm.C, n=Inf,sort.by="logFC",p=0.05)



IdsA <- str_split_fixed(rownames(tta), "\\.", n=2)[,1]
IdsB <- str_split_fixed(rownames(ttb), "\\.", n=2)[,1]
IdsC <- str_split_fixed(rownames(ttc), "\\.", n=2)[,1]


geneNamesa <- mapIds(EnsDb.Hsapiens.v86, keys=IdsA, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
geneNamesb <- mapIds(EnsDb.Hsapiens.v86, keys=IdsB, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
geneNamesc <- mapIds(EnsDb.Hsapiens.v86, keys=IdsC, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')


geneInfoA <- select(org.Hs.eg.db, keys=unlist(geneNamesa), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")
geneInfoB <- select(org.Hs.eg.db, keys=unlist(geneNamesb), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")
geneInfoC<-select(org.Hs.eg.db, keys=unlist(geneNamesc), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")


tablea <- cbind(geneInfoA$GENENAME, signif(tta$table[,c(1,4,5)], digits=2))
names(tablea)[[1]] <- "Gene Info"
tablea <- tablea[!is.na(tablea$`Gene Info`),]
rownames(tablea) <- NULL
formattable(tablea[1:15,])


tableb <- cbind(geneInfoB$GENENAME, signif(ttb$table[,c(1,4,5)], digits=2))
names(tableb)[[1]] <- "Gene Info"
tableb <- tableb[!is.na(tableb$`Gene Info`),]
rownames(tableb) <- NULL
formattable(tableb)

tablec <- cbind(geneInfoC$GENENAME, signif(ttc$table[,c(1,4,5)], digits=2))
names(tablec)[[1]] <- "Gene Info"
tablec <- tablec[!is.na(tablec$`Gene Info`),]
rownames(tablec) <- NULL
formattable(tablec)




