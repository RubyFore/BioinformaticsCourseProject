
library(org.Hs.eg.db)
library(limma)
source("~/School/PHP2620_Bioinformatics/FinalProject/loadData.R")
library("edgeR")
library("formattable")
library('knitr')
### DIFFERENTIAL GENE EXPRESSION ANALYSIS ######

# Make design matrix 
design <- model.matrix(~0+responseType, data=colDat)
design2 <- model.matrix(~responseType, data=colDat)
colnames(design) <- c('normal', 'resistance', 'sensitive')
colnames(design2) <- c('normal', 'resistance', 'sensitive')

contrasts <- makeContrasts(trtA = resistance-normal, trtB=resistance-sensitive, trtC = sensitive-normal ,
                           levels=c("normal", 'resistance', 'sensitive'))

# make a DGElist object from the counts matrix
#choosing rows that cqn did not remove 
dgeObj <- DGEList(counts=cqnObj$counts, lib.size=NULL, group=c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4)))
dgeObj2 <- DGEList(counts=cqnObj$counts, lib.size=NULL, group=c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4)))
dgeObj2$offset <- cqnObj$offset
# ? dgeObj$counts <- dgeObj$counts +1
dgeObj$offset <- cqnObj$glm.offset
dgeObj.Normal <- DGEList(counts=counts,group=c(rep('normal', 8), rep('resistance', 8), rep('sensitive', 4)))
dgeObj.Normal <- calcNormFactors(dgeObj.Normal)
dgeObj.S3 <- dgeObj

# estimate dispersion
dgeObj <- estimateGLMCommonDisp(dgeObj, design)
dgeObj <- estimateGLMTagwiseDisp(dgeObj, design)
dgeObj.Normal <- estimateDisp(dgeObj.Normal, design)
# for dge Obj2, with better offsets from cqn

dgeObj2 <- estimateGLMCommonDisp(dgeObj2, design)
dgeObj2 <- estimateGLMTagwiseDisp(dgeObj2, design)


dgeObj.S3$common.dispersion<- estimateGLMCommonDisp(dgeObj$counts, design, offset=dgeObj$offset)
dgeObj.S3$tagwise.dispersion <- estimateGLMTagwiseDisp(dgeObj.S3$counts, design, offset=dgeObj$offset, dispersion=dgeObj.S3$common.dispersion)
#fitting model
efit.cqn.s3 <- glmFit(dgeObj.S3$counts, design = design, offset=dgeObj.S3$offset, dispersion=dgeObj.S3$tagwise.dispersion, prior.count=1)
elrt.cqn.s3 <- glmLRT(efit.cqn.s3, contrast = contrasts[,'trtA'])

efit.norm <- glmFit(dgeObj.Normal, design = design, prior.count=1)
elrt.norm.A <- glmLRT(efit.norm, contrast = contrasts[,'trtA'])
elrt.norm.B <- glmLRT(efit.norm, contrast = contrasts[,'trtB'])
elrt.norm.C <- glmLRT(efit.norm, contrast = contrasts[,'trtC'])

efit.cqn <- glmFit(dgeObj, design = design, prior.count=1)
elrt.cqn <- glmLRT(efit.cqn, contrast = contrasts[,'trtA'])

efit.cqn.2 <-glmFit(dgeObj2, design = design, prior.count=1)
elrt.cqn2 <- glmLRT(efit.cqn.2, contrast = contrasts[,'trtA'])

# gives sensible results
et <- exactTest(dgeObj.Normal)

# using f test instead
#qlFit <- glmQLFit(dgeObj, design=design)
#eft.cqn <- glmQLFTest(qlFit, contrast=contrasts[,'trtA'], coef=1)
topTags(elrt.cqn.s3, n = 15)
tta <- topTags(elrt.norm, n=15)
ttb <- topTags(elrt.norm.B, n=15)
ttc <- topTags(elrt.norm.C, n=15)
topTags(elrt.cqn.2, n=15)
topTags(et, n=15)[rownames(topTags(et, n=15)) %in%rownames(topTags(elrt.cqn, n = 15)),]
topTags(elrt.cqn.s3, n = 15)[rownames(topTags(elrt.cqn.s3, n = 15))%in%rownames(topTags(elrt.cqn, n=15)),]
topTags(elrt.cqn, n = 15)[rownames(topTags(elrt.cqn, n = 15))%in%rownames(topTags(elrt.cqn.s3, n=15)),]
topTags(elrt.norm, n=15)[rownames(topTags(elrt.norm, n=15))%in%rownames(topTags(et, n=15)),]
big <- "ENSG00000228131.1"
# runing regressions - exactTest?
#normVres <- exactTest(dgeObj, pair= c('normal','resistance'), dispersion = 'common')
#resVsens <- exactTest(dgeObj, pair=c('sensitive', 'resistance'), dispersion='common')

# doing DE analysis another way using limma
#fit <- lmFit(cqnObj$y, design)
#cont.fit <- contrasts.fit(fit, contrasts)
#cont.fit <- eBayes(cont.fit)

#topA<-topTable(cont.fit,coef=1, number=Inf,sort.by='p', adjust.method='BH', p.value=0.05)


#Idss3 <- str_split_fixed(rownames(topTags(elrt.cqn.s3, n=20)), "\\.", n=2)[,1]
IdsA <- str_split_fixed(rownames(tta), "\\.", n=2)[,1]
IdsB <- str_split_fixed(rownames(ttb), "\\.", n=2)[,1]
IdsC <- str_split_fixed(rownames(ttc), "\\.", n=2)[,1]

#Ids.cqn2 <- str_split_fixed(rownames(topTags(elrt.cqn2, n=15)), "\\.", n=2)[,1]
rowIds <- str_split_fixed(rownames(efit.norm), "\\.", n=2)[,1]
#Ids <-str_split_fixed(rownames(topTags(normVres)$table), "\\.", n=2)[,1]
#Ids <-str_split_fixed(rownames(topA), "\\.", n=2)[,1]

geneNamesa <- mapIds(EnsDb.Hsapiens.v86, keys=IdsA, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
geneNamesb <- mapIds(EnsDb.Hsapiens.v86, keys=IdsB, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')
geneNamesc <- mapIds(EnsDb.Hsapiens.v86, keys=IdsC, column = c("GENENAME"), keytype = "GENEID", multiVals = 'list')


#rowEntrez<- mapIds(EnsDb.Hsapiens.v86, keys=rowIds, column = c("ENTREZID", "GENEID"), keytype = "GENEID", multiVals = 'list')
#rowEntrez <- select(org.Hs.eg.db, keys=rowIds, columns=c("ENTREZID", "ENSEMBL"), keytype="ENSEMBL", multiVals='first')
formattable(results4[1:20,])
# getting the rest of the info from org.HS.eg.db
geneInfoA <- select(org.Hs.eg.db, keys=unlist(geneNamesa), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")
geneInfoB <- select(org.Hs.eg.db, keys=unlist(geneNamesb), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")
geneInfoC<-select(org.Hs.eg.db, keys=unlist(geneNamesc), columns=c("SYMBOL", "GENENAME"), keytype="SYMBOL")


tablea <- cbind(geneInfoA$GENENAME, signif(tta$table[,c(1,4,5)], digits=2))
names(tablea)[[1]] <- "Gene Info"
rownames(tablea) <- NULL
formattable(tablea)


tableb <- cbind(geneInfoB$GENENAME, signif(ttb$table[,c(1,4,5)], digits=2))
names(tableb)[[1]] <- "Gene Info"
rownames(tableb) <- NULL
formattable(tableb)

tablec <- cbind(geneInfoC$GENENAME, signif(ttc$table[,c(1,4,5)], digits=2))
names(tablec)[[1]] <- "Gene Info"
rownames(tablec) <- NULL
formattable(tablec)


mapIds(EnsDb.Hsapiens.v86, keys="ENSG00000111087", column="GENENAME", keytype="GENEID")
# gene ontology analysis
go <- goana(elrt.norm)





