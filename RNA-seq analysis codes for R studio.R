
#data clearing
tmp <- read.csv("counts of gene expression.csv", header=T) 

rownames(tmp) <- tmp[,2]
tmp <- tmp[,-c(1,2)]

for (i in 1:9)
{
  colnames(tmp)[i] <- strsplit( substr(colnames(tmp),14,30)[i], ".", fixed=T)[[1]][1]
}
o <- order(colnames(tmp))

tmp <- tmp[,o]
cc <- tmp 
ccl <- apply(cc, c(1,2), function(x) log(x+1))

annot <- read.table("genes list", sep="\t", header=F) 
rownames(annot) <- annot[,1]


#PCA analysis
library(PCAtools)

mdata <- data.frame(c("T","T","T","T","T","T","WT","WT","WT"))
rownames(mdata) <- colnames(ccl)
colnames(mdata) <- "group"

p<- pca(ccl, metadata = mdata)
biplot(p, showLoadings = F )
biplot(p,   x = "PC1", y = "PC2",showLoadings = F , 
       colby="group", legendPosition = 'right')
biplot(p,   x = "PC1", y = "PC3",showLoadings = F , 
       colby="group", legendPosition = 'right')
biplot(p,   x = "PC2", y = "PC3",showLoadings = F , 
       colby="group", legendPosition = 'right')
biplot(p,   x = "PC1", y = "PC4",showLoadings = F , 
       colby="group", legendPosition = 'right')


#EdgeR
library(edgeR)
library(Vennerable)

targets <- data.frame(c("TA","TA","TA","TB","TB","TB","WT","WT","WT"))
colnames(targets)  <- "group"

Treat <- factor(targets$group)
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

y <- DGEList(counts=cc)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

cma <- makeContrasts( TA-WT, levels=design)
lrt <- glmLRT(fit, contrast=cma)
topOnes_TA_WT <- topTags(lrt, n=dim(cc)[1])$table
length(which(topOnes_TA_WT$FDR<0.01 ))

tmp <- cbind(topOnes_TA_WT, annot[rownames(topOnes_TA_WT),], cc[rownames(topOnes_TA_WT),])
write.csv(tmp, file="topOnes_TA_WT.csv")

cm <- makeContrasts( TA, levels=design)
lrt <- glmLRT(fit, contrast=cm)
topOnes <- topTags(lrt, n=dim(cc)[1])$table

cmb <- makeContrasts( TB-WT, levels=design)
lrt <- glmLRT(fit, contrast=cmb)
topOnes_TB_WT <- topTags(lrt, n=dim(cc)[1])$table
length(which(topOnes_TB_WT$FDR<0.01 ))
tmp <- cbind(topOnes_TB_WT, annot[rownames(topOnes_TB_WT),], cc[rownames(topOnes_TB_WT),])
write.csv(tmp, file="topOnes_TB_WT.csv")

cmab <- makeContrasts( TA-TB, levels=design)
cmab
lrt <- glmLRT(fit, contrast=cmab)
topOnes_TA_TB <- topTags(lrt, n=dim(cc)[1])$table
length(which(topOnes_TA_TB$FDR<0.01 ))

thr <- 0.01
TA_TB <- rownames(topOnes_TA_TB)[which(topOnes_TA_TB$FDR<thr )]
TA_WT <- rownames(topOnes_TA_WT)[which(topOnes_TA_WT$FDR<thr )]
TB_WT <- rownames(topOnes_TB_WT)[which(topOnes_TB_WT$FDR<thr )]

vennD <- Venn ( list(TA_TB, TA_WT, TB_WT), SetNames = c("TA_TB", "TA_NT", "TB_NT")) # NT:Non-Stimulated is the WT group
plot(vennD, doWeights = T, type = "circles")  # Figure S8 A Plotting

which(topOnes_TA_WT$FDR<0.01 )
toptaname <- rownames(topOnes_TA_WT)[which(topOnes_TA_WT$FDR<0.01 )]
length(toptaname)
toptaname

which(topOnes_TA_TB$FDR<0.01 )
toptabname <- rownames(topOnes_TA_TB)[which(topOnes_TA_TB$FDR<0.01 )]
length(toptabname)
toptabname

which(topOnes_TB_WT$FDR<0.01 )
toptbname <- rownames(topOnes_TB_WT)[which(topOnes_TB_WT$FDR<0.01 )]
length(toptbname)
toptbname


#cc mean
ccta <- cc[1:3]
cctb <- cc[4:6]
ccwt <- cc[7:9]

ccta$mean <- apply(cc[,1:3], 1, mean, na.rm=T)
cctamean <- ccta[,4]

cctb$mean <- apply(cctb[,1:3], 1, mean, na.rm=T)
cctbmean = cctb[,4]

ccwt$mean <- apply(ccwt[,1:3], 1, mean, na.rm=T)
ccwtmean = ccwt[,4]

ccam <- c(cctamean, cctbmean, ccwtmean)
ccallmean <- matrix (c(ccam), nrow = 61487, ncol=3, byrow = FALSE)
ccallmean

rownamenew = rownames(cc)
colnamenew = c("TAmean", "TBmean", "NSmean")
rownames(ccallmean) = rownamenew
colnames(ccallmean) = colnamenew
cclallmean <- apply(ccallmean, c(1,2), function(x) log(x+1))

# Figure S8 B-D Plotting
plot(cclallmean[,1], cclallmean[,3], cex=0.3)    #### TA_WT
points(cclallmean[toptaname,1],cclallmean[toptaname,3], cex=0.5, col="orange") # Figure S8B

plot(cclallmean[,2], cclallmean[,3], cex=0.3)      #### TB_WT
points(cclallmean[toptbname,2],cclallmean[toptbname,3], cex=0.3, col="green") # Figure S8C

plot(cclallmean[,1], cclallmean[,2], cex=0.3)      #### TA_TB
points(cclallmean[toptabname,1],cclallmean[toptabname,2], cex=0.3, col="red") # Figure S8D


# Figure S9A
library(pheatmap)

sds <- apply(ccl, 1, sd)
o <- order(sds, decreasing = TRUE)
ht <- ccl[o[1:500], ]

annot_col <- data.frame(c("TA","TA","TA","TB","TB","TB","WT","WT","WT"))
rownames(annot_col) <- colnames(ccl)
colnames(annot_col) <- "group"
annot_row <- annot[rownames(ht),2]

pheatmap(ht, annotation_col = annot_col, labels_row=annot_row, fontsize_row = 4, fontsize_col = 7, cluster_rows=T)
#N1, N2 and N3 in the Figure S9A are the WT1, WT2 and WT3, respectively. # Figure S9A
write.csv(ht, file = "List S2-overall heatmap list-top500.csv", quote=F, row.names = T)
# List S2

# Figure S9B
library(Vennerable)
library(pheatmap)

jl <- read.csv("List S3-Antioxidant related genes in human.csv", header=T) 
# List S3

ens_jl <- jl[,1]
ens_jl2 <- intersect(ens_jl, rownames(ccl))
ht2 <- ccl[ens_jl2, ]
annot_col <- data.frame(c("TA","TA","TA","TB","TB","TB","WT","WT","WT"))
rownames(annot_col) <- colnames(ccl)
colnames(annot_col) <- "group"
annot_row <-  data.frame(jl$notes)
rownames(annot_row) <- ens_jl

pheatmap(ht2, annotation_col = annot_col, annotation_row=annot_row,
         fontsize_row = 4, fontsize_col = 7, cluster_rows=T)    # Figure S9B
#N1, N2 and N3 in the Figure S9A are the WT1, WT2 and WT3, respectively.

write.csv(ht2, file = "List S4-Clustering list of expression levels among antioxidant related genes.csv", quote=F, row.names = T)
# List S4
