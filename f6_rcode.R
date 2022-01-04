library(limma)
library(edgeR)
library(dplyr)
library(gplots)
library(org.Mm.eg.db)
library(ggplot2)
library(Mus.musculus)
library(ggvenn)

DGEobj = read.csv("merged_FC.csv", header= TRUE)
group=as.factor(c('IL1','IL1','IL1','LPS','LPS','LPS','PAM','PAM','PAM','CON','CON','CON') )
DGEobj = tibble::column_to_rownames(DGEobj[-1], var="gene")
DGEobj = DGEList(counts = DGEobj, group=group)
keep <- filterByExpr(DGEobj, min.count = 25, group=group)
DGEobj <- DGEobj[keep,,keep.lib.sizes=FALSE]
DGEobj <- calcNormFactors(DGEobj)

logCPM <- cpm(DGEobj, log=TRUE, prior.count = 3)

design <- model.matrix(~0+group)
contr.matrix <- makeContrasts(IL1vsCon = groupIL1 - groupCON, 
                              LPSvsCon = groupLPS - groupCON,
                              PAMvsCon = groupPAM - groupCON,
                              LPSvsIL1 = groupLPS - groupIL1,
                              levels = colnames(design))

fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
dt <- decideTests(efit, adjust.method = "BH", p.value = 0.01,lfc = 1)
topTreat(efit, coef=4, n=Inf)
IL1vsCon <- topTreat(efit, coef=1, n=Inf)
LPSvsCon <- topTreat(efit, coef=2, n=Inf)
PAMvsCon <- topTreat(efit, coef=3, n=Inf)
LPSvsIL1 <- topTreat(efit, coef=4, n=Inf)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)

foo <- list(A = rownames( dt[ dt[,1]!=0, ] ),
            B = rownames( dt[ dt[,2]!=0, ] ),
            C = rownames( dt[ dt[,4]!=0, ] ),
            D = rownames( dt[ dt[,3]!=0, ] ) )
names(foo) <- c("IL1", "LPS", "PAM", "LPS vs IL1")

#F6A
win.graph()
ggvenn(foo[1:3], fill_color=c("seagreen", "royalblue", "indianred3"), stroke_size = 1, show_percentage = FALSE, text_size = 8, set_name_size = 10 )



#F6B
IL1vsCon.genes <- rownames( dt[ dt[,1]!=0, ] )
LPSvsCon.genes <- rownames( dt[ dt[,2]!=0, ] )
LPSvsIL1.genes <- rownames( dt[ dt[,4]!=0, ] )
PAMvsCon.genes <- rownames( dt[ dt[,3]!=0, ] )


i <- c(LPSvsCon.genes, IL1vsCon.genes, PAMvsCon.genes)
i <- unique(i)

mycol <- colorpanel(100,"blue","white","red")
win.graph()
lmat = rbind(c(0,3, 0),c(2,1, 0),c(0,0, 4))


heatmapdata <- logCPM[, c(10:12, 1:3, 7:9, 4:6)]
coolmap(heatmapdata[i, ], margins=c(7,7), lmat = lmat, lhei=c(1.5,4, 1), 
        lwid=c(1,4, 2), show.dendrogram="row", srtCol=45,  adjCol = c(1,1) , 
        labRow = "", col=mycol, labCol = c("Con","Con","Con", "IL-1","IL-1",
                                           "IL-1","PAM","PAM","PAM", "LPS",
                                           "LPS","LPS"),
        key = TRUE, keysize = 5, linkage.col = "none", linkage.row = "ward.D2") 

geneid <- rownames(logCPM)
genes <- select(org.Mm.eg.db, keys = geneid, column = c('SYMBOL'), keytype = "ENTREZID")
genes <- genes[!is.na(genes), ]


#F6C-D
LPSvsIL1$logP <- -log10( LPSvsIL1$adj.P.Val )
x <- org.Mm.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
xx <- as.list(org.Mm.egGO2EG)
cytokines <- xx$`GO:0005125`
nfkbregs <- xx$`GO:0032088`
nfkbregs <- c(nfkbregs, "18037", '80859') #adding IkBe and IkBz


# Fig 3D cytokines regs
foo <- complete.cases(match(rownames( LPSvsIL1 ), cytokines))
foo2 <- rownames(LPSvsIL1)[foo]
foo2 <- select(org.Mm.eg.db, keys = foo2, column = c('SYMBOL'), keytype = "ENTREZID")
foo2$logP <- LPSvsIL1[foo2$ENTREZID, ]$logP
foo2$logFC <- LPSvsIL1[foo2$ENTREZID, ]$logFC
foo2$DE <- foo2$logP > 2 & foo2$logFC > 2
col <- rep(FALSE, length(foo))
col[match( foo2$ENTREZID[foo2$DE == TRUE], rownames(LPSvsIL1) )] <- TRUE

ggplot(LPSvsIL1) + 
  geom_point(size = 2, aes(x = logFC, y = logP, color = col, shape = col, fill = col) ) + 
  scale_color_manual(LPSvsIL1, values=c("grey60","black" ))  +
  scale_shape_manual(LPSvsIL1, values = c(19, 21)) + 
  scale_fill_manual(LPSvsIL1, values = c("grey60","red" ), aesthetics = "fill", breaks = waiver(), labels = c("", "Cytokines")) +
  theme_classic() +
  xlab("Fold Change (Log2)") +
  ylab("p-value (Log 10)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(alpha = FALSE) + 
  theme(legend.title = element_blank(), legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold")) 


# Fig 3E NFKB regs
foo <- complete.cases(match(rownames( LPSvsIL1 ), nfkbregs))
foo2 <- rownames(LPSvsIL1)[foo]
foo2 <- select(org.Mm.eg.db, keys = foo2, column = c('SYMBOL'), keytype = "ENTREZID")
foo2$logP <- LPSvsIL1[foo2$ENTREZID, ]$logP
foo2$logFC <- LPSvsIL1[foo2$ENTREZID, ]$logFC
foo2$DE <- foo2$logP > 2 & foo2$logFC > 2


col <- rep(FALSE, length(foo))
col[match( foo2$ENTREZID[foo2$DE == TRUE], rownames(LPSvsIL1) )] <- TRUE

linepoints <- cbind(c(4.1, 4.1, 4.6, 5.1, 4, 4.1, 4.1, 4.4), c(6, 5.5, 4.25, 4, 3.5, 3, 2.1, 1.5), c(1, 2, 3, 4, 5, 6, 7, 8))
linepoints2 <- cbind(foo2$logFC[foo2$DE == TRUE], foo2$logP[foo2$DE == TRUE], c(1, 2, 3, 4, 5, 6, 7, 8))
linepoints <- rbind(linepoints, linepoints2)
linepoints <- as.data.frame(linepoints)

win.graph()

ggplot(LPSvsIL1) + 
  geom_point(size = 2, aes(x = logFC, y = logP, color = col, shape = col, fill = col) ) + 
  scale_color_manual(LPSvsIL1, values=c("grey60","black" ))  +
  scale_shape_manual(LPSvsIL1, values = c(19, 21)) + 
  scale_fill_manual(LPSvsIL1, values = c("grey60","red" ), aesthetics = "fill", breaks = waiver(), labels = c("", "NF-kB regulators")) +
  #annotate("text", x=c(5, 5, 5, 6, 5, 5, 5, 5), y=foo2$logP[foo2$DE == TRUE]+0.2, label=foo2$SYMBOL[foo2$DE == TRUE]) +
  annotate("text", x=c(5, 5, 5.5, 6, 5, 5, 5, 5), y=c(6, 5.5, 4.25, 4, 3.5, 3, 2, 1.5), label=foo2$SYMBOL[foo2$DE == TRUE], size= 5) +
  theme_classic() +
  xlab("Fold Change (Log2)") +
  ylab("p-value (Log 10)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(alpha = FALSE) + 
  theme(legend.title = element_blank(), legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"))  +
  geom_line(data = linepoints, aes(x = V1, y = V2, group = V3 ))


