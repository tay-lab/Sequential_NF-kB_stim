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

#Figure 3F qPCR

ligandorder <- read.csv("20201124_LOqPCR3_tidy.csv", header= TRUE)

ligandstats <- ligandorder %>% group_by(Condition, Gene) %>% 
  summarize(mean = mean(Cq), sd = sd(Cq), n = n(), se = sd / sqrt(n))

ligandstats$Condition <- factor(ligandstats$Condition,levels = c("Con", "IL-1_30", "IL-1_60", "IL-1_120",
                                                                 "IL-1_MH", "IL-1_H", "LPS"))

LPS_Nfkbia <- ligandorder %>% filter(Condition == "LPS" & Gene == "Nfkbia")
IL1_120_Nfkbia <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Nfkbia")
IL1_MH_Nfkbia <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Nfkbia")

LPS_Nfkbie <- ligandorder %>% filter(Condition == "LPS" & Gene == "Nfkbie")
IL1_120_Nfkbie <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Nfkbie")
IL1_MH_Nfkbie <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Nfkbie")

LPS_Tnfaip3 <- ligandorder %>% filter(Condition == "LPS" & Gene == "Tnfaip3")
IL1_120_Tnfaip3 <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Tnfaip3")
IL1_MH_Tnfaip3 <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Tnfaip3")

LPS_Il23a <- ligandorder %>% filter(Condition == "LPS" & Gene == "Il23a")
IL1_120_Il23a <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Il23a")
IL1_MH_Il23a <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Il23a")

LPS_Cxcl3 <- ligandorder %>% filter(Condition == "LPS" & Gene == "Cxcl3")
IL1_120_Cxcl3 <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Cxcl3")
IL1_MH_Cxcl3 <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Cxcl3")

LPS_Cxcl2 <- ligandorder %>% filter(Condition == "LPS" & Gene == "Cxcl2")
IL1_120_Cxcl2 <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Cxcl2")
IL1_MH_Cxcl2 <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Cxcl2")

LPS_Csf3 <- ligandorder %>% filter(Condition == "LPS" & Gene == "Csf3")
IL1_120_Csf3 <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Csf3")
IL1_MH_Csf3 <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Csf3")

LPS_Csf2 <- ligandorder %>% filter(Condition == "LPS" & Gene == "Csf2")
IL1_120_Csf2 <- ligandorder %>% filter(Condition == "IL-1_120" & Gene == "Csf2")
IL1_MH_Csf2 <- ligandorder %>% filter(Condition == "IL-1_MH" & Gene == "Csf2")

pvals <- c()
pvals[1] <- t.test(LPS_Nfkbia$Cq, IL1_120_Nfkbia$Cq)$p.value
pvals[2] <- t.test(LPS_Nfkbia$Cq, IL1_MH_Nfkbia$Cq)$p.value
pvals[3] <- t.test(LPS_Nfkbie$Cq, IL1_120_Nfkbie$Cq)$p.value
pvals[4] <- t.test(LPS_Nfkbie$Cq, IL1_MH_Nfkbie$Cq)$p.value
pvals[5] <- t.test(LPS_Tnfaip3$Cq, IL1_120_Tnfaip3$Cq)$p.value
pvals[6] <- t.test(LPS_Tnfaip3$Cq, IL1_MH_Tnfaip3$Cq)$p.value
pvals[7] <- t.test(LPS_Il23a$Cq, IL1_120_Il23a$Cq)$p.value
pvals[8] <- t.test(LPS_Il23a$Cq, IL1_MH_Il23a$Cq)$p.value
pvals[9] <- t.test(LPS_Cxcl3$Cq, IL1_120_Cxcl3$Cq)$p.value
pvals[10] <- t.test(LPS_Cxcl3$Cq, IL1_MH_Cxcl3$Cq)$p.value
pvals[11] <- t.test(LPS_Cxcl2$Cq, IL1_120_Cxcl2$Cq)$p.value
pvals[12] <- t.test(LPS_Cxcl2$Cq, IL1_MH_Cxcl2$Cq)$p.value
pvals[13] <- t.test(LPS_Csf3$Cq, IL1_120_Csf3$Cq)$p.value
pvals[14] <- t.test(LPS_Csf3$Cq, IL1_MH_Csf3$Cq)$p.value
pvals[15] <- t.test(LPS_Csf2$Cq, IL1_120_Csf2$Cq)$p.value
pvals[16] <- t.test(LPS_Csf2$Cq, IL1_MH_Csf2$Cq)$p.value
p.adjust(pvals, method = 'fdr', n = 16)

ligandstats %>% filter(Condition %in% c("IL-1_120", "IL-1_MH", "LPS")) %>%
  ggplot(aes(x=Gene, y=mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("steelblue1", "steelblue4", "red3"), 
                    labels = c("0.2 ng", "1 ng", "LPS")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3,
                position=position_dodge(.9), size = 0.75) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"), legend.position="bottom",
        legend.key.size = unit(0.4, 'cm'), axis.text=element_text(size=11), 
        axis.title=element_text(size=14), legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) +
  labs(x = "Gene ID",
       y = "Fold change over control",
       fill = "Treatment") 


