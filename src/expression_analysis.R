
# R package dependencies -------------------------------------------------------
# Comment this line out if your machine already has these packages installed
install.packages(c("matrixStats","ggplot2","MASS","scales"),repos='https://watson.nci.nih.gov/cran_mirror/')

## FPKM matrix input -----------------------------------------------------------
data1 <- read.table("../RESULTS/normout/genes.count_table",header=TRUE,sep="\t", stringsAsFactors = FALSE)
attr.table <- read.table("../RESULTS/normout/genes.attr_table",header=TRUE,sep="\t",stringsAsFactors = FALSE)
data1$gene_short_name <- attr.table$gene_short_name
# write out complete gene expression matrix
write.csv(data1, file="../RESULTS/ES&DE_dataset_cufflinksFPKM.csv")

# Generate row stats --------------------------------------------------  
# rename data1 as database to fork downstream additions
database <- data1
# generate group shortnames for selecting columns
ES <- grep("ES", colnames(database),ignore.case=F)
DE <- grep("DefEnd", colnames(database),ignore.case=F)
 library(matrixStats)
## calculate mean by treatment type
database$es_mean <- rowMeans(database[,ES], na.rm=T)
database$de_mean <- rowMeans(database[,DE], na.rm=T)
 # calculate log2 change
database$xlog2 <- log2(database$de_mean/database$es_mean)
 ## function to compare by row, returns t distribution
## The function is currently defined as
row.t <- function(mat1,mat2){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n1<-dim(mat1)[2]
  n2<-dim(mat2)[2] 
  n<-n1+n2 
  m1<-rowMeans(mat1,na.rm=TRUE) 
  m2<-rowMeans(mat2,na.rm=TRUE) 
  v1<-rowVars(mat1,na.rm=TRUE) 
  v2<-rowVars(mat2,na.rm=TRUE) 
  vpool<-(n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2 
  tstat<-sqrt(n1*n2/n)*(m1-m2)/sqrt(vpool) 
  return(tstat)}
 # calculate t-distribution for hES vs. Definitive endoderm
database$tstat <- row.t(database[,ES],database[,DE])
 # express t-dist as p-value
database$p <- 2*pt(-abs(database$tstat),df=2)
 # calculate Bonferroni correction
database$Bonf_p <- p.adjust(database$p, method = 'bonferroni', n = length(database$p))
database <- database[order(-database$xlog2),]
# output matrix with stats
write.csv(database,file="../RESULTS/ES&DE_dataset_cufflinksFPKM_stats.csv")
#--------------------------------------------------------------------------------
# SET CRITERIA FOR INCLUSION
# Fmpk cut-off 
fco <- 0.4
# P-value cut-off
pco <- 0.05
#--------------------------------------------------------------------------------
#Count the number of reads with Fmpk > 2 in each group 
database$es_count <- rowSums(database[,ES] >= fco)
database$de_count <- rowSums(database[,DE] >= fco)
# generate an output (present = 1, not present = 0) for each group
database$es_present <- ifelse(database$es_count > 1, 1, 0)
database$de_present <- ifelse(database$de_count > 1, 1, 0)

# Generate "Volcano" expression plot -------------------------------------------
# open png device
png(filename = "../RESULTS/volcano_plot.png",
    width = 1200,
    height = 1200,
    units = "px",
    pointsize =10,
    bg = "white")
par(lwd=6,
    mgp=c(9,3,0),
    # c(bottom, left, top, right)
    mar=c(12,15,1,1),
    font.lab=2,
    font.axis=1,
    cex.lab=5,
    font.main=2,
    cex.axis=5,
    bg="white")
# Axis Labels
ylab.name <- expression(paste("-log"[10],"(p-value)"))
main.lab <- ""
xlab.name <- expression(paste("Expression ratio DE/ES (log"[2]," FPKM)"))
# subset for plotting
# apply P-value cut-off
vplot <- subset(database, database$p <= pco)
# Gene expression must be present in at least one condition
vplot.1 <- subset(vplot[51:nrow(vplot),],
                  vplot$es_present ==1 | vplot$de_present == 1)
vplot.2 <- subset(database, database$p > pco)
# subset to genes expressed in both conditions
db <- subset(database, database$xlog2 != Inf|database$xlog2 != NA)
# subset top 20 up-regulated genes (formerly top 50)
db.50 <- db[1:20,]
# subset LGR5
lgr5 <- db[24,]
plot(vplot.1$xlog2,-log10(vplot.1$p),
     type="n",
     ylab=ylab.name,
     xlab=xlab.name,
     col=rgb(0,0,139,90,maxColorValue=255),
     ylim= c(0,3.5),cex=2,xlim= c(-10,12))
grid(lwd=3)
lines(vplot.1$xlog2,-log10(vplot.1$p),
      type="p",
      col=rgb(65,105,225,100,maxColorValue=255),
        ylim= c(0,5),
        cex=2,
        xlim= c(-10,12))
lines(vplot.2$xlog2,-log10(vplot.2$p),
        type="p",
        col=rgb(112,128,144,100,maxColorValue=255),
        cex=2)
# points for top 20 upregulated
lines(db.50$xlog2,-log10(db.50$p),
        type="p",
        col=rgb(178,34,35,150,maxColorValue=255),
        cex=2,
        pch=21,
        bg="white")
# pont for LGR5
lines(lgr5$xlog2,-log10(lgr5$p),
        type="p",
        col=rgb(178,34,35,150,maxColorValue=255),
        cex=4,
        pch=21,
        bg="red")
n <- nrow(vplot)
# add arrow for LGR5
arrows(lgr5$xlog2-0.5,-log10(lgr5$p)+0.1,
       x1 = lgr5$xlog2,y1 = -log10(lgr5$p),
       lwd =5, col= "blue",length =0.1)
#Add gene name labels to top 30 upregulated genes
text(lgr5$xlog2-0.5,-log10(lgr5$p)+0.1,
     labels = lgr5$gene_short_name,
     col = "black",cex = 4,font = 4,pos = 2)
text(db.50$xlog2[1:20],(-log10(db.50$p)[1:20]),
     labels = db.50$gene_short_name[1:20],
     col="black",cex=2,font=2)
# Close plotting device
dev.off()

# Supplemental Table 1 ---------------------------------------------------------
db <- subset(database, database$xlog2 != Inf|database$xlog2 != NA)
db.50 <- db[1:50,]
write.csv(db.50,file="../RESULTS/top50upregulated_ES&DE.csv")

# Import Curated list of Wnt targets -------------------------------------------
list <- read.csv("../DATA/Wnt_Gene_List.csv",header=FALSE)
wnt.genes <- list[,1]
# retrieve rows with gene names matching wnt,genes
genes <- grep(paste(wnt.genes,collapse="|"),database$gene_short_name)
wnt.expression <- as.data.frame(database[genes,])
wnt.expression <- wnt.expression[order(wnt.expression$xlog2),]
write.csv(wnt.expression, file = "../RESULTS/wnt_gene_expression.csv") # Supplemental Table 2
wnt.expression <- subset(wnt.expression, wnt.expression$xlog2 != abs(Inf))
# sort by absolute value of log2 expression ratio (biggest differences first)
wnt.expression <- wnt.expression[order(-abs(wnt.expression$xlog2)),]
# Plot expression in ggplot2 ---------------------------------------------------
library(ggplot2)
library(MASS)
library(scales)
library(grid)
png(filename = "../RESULTS/wnt_scatter.png",
      width = 1200,
      height = 1200,
      units = "px",
      pointsize =10,
      bg = "white") 
print(ggplot(data = wnt.expression) + 
          geom_abline(intercept = 0 , slope = 1,
                      color = "grey", size = 2, linetype = "dashed") +
          geom_point(shape = 21,alpha = I(0.95),
                     aes(de_mean, es_mean, size = -log(p,10),fill = xlog2)) +
          guides(fill = guide_colorbar(title = expression(paste("DE/ES (log"[2]," FPKM)"))),
                 size = guide_legend(title = expression(paste("-log"[10],"(p-value)")))) +
          scale_fill_gradient2(low = "blue", high = "red") +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          annotation_logticks(size = 3) +
          geom_text(data = wnt.expression[c(1:(grep("RSPO2",wnt.expression$gene_short_name)-1),(grep("RSPO2",wnt.expression$gene_short_name)+2):23),],
                    aes(de_mean,es_mean,label=gene_short_name),
                    size = 10, fontface = 2,
                    # nudge labels to improve visibility
                    vjust = ifelse(wnt.expression[c(1:(grep("RSPO2",wnt.expression$gene_short_name)-1),(grep("RSPO2",wnt.expression$gene_short_name)+2):23),]$gene_short_name == "TP53TG5" | wnt.expression[c(1:(grep("RSPO2",wnt.expression$gene_short_name)-1),(grep("RSPO2",wnt.expression$gene_short_name)+2):23),]$gene_short_name == "WIF1", -0.45, 0.5),
                    # move labels for 3 genes with 0 expression in DE out of view rather than plot names underneath y-axis (unreadable). Expression data available in Supplemental Table 2
                    hjust = ifelse(wnt.expression[c(1:(grep("RSPO2",wnt.expression$gene_short_name)-1),(grep("RSPO2",wnt.expression$gene_short_name)+2):23),]$de_mean < 1, 2, 0.5)) +
          # Offset RSPO2 to improve visibility
          geom_text(data = wnt.expression[grep("RSPO2",wnt.expression$gene_short_name),],
                    aes(de_mean,es_mean+7,label = gene_short_name),
                    size = 10,fontface=2) +
          # Offset RSPO1 to improve visibility
          geom_text(data = wnt.expression[grep("RSPO1",wnt.expression$gene_short_name),],
                    aes(de_mean+0.1,es_mean+2,label = gene_short_name),
                    size = 10,fontface=2) +
          scale_size(range=c(0.25,25)) +
          labs(x = "\nDefinitive Endoderm (Mean FPKM)",
               y = "ES cells (Mean FPKM)\n") +
          coord_fixed(ratio = 1) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position = "bottom",
                legend.key.size = unit(2, "cm"),
                legend.text = element_text(size = 20),
                legend.title = element_text(size = 20, face =2),
                axis.text = element_text(size=30, color = "black", face = 2),
                axis.title = element_text(size = 35, face = 2),
                panel.background = element_rect(fill = "grey85", color = "black"))
      ) # close print window
dev.off()
