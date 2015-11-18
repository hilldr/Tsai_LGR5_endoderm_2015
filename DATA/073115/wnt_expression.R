database <- read.csv(file="./ES&DE_FPKM_database.csv",header = T)
list <- read.csv("./Wnt_Gene_List.csv",header=FALSE)
wnt.genes <- list[,1]
genes <- grep(paste(wnt.genes,collapse="|"),database$gene_short_name)
wnt.expression <- as.data.frame(database[genes,])
library(ggplot2)
library(MASS)
library(scales)
plot <- ggplot(data=wnt.expression) + 
geom_point(shape=21,alpha=I(0.8),aes(de_mean,es_mean,size=-log(p,10),fill=-log(p,10))) +
scale_fill_gradientn(colours = c("white","red"),limits = c(min(wnt.expression$xlog2), max(wnt.expression$xlog2)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
annotation_logticks(size = 1.5) +
geom_text(aes(de_mean,es_mean,label=gene_short_name,size=-log(p,10)),fontface=2) +
scale_size(range=c(0.25,8)) +
geom_abline(intercept=0,slope=1, color = "grey",size = 3, alpha=I(0.25)) +
labs(x="\nDefinitive Endoderm (Mean FPKM)", y="ES cells (Mean FPKM)\n") +
coord_fixed(ratio = 1) +
theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size=20, color = "black"),
      axis.title = element_text(size = 25),
      panel.background = element_rect(fill = "white", color = "black"))
pdf(file= "wnt_expression.pdf", width = 12, height = 12)
print(plot)
dev.off()
