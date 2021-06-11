###########################################
# fischerlab.org
# Explore SILAC data
# ESF 2020-09-24
# Explore SILAC data
# Tmp script - no repository yet
###########################################


# Load dependencies:
#source("~/R-helper-functions/mypairs.R")
library("gplots")
library("stringr")
library("limma")
library("dplyr")
library("readxl")
library(data.table)
suppressPackageStartupMessages(library("gridExtra"))
library(ggpubr)

# Get and set working directory
# list PD protein file names from "proteins" subfolder
upgenes <- list()
downgenes <- list()
plotlist <- list()
results <- "results_no_less_3"
for(i in c("setA","setB","setC","setD")) {
	message(i)
	# Load data
	D.exp_A <- as.data.frame(read.delim(paste0("./data/", i, ".txt"), sep="\t", na = "NaN"))

	#### Explore dataset A ####
	#Filter data to include only data with a minimum of 2 peptides (3 would be better)
	sel <- D.exp_A$Rep1.Number.of.peptides >= 3 & D.exp_A$Rep2.Number.of.peptides >= 3 & D.exp_A$Rep3.Number.of.peptides >= 3
	### Total of 3741 proteins with min of 2 peptides
	D.exp_A <- D.exp_A[sel,]

	# Explore data quality
	# create matrix with heavy to light ratios D - include mean/median at this point
	D <- as.matrix(D.exp_A[,c(5,6,8,9,11,12)])
	rownames(D) <- D.exp_A$Protein.Id

	#Create log2 transformed matrix
	A <- log2(D)

	#Check cross-correlation
	# log2 transformed ratios
	#mypairs(A)
	#Modest correlation
	
	AA <- as.data.table(A)
	rp1 = ggscatter(AA, x = "Rep1.sum.H_L.ratio", y = "Rep2.sum.H_L.ratio", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + labs(x=paste0("Rep1.sum.H_L.ratio"),y=paste0("Rep2.sum.H_L.ratio")) +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
	rp2 = ggscatter(AA, x = "Rep2.sum.H_L.ratio", y = "Rep3.sum.H_L.ratio", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + labs(x=paste0("Rep2.sum.H_L.ratio"),y=paste0("Rep3.sum.H_L.ratio")) +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
	rp3 = ggscatter(AA, x = "Rep1.sum.H_L.ratio", y = "Rep3.sum.H_L.ratio", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + labs(x=paste0("Rep1.sum.H_L.ratio"),y=paste0("Rep3.sum.H_L.ratio")) +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
    plotlist[["a"]] = rp1
	plotlist[['b']] = rp2
	plotlist[['c']] = rp3
	glist <- lapply(plotlist, ggplotGrob)
    ggsave(paste0(results,"/scatter_",i,".pdf"), marrangeGrob(grobs = glist, layout_matrix =matrix(1:3,  nrow = 1,ncol=3, byrow=TRUE)),width=15,height=6)		
			
	## Data does not look terrible, but after applying filter for 2 peptides only left with few proteins. 
	## Maybe filter later for outliers
	## Correlation is modest. With 2 peptides no advantage for median so go with sum

	#subset to only sum
	A <- A[,c(1,3,5)]
	#normalize data
	plotDensities(A)
	Anorm <- normalizeBetweenArrays(A)
	plotDensities(Anorm)
	#mypairs(Anorm)

	### Apply moderated t test
	fit <- lmFit(Anorm)
	fit <- eBayes(fit)
	tt <- topTable(fit, number=nrow(Anorm))
	
	D.exp_AA <- as.data.table(D.exp_A,keep.rownames=T)
	D.exp_AA <- D.exp_AA[,.(ID=Protein.Id,Gene.Symbol,Description)]
	ttt <- as.data.table(tt,keep.rownames=T)
	setnames(ttt,"rn","ID")
	ttt <- merge(D.exp_AA,ttt,by="ID",all.y=T)
	write.csv(ttt, file = paste0(results,"/hitlist_heavy_to_light_",i,".csv"), row.names = T)

	upgenes[[i]] <- ttt[logFC >= log2(1.2) & adj.P.Val < 0.05,ID]
	downgenes[[i]] <- ttt[logFC <= -log2(1.2) & adj.P.Val < 0.05,ID]
	
	# Set thresholds for graphs
	pval <- 0.05
	lfc <- 1
	
	ttt <- ttt[!is.na(logFC),]
	hits <- ttt$logFC >= lfc & ttt$P.Value <= pval | ttt$logFC <= -lfc & ttt$P.Value <= pval
	maxFC <- abs(max(ttt$logFC))+2
	maxP <- -log10(min(ttt$P.Value))+2

	#  volcano plot
	pdf(file = paste0(results,"/volcano_heavy_to_light_",i,".pdf"), width = 6, height = 6)
	plot(x=ttt$logFC,
		 y=-log10(ttt$P.Value), # data in the limma table is already log2 transformed
		 ylab="-log10 P value", # x-axis label
		 xlab=paste("log2FC Heavy / Light"), # y-axis label
		 cex=0.8, # point size
		 xlim=c(-maxFC,maxFC),
		 ylim=c(0,maxP),
		 pch=20, # point style
		 col=ifelse(hits, "sienna1", "#22222222"))
	abline(h=-log10(pval), lty=2)  # add vertical dashed line for p-value cutoff     
	abline(v=c(-1,1)*(lfc), lty=2)  # add two horizontal dashed lines for log fold change cutoffs
	if (sum(hits)>0){  # adds gene symbol as a text label to points if they are beyond the cutoffs
	  text(y = -log10(ttt[hits,P.Value]),   # x coordinate of label
		   x = ttt[hits,logFC],   # y coordinate of label
		   adj = c(-0.25,0.55),   # adjustment to offset label from point
		   labels=ttt[hits,Gene.Symbol],  # character strings to be used as labels
		   col="black",   # text color
		   cex=0.6)}   # text size
	dev.off()


	#### Create a heatmap plot summarizing the data
	#create color palette from green to red

	hits2 <- rownames(A) %in% ttt[hits, ID] 

	#Create matrix for heatmap
	m <- A[hits2,]
	sel <- rownames(m) != ""
	m <- m[sel,]
	m <- m[!(is.infinite(m[,1]) | is.infinite(m[,2]) | is.infinite(m[,3])),]
	my_palette <- colorRampPalette(c("royalblue1", "black", "sienna1"))(n=299)
	### Heatmap for % RA
	col_breaks = c(seq(-0.75,-0.25, length=100),               # for green
				   seq(-0.24, 0.24, length=100),            # for black
				   seq(0.25, 0.75, length=100))             # for red
	pdf(file = paste0(results,"/heatmap_heavy_to_light_",i,".pdf"), width = 6, height = 6)
	library(ggplot2)
	library(gplots)
	heatmap.2(m,
			  na.rm = T,
			  #cellnote = format(round(a, 2), nsmall = 2),
			  main = "Relative abundance",
			  notecol = "black",
			  #density.info = "none",
			  trace = "none",
			  margins = c(12,9),
			  col=my_palette,
			  # breaks = col_breaks,
			  #hclustfun = function(x) hclust(a, method ="complete"),
			  dendrogram = c("both"))

	dev.off()
	rm(my_palette,col_breaks)
}


library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)
temp <- venn.diagram(list(setA_up=upgenes[["setA"]],setB_up=upgenes[["setB"]],setC_down=downgenes[["setC"]],setD_down=downgenes[["setD"]]),  fill = c("green","yellow","red","blue"), alpha = c(0.5, 0.5, 0.5, 0.5), cex = 1,cat.fontface = 2,lty =2, filename = NULL,category.names=c('setA_up','setB_up','setC_down','setD_down'))
plot.new()
pdf(paste0(results,"/updown.veen.pdf"), width = 6, height = 6)
grid.draw(temp)
dev.off()


