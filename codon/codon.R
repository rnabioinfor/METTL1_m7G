library(seqinr)
library(GenomicFeatures)
library(ggpubr)
library(data.table)
library(Biostrings)
library(cowplot)
library(gridExtra)
library(gtools)

txdb <- loadDb(paste0("mm10.gencode.sqlite"))
load(paste0("mm10.gff.rda"))
load(paste0("mm10.txlens.rda"))
txlens <- txlens[cds_len>0]
txlens[, maxlen := max(cds_len), by = gene_id]
txlensMax <- txlens[cds_len==maxlen]
txlensMax <- unique(txlensMax,by="gene_id")
txlensMax[startsWith(gene_id,"ENS"),gene_id := sub("\\.\\d+$","",gene_id)]
geneinfo <- setDT(as.data.frame(gff))
geneinfo <- geneinfo[ type=="transcript" & transcript_type == "protein_coding"]
geneinfo <- geneinfo[,.(transcript_id,gene_id,gene_name)]
geneinfo[startsWith(gene_id,"ENS"),gene_id := sub("\\.\\d+$","",gene_id)]
txlens[,c('tx_id','gene_id','nexon') := NULL]
setnames(txlens,c('tx_name'),c('transcript_id'))
transcript_seqs <- read.fasta(paste0('mm10.txdb.fa'), seqtype = 'DNA', as.string = T)
transcript_seqs <- data.table(transcript_id=names(transcript_seqs), seq=as.character(transcript_seqs))
transcript_seqs <- merge(x=transcript_seqs,y=txlens[,c('transcript_id','cds_len','utr5_len','tx_len')], by='transcript_id', all.x=T)
transcript_seqs <- transcript_seqs[geneinfo, on="transcript_id"]
transcript_seqs <- transcript_seqs[cds_len%%3==0]
transcript_seqs[,'cds' := mapply(function(seq,i,j) {substr(seq,i+1,i+j)}, seq, utr5_len, cds_len)]
transcript_seqs[,'codon_seqs' := mapply(function(seq,i,j) {strsplit(gsub("([[:alnum:]]{3})", "\\1 ", substr(seq,i+1,i+j)), ' ')[[1]]}, seq, utr5_len, cds_len)]
transcript_seqs[,c('seq','cds') := NULL]
transcript_seqs[,transcript_id2 := sub("\\.\\d+$","",transcript_id)]
#transcript_seqs <- transcript_seqs[transcript_id %in% txlensMax$tx_name]
cano <- fread("vm23_canonical")
cano[,transcript2 := sub("\\.\\d+$","",transcript)]
transcript_seqs <- transcript_seqs[transcript_id2 %in% cano$transcript2]
degA <- fread(paste0("hitlist_heavy_to_light_setA.csv"),sep=",")
degB <- fread(paste0("hitlist_heavy_to_light_setB.csv"),sep=",")
degA <- degA[,.(gene_name=Gene.Symbol,logFC,pvalue_deg=P.Value,qvalue_deg=adj.P.Val)]
degB <- degB[,.(gene_name=Gene.Symbol,logFC,pvalue_deg=P.Value,qvalue_deg=adj.P.Val)]
degA[, updown := "none"]
#degA[logFC > -log2(1.2) & logFC < log2(1.2) & pvalue_deg > 0.05, updown := "none"]
degA[logFC <= -log2(1.2) & qvalue_deg < 0.05, updown := "down"]
degA[logFC >= log2(1.2) & qvalue_deg < 0.05, updown := "up"]
degB[, updown := "none"]
#degB[logFC > -log2(1.2) & logFC < log2(1.2) & pvalue_deg > 0.05, updown := "none"]
degB[logFC <= -log2(1.2) & qvalue_deg < 0.05, updown := "down"]
degB[logFC >= log2(1.2) & qvalue_deg < 0.05, updown := "up"]
upgenes <- intersect(degA[updown=="up",gene_name],degB[updown=="up",gene_name])
downgenes <- intersect(degA[updown=="down",gene_name],degB[updown=="down",gene_name])
#nongenes <- intersect(degA[updown=="none",gene_name],degB[updown=="none",gene_name])
nongenes <- degA[updown=="none",gene_name]
deg <- fread("proteomics.txt")
deg <- deg[,.(gene_name=ID,logFC,pvalue_deg=P.Value,qvalue_deg=adj.P.Val)]
deg[,updown := "none"]
deg[logFC >= log2(1.2) & qvalue_deg < 0.05, updown := "up"]
deg[logFC <= -log2(1.2) & qvalue_deg < 0.05, updown := "down"]
deg[,log2pvalue:=-log2(qvalue_deg)]
deg <- deg[gene_name != ""]
data_summary <- function(codon_seqs, transcript_id, gene_id, gene_name) {
    temp <- table(codon_seqs)
    res = data.table(transcript_id=rep(transcript_id,length(temp)), gene_id=rep(gene_id,length(temp)), gene_name=rep(gene_name,length(temp)), codon=toupper(names(temp)),freq=as.vector(temp))
    return(res)
}
transcript_seqs_obj <- transcript_seqs[gene_name %in% deg$gene_name]
cct <- mapply(data_summary, transcript_seqs_obj$codon_seqs, transcript_seqs_obj$transcript_id, transcript_seqs_obj$gene_id, transcript_seqs_obj$gene_name, SIMPLIFY = FALSE)
rl <- rbindlist(cct)
rl <- rl[!codon %in% c("TAG","TAA","TGA")]
#rl <- rl[nchar(codon) == 3]
rl[,sumv:=sum(freq),by='transcript_id']
rl[,perc:=round(freq/sumv,5),by='transcript_id']
rl[,freq2:=sum(freq),by='codon']
rl[,total:=sum(freq)]
rl[,per1k:=(freq*1000)/sumv]
msd <- rl[,.(meanvalue=mean(per1k),std=sd(per1k)),by=c("codon")]
rl <- merge(rl,msd,by="codon")
#rl[,zscore := mapply(function(q,m,n,k) {zscoreHyper(q,m,n,k)}, freq, freq2, total, sumv)]
rl[,zscore := mapply(function(q,m,n) {(q-m)/n}, per1k, meanvalue, std)]
#rl[,pvalue := mapply(function(q,m,n,k) {phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)}, freq, freq2, total, sumv)]
rl[, pvalue := mapply(function(freq,sumv,freq2,total) phyper(freq, freq2, total-freq2, sumv, lower.tail=F), freq, sumv, freq2, total)]
#rl[,padj := p.adjust(pvalue)]
#a<-fread("proteomic.txt",head=F)
#rl <- rl[gene_name %in% a$V1]
#rl[,codonType:="Non_m7G"]
m7Gcodon = c("GCT", "GCG", "GCA", "AGA", "AAC", "TGC", "GGT", "ATT", "AAG", "AAA", "ATG", "TTC", "ACA", "TGG", "TAC", "CCT", "CCG", "CCA", "GTT", "GTG", "GTA","ATG", "GGG", "ATA", "TCC", "ACT")
#rl[codon %in% m7Gcodon,codonType:="m7G"]
#fwrite(rl,file="proteomics.genes.codon.stat.xls",sep="\t")
rl_ud <- merge(rl,deg,all.y=T,by="gene_name")
rl_ud <- rl_ud[!is.na(codon)]
objcodon <- c("AGA")
objcodon <- objcodon
plotlist <- list()
objCodonStat <- rl[codon %in% objcodon]
objCodonStat2 <- objCodonStat[,sum(zscore),by=c("gene_name")]
setnames(objCodonStat2,"V1","per1k")
objCodonStat2 <- merge(objCodonStat2,deg,all.y=T,by="gene_name")
objCodonStat2 <- objCodonStat2[!is.na(per1k)]
objCodonStat2[,fcAndPvalue:=logFC*(-log2(pvalue_deg))]
rp2 <- ggscatter(objCodonStat2[per1k > -0.75], x = "per1k", y = "fcAndPvalue", title="AGA", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + labs(x="z-scores for AGA/AGG codons",y="Expression changes (-log2Pvalue x log2FC)") +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
objCodonStat2 <- objCodonStat2[per1k>0]
objCodonStat2[,per1k:=log2(per1k)]
my_comparisons <- list(c("up", "none"),c("down", "none"),c("down", "up"))
objCodonStat2$updown = factor(objCodonStat2$updown,c("down","none","up"))
rp3 <- ggboxplot(objCodonStat2, x = 'updown', y = 'per1k', color = 'updown', palette = "jco") +labs(y="Z-scores of AGA codon (log2)",xlab="") + stat_compare_means(comparisons = my_comparisons,method="t.test",label =  "p.signif")+theme(legend.position = "none")
ggsave(rp3,file=paste0("boxplot_zscore_per1k.pdf"),width=3,height=5)   
rp5 <- ggecdf(objCodonStat2,x='per1k', color='updown', linetype = "updown", palette = "jco", legend.title = "") + labs(x="", y="Cumulative changes")
plotlist[[paste0('a')]] = rp2
plotlist[[paste0('b')]] = rp3
plotlist[[paste0('c')]] = rp5
glist <- lapply(plotlist, ggplotGrob)
ggsave(paste0("scatter_plot_zscore_per1k.pdf"), marrangeGrob(grobs = glist, layout_matrix =matrix(1:3,  nrow = 1,ncol=3, byrow=TRUE)),width=12,height=8)   
#scatter
res <- data.table()
allcodons <- unique(rl_ud$codon)
allcodons <- allcodons[!allcodons %in% c("TAG","TAA","TGA")]
for(i in allcodons) {
    objrl <- rl_ud[pvalue<0.01 & codon == i]
    nd <- nrow(objrl[updown=="down"])
    nu <- nrow(objrl[updown=="up"])
    if(nu>0) {
        duRatio=nu/nd
        downratio=nu/nrow(objrl)
        res <- rbind(res,data.table(codon=i,duRatio=duRatio,downratio=downratio,nu=nu))
    }
}
res[,type:="other"]
res[codon %in% m7Gcodon,type:="m7G"]
genomeDuRatio <- nrow(deg[updown=="up"])/nrow(deg[updown=="down"])
genomeDownratio <- nrow(deg[updown=="up"])/nrow(deg)
genomeup <- nrow(deg[updown=="up"])
fwrite(res,file="upRatio_downRatio.txt",sep="\t")
ggscatter(res, x = "downratio", y = "duRatio", color = "type", palette = "jco", size = 2, alpha= 0.6, label = "codon", repel = TRUE, label.select = m7Gcodon) + geom_hline(yintercept=genomeDuRatio, linetype="dashed", color = "red") + geom_vline(xintercept=genomeDownratio, linetype="dashed", color = "red")+xlab("% up-regulated proteins")+ylab("U/D ratio")
ggsave("codon_ratio_stat_percent.pdf")
ggscatter(res, x = "nu", y = "duRatio", color = "type", palette = "jco", size = 2, alpha= 0.6, label = "codon", repel = TRUE, label.select = m7Gcodon) + geom_hline(yintercept=genomeDuRatio, linetype="dashed", color = "red")+xlab("#up-regulated proteins")+ylab("U/D ratio")
ggsave("codon_ratio_stat_numberOfUp.pdf")  
#heatmap
testt <- as.data.frame(dcast(rl, gene_name~codon, value.var="zscore",fun.aggregate=sum))
rownames(testt) <- testt[,1]
testt <- testt[,-1]
library("pheatmap")
pheatmap(testt, scale="row",cluster_cols=T, color = colorRampPalette(colors = c("cyan","black","yellow"))(10),cluster_rows=T, show_rownames=F, fontsize_col = 5, filename="codon_zscore_heatmap_1.pdf", width=6, height=8)
#dendrogram
library(dendextend)
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
dist_method <- c("maximum")
dist_method <- c("euclidean","maximum","manhattan","canberra","binary","minkowski")
testt <- as.data.frame(dcast(rl_ud, gene_id~codon,value.var="zscore",fun.aggregate=sum))
rownames(testt) <- testt[,1]
testt <- testt[,-1]
pdf("test.pdf",width=10,height=6)
#par(mfrow=c(1, 2))
#highlight <- c("AGA","CAA","GAA","GGG","AGG","GAA")
highlight <- c("ATG","AGA","TGC","AAA","AAG","GCT","ACT","CCA")
for(i in hclust_methods) {
    for(j in dist_method) {
        message(i)
        dend <- as.dendrogram(hclust(dist(as.matrix(t(testt)),method=j), method=i))
        dend2 <- color_labels(dend, labels = highlight , col = 2)
        plot(dend2, main=paste0(j,"-",i), xlab=paste0(j,"-",i))
    }
}
dev.off()
#heatmap2
library("gplots")
library("heatmap3")
gAnnoData = data.table(codon=colnames(testt),tRNAseq="green")
gAnnoData[codon %in% highlight,tRNAseq:="red"]
gAnnoData[,m7G:="blue"]
gAnnoData[codon %in% m7Gcodon,m7G:="yellow"]
rAnnoData <- data.table(gene_id=rownames(testt),type="green")
rAnnoData <- merge(rAnnoData,unique(rl_ud[,.(gene_id,updown)]),all.y=T,by="gene_id")
rAnnoData <- unique(rAnnoData,by="gene_id")
rAnnoData[,proteomics:="gray"]
rAnnoData[updown=="up",proteomics:="red"]
rAnnoData[updown=="down",proteomics:="green"]
mycol <- colorRampPalette(colors = c("cyan","black","yellow"))(10)
pdf("test_heatmap.pdf",width=10,height=8)
for(i in c("ward.D","ward.D2")) {
    for(j in c("euclidean","maximum")) {
        message(i)
        hc <- hclust(dist(as.matrix(t(testt)),method=j), method=i)
        hr <- hclust(dist(as.matrix(testt),method=j), method=i)
        heatmap.3(as.matrix(testt), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram = "col", col=mycol, scale="row", density.info="none", trace="none",ColSideColors=as.matrix(gAnnoData[,.(tRNAseq,m7G)]),RowSideColors=as.matrix(t(rAnnoData[,.(proteomics)])), main=paste0(j,"-",i), xlab=paste0(j,"-",i),key=F,labRow=F)
        #heatmap3(as.matrix(testt), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), showColDendro=T, showRowDendro = F, col=mycol, scale="row", density.info="none", trace="none",ColSideColors=gAnnoData[,.(type,color)], main=paste0(j,"-",i), xlab=paste0(j,"-",i),key=F,RowSideLabs=F)
    }
}
dev.off()
objcodon <- 'AGA'
objCodonStat <- rl_ud[pvalue<0.01 & codon==objcodon]
nd <- nrow(objCodonStat[updown=="down"])
nu <- nrow(objCodonStat[updown=="up"])
objupratio <-nu/nd
objudratio <-nu/nrow(objCodonStat)
nonObjCodonStat <- rl_ud[pvalue>=0.01 & codon==objcodon]
res <- data.table()
for(i in seq(10000)) {
    temp <- nonObjCodonStat[sample(.N, nrow(objCodonStat), FALSE)]
    nd <- nrow(temp[updown=="down"])
    nu <- nrow(temp[updown=="up"])
    if(nd>0) {
        res <- rbind(res,data.table(rep=i,duRatio=nu/nd,udratio=nu/nrow(temp)))
    }
}
obs_pval1 <- 2*mean(res$udratio>=objudratio)
obs_pval2 <- 2*mean(res$duRatio>=objupratio)
objCodonStat2 <- rl_ud[codon==objcodon]
objCodonStat2[,type:="other"]
objCodonStat2[pvalue<0.01,type:="enriched"]
res2 <- data.table()
for(i in seq(0,max(objCodonStat2$per1k-5),5)) {
    temp <- objCodonStat2[per1k>=i & per1k <i+5]
    nd <- nrow(temp[updown=="down"])
    nu <- nrow(temp[updown=="up"])
    if(nd>0) {
        res2 <- rbind(res2,data.table(rep=i,duRatio=nu/nd,udratio=nu/nrow(temp)))
    }
}
pdf("udratio_distribution.pdf")
par(mfrow = c(2,2))
hist(res$udratio,col = "black", freq=FALSE, breaks = 100,xlim=c(0,0.8),xlab="% up-regulated proteins",ylab="Percentage",main=paste0("pvalue",":",obs_pval1))
abline(v=objudratio, col = "blue", lwd = 2)
hist(res$duRatio,col = "black", freq=FALSE, breaks = 100,xlab="U/D ratio",ylab="Percentage",main=paste0("pvalue",":",obs_pval2))
abline(v=objupratio, col = "blue", lwd = 2)
objCodonStat2[,fcAndPvalue:=logFC*(-log2(qvalue_deg))]
ggscatter(objCodonStat2, x = "per1k", y = "fcAndPvalue", color="type", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + labs(x="#AGA codon per 1000 codons",y="Expression changes (-log2Pvalue x log2FC)") +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
#objCodonStat2[,type2:="other"]
#objCodonStat2[updown=="up",type2:="up"]
ggline(res2[rep<=20], x = "rep", y = "duRatio",linetype="dotted",shape=18)+ylab("U/D ratio")+xlab("# AGA codon per 1000 codons")
dev.off()
data_summary1 <- function(codon_seqs, transcript_id, gene_id, gene_name) {
    #message(transcript_id)
    codon_seqs <- unlist(codon_seqs)
    res <- data.table()
    nl <- length(codon_seqs)
    nw <- floor(nl/15)*15
    for(i in seq(15)) {
        ui <- floor((nl-i+1)/15)*15
        temp <-data.table(codon=codon_seqs[i:ui])
        temp[,index:=.I]
        temp[,window := floor((index-1)/15)+1]
        temp <- temp[,.N,by=c("codon","window")]
        temp[,window:=window*15-15+i]
        res <- rbind(res,temp)
    }
    res[,transcript_id:=transcript_id]
    res[,gene_id:=gene_id]
    res[,gene_name:=gene_name]
    res[,nw:=nl]
    return(res)
}
#codon_seqs <- transcript_seqs$codon_seqs[1]; transcript_id<-transcript_seqs$transcript_id[1]; gene_id<-transcript_seqs$gene_id[1]; gene_name<-transcript_seqs$gene_name[1]
transcript_seqs_obj <- transcript_seqs[gene_name %in% deg$gene_name]
cct3 <- mapply(data_summary1, transcript_seqs_obj$codon_seqs, transcript_seqs_obj$transcript_id, transcript_seqs_obj$gene_id, transcript_seqs_obj$gene_name, SIMPLIFY = FALSE)
rl3 <- rbindlist(cct3)
objcodons <- c("aga")
for(gene in c("Hmga2","Ash2l","Cdk6")) {
    tempForPlot <- rl3[gene_name==gene & codon %in% objcodons]
    nw <- unique(tempForPlot$nw)
    tempForPlot <- tempForPlot[,.(value=sum(N)),by=c("gene_name","window")]
    tempForPlot[,mean:=mean(value),by="gene_name"]
    tempForPlot[,value:=value-mean]
    tempForPlot[,type:="non"]
    tempForPlot[value>0,type:="up"]
    tempForPlot[value<0,type:="down"]
    ggplot(data = tempForPlot,aes(x = window, y = value, fill = type))+ xlim(1,nw)+geom_area()+theme_classic()+guides(fill = FALSE)
    ggsave(paste0(gene,"_objcodon_plot.pdf"),width=8,heigh=6)
    #pdf(paste0(gene,"_objcodon_plot.pdf"))
    #plot(tempForPlot$window, tempForPlot$value, type="h",xlim=c(1,nw))
    #dev.off()
}
data_summary2 <- function(codon_seqs, transcript_id, gene_id, gene_name, runLen) {
    #message(transcript_id)
    codon_seqs <- unlist(codon_seqs)
    nl <- length(codon_seqs)
    temp <- data.table()
    for(i in seq(runLen)) {
        temp <- cbind(temp,codon_seqs[seq(i,nl-runLen+i,1)])
    }
    colheads <- paste0("v",seq(runLen))
    setnames(temp,colheads)
    temp[, codonRun := Reduce(function(...) paste(..., sep = ""), .SD), .SDcols = colheads]
    temp <- table(temp$codonRun)
    res = data.table(transcript_id=rep(transcript_id,length(temp)), gene_id=rep(gene_id,length(temp)), gene_name=rep(gene_name,length(temp)), codon=toupper(names(temp)),freq=as.vector(temp))
    return(res)
}
data_summary3 <- function(codon_seqs, transcript_id, gene_id, gene_name, runLen) {
    codon_seqs <- unlist(codon_seqs)
    nl <- length(codon_seqs)
    colheads <- paste0("v",seq(runLen))
    temp <- data.table()
    for(r in 1:100) {
        codon_seqs <- codon_seqs[sample(seq(length(codon_seqs)),length(codon_seqs),replace=F)]
        temptemp <- data.table()
        for(i in seq(runLen)) {
            temptemp <- cbind(temptemp,codon_seqs[seq(i,nl-runLen+i,1)])
        }
        setnames(temptemp,colheads)
        temp <- rbind(temp,temptemp)
    }
    temp[, codonRun := Reduce(function(...) paste(..., sep = ""), .SD), .SDcols = colheads]
	gene_id<-transcript_seqs$gene_id[1]; gene_name<-transcript_seqs$gene_name[1]; runLen <- 2
    temp <- table(temp$codonRun)
    temp <- temp/100
    res = data.table(transcript_id=rep(transcript_id,length(temp)), gene_id=rep(gene_id,length(temp)), gene_name=rep(gene_name,length(temp)), codon=toupper(names(temp)),freq=as.vector(temp))
    return(res)
}
objcodon <- list()
objcodon[[1]] <- c('AGA')
codonObjRun <- as.data.table(permutations(n=2,r=2,v=c("AGA","AGG"), repeats.allowed=T))
codonObjRun[,run:=paste0(V1,V2)]
#objcodon[[2]] <- codonObjRun$run
objcodon[[2]] <- 'AGAAGA'
codonObjRun <- as.data.table(permutations(n=2,r=3,v=c("AGA","AGG"), repeats.allowed=T))
codonObjRun[,run:=paste0(V1,V2,V3)]
#objcodon[[3]] <- codonObjRun$run
objcodon[[3]] <- paste0("AGA",allcodons,"AGA")
codonnames <- c('AGA','AGAAGA','AGANNNAGA')
numcodon <- c(1,2,3)
for(i in seq(length(numcodon))) {
    #codonObjRun=c(paste0("AGA"))
    codonObjRun <- objcodon[[i]]
    codonname <- codonnames[i]
    message(codonname)
    message(i)
    transcript_seqs_obj <- transcript_seqs[gene_name %in% deg$gene_name]
    cct5 <- mapply(data_summary2, transcript_seqs$codon_seqs, transcript_seqs$transcript_id, transcript_seqs$gene_id, transcript_seqs$gene_name, numcodon[i], SIMPLIFY = FALSE)
    rl5 <- rbindlist(cct5)
    rl5 <- merge(rl5,transcript_seqs[,.(transcript_id,cds_len)],all.x=T,by="transcript_id")
    rl5[,codon_len:=cds_len/3]
    rl5[,per1k:=(freq*1000)/codon_len]
    msd <- rl5[,.(meanvalue=mean(per1k),std=sd(per1k)),by=c("codon")]
    rl5 <- merge(rl5,msd,by="codon")
    rl5[,zscore := mapply(function(q,m,n) {(q-m)/n}, per1k, meanvalue, std)]
    forDotPlot <- merge(rl5,deg,all.x=T,by="gene_name",allow.cartesian=TRUE)
    forDotPlot <- forDotPlot[updown=="up"]
    forDotPlot[,type:="other"]
    forDotPlot[codon %in% codonObjRun,type:="AGA"]
    setorder(forDotPlot,by=-zscore)
    fwrite(forDotPlot[type=="AGA"],file=paste0("codonRun_list.",codonname,".txt"),sep="\t")
    #fwrite(forDotPlot,file=paste0("codonRun_list.",codonname,".txt"),sep="\t")
    #ggscatter(rl5, x = "per1k", y = "zscore", color="type", size = 2, alpha= 0.6) + labs(x="#AGA codon per 1000 codons",y="Expression changes log2(fold change)") +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
    #ggsave(paste0("test_all.pdf"), width=10, height=6)
    if(TRUE) {
        tempForPlot <- rl5[codon %in% codonObjRun]
        tempForPlot <- tempForPlot[,sum(per1k),by=c("gene_id")]
        tempForPlot[,type:="Actual"]
        cct6 <- mapply(data_summary3, transcript_seqs_obj$codon_seqs, transcript_seqs_obj$transcript_id, transcript_seqs_obj$gene_id, transcript_seqs_obj$gene_name, numcodon[i], SIMPLIFY = FALSE)
        rl6 <- rbindlist(cct6)
        rl6 <- merge(rl6,transcript_seqs_obj[,.(transcript_id,cds_len)],all.x=T,by="transcript_id")
        rl6[,codon_len:=cds_len/3]
        rl6[,per1k:=(freq*1000)/codon_len]
        msd <- rl6[,.(meanvalue=mean(per1k),std=sd(per1k)),by=c("codon")]
        rl6 <- merge(rl6,msd,by="codon")
        rl6[,zscore := mapply(function(q,m,n) {(q-m)/n}, per1k, meanvalue, std)]
    
        tempForPlot2 <- rl6[codon %in% codonObjRun]
        tempForPlot2 <- tempForPlot2[,sum(per1k),by=c("gene_id")]
        tempForPlot2[,type:="Random"]
        tempForPlot <- rbind(tempForPlot,tempForPlot2)
        tempForPlot[,V1:=log2(V1)]
        my_comparisons <- list(c("up", "none"),c("down", "none"),c("down", "up"))
        ggboxplot(tempForPlot, x = 'type', y = 'V1', color = 'type', palette = "jco") +labs(y="Z-score for codon run") + stat_compare_means(comparisons = my_comparisons,method = "t.test",label =  "p.signif")
        ggsave(paste0("all.upOther_AGAenriched.",codonname,".pdf"), width=10, height=6)
    }
    #rl5 <- merge(rl5,deg,all.y=T,by="gene_name",allow.cartesian=TRUE)
    #rl5 <- rl5[!is.na(codon)]
    tempForPlot <- rl5[codon %in% codonObjRun]
    tempForPlot <- tempForPlot[,sum(per1k),by=c("gene_name")]
    tempForPlot <- merge(tempForPlot,deg,all.y=T,by="gene_name")
    tempForPlot <- tempForPlot[!is.na(V1)]
    tempForPlot[,type:="Other"]
    tempForPlot[updown=="up",type:="up"]   
    tempForPlot[,V1:=log2(V1)]
    my_comparisons <- list(c("up", "none"),c("down", "none"),c("down", "up"))
    ggboxplot(tempForPlot, x = 'updown', y = 'V1', color = 'updown', palette = "jco") + labs(y="Z-score for codon run") + stat_compare_means(comparisons = my_comparisons,method = "t.test",label =  "p.signif")
    ggsave(paste0("enriched.upOther_all.",codonname,".pdf"), width=10, height=6)
    #tempForPlot <- tempForPlot[gene_name %in% objCodonStat$gene_name]
    #ggboxplot(tempForPlot, x = 'type', y = 'V1', color = 'type', palette = "jco") + labs(y="#codon run per 1k codons") + stat_compare_means(aes(group = type),label =  "p.signif", hide.ns = TRUE)
    #ggsave(paste0("enriched.upOther_AGAenriched.",codonname,".pdf"), width=10, height=6)
}