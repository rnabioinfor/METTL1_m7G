#TCGA tRNA
library("TCGAbiolinks")
library(ggpubr)
library("data.table")
library("gridExtra")

misodifffiles = list.files("/n/data2/bch/hemonc/gregory/qiliu/projects/data_junho/GDCdata/", pattern = "*.mRNA_readcount_data.txt$", full.names = TRUE, ignore.case = TRUE)
objgene="METTL1"
dat = data.frame(ctype=NULL,exp=NULL,sample=NULL)
alldat = data.frame(ctype=NULL,exp=NULL,sample=NULL)
#statistic = data.frame()
for(infile in misodifffiles){
       aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
       name = sub(".*TCGA-(\\w+).mRNA_readcount_data.txt","\\1",infile)
       sampleBarcodes = colnames(aa)
       dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")
       dataSmNT <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "NT")
       expTP = aa[aa[,1]==objgene,dataSmTP]
       expNT = aa[aa[,1]==objgene,dataSmNT]
       a1 = data.frame(ctype=name,exp=median(as.numeric(expTP)))
       alldat = rbind(alldat,a1)
}
alldat <- alldat[!is.na(alldat$exp),]
tumorHigh <- as.vector(alldat[alldat[,2] > quantile(alldat$exp,.75,na.rm=T),1])
tumorLow <- as.vector(alldat[alldat[,2] < quantile(alldat$exp,.25,na.rm=T),1])

objtypes <- c("LGG","MESO","THCA","SARC","READ","ACC","BRCA", "KIRC","CHOL","GBM","SKCM","OV", "DLBC","UVM","TGCT","LUSC","THYM", "KICH","LAML","PRAD","UCEC", "STAD","ESCA","HNSC","LIHC","COAD","BLCA","KIRP","UCS", "PCPG","LUAD","CESC","PAAD")
m7GtRNA = c("Ala-AGC", "Ala-CGC", "Ala-TGC", "Arg-TCT", "Asn-GTT", "Cys-GCA", "Gly-ACC", "Ile-AAT", "Lys-CTT", "Lys-TTT", "Met-CAT", "Phe-GAA", "Thr-TGT", "Trp-CCA", "Tyr-GTA", "Pro-AGG", "Pro-CGG", "Pro-TGG", "Val-AAC", "Val-CAC", "Val-TAC","iMet-CAT")

highlist = list()
lowlist = list()
for (project in objtypes) {
       message(project)
       infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
         next
       }
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes
       dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")

       for(tRNA in m7GtRNA) {  
           expTP = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),dataSmTP]
           expTP <- colSums(expTP)
           if(project %in% tumorHigh) {
                highlist[[tRNA]] <- c(highlist[[tRNA]], as.numeric(expTP))
           } else if (project %in% tumorLow){
                lowlist[[tRNA]] <- c(lowlist[[tRNA]], as.numeric(expTP))
           }
       }
}

plotlist = list()
for(tRNA in m7GtRNA) {
    htl <- highlist[[tRNA]]
    a1 = data.frame(exp=as.numeric(htl),sample=rep("High",length(htl)))
    ltl <- lowlist[[tRNA]]
    a2 = data.frame(exp=as.numeric(ltl),sample=rep("Low",length(ltl)))
    dat = rbind(a1,a2)
    p <- ggboxplot(dat, x = 'sample', y = 'exp', color = 'sample', palette = "jco")+yscale("log2", .format = TRUE)+ labs(x=tRNA, y="log2(Normalized expression)") + stat_compare_means(aes(group = sample),label =  "p.signif",hide.ns = TRUE)
    #ggsave(paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",tRNA,".tumorHL.pdf"), width=10, height=6)
    plotlist[[tRNA]] <- p
}   
glist <- lapply(plotlist, ggplotGrob)
ggsave("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/m7G.tRNA.tumorHL.pdf", width = 16, height = 16, marrangeGrob(grobs = glist, nrow=4, ncol=6, top=NULL))

#correlation and highlow separately
library("TCGAbiolinks")
library(ggpubr)
library("data.table")
library("gridExtra")

misodifffiles = list.files("/n/data2/bch/hemonc/gregory/qiliu/projects/data_junho/GDCdata/", pattern = "*.mRNA_readcount_data.txt$", full.names = TRUE, ignore.case = TRUE)
objtypes <- c("LGG","MESO","THCA","SARC","READ","ACC","BRCA", "KIRC","CHOL","GBM","SKCM","OV", "DLBC","UVM","TGCT","LUSC","THYM", "KICH","LAML","PRAD","UCEC", "STAD","ESCA","HNSC","LIHC","COAD","BLCA","KIRP","UCS", "PCPG","LUAD","CESC","PAAD")
m7GtRNA = c("Ala-AGC", "Ala-CGC", "Ala-TGC", "Arg-TCT", "Asn-GTT", "Cys-GCA", "Gly-ACC", "Ile-AAT", "Lys-CTT", "Lys-TTT", "Met-CAT", "Phe-GAA", "Thr-TGT", "Trp-CCA", "Tyr-GTA", "Pro-AGG", "Pro-CGG", "Pro-TGG", "Val-AAC", "Val-CAC", "Val-TAC","iMet-CAT")

for(infile in misodifffiles){
       aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
       project = sub(".*TCGA-(\\w+).mRNA_readcount_data.txt","\\1",infile)
       message(project)
       gene1 = "METTL1"
       exp1 = aa[aa[,1]==gene1,]
       cct <- as.data.table(t(exp1),keep.rownames=T)       
       infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
         next
       }
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes       
       plotlist = list()
       for(tRNA in m7GtRNA) {
           expTP = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),]
           expTP <- colSums(expTP)
           expTP <- as.data.table(expTP,keep.rownames=T)
           df <- merge(expTP,cct,by="rn",all.x=T)
           setnames(df,c("sample","tRNA","METTL1"))
            rp = ggscatter(df, x = "METTL1", y = "tRNA", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + geom_density2d() + labs(x=paste0("log2 TPM (",gene1,")"),y=paste0("log2 TPM (",tRNA,")"),title=tRNA) +stat_cor(method = "pearson", label.x.npc = 0.2, label.y.npc = 0.8)
            plotlist[[tRNA]] = rp
       }
       glist <- lapply(plotlist, ggplotGrob)
       ggsave(paste0(project,".mett1_tRNA_correlation_.pdf"), marrangeGrob(grobs = glist, layout_matrix =matrix(1:25,  nrow = 5, ncol=5, byrow=TRUE)),width=20,height=20)
}

#high lower classification
for(infile in misodifffiles){
       aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
       project = sub(".*TCGA-(\\w+).mRNA_readcount_data.txt","\\1",infile)
       message(project)
       sampleBarcodes = colnames(aa)
       dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")
       gene1 = "METTL1"
       expTP = aa[aa[,1]==gene1,dataSmTP]       
       tumorHigh <- names(expTP)[expTP > quantile(expTP,.75,na.rm=T)]
       tumorLow <- names(expTP)[expTP < quantile(expTP,.25,na.rm=T)]
       
       infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
         next
       }
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes
       plotlist = list()
       for(tRNA in m7GtRNA) {
            tumorHigh <- tumorHigh[tumorHigh %in% colnames(aa)]
            expH = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),tumorHigh]
            expH <- colSums(expH)
            expH <- as.numeric(expH)
            tumorLow <- tumorLow[tumorLow %in% colnames(aa)]
            expL = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),tumorLow]
            expL <- colSums(expL)
            expL <- as.numeric(expL)
            a1 = data.frame(exp=as.numeric(expH),sample=rep("High",length(expH)))
            a2 = data.frame(exp=as.numeric(expL),sample=rep("Low",length(expL)))
            dat = rbind(a1,a2)            
            p <- ggboxplot(dat, x = 'sample', y = 'exp', color = 'sample', palette = "jco")+yscale("log2", .format = TRUE)+ labs(x=tRNA, y="log2(Normalized expression)") + stat_compare_means(aes(group = sample),label =  "p.signif",hide.ns = TRUE)
            plotlist[[tRNA]] <- p
        }
        glist <- lapply(plotlist, ggplotGrob)
        ggsave(paste0(project,".m7G.tRNA.tumorHL.pdf"), width = 16, height = 16, marrangeGrob(grobs = glist, nrow=4, ncol=6, top=NULL))
}

##high lower classification group
objcodons <- c("Pro-TGG","Met-CAT","Lys-TTT","Lys-CTT","Ile-AAT","Cys-GCA","Arg-TCT","Ala-AGC")
for(infile in misodifffiles){
       aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
       project = sub(".*TCGA-(\\w+).mRNA_readcount_data.txt","\\1",infile)
       message(project)
       sampleBarcodes = colnames(aa)
       dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")
       gene1 = "METTL1"
       expTP = aa[aa[,1]==gene1,dataSmTP]       
       tumorHigh <- names(expTP)[expTP > quantile(expTP,.75,na.rm=T)]
       tumorLow <- names(expTP)[expTP < quantile(expTP,.25,na.rm=T)]
       
       infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
         next
       }
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes
       aa <- aa[!grepl("-up",rownames(aa)),]
       aa <- aa[!grepl("-down",rownames(aa)),]
       
       codons <- gsub("-chr.*|-\\d+.*","",gsub(".*tRNA-","",rownames(aa)))
       
       codons<- codons[!duplicated(codons)]
       #codons <- gsub("SerTGA","Ser",codons)

       plotlist = list()
       datall = data.frame()
       for(tRNA in codons) {
            tumorHigh <- tumorHigh[tumorHigh %in% colnames(aa)]
            expH = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),tumorHigh]
            expH <- colSums(expH)
            expH <- as.numeric(expH)
            tumorLow <- tumorLow[tumorLow %in% colnames(aa)]
            expL = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),tumorLow]
            expL <- colSums(expL)
            expL <- as.numeric(expL)
            
            type=ifelse(tRNA %in% m7GtRNA, "M","O")
            
            a1 = data.frame(tRNA=tRNA, type=type, exp=as.numeric(expH),sample=rep("High",length(expH)))
            a2 = data.frame(tRNA=tRNA, type=type, exp=as.numeric(expL),sample=rep("Low",length(expL)))
            dat = rbind(a1,a2)
            datall <- rbind(datall,dat)
        }
        datall <- datall[rev(order(datall$type)),]
        datall <- datall[datall$type=="O"|datall$tRNA %in%objcodons, ]
        datall$tRNA =factor(datall$tRNA,level=unique(datall$tRNA))
        
        ggboxplot(datall, x = 'tRNA', y = 'exp', group='sample', color = 'sample', palette = "jco")+yscale("log2", .format = TRUE)+ labs(x="", y="log2(Normalized expression)") + theme(axis.text.x = element_text(angle = 60, hjust = 1))+ stat_compare_means(aes(group = sample),label =  "p.signif",hide.ns = TRUE)
        ggsave(paste0(project,".m7G.tRNA.tumorHL.pdf"), width = 15, height = 6)  

        aat <- setDT(datall)
        aat[,type:=NULL]
        aat <- aat[,median(exp),by=c("tRNA","sample")]
        aat[,type:=ifelse(tRNA%in%m7GtRNA, "M","O")]
        
        ggboxplot(aat, x = 'type', y = 'V1', group='sample', color = 'sample', palette = "jco")+yscale("log2", .format = TRUE)+ labs(x="", y="log2(Normalized expression)") + stat_compare_means(aes(group = sample),label =  "p.signif",hide.ns = TRUE)
        ggsave(paste0(project,".m7G.group.tumorHL.pdf"), width = 8, height = 6)
}


#correlation between Arg-TCT and METTL1
misodifffiles = list.files("/n/data2/bch/hemonc/gregory/qiliu/projects/data_junho/GDCdata/", pattern = "*.mRNA_readcount_data.txt$", full.names = TRUE, ignore.case = TRUE)
objtypes <- c("LGG","MESO","THCA","SARC","READ","ACC","BRCA", "KIRC","CHOL","GBM","SKCM","OV", "DLBC","UVM","TGCT","LUSC","THYM", "KICH","LAML","PRAD","UCEC", "STAD","ESCA","HNSC","LIHC","COAD","BLCA","KIRP","UCS", "PCPG","LUAD","CESC","PAAD")
plotlist = list()
for(infile in misodifffiles){
       aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
       project = sub(".*TCGA-(\\w+).mRNA_readcount_data.txt","\\1",infile)
       message(project)
       gene1 = "METTL1"
       exp1 = aa[aa[,1]==gene1,]
       cct <- as.data.table(t(exp1),keep.rownames=T)       
       infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/",project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
         next
       }
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes       
       tRNA ="Arg-TCT"
       expTP = aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),]
       expTP <- colSums(expTP)
       expTP <- as.data.table(expTP,keep.rownames=T)
       df <- merge(expTP,cct,by="rn",all.x=T)
       setnames(df,c("sample","tRNA","METTL1"))
       df[,tRNA:=log2(tRNA)]
       df <- df[!is.na(tRNA) & !is.na(METTL1)]
       df[,METTL1:=log2(as.numeric(METTL1))]
       rp = ggscatter(df, x = "METTL1", y = "tRNA", size = 2, alpha= 0.6, add = "reg.line", add.params = list(color = "#00AFBB", fill = "lightgray"), conf.int = TRUE) + geom_density2d() + labs(x=paste0("log2 TPM (",gene1,")"),y=paste0("log2 TPM (",tRNA,")"),title=project) +stat_cor(method = "pearson")
       plotlist[[project]] = rp
}

glist <- lapply(plotlist, ggplotGrob)
ggsave(paste0("mett1_tRNA_correlation_.pdf"), marrangeGrob(grobs = glist, layout_matrix =matrix(1:35,  nrow = 7, ncol=5, byrow=TRUE)),width=30,height=30)



#UGA frequency for top 1000 differential expressed genes in TCGA
#/n/data2/bch/hemonc/gregory/qiliu/O2/estb/tcga/codon_freq
library("TCGAbiolinks")
library(ggpubr)
library("data.table")
library(seqinr)
library(GenomicFeatures)
library(ggpubr)
library(data.table)
library(Biostrings)
library(cowplot)
library(gridExtra)

txdb <- loadDb(paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/ribo/db/annotation/hg38.gencode.sqlite"))
load(paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/ribo/db/annotation/hg38.gff.rda"))
load(paste0("/n/data2/bch/hemonc/gregory/qiliu/O2/ribo/db/annotation/hg38.txlens.rda"))

txlens <- txlens[cds_len>0]
txlens[, maxlen := max(cds_len), by = gene_id]
txlensMax <- txlens[cds_len==maxlen]
txlensMax[startsWith(gene_id,"ENS"),gene_id := sub("\\.\\d+$","",gene_id)]
geneinfo <- setDT(as.data.frame(gff))
geneinfo <- geneinfo[ type=="transcript" & transcript_type == "protein_coding"]
geneinfo <- geneinfo[,.(transcript_id,gene_id,gene_name)]
geneinfo[startsWith(gene_id,"ENS"),gene_id := sub("\\.\\d+$","",gene_id)]

txlens[,c('tx_id','gene_id','nexon') := NULL]
setnames(txlens,c('tx_name'),c('transcript_id'))
transcript_seqs <- read.fasta(paste0('/n/data2/bch/hemonc/gregory/qiliu/O2/ribo/db/mRNA/hg38.txdb.fa'), seqtype = 'DNA', as.string = T)
transcript_seqs <- data.table(transcript_id=names(transcript_seqs), seq=as.character(transcript_seqs))
transcript_seqs <- merge(x=transcript_seqs,y=txlens[,c('transcript_id','cds_len','utr5_len','tx_len')], by='transcript_id', all.x=T)
transcript_seqs <- transcript_seqs[geneinfo, on="transcript_id"]
transcript_seqs[,'cds' := mapply(function(seq,i,j) {substr(seq,i+1,i+j)}, seq, utr5_len, cds_len)]

transcript_seqs[,'codon_seqs' := mapply(function(seq,i,j) {strsplit(gsub("([[:alnum:]]{3})", "\\1 ", substr(seq,i+1,i+j)), ' ')[[1]]}, seq, utr5_len, cds_len)]
transcript_seqs[,c('seq','cds') := NULL]
transcript_seqs <- transcript_seqs[gene_id %in% txlensMax$gene_id]

data_summary <- function(codon_seqs, transcript_id, gene_id, gene_name) {
    temp <- table(codon_seqs)
    res = data.table(transcript_id=rep(transcript_id,length(temp)), gene_id=rep(gene_id,length(temp)), gene_name=rep(gene_name,length(temp)), codon=toupper(names(temp)),freq=as.vector(temp))
    return(res)
}
cct <- mapply(data_summary, transcript_seqs$codon_seqs, transcript_seqs$transcript_id, transcript_seqs$gene_id, transcript_seqs$gene_name, SIMPLIFY = FALSE)
rl <- rbindlist(cct)
rl[,sumv:=sum(freq),by='transcript_id']
rl[,perc:=round(freq/sumv,5),by='transcript_id']
rl[,freq2:=sum(freq),by='codon']
rl[,total:=sum(freq)]
#rl[, pvalue := mapply(function(freq,sumv,freq2,total) phyper(freq, freq2, total-freq2, sumv, lower.tail=F), freq, sumv, freq2, total)]
#rl[,padj := p.adjust(pvalue)]

objcodons <- "AGA"
objtypes <- c("LGG","MESO","THCA","SARC","READ","ACC","BRCA", "KIRC","CHOL","GBM","SKCM","OV", "DLBC","UVM","TGCT","LUSC","THYM", "KICH","LAML","PRAD","UCEC", "STAD","ESCA","HNSC","LIHC","COAD","BLCA","KIRP","UCS", "PCPG","LUAD","CESC","PAAD")
plotlist = list()
for(project in objtypes){
    infile <- paste0("/n/data2/bch/hemonc/gregory/qiliu/projects/data_junho/GDCdata/TCGA-",project,".mRNA_readcount_data.txt")
    aa = read.table(infile,sep="\t",header = T,as.is=T,check.names=FALSE)
    sampleBarcodes = colnames(aa)
    message(project)
    dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")
    dataSmNT <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "NT")
    if(length(dataSmNT) < 3) {
        next
    }
    genes <- unique(aa[,1])
    siggenes <- data.frame()
    othergenes <- data.frame()
    for(i in 1:length(genes)) {
        objgene <- genes[i]
        expTP = aa[aa[,1]==objgene,dataSmTP]
        expNT = aa[aa[,1]==objgene,dataSmNT]
        if(length(expTP) >= 3 & length(expNT) >= 3) {
            pvalue <- wilcox.test(as.numeric(expTP), as.numeric(expNT))$p.value
            if(is.na(pvalue)){
                next
            }
            if(pvalue < 0.001) {
                siggenes <- rbind(siggenes,data.frame(genes[i],pvalue))
            } else if (pvalue >= 0.05) {
                othergenes <- rbind(othergenes,data.frame(genes[i],pvalue))
            }
        }
    }
    siggenes <- as.data.table(siggenes)
    siggenes <- siggenes[order(pvalue)]
    othergenes <- as.data.table(othergenes)
    othergenes <- othergenes[rev(order(pvalue))]
    if(nrow(siggenes) == 0 | nrow(othergenes) == 0) {
        next
    }
    objgenes <- as.character(siggenes[1:1000,genes.i.])
    bggenes <- as.character(othergenes[1:1000,genes.i.])

    rl2 <- copy(rl)
    rl2[,objcodon:='F']
    if(length(objcodons) > 0) {
        rl2[codon %in% objcodons, objcodon := 'T']
    } else {
        rl2[, objcodon := 'T']
    }
    rl2 <- rl2[objcodon=='T',]
    rl2[gene_name %in% objgenes, objgene:='Input']
    rl2[gene_name %in% bggenes, objgene:='Background']
    rl2 <- rl2[!is.na(objgene)]
    #caa = rl2[,sum(perc),by=c('gene_name','objcodon')]
    #setnames(caa,'V1','freq')
    #caa = caa[objcodon=='T',]
    #caa <- merge(caa,unique(rl2[,list(gene_name, objgene)]),by="gene_name")
    p <- ggboxplot(rl2,x='objgene',y='perc', color='objgene', palette = "jco",title=project) + yscale("log2") + labs(x="", y="Codon frequency") + theme(axis.text.x = element_text(angle = 60, hjust = 1),legend.title = element_blank()) + stat_compare_means(aes(group = objgene),label =  "p.signif",hide.ns = TRUE)
    plotlist[[project]] = p
}

glist <- lapply(plotlist, ggplotGrob)
ggsave(paste0("tRNA.density.list.pdf"), width = 30, height = 30, marrangeGrob(grobs = glist, nrow=6, ncol=ceiling(length(plotlist)/6), top=NULL))

#TCGA tRNA dotplot
library("TCGAbiolinks")
library(ggpubr)
library("data.table")
library("gridExtra")
#objtypes <- c("BLCA", "BRCA", "ESCA","LUAD", "LUSC", "READ", "STAD")
objtypes <- c("LGG","MESO","THCA","SARC","READ","ACC","BRCA", "KIRC","CHOL","GBM","SKCM","OV", "DLBC","UVM","TGCT","LUSC","THYM", "KICH","LAML","PRAD","UCEC", "STAD","ESCA","HNSC","LIHC","COAD","BLCA","KIRP","UCS", "PCPG","LUAD","CESC","PAAD")
m7GtRNA = c("Ala-AGC", "Ala-CGC", "Ala-TGC", "Arg-TCT", "Asn-GTT", "Cys-GCA", "Gly-ACC", "Ile-AAT", "Lys-CTT", "Lys-TTT", "Met-CAT", "Phe-GAA", "Thr-TGT", "Trp-CCA", "Tyr-GTA", "Pro-AGG", "Pro-CGG", "Pro-TGG", "Val-AAC", "Val-CAC", "Val-TAC","iMet-CAT", "Gly-CCC", "Ile-TAT", "Ser-GGA", "Thr-AGT")
reslist = data.table()
for (project in objtypes) {
       message(project)
       infile <- paste0(project,'.expression.matrix.TMM')
       if(!file.exists(infile)) {
            next
       }
       dat = data.frame(ctype=NULL,exp=NULL,sample=NULL)
       alldat = data.frame(ctype=NULL,exp=NULL,sample=NULL)
       aa = read.table(infile,sep="\t",header = T,as.is=T, check.names=FALSE)
       sampleBarcodes = colnames(aa)
       sampleBarcodes = sub("^(\\w+\\.\\w+\\.\\w+\\.\\w+).*","\\1",sampleBarcodes)
       sampleBarcodes <- gsub("\\.","-",sampleBarcodes,perl=T)
       colnames(aa) <- sampleBarcodes
       dataSmTP <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "TP")
       dataSmNT <- TCGAquery_SampleTypes(barcode = sampleBarcodes, typesample = "NT")
       for(i in 1:length(m7GtRNA)) {
           #message(m7GtRNA[i])
           tRNA <- m7GtRNA[i]
           expTP = as.data.frame(aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),dataSmTP])
           expTP <- colMeans(expTP)
           a1 = data.frame(ctype=rep(tRNA,length(expTP)),exp=as.numeric(expTP),sample=rep("TP",length(expTP)))
           expNT = as.data.frame(aa[grep(paste0(tRNA,".*\\d$"),rownames(aa)),dataSmNT])
           expNT <- colMeans(expNT)
           if(length(expTP[!is.nan(expTP)]) > 1 & length(expNT[!is.nan(expNT)]) > 1){
                a2 = data.frame(ctype=tRNA,exp=as.numeric(expNT),sample=rep("NT",length(expNT)))
                a1 = rbind(a1,a2)
                dat = rbind(dat,a1)
                dat <- as.data.table(dat)
                pvalue <- wilcox.test(dat[sample=="TP",exp],dat[sample=="NT",exp])$p.value
                fc <- log2(mean(dat[sample=="TP",exp])/mean(dat[sample=="NT",exp]))
                temp <- data.table(project=project,tRNA=tRNA,pvalue=pvalue,log2fc=fc)
                reslist <- rbind(reslist,temp)
           }
       }
}
reslist[,pvalue := abs(log10(pvalue))]
pl1 <-  ggplot(reslist, aes(y = factor(tRNA),  x = factor(project))) +  geom_tile(aes(fill = log2fc)) +  scale_fill_continuous(low = "white", high = "white")
pl1 + geom_point(aes(colour = log2fc,  size =pvalue)) + scale_color_gradient(low = "red",high = "green")+ scale_size(range = c(0, 10))+ theme_linedraw() + labs(x = "TCGA cancer types", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.text.x = element_text(size = 9, angle = 60, hjust = 1, colour = "black"))
ggsave("heatmap_dotplot_colMean.pdf",width=8,height=8)