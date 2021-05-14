# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# BiocManager::install(c("edgeR", "limma"))
library("edgeR")
library("limma")
library('stringr')

setwd('D:\\Transformer_switchgrass\\02_RNA-seq_results')
counts = read.table('Other_samples_read_counts_2.txt',head=T,sep='\t',stringsAsFactors=F,row.names=1)
counts = counts[,order(colnames(counts))]
study = unique(str_remove_all(as.character(colnames(counts)), "__(.)*$"))
res <- c()
# for Alkali
subdat <- counts[,grep(colnames(counts),pattern='Alkali',fixed = TRUE)]
group <- factor(str_remove_all(as.character(colnames(subdat)), "_\\d$"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("Alkali__Control", "Alkali__alkali_6h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Alkali__alkali_6h_vs_Control',colnames(res1),sep='')
res <- res1

FG_de.com <- exactTest(d, pair=c("Alkali__Control", "Alkali__alkali_24h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Alkali__alkali_24h_vs_Control',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

# for Dehydration
subdat <- counts[,grep(colnames(counts),pattern='Dehydration',fixed = TRUE)]
group <- factor(str_remove_all(as.character(colnames(subdat)), "_\\d$"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("Dehydration__C1", "Dehydration__D1"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Dehydration__D1_vs_C1',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Dehydration__C1", "Dehydration__R1"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Dehydration__R1_vs_C1',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Dehydration__C1", "Dehydration__D2"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Dehydration__D2_vs_C1',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

# for Drought
subdat <- counts[,grep(colnames(counts),pattern='Drought',fixed = TRUE)]
group <- factor(str_remove_all(as.character(colnames(subdat)), "_\\d$"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("Drought__ADC_0d", "Drought__ADT_0d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_0d_vs_ADC_0d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Drought__ADC_6d", "Drought__ADT_6d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_6d_vs_ADC_6d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Drought__ADC_12d", "Drought__ADT_12d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_12d_vs_ADC_12d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Drought__ADC_18d", "Drought__ADT_18d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_18d_vs_ADC_18d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Drought__ADC_24d", "Drought__ADT_24d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_24d_vs_ADC_24d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Drought__ADC_30d", "Drought__ADT_30d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Drought__ADT_30d_vs_ADC_30d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

# for Salt
subdat <- counts[,grep(colnames(counts),pattern='Salt',fixed = TRUE)]
group <- factor(str_remove_all(as.character(colnames(subdat)), "_\\d$"))
d <- DGEList(counts=subdat,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("Salt__ASC_0h", "Salt__AST_0h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_0h_vs_ASC_0h',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_12h", "Salt__AST_12h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_12h_vs_ASC_12h',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_24h", "Salt__AST_24h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_24h_vs_ASC_24h',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_48h", "Salt__AST_48h"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_48h_vs_ASC_48h',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_6d", "Salt__AST_6d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_6d_vs_ASC_6d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_12d", "Salt__AST_12d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_12d_vs_ASC_12d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_18d", "Salt__AST_18d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_18d_vs_ASC_18d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])

FG_de.com <- exactTest(d, pair=c("Salt__ASC_24d", "Salt__AST_24d"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Salt__AST_24d_vs_ASC_24d',colnames(res1),sep='')
res <- cbind(res[order(rownames(res)),],res1[order(rownames(res1)),])
write.table(res,'FC_switchgrass_other_stress.txt',quote=F,sep='\t')

