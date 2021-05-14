library(scater)
library('stringr')
library(devtools)
library(Biobase)
library(preprocessCore)

args = commandArgs(TRUE)
RC_file <- args[1] # your read counts file
TL_file <- args[2] # transcript length file

Exp <- read.table(RC_file,head=T,sep='\t',stringsAsFactors=F,row.names=1)
dat <- read.table(TL_file,head=T,sep='\t',stringsAsFactors=F,row.names=1)
Expression <- cbind(Exp[order(row.names(Exp)),],'Length'=dat[order(row.names(dat)),])
Len <- Expression$Length
tpm_table <- calculateTPM(as.matrix(Expression[,1:(ncol(Expression)-1)]), Len)
write.table(tpm_table,'TPM.txt',row.names=T,sep='\t',quote=F)

tissue_list <- unique(str_remove_all(as.character(colnames(Expression)), "_\\d$"))
tbl_median_counts <- data.frame(row.names = rownames(tpm_table))
for (i in 1:length(tissue_list)) {
  subdat <- tpm_table[,grep(colnames(tpm_table),pattern=tissue_list[i],fixed = TRUE)]
  if(!is.null(nrow(subdat))){
	  new_dat <- as.data.frame(rowMedians(subdat), optional = T)
	  colnames(new_dat) = paste(tissue_list[i],sep='')
	  tbl_median_counts <- cbind(tbl_median_counts, new_dat)
	}
}
write.table(tbl_median_counts,'Median_TPM.txt',sep='\t',quote=F)
