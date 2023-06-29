sample="WT_60"

#path=paste("C:/Users/vlakuz/OneDrive - KI.SE/Dokument/Documents/PIC/RNAPII_",sample,"/",sep="")
extract_ordered_gene_list<-function (sample)
{
setwd(dir = paste("C:/Users/vlakuz/OneDrive - KI.SE/Dokument/Documents/PIC/RNAPII_",sample,"/",sep=""))
#library(rtracklayer)

#NO INPUT NORMALIZATION

PI<-read.table("strongest_PI.txt", header =TRUE)

prot_coding<-read.table("C:/Users/vlakuz/OneDrive - KI.SE/Dokument/Documents/genes_etc/gene_lists/ngsplot/ngsplot_genes_hgnc.txt")

#PI<-PI[PI$gene_symbol %in% prot_coding$V1,]

#PI$alt<-(PI$ChIPTSScount/(PI$TSSend-PI$TSSstart))/(PI$ChIPGBcount/(PI$GBend-PI$GBstart)) thats how PI is calculated
PI$Flag="Other"
PI$Flag[(PI$PI>1) & (PI$PI<3) & (PI$ChIPGBcount/((PI$GBend-PI$GBstart)/1000) > 3)]="Elongating"
#PI$Flag[(PI$PI>10)& (PI$ChIPGBcount/((PI$GBend-PI$GBstart)/1000) < 2) & (PI$ChIPTSScount) > 10]="Stalled"
PI$Flag[(PI$PI>10)& (PI$ChIPGBcount/((PI$GBend-PI$GBstart)/1000) < 1) & (PI$ChIPTSScount) >8]="Stalled"
PI$Flag[PI$ChIPTSScount==0]="Silent"


length(which(PI$Flag=="Elongating"))
length(which(PI$Flag=="Stalled"))
length(which(PI$Flag=="Silent"))
length(which(PI$Flag=="Other"))

PI[c("GB1","GB2")]=NA
PI$TSS1<-(PI$TSSstart+PI$TSSend)/2
PI$TSS2<-((PI$TSSstart+PI$TSSend)/2)+1 # chr TSS TSS+1

PI$GB1[PI$strand=="+"]<-PI$GBstart[PI$strand=="+"]-999  #coordinates in the file are shifted +1000 of GB to calculate PI, so to create gene beds we need to shift back
PI$GB1[PI$strand=="-"]<-PI$GBstart[PI$strand=="-"] #999 determined by comparing to gtfs

PI$GB2[PI$strand=="+"]<-PI$GBend[PI$strand=="+"]
PI$GB2[PI$strand=="-"]<-PI$GBend[PI$strand=="-"]+999


#PI$Length<-PI$GBend-PI$GBstart



#PI[which(PI$Flag=="Elongating"),c("chr","TSSstart","TSSend")]


##gene list
write.table(PI[which(PI$Flag=="Elongating"),]$gene_id, file = paste("Elongating_",sample,".txt",sep="") ,quote = FALSE, col.names = FALSE, row.names = FALSE)
#close(("Elongating_WT_60.txt"))
write.table(PI[which(PI$Flag=="Stalled"),]$gene_id,file = paste("Stalled_",sample,".txt",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(PI[which(PI$Flag=="Silent"),]$gene_id, file=paste("no_RNAPII_",sample,".txt",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE)
## TSS bed files for each type of PI
write.table(PI[which(PI$Flag=="Elongating"),c("chr","TSS1","TSS2","gene_symbol","transcript_id","strand")], file = paste("Elongating_",sample,"_tss.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
write.table(PI[which(PI$Flag=="Stalled"),c("chr","TSS1","TSS2","gene_symbol","transcript_id","strand")],file = paste("Stalled_",sample,"_tss.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
write.table(PI[which(PI$Flag=="Silent"),c("chr","TSS1","TSS2","gene_symbol","transcript_id","strand")], file= paste("no_RNAPII_",sample,"_tss.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
## Genebody bed files for each type of PI
write.table(PI[which(PI$Flag=="Elongating"),c("chr.1","GB1","GB2","gene_symbol","transcript_id","strand")], file = paste("Elongating_",sample,"_gb.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
write.table(PI[which(PI$Flag=="Stalled"),c("chr.1","GB1","GB2","gene_symbol","transcript_id","strand")],file = paste("Stalled_",sample,"_gb.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
write.table(PI[which(PI$Flag=="Silent"),c("chr.1","GB1","GB2","gene_symbol","transcript_id","strand")], file = paste("no_RNAPII_",sample,"_gb.bed",sep=""),quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
#close()





print("We are done for now")
}

####NORMALIZATION WITH INPUT

PI<-read.table("strongest_PI.txt", header =TRUE)
PI$Flag[(PI$PI>1) & (PI$PI<3) & (((PI$ChIPGBcount-PI$InputGBcount)/((PI$GBend-PI$GBstart)/1000)) > 5)]="Elongating"
PI$Flag[(PI$PI>10)& ((PI$ChIPGBcount-PI$InputGBcount)/((PI$GBend-PI$GBstart)/1000)<0.01)]="Stalled"
#PI$Flag[(PI$TSSAvg==0)]="No_RNAPII"
length(which(PI$Flag=="Stalled"))
length(which(PI$Flag=="Elongating"))
PI[which(PI$Flag=="Elongating"),]$gene_id

RPM<-read.table("tag_count_normalize.txt",header = TRUE)
gtf <- rtracklayer::import('C:/Users/vlakuz/OneDrive - KI.SE/Dokument/Documents/gene_lists/gtf/hg38.ncbiRefSeq.sorted.gtf')
gtf_df=as.data.frame(gtf)
gtf_df<-gtf_df[which(gtf_df$type=="transcript"),]
RPM_names<-merge(RPM,gtf_df,by.x="gene_id",by.y="transcript_id")
RPM_names<-RPM_names[which(RPM_names$ChIP_tss==0),]
no_RNAPI_names<-unique(RPM_names$gene_name)
names<-unique(RPM_names$gene_name)

PI$GB_norm<-(PI$ChIPGBcount-PI$InputGBcount)/((PI$GBend-PI$GBstart)/1000)
PI$TSS_norm<-(PI$ChIPTSScount-PI$InputTSScount)

length(intersect(no_RNAPI_names,PI[which(PI$Flag=="Silent"),]$gene_id))

write.table(PI[which(PI$Flag=="Elongating"),]$gene_id,file = "Elongating_WT_60.txt",quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(PI[which(PI$Flag=="Stalled"),]$gene_id,file = "Stalled_WT_60.txt",quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_RNAPI_names, "no_RNAPII.txt",quote = FALSE, col.names = FALSE, row.names = FALSE)

#PI$length<-(PI$GBend-PI$GBstart)/1000
#PI$ChIPGBper1000<-PI$ChIPGBcount/PI$length
#PI$InputGBper1000<-PI$InputGBcount/PI$length
#PI$PI_alt<-(PI$ChIPTSScount-PI$InputTSScount)/(PI$ChIPGBper1000-PI$InputGBper1000)
