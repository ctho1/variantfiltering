# Filter CancerVar annotated Variant Calls
library(openxlsx); library(GenomicRanges)

cancervarpath<-"./Data/CancerVar/"
files<-list.files(cancervarpath)

variants<-read.delim(paste0(cancervarpath,files[1]),header = T,sep = "\t", stringsAsFactors = F)
variants<-cbind(sample=sub(".hg38_multianno.txt.cancervar","",files[1]),variants[,1:14])

for(i in 2:length(files)){
  tmp<-read.delim(paste0(cancervarpath,files[i]),skip = 1,header=F,sep = "\t",stringsAsFactors = F)
  tmp<-cbind(sample=sub(".hg38_multianno.txt.cancervar","",files[i]),tmp[,1:14])
  names(tmp)<-names(variants)
  variants<-rbind(variants,tmp)
}

# Tidy input
for(i in 1:length(variants[,1])){
  variants[i,15]<-strsplit(variants[i,15],split="\\#|\\ EVS")[[1]][2]
}

variants<-variants[,-c(10,12)]
colnames(variants)<-c("Sample","Chr","Start","End","Ref","Alt","Gene","Gene_Loc","Exonic_Consequence","avsnp147",
                      "Mutation","ClinVar","CancerVar")

# Filter intronic mutations
variants_filtered<-variants[variants$Gene_Loc=="exonic",]

# Filter CancerVar "Benign/Likely_benign"
variants_filtered<-variants_filtered[variants_filtered$CancerVar!="Benign/Likely_benign",]

# Filter Mutations that occur in > 2 Samples
mutation<-data.frame(Mutation=unique(variants_filtered$Mutation),
                     Frequency=NA)

for(i in 1:length(mutation[,1])){
  mutation[i,2]<-length(variants_filtered[as.character(mutation[i,1])==variants_filtered$Mutation,1])
}

keep<-as.character(mutation[mutation[,2]<=2,1])

variants_filtered<-variants_filtered[variants_filtered$Mutation%in%keep,]

# Filter Variants within repetetive Regions
repeats<-read.table("./Data/simple_repeats.txt",header=F, stringsAsFactors = F)
repeats<-repeats[,-c(4,5)]
colnames(repeats)<-c("chrom","start","end")

ranges1<-GRanges(seqnames = repeats$chrom,
                 ranges = IRanges(start = repeats$start,
                                  end = repeats$end))
ranges2<-GRanges(seqnames = paste0("chr",variants_filtered$Chr),
                 ranges = IRanges(start = variants_filtered$Start,
                                  end = variants_filtered$End))

intersection<-as.data.frame(granges(intersect(ranges1,ranges2)))
filter<-as.character(paste0(sub("chr","",intersection$seqnames),":",intersection$start,":",intersection$end))
keep<-rep(TRUE,length(variants_filtered[,1]))

for(i in 1:length(variants_filtered[,1])){
  tmp<-paste0(variants_filtered$Chr[i],":",variants_filtered$Start[i],":",variants_filtered$End[i])
  keep[i]<-ifelse(tmp%in%filter,FALSE,TRUE)
}
variants_filtered<-variants_filtered[keep,]

# Filter polymorphic genes / genes with common artifacts
polymorphic_genes<-c("HLA-A","HLA-B","HLA-C","TTN","AR","HLA-DQA2", "MUC20", "FAM136A", "ANKRD36B", 
                     "MADCAM1", "PABPC1", "RPL21", "ZNF273", "AHNAK2", "CEL", "DUX4L3", "MTCH2", 
                     "PSPH", "TUBBP5", "ADAM29", "HLA-DQA1", "PRB4", "APOE", "GALNT9", "PRB2",
                     "ZNF141", "ANKRD30A", "HNRNPCL3", "TIGD1","NCOR2")

keep<-rep(TRUE,length(variants_filtered[,1]))
for(i in 1:length(variants_filtered[,1])){
  if(variants_filtered$Gene[i]%in%polymorphic_genes==TRUE){
    keep[i]<-FALSE
  }
}
variants_filtered<-variants_filtered[keep,]

write.xlsx(variants_filtered, file="./Supplementary Table 2.xlsx")