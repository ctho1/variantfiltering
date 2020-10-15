#-----------------------------------------------------------------------------------
# Filter Variants
#                                                                     
# Christian Thomas
# c.thomas@uni-muenster.de                                                                 
# 
# 2020-09-11
#------------------------------------------------------------------------------------   

# Load libraries
library(openxlsx)
library(UpSetR)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
options(scipen=999)

# Read input from appreci8
input<-read.table("./appreci8/results_probably_true.txt",
                  header=T,sep="\t",stringsAsFactors=F) 

### 18148 Variants

# Tidy and add Sample Names
input$ExAC<-as.numeric(input$ExAC)
input$G1000<-as.numeric(input$G1000)
input$ESP6500<-as.numeric(input$ESP6500)
input$Called<-as.numeric(input$Called)

for(i in 1:length(input[,1])){
  input$Gene[i]<-strsplit(input$Gene[i],split=",",fixed=T)[[1]][1] 
}

sample_names<-read.xlsx("./appreci8/Sample_IDs.xlsx")
for(i in 1:length(input[,1])){
  input$Sample[i]<-sample_names[sample_names$Humangenetik%in%input$Sample[i],2]
}

# Prepare Appreci8 Calls for Cancervar annotation
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
df<-data.frame(CHROM=input$Chr,POS=input$Pos,ID=input$Sample,
               REF=input$Ref,ALT=input$Alt,QUAL=".",
               FILTER=".",INFO=input$Sample)
write.table(df,file="./appreci8/appreci8_variants_export.vcf",
            quote=FALSE,sep="\t",row.names=F,col.names=F)
# appreci8_variants_export.vcf is ready to be annotated with CancerVar using the following command:
# python path/to/CancerVar/CancerVar.py -b hg19 -i ./appreci8_variants_export.vcf --input_type=VCF -o ./output
# This command produces the output output.hg19_multianno.txt.cancervar

# Read CancerVar Output
cancervar<-read.delim("./appreci8/output.hg19_multianno.txt.cancervar",
                      header = T,sep = "\t", stringsAsFactors = F)

cancervar$Score<-NA
cancervar$Verdict<-NA

for(i in 1:length(cancervar[,1])){
  tmp<-strsplit(cancervar$CancerVar..CancerVar.and.Evidence[i],split="#")[[1]]
  tmp<-c(strsplit(tmp[1],split=": ")[[1]][2],strsplit(tmp[2],split=" ")[[1]][1])
  cancervar$Score[i]<-tmp[1]
  cancervar$Verdict[i]<-tmp[2]
}

cancervar$Score<-as.numeric(cancervar$Score)

# Add CancerVar Annotation to filtered Appreci8 Output
input$CancerVar_Verdict<-cancervar$Verdict
input$CancerVar_Score<-cancervar$Score

plot(density(input_filtered$CancerVar_Score))

# Variant Caller Statistics
caller<-input[,c(14:21)]
caller[is.na(caller)]<-0

tmp<-data.frame(Caller=colnames(caller),Variants=NA)

for(i in 1:length(colnames(caller))){
  tmp[i,2]<-sum(caller[,i])
}

mean(tmp$Variants)
min(tmp$Variants)
max(tmp$Variants)

ggplot(tmp, aes(x=Caller,y=Variants))+
  geom_bar(stat="identity",fill="cornflowerblue")

upset(caller,nsets=8,nintersects=20,order.by = "freq",
      sets.bar.color = "cornflowerblue")

# Filter Variants

# 1. Filter for Variants with AF<0.01 or NA in ExAC, ESP6500
input_filtered<-input[input$G1000<0.01|is.na(input$G1000)==TRUE,]
input_filtered<-input_filtered[input_filtered$ExAC<0.01|is.na(input_filtered$ExAC)==TRUE,]
input_filtered<-input_filtered[input_filtered$ESP6500<0.01|is.na(input_filtered$ESP6500)==TRUE,]

### 9139 Variants

# 2. Filter synonymous variants
input_filtered<-input_filtered[input_filtered$Type!="synonymous_variant",] #n=1
for(i in 1:20){
  input_filtered<-input_filtered[input_filtered$Type!=paste(c(rep("synonymous_variant,",i),
                                                              "synonymous_variant"),collapse=""),] 
}

### 7773 Variants

# 3. Filter Variants within repetetive Regions
repeats<-read.table("./Data/UCSC_RepeatMasker_hg19",header=F, stringsAsFactors = F)
repeats<-repeats[,-c(4,5,6)]
colnames(repeats)<-c("chrom","start","end")

ranges1<-GRanges(seqnames=repeats$chrom,
                 ranges=IRanges(start=repeats$start-10,end=repeats$end+10))
ranges2<-GRanges(seqnames=paste0("chr",input_filtered$Chr),
                 ranges=IRanges(start=input_filtered$Pos,end=input_filtered$Pos))

intersection<-as.data.frame(granges(intersect(ranges1,ranges2)))
filter<-as.character(paste0(sub("chr","",intersection$seqnames),":",
                            intersection$start,":",intersection$end))
keep<-rep(TRUE,length(input_filtered[,1]))

for(i in 1:length(input_filtered[,1])){
  tmp<-paste0(input_filtered$Chr[i],":",input_filtered$Pos[i],":",
              input_filtered$Pos[i])
  keep[i]<-ifelse(tmp%in%filter,FALSE,TRUE)
}
input_filtered<-input_filtered[keep,]

### 7469 Variants

# 4. Filter Variants occuring in > 4 Samples
mutation<-data.frame(Mutation=unique(input_filtered$Mutation),
                     Frequency=NA)

for(i in 1:length(mutation[,1])){
  mutation[i,2]<-length(input_filtered[as.character(mutation[i,1])==input_filtered$Mutation,1])
}

plot(density(mutation$Frequency),xaxt="n")
xtick<-seq(0, 25, by=1)
axis(side=1, at=xtick, labels = TRUE)
abline(v=5, col="blue")

keep<-as.character(mutation[mutation[,2]<=4,1])
input_filtered<-input_filtered[input_filtered$Mutation%in%keep,]

### 4134 Variants

# 5. Filter CancerVar Score <=3

plot(density(input_filtered$CancerVar_Score))
abline(v=4, col="blue")

input_cv_filtered<-input_filtered[input_filtered$CancerVar_Score>3,]

# Statistics
pathogenic<-input_cv_filtered[grepl("athogenic",input_cv_filtered$CancerVar_Verdict),]
variant_stats<-(input_cv_filtered %>% dplyr::group_by(Sample) %>% dplyr::tally())
median(variant_stats$n)
min(variant_stats$n)
max(variant_stats$n)
sd(variant_stats$n)
### 9 Variants classified as likely pathogenic / pathogenic

write.xlsx(input_cv_filtered,"./final_variants.xslx")

### 526 Variants

# 6. Filter Genes with mutations classified as uncertain significance that are mutated in <2 Samples
input_cv_filtered_uc<-input_cv_filtered[input_cv_filtered$CancerVar_Score<8,]
genes<-data.frame(Gene=unique(input_cv_filtered_uc$Gene),
                  Frequency=NA)

for(i in 1:length(genes[,1])){
  genes[i,2]<-length(input_cv_filtered_uc[as.character(genes[i,1])==input_cv_filtered_uc$Gene,1])
}

plot(density(genes$Frequency),xaxt="n")
xtick<-seq(0, 75, by=5)
axis(side=1, at=xtick, labels = TRUE)
abline(v=1, col="blue")

keep<-as.character(genes[genes[,2]>1,1])
input_cv_filtered1<-input_cv_filtered_uc[input_cv_filtered_uc$Gene%in%keep,]
length(unique(input_cv_filtered1$Gene))

### 152 Variants
