library("stringr")
library("vioplot")
library("ComplexHeatmap")
library("circlize")
library("vioplot")

MS <- read.delim("2.humanV3.intronic.txt", header=F)

colnames(MS) <- c("id","start","stop","bases.strand","repeat.length","flanking5","flanking3","genomic.position","refseq","ensembl","strand","gene.name","total.exons","exon.5","position.exon.5","exon.3","position.exon.3")
MS$position.exon.5 <- as.numeric(as.vector(MS$position.exon.5))
MS$position.exon.3 <- as.numeric(as.vector(MS$position.exon.3))
MS$exon.5 <- as.numeric(as.vector(MS$exon.5))
MS$exon.3 <- as.numeric(as.vector(MS$exon.3))

MS <-  MS[which(MS$position.exon.5<=50 | MS$position.exon.3>=-50),]

MS$id <- paste("chr", MS$id, sep=".")

MS$bases <- ""
MS[which(MS$strand == "+"),"bases" ] <- as.vector(MS[which(MS$strand == "+"), "bases.strand"])
MS[which(MS$strand == "-" & MS$bases.strand == "A"),"bases" ] <- "T"
MS[which(MS$strand == "-" & MS$bases.strand == "T"),"bases" ] <- "A"
MS[which(MS$strand == "-" & MS$bases.strand == "G"),"bases" ] <- "C"
MS[which(MS$strand == "-" & MS$bases.strand == "C"),"bases" ] <- "G"
MS <- MS[,-c(6,7,8)]
MS$type <- "none"
MS[which(MS$position.exon.5>50 & MS$position.exon.3>=-50),"type"] <- "acceptor"
MS[which(MS$position.exon.3<-50 & MS$position.exon.5<=50),"type"] <- "donor"
MS[which(MS$position.exon.5<=50 & MS$position.exon.3>=-50),"type"] <- "-/+"
MS.all <- MS

#Strange cases
MS <- MS[-which(MS$type=="none"),]
table(MS$type)
# -/+ donor acceptor	
# 2549 124522 165245

length(unique(MS[which(MS$type=="donor"),"refseq"]))
length(unique(MS[which(MS$type=="acceptor"),"refseq"]))

tp <- MS[which(MS$type=="-/+"),]

tp$type <- "acceptor"
MS[which(MS$type=="-/+"),"type"] <- "donor"

MS <- rbind(MS, tp)
table(MS$type)
# acceptor    donor 
# 167794 127071   

MS[which(MS$type=="acceptor"),"short.MS"] <- paste(MS[which(MS$type=="acceptor"),"position.exon.3"],
						   MS[which(MS$type=="acceptor"),"bases"],
						   MS[which(MS$type=="acceptor"),"repeat.length"],
						   sep="")
MS[which(MS$type=="donor"),"short.MS"] <- paste("+",
						MS[which(MS$type=="donor"),"position.exon.5"],
						MS[which(MS$type=="donor"),"bases"],
						MS[which(MS$type=="donor"),"repeat.length"],
						sep="")
write.table(MS, "MS.Candidats.Splicing.txt", sep="\t", row.names=F)

tp <- MS[which(MS$type=="acceptor"),]
table(tp$bases)
#   A      C      G      T 
#  20142  29449  10195 108008 

tp <- data.frame(table(MS$repeat.length))
write.table(tp,"3.repeat.length.txt", sep="\t")

############################################################
MS <- read.delim("MS.Candidats.Splicing.txt", sep="\t", header=T)

#only keep exon n and not n-1 
colnames(MS)[13] <- "nb.exon"
# $total.exons
colnames(MS)[10] <- "total.exon"

MS$dist<- "NA"
MS[which(MS$type=="donor"),"dist"] <- as.numeric(MS[which(MS$type=="donor"),"position.exon.5"])
MS[which(MS$type=="acceptor"),"dist"] <- as.numeric(-MS[which(MS$type=="acceptor"),"position.exon.3"] )
MS$dist <- as.numeric(MS$dist)
MS$type <- as.character(MS$type)
MS[which(MS$type=="donor"),"type"] <- "+"
MS[which(MS$type=="acceptor"),"type"] <- "-"

MS <- MS[,c("id","short.MS","bases.strand","repeat.length","type","dist","gene.name","strand","refseq","ensembl","nb.exon","total.exon","bases")]
write.table(MS, "MS.Candidats.Slicing.simple.version.txt", sep="\t", row.names=F)

MS. <- read.delim("MS.Candidats.Slicing.simple.version.txt", sep="\t", header=T)
nrow(MS.)
#294865
length(unique(MS.$id))
# 182803
length(unique(MS.$id))*100/nrow(MS.)
# 62%

############################################################
MS$exon.id <- paste(MS$refseq, MS$nb.exon, sep=".")
MS.a <- MS[which(MS$type=="-"),]

tp <- tp. <- NULL

for (i in unique(MS.a$exon.id)) {
  min. <- NULL
  min. <- min(MS.a[which(MS.a$exon.id==i),"dist"])
  tp <- MS.a[which(MS.a$exon.id==i & MS.a$dist== min.),]
  tp. <- rbind(tp., tp)
}

acceptor <- tp.
write.table(acceptor, "acceptorsave.txt", sep="\t", row.names=F)

MS.d <- MS[which(MS$type=="+"),]
tp <- tp. <- NULL

for (i in unique(MS.d$exon.id)) {
  min. <- NULL
  min. <- min(MS.d[which(MS.d$exon.id==i),"dist"])
  tp <- MS.d[which(MS.d$exon.id==i & MS.d$dist== min.),]
  tp. <- rbind(tp., tp)
}

donor <- tp.
write.table(donor, "donorsave.txt", sep="\t", row.names=F)

MS.closer <- rbind(donor, acceptor)
write.table(MS.closer, "MS.Candidats.Slicing.simple.version.CLOSER.txt", sep="\t", row.names=F)

nrow(MS.c)
# 229110

length(unique(MS.c$id))
#142400

######################################### 
##############  genomic MS   ############ 
#########################################

length(unique(MS.c$gene.name))
# 16834
length(unique(MS.c$refseq))
# 27267 
length(unique(MS.c$exon.id))
# 181755

length(unique(MS.c[which(MS.c$type=="-"),"gene.name"]))
# 15634
length(unique(MS.c[which(MS.c$type=="-"),"refseq"]))
# 25411
length(unique(MS.c[which(MS.c$type=="-"),"exon.id"]))
# 127016

length(unique(MS.c[which(MS.c$type=="+"),"gene.name"]))
# 15294
length(unique(MS.c[which(MS.c$type=="+"),"refseq"]))
# 24715 
length(unique(MS.c[which(MS.c$type=="+"),"exon.id"]))
# 102094

geno2 <-cbind(c(2052*100/18886, 1200*100/18886, 1504*100/18886, 14094*100/18886),
              c(142161*100/323916, 54739*100/323916, 79661*100/323916, 47355*100/323916))
colnames(geno2)=c("Genes \n n=18886", "Exons\n n=323916")
rownames(geno2) <- c("WithOut MSs","MS in flanking 5 seq","MS in flanking 3 seq","MS in flanking 5 & 3 seq")

svg("Genomics.Features.MS.2.svg")
barplot(geno2, col=c("#edf8fb","#b2e2e2","#66c2a4","#238b45"), ylab="Percentage of genomic features")
legend(2,100,legend=c("WithOut MSs","MS in flanking 5 seq","MS in flanking 3 seq","MS in flanking 5 & 3 seq"), pch=15,col=c("#edf8fb","#b2e2e2","#66c2a4","#238b45"),box.lty=0)
dev.off()

#######################################################################
#####################    MS.PER.Gene   ################################
#######################################################################

### Total GENE WITH/WITHOUT MS
EI <- read.delim("0.EXON.INTRO.LENGTH.txt", sep="\t", header=T)
tp <- data.frame(matrix(unlist(str_split(EI$gene,';')), nrow(EI), byrow=T))
tp$X1 <- gsub("gene_id ", "", tp$X1)
EI$gene.name <- tp$X1
### Gene +/-
tp <- tp. <- NULL

for (i in unique(MS.c$gene.name))i {
  tp <-c(i,as.numeric(length(unique(as.character(MS.c[which(MS.c$gene.name==i),"id"])))),
	 unique(MS.c[which(MS.c$gene.name==i),"total.exon"]))
  tp. <- rbind(tp., tp)
}

nbg <- data.frame(tp.)
colnames(nbg) <- c("gene.name", "nb.MS", "total.exon")
tp <- cbind(setdiff(EI$gene.name, as.vector(nbg$gene.name)),0,0)
colnames(tp) <- colnames(nbg)
tp <- data.frame(tp)
nbg$nb.MS <- as.numeric(as.vector(nbg$nb.MS))
nbg$gene.name <- as.vector(nbg$gene.name)
nbg$total.exon <- as.vector(nbg$total.exon)
tp$nb.MS <- as.numeric(as.vector(tp$nb.MS))
tp$gene.name <- as.vector(tp$gene.name)
tp$total.exon <- as.vector(tp$total.exon)
nbg <- rbind(tp,nbg)
par(mfrow = c(1,1))

svg("MS.PER.Gene.svg")
vioplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", ylab="Number of Microsatellites", main="gene")
boxplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", add=T, outline=F)
  axis(2, at=c(0,0.3010299957,0.4771212547,0.7781512504,1.0413926852,1.414973348,1.7075701761,2.0043213738), labels=c(0,1,2,5,10,25,50,100), las=2)
dev.off()

### Gene +
MS.cPlus<- MS.c[which(MS.c$type=="+"),]
tp <- tp. <- NULL

for (i in unique(MS.cPlus$gene.name)) {
  tp <-c(i,as.numeric(length(unique(as.character(MS.cPlus[which(MS.cPlus$gene.name==i),"id"])))),
	 unique(MS.cPlus[which(MS.cPlus$gene.name==i),"total.exon"]))
  tp. <- rbind(tp., tp)
}

nbg <- data.frame(tp.)
colnames(nbg) <- c("gene.name", "nb.MS", "total.exon")
tp <- cbind(setdiff(EI$gene.name, as.vector(nbg$gene.name)),0,0)
colnames(tp) <- colnames(nbg)
tp <- data.frame(tp)
nbg$nb.MS <- as.numeric(as.vector(nbg$nb.MS))
nbg$gene.name <- as.vector(nbg$gene.name)
nbg$total.exon <- as.vector(nbg$total.exon)
tp$nb.MS <- as.numeric(as.vector(tp$nb.MS))
tp$gene.name <- as.vector(tp$gene.name)
tp$total.exon <- as.vector(tp$total.exon)
nbg <- rbind(tp,nbg)
par(mfrow = c(1,1))

svg("MS.PER.Gene.donor.svg")
vioplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", ylab="Number of Microsatellites", main="gene.donor")
boxplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", add=T, outline=F)
axis(2, at=c(0,0.3010299957,0.4771212547,0.7781512504,1.0413926852,1.414973348,1.7075701761,2.0043213738), labels=c(0,1,2,5,10,25,50,100), las=2)
dev.off()

### Gene -
MS.cless<- MS.c[which(MS.c$type=="-"),]
tp <- tp. <- NULL

for (i in unique(MS.cless$gene.name)){
  tp <-c(i,as.numeric(length(unique(as.character(MS.cless[which(MS.cless$gene.name==i),"id"])))),
	 unique(MS.cless[which(MS.cless$gene.name==i),"total.exon"]))
  tp. <- rbind(tp., tp)
}

nbg <- data.frame(tp.)
colnames(nbg) <- c("gene.name", "nb.MS", "total.exon")
tp <- cbind(setdiff(EI$gene.name, as.vector(nbg$gene.name)),0,0)
colnames(tp) <- colnames(nbg)
tp <- data.frame(tp)
nbg$nb.MS <- as.numeric(as.vector(nbg$nb.MS))
nbg$gene.name <- as.vector(nbg$gene.name)
nbg$total.exon <- as.vector(nbg$total.exon)
tp$nb.MS <- as.numeric(as.vector(tp$nb.MS))
tp$gene.name <- as.vector(tp$gene.name)
tp$total.exon <- as.vector(tp$total.exon)
nbg <- rbind(tp,nbg)
par(mfrow = c(1,1))

svg("MS.PER.Gene.accep.svg")
vioplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", ylab="Number of Microsatellites", main="gene.accep")
boxplot(log(nbg$nb.MS+1,10), xaxt="n", yaxt="n", add=T, outline=F)
axis(2, at=c(0,0.3010299957,0.4771212547,0.7781512504,1.0413926852,1.414973348,1.7075701761,2.0043213738), labels=c(0,1,2,5,10,25,50,100), las=2)
dev.off()

nbg. <- nbg[-which(as.numeric(nbg$total.exon)==0),]
svg("Cor.nb.exon.nb.MS.svg")
plot(log(nbg.$nb.MS+1,10), log(as.numeric(nbg.$total.exon)+1,10), pch=19, col="grey", xlab="Number of Microsatellites", ylab="Number of Exons", xaxt="n", yaxt="n", xlim=c(0,2.5), ylim=c(0,2.5))
axis(1, at=c(0,0.3010299957,0.4771212547,0.7781512504,1.0413926852,1.414973348,1.7075701761,2.0043213738,2.3996737215), labels=c(0,1,2,5,10,25,50,100,250), las=1)
axis(2, at=c(0,0.3010299957,0.4771212547,0.7781512504,1.0413926852,1.414973348,1.7075701761,2.0043213738,2.3996737215), labels=c(0,1,2,5,10,25,50,100,250), las=2)

Reg <- lm(log(as.numeric(nbg.$total.exon)+1,10)~log(nbg.$nb.MS+1,10))
abline(Reg,lty=3 )
cor.test(log(as.numeric(nbg.$total.exon)+1,10),log(nbg.$nb.MS+1,10))
legend(1.7,0.2, legend = "p-val<2.2e-16; R:0.83", cex=0.8)
dev.off()

### Transcript
tp <- tp. <- NULL

for (i in unique(MS.c$refseq)){
  tp <-c(i,as.numeric(length(which(MS.c$refseq==i))),unique(MS.c[which(MS.c$refseq==i),"total.exon"]))
  tp. <- rbind(tp., tp)
}

nbt <- data.frame(tp.)
colnames(nbt) <- c("transcrit", "nb.MS", "total.exon")
nbt$nb.MS <- as.numeric(as.vector(nbt$nb.MS))

vioplot(log(nbt$nb.MS,10), xaxt="n", yaxt="n", ylab="Number of Microsatellites")
axis(2, at=c(0,0.30103,0.69897,1,1.39794,1.69897,2), labels=c(1,2,5,10,25,50,100), las=2)

MS.c$exon.id <- as.vector(MS.c$exon.id)
MS.c$id <- as.vector(MS.c$id)
tp <- table(MS.c$exon.id, MS.c$id)

MS.c$exon.. <- paste(MS.c$refseq, MS.c$nb.exon, sep=".")
length(unique(MS.c[which(MS.c$type=="-" & MS.c$dist<20),"exon.id"]))
#82526
length(unique(MS.c[which(MS.c$type=="-" & MS.c$dist<20),"gene.name"]))
#13865
length(unique(MS.c[which(MS.c$type=="-" & MS.c$dist<20),"refseq"]))
#22577

load("microsatellites.listred.RData")
table(microsred$repeat_unit_length)
# 1       2       3       4       5 
# 9407309  226212   31094   35517    5934 

microsred <- microsred[which(microsred$repeat_unit_length==1),]
length(microsred$repeat_unit_length)
# 9,407,309

#Table to annotate with Annovar
all.MS <-  microsred[,c(1,2,2)] 

all.MS[,4] <- "-"
all.MS[,5] <-  microsred[,c(8)] 
write.table(all.MS,"all.MS.txt", sep="\t", row.names = F, col.names = F, quote=F)

all.MS.annot <- read.delim("all.MS.annot.txt.hg19_multianno.txt", sep="\t", header=T)
all.MS.annot$id <- paste(all.MS.annot$Chr,all.MS.annot$Start, sep=".")
microsred$id <- paste(microsred$chromosome, microsred$location, sep=".")
which(all.MS.annot$id!=microsred$id)

all.MS <- cbind(all.MS.annot,microsred)

all.MS <- all.MS[,c(11,1,2,24,5,6,7,16,19)]
all.MS$genomic.type <- "NA"
# UTRs
all.MS[which(all.MS$Func.refGene=="UTR5;UTR3" | all.MS$Func.refGene== "UTR3" | all.MS$Func.refGene== "UTR5"),"genomic.type"] <- "UTRs"
# Not Use
all.MS[which(all.MS$Func.refGene== "ncRNA_splicing" | all.MS$Func.refGene== "ncRNA_UTR5" |all.MS$Func.refGene== "exonic;splicing" |  all.MS$Func.refGene== "splicing" | all.MS$Func.refGene=="ncRNA_intronic" | all.MS$Func.refGene== "ncRNA_exonic" | all.MS$Func.refGene== "ncRNA_exonic;splicing"| all.MS$Func.refGene== "upstream;downstream"),"genomic.type"] <- "notuse"
# Intronic
all.MS[which(all.MS$Func.refGene=="intronic"),"genomic.type"] <- "intronic"
#intergenic
all.MS[which(all.MS$Func.refGene== "downstream" | all.MS$Func.refGene== "upstream"| all.MS$Func.refGene== "intergenic"),"genomic.type"] <- "intergenic"
all.MS[which(all.MS$Func.refGene== "exonic"),"genomic.type"] <- "exonic"                           
 
barplot(table(all.MS[which(all.MS$genomic.type=="exonic"),"repeat_times"]), las=2, log="y", xlim=c(0,50))   
barplot(table(all.MS[which(all.MS$genomic.type=="intronic"),"repeat_times"]), las=2, log="y", xlim=c(0,50))   
barplot(table(all.MS[which(all.MS$genomic.type=="UTRs"),"repeat_times"]), las=2, log="y", xlim=c(0,50))   
barplot(table(all.MS[which(all.MS$genomic.type=="intergenic"),"repeat_times"]), las=2, log="y", xlim=c(0,50))   
barplot(table(all.MS[which(all.MS$Func.refGene=="intronic"),"repeat_times"]), las=2, log="y", xlim=c(0,50))   

#Exonic
exonic <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(exonic) <- 5:25
exonic[1,] <- 0                       
tp <- table(all.MS[which(all.MS$genomic.type=="exonic"),"repeat_times"])

for (i in 1:length(as.integer(names(tp)))){
  exonic[1,which(colnames(exonic[1,])==as.integer(names(tp))[i])] <- tp[i]
}

#intronic
intronic <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(intronic) <- 5:25
intronic[1,] <- 0                       
tp <- table(all.MS[which(all.MS$genomic.type=="intronic"),"repeat_times"])

for (i in 1:length(as.integer(names(tp)))){
  intronic[1,which(colnames(intronic[1,])==as.integer(names(tp))[i])] <- tp[i]
}

#UTR
UTR <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(UTR) <- 5:25
UTR[1,] <- 0                       
tp <- table(all.MS[which(all.MS$genomic.type=="UTRs"),"repeat_times"])

for (i in 1:length(as.integer(names(tp)))){
  UTR[1,which(colnames(UTR[1,])==as.integer(names(tp))[i])] <- tp[i]
}

#intergenic
intergenic <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(intergenic) <- 5:25
intergenic[1,] <- 0                       
tp <- table(all.MS[which(all.MS$genomic.type=="intergenic"),"repeat_times"])

for (i in 1:length(as.integer(names(tp)))){
  intergenic[1,which(colnames(intergenic[1,])==as.integer(names(tp))[i])] <- tp[i]
}

svg("DistributionMS.length.svg")
par(mfrow = c(2,2))
barplot(log(as.numeric(exonic[1,])+1), las=2, main="Exonic", col="#7fc97f", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)

barplot(log(as.numeric(UTR[1,])+1), las=2, main="UTRs", col="#beaed4", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)

barplot(log(as.numeric(intronic[1,])+1), las=2, main="Intronic", col="#fdc086", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)

barplot(log(as.numeric(intergenic[1,])+1), las=2, main="intergenic", col="#ffff99", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)
dev.off()

##### MS Candidates

#flanking MS donor
dMS <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(dMS) <- 5:25
dMS[1,] <- 0                       
tp <- table(MS.c[which(MS.c$type=="+"),"repeat.length"])

for (i in 1:length(as.integer(names(tp)))) {
  dMS[1,which(colnames(dMS[1,])==as.integer(names(tp))[i])] <- tp[i]
}

#flanking MS accept
aMS. <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(aMS.) <- 5:25
aMS.[1,] <- 0                       
tp <- table(MS.c[which(MS.c$type=="-" ),"repeat.length"])

for (i in 1:length(as.integer(names(tp)))) {
  aMS.[1,which(colnames(aMS.[1,])==as.integer(names(tp))[i])] <- tp[i]
}



#flanking MS accept all 20
pptMS <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(pptMS) <- 5:25
pptMS[1,] <- 0                       
tp <- table(MS.c[which(MS.c$type=="-" & MS.c$dist <=20),"repeat.length"])

for (i in 1:length(as.integer(names(tp)))) {
  pptMS[1,which(colnames(pptMS[1,])==as.integer(names(tp))[i])] <- tp[i]
}

#flanking MS accept
aMS <- data.frame(matrix(ncol = 21, nrow = 1))
colnames(aMS) <- 5:25
aMS[1,] <- 0                       
tp <- table(MS.c[which(MS.c$type=="-" & MS.c$dist >20),"repeat.length"])

for (i in 1:length(as.integer(names(tp)))) {
  aMS[1,which(colnames(aMS[1,])==as.integer(names(tp))[i])] <- tp[i]
}

svg("Distribution.MS.candidat.length.svg")
par(mfrow = c(2,2))

barplot(log(as.numeric(dMS[1,])+1), las=2, main="Flanking MS donor", col="#b2e2e2", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)
barplot(log(as.numeric(aMS[1,])+1), las=2, main="Flanking MS acceptor distant (d>20bp)", col="#66c2a4", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)
barplot(log(as.numeric(pptMS[1,])+1), las=2, main="Flanking MS ppt (d<20bp)", col="#238b45", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)
barplot(log(as.numeric(aMS.[1,])+1), las=2, main="Flanking MS acceptor distant ", col="#66c2a4", ylab="log10(number of MS)")
axis(1, at= seq(1,24.8, by=1.19), labels=5:25,cex.axis=0.6)
dev.off()


ad <- "#d7191c"
th <- "#fdae61"
cy <- "#a6d96a"
gu <- "#1a9641"

# 
MS <- read.delim("MS.Candidats.Slicing.simple.version.txt", sep="\t", header=T, as.is=T)
MS$exon.percent <- MS$nb.exon/MS$total.exon*100

###Position of MS in transcrit
svg("Localization.of.MS.in.gene.5.svg")
par(mfrow = c(1,1))
plot(1,5, xlim=c(0,12), col="white", ylim=c(-10,110), xaxt="n", las=2, ylab="position in gene (Percent)", main="Localization of MS in gene (MS length≥5)")
vioplot(MS[which(MS$bases.strand == "A" & MS$type =="-"),"exon.percent"], at=1, add=T, col=ad)
vioplot(MS[which(MS$bases.strand == "T" & MS$type =="-"),"exon.percent"], at=2, add=T, col=th)
vioplot(MS[which(MS$bases.strand == "C" & MS$type =="-"),"exon.percent"], at=3, add=T,col=cy)
vioplot(MS[which(MS$bases.strand == "G" & MS$type =="-"),"exon.percent"], at=4, add=T, col=gu)
boxplot(MS[which(MS$bases.strand == "A" & MS$type =="-"),"exon.percent"], at=1, add=T, col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "T" & MS$type =="-"),"exon.percent"], at=2, add=T, col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "C" & MS$type =="-"),"exon.percent"], at=3, add=T,col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "G" & MS$type =="-"),"exon.percent"], at=4, add=T, col="grey",boxwex=0.5, yaxt="n")
vioplot(MS[which(MS$bases.strand == "A" & MS$type =="+"),"exon.percent"], at=6, add=T, col=ad)
vioplot(MS[which(MS$bases.strand == "T" & MS$type =="+"),"exon.percent"], at=7, add=T, col=th)
vioplot(MS[which(MS$bases.strand == "C" & MS$type =="+"),"exon.percent"], at=8, add=T, col=cy)
vioplot(MS[which(MS$bases.strand == "G" & MS$type =="+"),"exon.percent"], at=9, add=T, col=gu)
boxplot(MS[which(MS$bases.strand == "A" & MS$type =="+"),"exon.percent"], at=6, add=T, col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "T" & MS$type =="+"),"exon.percent"], at=7, add=T, col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "C" & MS$type =="+"),"exon.percent"], at=8, add=T,col="grey",boxwex=0.5, yaxt="n")
boxplot(MS[which(MS$bases.strand == "G" & MS$type =="+"),"exon.percent"], at=9, add=T, col="grey",boxwex=0.5, yaxt="n")
legend(10,80, legend=c("Adenosine","Uracile","Cytosine","Guanine"), pch=15, col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
axis(1, at=c(2.5,7.5), labels=c("Acceptor","Donor"))
axis(2, at=c(-10,110), labels=c("Start","End"), las=2)
dev.off()

########################################################################

ad <- "#d7191c"
th <- "#fdae61"
cy <- "#a6d96a"
gu <- "#1a9641"

#### Acceptor site
Ac.A <- MS[which(MS$bases == "A" & MS$type =="-"),c(6,4)]
Ac.T <- MS[which(MS$bases == "T" & MS$type =="-"),c(6,4)]
Ac.C <- MS[which(MS$bases == "C" & MS$type =="-"),c(6,4)]
Ac.G <- MS[which(MS$bases == "G" & MS$type =="-"),c(6,4)]
Ac.All <- rbind(Ac.A,Ac.T, Ac.C, Ac.G)
#### Donor site
do.A <- MS[which(MS$bases == "A" & MS$type =="+"),c(6,4)]
do.T <- MS[which(MS$bases == "T" & MS$type =="+"),c(6,4)]
do.C <- MS[which(MS$bases == "C" & MS$type =="+"),c(6,4)]
do.G <- MS[which(MS$bases == "G" & MS$type =="+"),c(6,4)]

do.All <- rbind(do.A,do.T, do.C, do.G)

svg("1.total.events.svg")
total.events <- c(nrow(Ac.A), nrow(Ac.T),nrow(Ac.C), nrow(Ac.G),nrow(do.A), nrow(do.T),nrow(do.C), nrow(do.G))
names(total.events) <- c("Acceptor.A","Acceptor.T","Acceptor.C","Acceptor.G","Donor.A","Donor.T","Donor.C","Donor.G")
barplot(total.events, las=2, col=c(rep("chartreuse4",4), rep("indianred4",4)), ylab="nb of MS", main="Number of MS in 50 bp surrounding exons", log="y")
dev.off()

tab.Ac.A <- as.data.frame.matrix(table(factor(Ac.A$repeat.length,levels=5:40),factor(Ac.A$dist,levels=0:50)))
tab.Ac.T <- as.data.frame.matrix(table(factor(Ac.T$repeat.length,levels=5:40),factor(Ac.T$dist,levels=0:50)))
tab.Ac.C <- as.data.frame.matrix(table(factor(Ac.C$repeat.length,levels=5:40),factor(Ac.C$dist,levels=0:50)))
tab.Ac.G <- as.data.frame.matrix(table(factor(Ac.G$repeat.length,levels=5:40),factor(Ac.G$dist,levels=0:50)))
tab.do.A <- as.data.frame.matrix(table(factor(do.A$repeat.length,levels=5:40),factor(do.A$dist,levels=0:50)))
tab.do.T <- as.data.frame.matrix(table(factor(do.T$repeat.length,levels=5:40),factor(do.T$dist,levels=0:50)))
tab.do.C <- as.data.frame.matrix(table(factor(do.C$repeat.length,levels=5:40),factor(do.C$dist,levels=0:50)))
tab.do.G <- as.data.frame.matrix(table(factor(do.G$repeat.length,levels=5:40),factor(do.G$dist,levels=0:50)))
tab.do.all <- (tab.do.G + tab.do.C + tab.do.A + tab.do.T)
tab.Ac.all <- (tab.Ac.G + tab.Ac.C + tab.Ac.A + tab.Ac.T)

########### PLOT ###############

par(mfrow = c(10,10))

#### Nb distance Acceptor

par(mfrow = c(2,2))
barplot(apply(tab.Ac.A, 2, sum),col=ad, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Adenosine - Acceptor")
barplot(apply(tab.Ac.T, 2, sum),col=th, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Uracile - Acceptor (ppt)")
barplot(apply(tab.Ac.G, 2, sum),col=gu, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Guanine - Acceptor")
barplot(apply(tab.Ac.C, 2, sum),col=cy, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Cytosine - Acceptor (ppt)")

#### Nb distance Donor
barplot(apply(tab.do.A, 2, sum),col=ad, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Adenosine - donor")
barplot(apply(tab.do.T, 2, sum),col=th, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Urdoile - donor (ppt)")
barplot(apply(tab.do.G, 2, sum),col=gu, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Guanine - donor")
barplot(apply(tab.do.C, 2, sum),col=cy, xlab="distance", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Cytosine - donor (ppt)")

#### Nb repeat length Acceptor
temp  <- apply(tab.Ac.A, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=ad,xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Adenosine - Acceptor")

temp  <- apply(tab.Ac.T, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=th, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Uracile - Acceptor (ppt)")

temp  <- apply(tab.Ac.G, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=gu, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Guanine - Acceptor")

temp  <- apply(tab.Ac.C, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=cy, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Cytosine - Acceptor (ppt)")

#### Nb repeat length donor
temp  <- apply(tab.do.A, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=ad,xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Adenosine - donor")

temp  <- apply(tab.do.T, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=th, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Uracile - donor (ppt)")

temp  <- apply(tab.do.G, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=gu, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Guanine - donor")

temp  <- apply(tab.do.C, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",col=cy, xlab="repeat length", ylab="number of MS", las=2, cex.axis=0.8, cex.names=0.65, main="Cytosine - donor (ppt)")

#############################
###         HeatMap        ##
#############################
couleur= colorRamp2(c(0,1,3,20,50,300,1000), col=c("white","#fef0d9","#fdd49e","#fdbb84", "#fc8d59","#e34a33","#b30000"))

###############ACCEPTOR#################

jpeg("2.heatmap.Acceptor.A.jpeg")
Heatmap(tab.Ac.A,
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (adenosine)")
dev.off()

jpeg("2.heatmap.Acceptor.T.jpeg")
Heatmap(tab.Ac.T,
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (uracile)")
dev.off()

couleur= colorRamp2(c(0,0.1,0.3,2,5,30,100), col=c("white","#fef0d9","#fdd49e","#fdbb84", "#fc8d59","#e34a33","#b30000"))

jpeg("2.heatmap.Acceptor.C.jpeg")
Heatmap(tab.Ac.C,
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (cytosine)")
dev.off()

jpeg("2.heatmap.Acceptor.G.jpeg")
Heatmap(tab.Ac.G, 
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (guanine)")
dev.off()

################DONOR##################
jpeg("2.heatmap.donor.A.jpeg")
Heatmap(tab.do.A, 
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (adenosine)")
dev.off()

jpeg("2.heatmap.donor.T.jpeg")
Heatmap(tab.do.T, 
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (uracile)")
dev.off()

jpeg("2.heatmap.donor.C.jpeg")
Heatmap(tab.do.C, 
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (cytosine)")
dev.off()

jpeg("2.heatmap.donor.G.jpeg")
Heatmap(tab.do.G, cluster_rows = F, cluster_columns = F, col= couleur, 
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (guanine)")
dev.off()

### boxplot.Acceptor-length.vs.distance ###
par(mfrow = c(2,2))
boxplot(Ac.A$dist~factor(Ac.A$repeat.length, levels=5:40), col=ad,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Adenosine : Acceptor site')
boxplot(Ac.T$dist~factor(Ac.T$repeat.length, levels=5:40), col=th,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Uracile : Acceptor site (ppt)')
boxplot(Ac.G$dist~factor(Ac.G$repeat.length, levels=5:40), col=gu,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Guanine : Acceptor site')
boxplot(Ac.C$dist~factor(Ac.C$repeat.length, levels=5:40), col=cy,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='cytosine : Acceptor site (ppt)')

### boxplot.donor-length.vs.distance ####
boxplot(do.A$dist~factor(do.A$repeat.length, levels=5:40), col=ad,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Adenosine : donor site')
boxplot(do.T$dist~factor(do.T$repeat.length, levels=5:40), col=th,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Uracile : donor site (ppt)')
boxplot(do.G$dist~factor(do.G$repeat.length, levels=5:40), col=gu,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='Guanine : donor site')
boxplot(do.C$dist~factor(do.C$repeat.length, levels=5:40), col=cy,las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='cytosine : donor site (ppt)')

### boxplot.Acceptor-distance.vs.length ###
boxplot(Ac.A$repeat.length~Ac.A$dist , col=ad, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Adenosine : Acceptor site')
boxplot(Ac.T$repeat.length~Ac.T$dist , col=th, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Uracile : Acceptor site (ppt)')
boxplot(Ac.G$repeat.length~Ac.G$dist , col=gu, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Guanine : Acceptor site')
boxplot(Ac.C$repeat.length~Ac.C$dist , col=cy, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Cytosine: Acceptor site (ppt)')

### boxplot.donor-distance.vs.length ####
boxplot(do.A$repeat.length~do.A$dist , col=ad, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Adenosine : donor site')
boxplot(do.T$repeat.length~do.T$dist , col=th, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Uracile : donor site (ppt)')
boxplot(do.G$repeat.length~do.G$dist , col=gu, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Guanine : donor site')
boxplot(do.C$repeat.length~do.C$dist , col=cy, las=2, ylim=c(0,50), outline=F, xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Cytosine: donor site (ppt)')

### plot.Acceptor-distance.vs.length ###
plot(Ac.A$repeat.length~Ac.A$dist , col=ad, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Adenosine : Acceptor site')
plot(Ac.T$repeat.length~Ac.T$dist , col=th, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Uracile : Acceptor site (ppt)')
plot(Ac.G$repeat.length~Ac.G$dist , col=gu, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Guanine : Acceptor site')
plot(Ac.C$repeat.length~Ac.C$dist , col=cy, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Cytosine: Acceptor site (ppt)')

### plot.donor-distance.vs.length ####
plot(do.A$repeat.length~do.A$dist , col=ad, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Adenosine : donor site')
plot(do.T$repeat.length~do.T$dist , col=th, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Uracile : donor site (ppt)')
plot(do.G$repeat.length~do.G$dist , col=gu, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Guanine : donor site')
plot(do.C$repeat.length~do.C$dist , col=cy, las=2, ylim=c(0,50), xlim=c(5,40), ylab="repeat length", xlab="distance from exon", main='Cytosine: donor site (ppt)')

###################################
###############ACCEPTOR############
jpeg("2.heatmap.Acceptor.All.jpeg")
Heatmap(tab.Ac.all,
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (A/U/C/G)")
dev.off()

################DONOR##################
jpeg("2.heatmap.donor.All.jpeg")
Heatmap(tab.do.all,
        cluster_rows = F, cluster_columns = F, col= couleur,
        heatmap_legend_param = list(at=c(0,1,3,20,50,300,1000)),
        column_title="Distance (bp)", column_title_side = "bottom",
        row_title="Repeat length", row_title_side = "right",
        name="nb of MS (A/U/C/G)")
dev.off()

###Acceptor
boxplot(Ac.All$dist~factor(Ac.All$repeat.length, levels=5:40), las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='A/T/C/G : Acceptor site')
barplot(apply(tab.Ac.all, 2, sum), xlab="distance", ylab="number", las=2, cex.axis=0.8, cex.names=0.85, main='A/T/C/G : Acceptor site')
temp  <- apply(tab.Ac.all, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y",xlab="repeat length", ylab="number", las=2, cex.axis=0.8, cex.names=0.85, main='A/T/C/G : Acceptor site')
boxplot(Ac.All$repeat.length~Ac.All$dist , las=2, ylim=c(0,50), xlim=c(5,40), xlab="repeat length", ylab="distance from exon", main='A/T/C/G : Acceptor site')

###Donor
boxplot(do.All$dist~factor(do.All$repeat.length, levels=5:40), las=2, ylim=c(0,50), outline=F, xlim=c(1,35), xlab="repeat length", ylab="distance from exon", main='A/T/C/G  : donor site')
barplot(apply(tab.do.all, 2, sum), xlab="distance", ylab="number", las=2, cex.axis=0.8, cex.names=0.85, main='A/T/C/G  : donor site')
temp  <- apply(tab.do.all, 1, sum)
temp[which(temp == 0)] <- 0.1
barplot(temp, log="y", xlab="repeat length", ylab="number", las=2, cex.axis=0.8, cex.names=0.85, main='A/T/C/G  : donor site')
boxplot(do.All$repeat.length~do.All$dist , las=2, ylim=c(0,50), xlim=c(5,40), xlab="repeat length", ylab="distance from exon", main='A/T/C/G  : donor site')

################################################
########           SOMME                ########
################################################
########ACCEPTOR
Acc.All <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.all)) {
  for (j in 1:length(tab.Ac.all[k,])) {
    for (i in 1:rownames(tab.Ac.all[k,])) {
      Acc.All[k,i+j-1] <- Acc.All[k,i+j-1] + tab.Ac.all[k,j]
    }
  }
}


Acc.A <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.A)) {
  for (j in 1:length(tab.Ac.A[k,])) {
    for (i in 1:rownames(tab.Ac.A[k,])) {
      Acc.A[k,i+j-1] <- Acc.A[k,i+j-1] + tab.Ac.A[k,j]
    }
  }
}

Acc.T <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.T)) {
  for (j in 1:length(tab.Ac.T[k,])) {
    for (i in 1:rownames(tab.Ac.T[k,])) {
      Acc.T[k,i+j-1] <- Acc.T[k,i+j-1] + tab.Ac.T[k,j]
    }
  }
}

Acc.C <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.C)) {
  for (j in 1:length(tab.Ac.C[k,])) {
    for (i in 1:rownames(tab.Ac.C[k,])) {
      Acc.C[k,i+j-1] <- Acc.C[k,i+j-1] + tab.Ac.C[k,j]
    }
  }
}

Acc.G <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.G)) {
  for (j in 1:length(tab.Ac.G[k,])) {
    for (i in 1:rownames(tab.Ac.G[k,])) {
      Acc.G[k,i+j-1] <- Acc.G[k,i+j-1] + tab.Ac.G[k,j]
    }
  }
}

Acc.All <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.all)) {
  for (j in 1:length(tab.Ac.all[k,])) {
    for (i in 1:rownames(tab.Ac.all[k,])) {
      Acc.All[k,i+j-1] <- Acc.All[k,i+j-1] + tab.Ac.all[k,j]
    }
  }
}

##### Barplot
svg("3.barplot.Acceptor.somme.ALL.svg")
par(mfrow = c(1,1))
tp.all <- rev(apply(Acc.All[,1:51], 2, sum))
names(tp.all) <- -50:0
tp.A <-  rev(apply(Acc.A[,1:51], 2, sum))
names(tp.A) <- -50:0
tp.T <-  rev(apply(Acc.T[,1:51], 2, sum))
names(tp.T) <- -50:0
tp.G <-  rev(apply(Acc.G[,1:51], 2, sum))
names(tp.G) <- -50:0
tp.C <-  rev(apply(Acc.C[,1:51], 2, sum))
names(tp.C) <- -50:0
df.bar <- barplot(tp.all, las=2, col="grey90", ylim=c(0,35000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3)
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)
legend(2,35000, legend=c("Adenosine","Thymine","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

############## Donor
doo.A <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.A)) {
  for (j in 1:length(tab.do.A[k,])) {
    for (i in 1:rownames(tab.do.A[k,])) {
      doo.A[k,i+j-1] <- doo.A[k,i+j-1] + tab.do.A[k,j]
    }
  }
}

doo.T <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.T)) {
  for (j in 1:length(tab.do.T[k,])) {
    for (i in 1:rownames(tab.do.T[k,])) {
      doo.T[k,i+j-1] <- doo.T[k,i+j-1] + tab.do.T[k,j]
    }
  }
}

doo.C <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.C)) {
  for (j in 1:length(tab.do.C[k,])) {
    for (i in 1:rownames(tab.do.C[k,])) {
      doo.C[k,i+j-1] <- doo.C[k,i+j-1] + tab.do.C[k,j]
    }
  }
}

doo.G <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.G)) {
  for (j in 1:length(tab.do.G[k,])) {
    for (i in 1:rownames(tab.do.G[k,])) {
      doo.G[k,i+j-1] <- doo.G[k,i+j-1] + tab.do.G[k,j]
    }
  }
}

doo.All <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.all)) {
  for (j in 1:length(tab.do.all[k,])) {
    for (i in 1:rownames(tab.do.all[k,])) {
      doo.All[k,i+j-1] <- doo.All[k,i+j-1] + tab.do.all[k,j]
    }
  }
}

svg("3.barplot.Donor.somme.ALL.svg")
par(mfrow = c(1,1))
tp.all <- apply(doo.All[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(doo.A[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(doo.T[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(doo.G[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(doo.C[,1:51], 2, sum)
names(tp.C) <- 0:50

df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,35000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)

legend(2,35000, legend=c("Adenosine","Uracile","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

svg("3.barplot.Donor.somme.ALL.v2.svg")
par(mfrow = c(1,1))
tp.all <- apply(doo.All[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(doo.A[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(doo.T[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(doo.G[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(doo.C[,1:51], 2, sum)
names(tp.C) <- 0:50

df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,35000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)

legend(2,35000, legend=c("Adenosine","Uracile","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

Acc.All <- Acc.A <- Acc.T <- Acc.C <- Acc.G <- NULL
do.All <- do.A <- do.T <- do.C <- do.G <- NULL

################################################
########          Start / Stop          ########
################################################

svg("3.barplot.Acceptor.Stop.ALL.svg")
par(mfrow = c(1,1))
tp.all <- rev(apply(tab.Ac.all[,1:51], 2, sum))
names(tp.all) <- -50:0
tp.A <-  rev(apply(tab.Ac.A[,1:51], 2, sum))
names(tp.A) <- -50:0
tp.T <-  rev(apply(tab.Ac.T[,1:51], 2, sum))
names(tp.T) <- -50:0
tp.G <-  rev(apply(tab.Ac.G[,1:51], 2, sum))
names(tp.G) <- -50:0
tp.C <-  rev(apply(tab.Ac.C[,1:51], 2, sum))
names(tp.C) <- -50:0
df.bar <- barplot(tp.all, las=2, col="grey90", ylim=c(0,10000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3)
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)
legend(2,10000, legend=c("Adenosine","Thymine","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

export <- cbind(tp.all,tp.A, tp.T, tp.G, tp.C)
write.table(export,"3.Acceptor.Stop.All.txt", row.names=T)

svg("3.barplot.Donor.Start.ALL.svg")
par(mfrow = c(1,1))
tp.all <- apply(tab.do.all[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(tab.do.A[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(tab.do.T[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(tab.do.G[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(tab.do.C[,1:51], 2, sum)
names(tp.C) <- 0:50
df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,1500))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)
legend(2,1500, legend=c("Adenosine","Thymine","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

svg("3.barplot.Donor.Start.ALL.v2.svg")
par(mfrow = c(1,1))
tp.all <- apply(tab.do.all[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(tab.do.A[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(tab.do.T[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(tab.do.G[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(tab.do.C[,1:51], 2, sum)
names(tp.C) <- 0:50
df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,10000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)
legend(2,10000, legend=c("Adenosine","Thymine","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

#######################################################

Acc.All. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.all)) {
  for (i in 1:ncol(tab.Ac.all)) {
    Acc.All.[k,i + as.numeric(rownames(tab.Ac.all)[k])] <- tab.Ac.all[k,i]
  }
}


Acc.A. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.A)) {
  for (i in 1:ncol(tab.Ac.A)) {
    Acc.A.[k,i + as.numeric(rownames(tab.Ac.A)[k])] <- tab.Ac.A[k,i]
  }
}

Acc.T. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.T)) {
  for (i in 1:ncol(tab.Ac.T)) {
    Acc.T.[k,i + as.numeric(rownames(tab.Ac.T)[k])] <- tab.Ac.T[k,i]
  }
}

Acc.C. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.C)) {
  for (i in 1:ncol(tab.Ac.C)) {
    Acc.C.[k,i + as.numeric(rownames(tab.Ac.C)[k])] <- tab.Ac.C[k,i]
  }
}

Acc.G. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.Ac.G)) {
  for (i in 1:ncol(tab.Ac.G)) {
    Acc.G.[k,i + as.numeric(rownames(tab.Ac.G)[k])] <- tab.Ac.G[k,i]
  }
}

########DONOR
do.All. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.all)) {
  for (i in 1:ncol(tab.do.all)) {
    do.All.[k,i + as.numeric(rownames(tab.do.all)[k])] <- tab.do.all[k,i]
  }
}

do.A. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.A)) {
  for (i in 1:ncol(tab.do.A)) {
    do.A.[k,i + as.numeric(rownames(tab.do.A)[k])] <- tab.do.A[k,i]
  }
}

do.T. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.T)) {
  for (i in 1:ncol(tab.do.T)) {
    do.T.[k,i + as.numeric(rownames(tab.do.T)[k])] <- tab.do.T[k,i]
  }
}

do.C. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.C)) {
  for (i in 1:ncol(tab.do.C)) {
    do.C.[k,i + as.numeric(rownames(tab.do.C)[k])] <- tab.do.C[k,i]
  }
}

do.G. <-data.frame(matrix(0,ncol=90,nrow=36))
for(k in 1:nrow(tab.do.G)) {
  for (i in 1:ncol(tab.do.G)) {
    do.G.[k,i + as.numeric(rownames(tab.do.G)[k])] <- tab.do.G[k,i]
  }
}

svg("3.barplot.Acceptor.Start.ALL.svg")
par(mfrow = c(1,1))
tp.all <- rev(apply(Acc.All.[,1:51], 2, sum))
names(tp.all) <- -50:0
tp.A <-  rev(apply(Acc.A.[,1:51], 2, sum))
names(tp.A) <- -50:0
tp.T <-  rev(apply(Acc.T.[,1:51], 2, sum))
names(tp.T) <- -50:0
tp.G <-  rev(apply(Acc.G.[,1:51], 2, sum))
names(tp.G) <- -50:0
tp.C <-  rev(apply(Acc.C.[,1:51], 2, sum))
names(tp.C) <- -50:0
df.bar <- barplot(tp.all, las=2, col="grey90", ylim=c(0,7000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3)
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)
legend(2,7000, legend=c("Adenosine","Thymine","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

svg("3.barplot.Donor.Stop.ALL.svg")
par(mfrow = c(1,1))
tp.all <- apply(do.All.[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(do.A.[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(do.T.[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(do.G.[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(do.C.[,1:51], 2, sum)
names(tp.C) <- 0:50

df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,1500))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)

legend(2,1500, legend=c("Adenosine","Uracile","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

svg("3.barplot.Donor.Stop.ALL.svg")
par(mfrow = c(1,1))
tp.all <- apply(do.All.[,1:51], 2, sum)
names(tp.all) <- 0:50
tp.A <- apply(do.A.[,1:51], 2, sum)
names(tp.A) <- 0:50
tp.T <- apply(do.T.[,1:51], 2, sum)
names(tp.T) <- 0:50
tp.G <- apply(do.G.[,1:51], 2, sum)
names(tp.G) <- 0:50
tp.C <- apply(do.C.[,1:51], 2, sum)
names(tp.C) <- 0:50

df.bar <- barplot(tp.all, las=2, col="grey90",ylim=c(0,10000))
lines(x = df.bar, y = tp.A, col="#e66101", lty=2, lwd=3 )
lines(x = df.bar, y = tp.T, col="#d7191c", lty=2, lwd=3)
lines(x = df.bar, y = tp.C, col="#2c7bb6", lty=2, lwd=3)
lines(x = df.bar, y = tp.G, col="#1a9641", lty=2, lwd=3)

legend(2,10000, legend=c("Adenosine","Uracile","Cytosine","Guanine"), lty=2, lwd=3,col=c("#e66101","#d7191c","#2c7bb6", "#1a9641"),box.lty=0)
dev.off()

######## Size MS ###############
#######################
# Acceptor Site # All #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="-"),]
tp <- table(tp$repeat.length, tp$dist)

vec <- vec. <- NULL

svg("3.Acceptor.Site.all.distance.svg")
par(mfrow = c(1,1))
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,20), las=2, xlab="Distance from Exon", ylab="repeat length", main="All MS in acceptor Site")
for (j in 1:50) {
  for (i in 1:36) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n")
  vec. <- NULL
} 
dev.off()


svg("3.Acceptor.Site.Nucleos.distance.svg")
par(mfrow = c(2,2))
#######################
# Acceptor Site # thymine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="-" & tp$bases=="T"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Thymine in acceptor Site")
for (j in 1:50) {
  for (i in 1:20) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#d71719ff")
  vec. <- NULL
}

#######################
# Acceptor Site # Adenosine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="-" & tp$bases=="A"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Adénosine in acceptor Site")
for (j in 1:50) {
  for (i in 1:20) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#e66100ff")
  vec. <- NULL
} 

#######################
# Acceptor Site # Cytosine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="-" & tp$bases=="C"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Cytosine in acceptor Site")
for (j in 1:50) {
  for (i in 1:10) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#a6d968ff")
  vec. <- NULL
} 
#######################
# Acceptor Site # Guanine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="-" & tp$bases=="G"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Guanine in acceptor Site")
for (j in 1:50) {
  for (i in 1:10) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#19943fff")
  vec. <- NULL
} 

dev.off()

######## Size MS ###############
#######################
# Donor Site # All #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="+"),]
tp <- table(tp$repeat.length, tp$dist)

vec <- vec. <- NULL

svg("3.Donor.Site.all.distance.svg")
par(mfrow = c(1,1))
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,20), las=2, xlab="Distance from Exon", ylab="repeat length", main="All MS in Donor Site")
for (j in 1:50) {
  for (i in 1:10) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n")
  vec. <- NULL
} 
dev.off()

svg("3.Donor.Site.Nucleos.distance.svg")
par(mfrow = c(2,2))
#######################
# Donor Site # thymine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="+" & tp$bases=="T"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Thymine in Donor Site")
for (j in 1:50) {
  for (i in 1:20) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#d71719ff")
  vec. <- NULL
}

#######################
# Donor Site # Adenosine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="+" & tp$bases=="A"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Adénosine in Donor Site")
for (j in 1:50) {
  for (i in 1:20) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#e66100ff")
  vec. <- NULL
} 

#######################
# Donor Site # Cytosine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="+" & tp$bases=="C"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Cytosine in Donor Site")
for (j in 1:50) {
  for (i in 1:10) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#a6d968ff")
  vec. <- NULL
}

#######################
# Donor Site # Guanine #
#######################
tp <- MS.c
tp <- tp[which(tp$type=="+" & tp$bases=="G"),]
tp <- table(tp$repeat.length, tp$dist)
vec <- vec. <- NULL
plot(1,1, col="white", xlim=c(0,50)
     ,ylim=c(5,15), las=2, xlab="Distance from Exon", ylab="repeat length", main=" MS -Guanine in Donor Site")
for (j in 1:50) {
  for (i in 1:10) {
    vec <- rep(i+4,tp[i,j])
    vec. <- c(vec.,vec)
  }
  boxplot(vec., add=T, at=j, outline=F,xaxt="n", yaxt="n", col="#19943fff")
  vec. <- NULL
} 

dev.off()

#######################
tp <- MS.c

EI <- read.delim("0.EXON.INTRO.LENGTH.txt", sep="\t", header=T)
ffo  <- data.frame(matrix(unlist(str_split(as.character(EI$gene),' ')), nrow=nrow(EI), byrow=T))
colnames(ffo)[c(2,4)] <- c("gene.name.check", "transcript.id")

EI <- cbind(EI, ffo[,c(2,4)])      

EI$transcript.id <- gsub(";","",EI$transcript.id)
EI$exon.id <- paste(EI$transcript.id , EI$nb.exons, sep=".")
length(unique(EI$gene))
#29810

length(unique(EI$gene.name.check))
#18886

MS.c.e <- merge(EI, MS.c,  by.x="exon.id",by.y="exon.id",all=T)
MS.c.e$type <- as.character(MS.c.e$type)
MS.c.e[which(is.na(MS.c.e$type)),"type"] <- "No MS"

svg("3.Exon.Length.svg")
boxplot(log10(MS.c.e$length.exon)~MS.c.e$type, ylab="log10(Exon length) in bp", yaxt="n", xaxt="n",col=c("#edf8fb","#b2e2e2","#66c2a4"),ylim=c(0,3), outline=F)
axis(2, at= c(0,1 ,2,3,4), labels=c(0,10,100,1000,10000), las=2)    
axis(1, at= c(1 ,2,3), labels=c("MS in Intronic Flanking 5"," MS in Intronic Flanking 3","WithOut MS"), las=2)
dev.off()

svg("3.intro.n.plus.1.Length.svg")
boxplot(log10(MS.c.e$intro.n.plus.1)~MS.c.e$type, ylab="log10(Intron n+1 length) in bp", yaxt="n",ylim=c(0,9), xaxt="n",col=c("#edf8fb","#b2e2e2","#66c2a4"), outline=F)
axis(2, at= c(0,1 ,2,3,4,5,6,7,8,9), labels=c(0,10,100,1000,10000, 100000,1000000, 100000000, 100000000,1000000000), las=2)    
axis(1, at= c(1 ,2,3), labels=c("MS in Intronic Flanking 5"," MS in Intronic Flanking 3","WithOut MS"), las=2)
dev.off()

svg("3.intro.n.less.svg")
boxplot(log10(MS.c.e$intro.n.less.1)~MS.c.e$type, ylab="log10(Intron n-1 length) in bp", yaxt="n", xaxt="n",ylim=c(0,9),col=c("#edf8fb","#b2e2e2","#66c2a4"), outline=F)
axis(2, at= c(0,1 ,2,3,4,5,6,7,8,9), labels=c(0,10,100,1000,10000, 100000,1000000, 100000000, 100000000,1000000000), las=2)   
axis(1, at= c(1 ,2,3), labels=c("MS in Intronic Flanking 5"," MS in Intronic Flanking 3","WithOut MS"), las=2)
dev.off()

