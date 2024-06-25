library(vioplot)
library(party)
library(stringr)

# read output file Exome RNASeq per patient simplified
ER <- read.delim("10.Exome.RNA-seq.Output.MS.MSI.txt", sep="\t", header=T)
freq <- read.delim("Frequence.global.10.2019.txt", sep="\t", header=T)
colnames(freq)[c(2:5)] <- c("freq.T.N.High.I","freq.T.N.low.I","freq.T.N.High.S","freq.T.N.low.S")

#filter values under 5%
freq[which(freq$freq.T.N.High.I<5),"freq.T.N.High.I"] <- 0
freq[which(freq$freq.T.N.low.I.I<5),"freq.T.N.low.I"] <- 0
freq[which(freq$freq.T.N.High.S<3),"freq.T.N.High.S"] <- 0
freq[which(freq$freq.T.N.low.S<3),"freq.T.N.low.S"] <- 0
ER <- merge(ER, freq, by="exon.id")
ER$id.exon.MS <- paste(ER$exon.id, ER$id.MS, sep=".")
ER$id.unique <- paste(ER$exon.id, ER$sample.MS.id, ER$shortMS,sep=".")

# Add Min Del to ER table
ER.mut <- read.delim("12.ER.min.del.annot.mut.txt", sep="\t", header=T)
ER.mut <- merge(ER.mut, freq, by="exon.id")
ER.mut$id.unique <- paste(ER.mut$exon.id, ER.mut$sample.MS.id, ER.mut$short.MS,sep=".")
smER <- data.frame(ER[,c("id.unique")])
colnames(smER) <- "id.unique"
smER.mut <- ER.mut[,c("id.unique","min.del")]
smER <- merge(smER, smER.mut,by="id.unique", all=T)
ER <- merge(ER, smER,by="id.unique", all=T)
ER$id.MS <- as.character(ER$id.MS)
#PSE devient PSI
ER$ratio.exclu.inclusion.tumoral <- 1-ER$ratio.exclu.inclusion.tumoral
ER$ratio.exclu.inclusion.normal <- 1-ER$ratio.exclu.inclusion.normal

ER$exon.id <- as.character(ER$exon.id)
library(stringr)
ER$exon.position <- data.frame(matrix(unlist(str_split(ER$exon.id[],'\\.')),
				      nrow=nrow(ER), byrow=T))[3]

#MS candidates
mut.dyn <- read.delim("96.fonctinnal.C.txt", sep='\t')
mut.dyn$id.MS <- as.character(mut.dyn$id.MS)
mut.dyn$couleur <- as.character(mut.dyn$couleur)

pdf("Supplementary.fig.96.MS.pdf")
par(mfrow = c(2,2))

for (i in mut.dyn$id.MS) {
  soustract <- ER[which(ER$id.MS==i),]
  titre1 <- paste("gene:",unique(soustract$gene.name),",exon:",unique(soustract$exon.position), ",type:",mut.dyn[which(mut.dyn$id.MS==i),"Type.of.Exon"])
  titre2 <- paste("MS:",unique(soustract$id.MS),",length:",unique(soustract$length),",distance:",unique(soustract$dist))
  boxplot(soustract$ratio.exclu.inclusion.tumoral~soustract$min.del, pch=4, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"], ylab="Tumoral PSI", xlab="Deletion size", main=titre1,sub=titre2, las=2,cex.main=0.8, outline=F, frame=FALSE)
  stripchart(soustract$ratio.exclu.inclusion.tumoral~soustract$min.del, cex=1, pch=19,vertical=T,method="jitter", add=T,col=c("grey"))
}

dev.off()

tab <- data.frame(matrix(ncol = 2, nrow = 7))
colnames(tab) <- c("deletion size", "PSI")
mut.dyn[which(mut.dyn$couleur=="#ffaa00ff"),"couleur"] <- "#fcbb36ff"
mut.dyn[which(mut.dyn$couleur=="#ffeeaaff"),"couleur"] <- "#95b9adff"

svg("All.data.V1.svg")
par(mfrow = c(1,1))
soustract <- ER[which(ER$id.MS== mut.dyn$id.MS[1]),]

#First Graph
for (j in -6:0) {
  tab[j+7,] <- c(j,mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==j)],
			na.rm=T))
}

plot(tab, pch=20, ylim=c(0,1), las=2, cex=0.7, frame=F, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"])
lines(tab, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"])

for (i in mut.dyn$id.MS)i {
  soustract <- ER[which(ER$id.MS==i),]
  for (j in -6:0){
    tab[j+7,] <- c(j,mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==j)],
			    na.rm=T))
  }
  points(tab, pch=20, ylim=c(0,1), las=2, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"], cex=0.7)
  lines(tab, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"])
}

####Mean of Means

for (j in -6:0){
  tab[j+7,] <- c(j,mean(ER[which(ER$min.del==j),"ratio.exclu.inclusion.tumoral"], na.rm=T)) 
}

points(tab, pch=20, ylim=c(0,1.5), las=2, cex=1.5, col="black")
lines(tab, col="black")

dev.off()

svg("All.data.V2.svg")
par(mfrow = c(1,1))
soustract <- ER[which(ER$id.MS== mut.dyn$id.MS[1]),]

#First Graph
for (j in -6:0) {
  tab[j+7,] <- c(j, mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==j)],
			 na.rm=T) / mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==0)], na.rm=T))
}

plot(tab, pch=20, ylim=c(0,1.5), las=2, cex=0.7, frame=F, col=mut.dyn[which(mut.dyn$id.MS==mut.dyn$id.MS[1]),"couleur"])
lines(tab, col=mut.dyn[which(mut.dyn$id.MS==mut.dyn$id.MS[1]),"couleur"])

for (i in mut.dyn$id.MS) {
  soustract <- ER[which(ER$id.MS==i),]
  for (j in -6:0) {
    tab[j+7,] <- c(j,mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==j)],
			  na.rm=T) / mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==0)], na.rm=T))
  }
  points(tab, pch=20, ylim=c(0,2), las=2, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"], cex=0.7, frame=F)
  lines(tab, col=mut.dyn[which(mut.dyn$id.MS==i),"couleur"])
}

####Mean des Means
all. <- all <- NULL
for (i in mut.dyn$id.MS) {
  all <- which(ER$id.MS==i)
  all. <- c(all., all)
}

soustract <- ER[all.,]

for (j in -6:0) {
  tab[j+7,] <- c(j,mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==j)],
			na.rm=T) / mean(soustract$ratio.exclu.inclusion.tumoral[which(soustract$min.del==0)], na.rm=T))
}

points(tab, pch=20, ylim=c(0,1.5), las=2, cex=1.5, col="black")
lines(tab, col="black")
dev.off()

alt <- mut.dyn[which(mut.dyn$Type.of.Exon=="Alternative"),"id.MS"]
con <- mut.dyn[which(mut.dyn$Type.of.Exon=="Constitutive"),"id.MS"]

alt.l. <- alt.l <- NULL
for (i in alt) {
  alt.l <- which(ER$id.MS==i)
  alt.l. <-c(alt.l.,alt.l )
}

con.l. <- con.l <- NULL
for (i in con) {
  con.l <- which(ER$id.MS==i)
  con.l. <-c(con.l.,con.l )
}

svg("20.boxplot.AlterVsConstit.svg")
par(mfrow = c(2,1))
boxplot(ER[alt.l.,"ratio.exclu.inclusion.tumoral"]~ER[alt.l.,"min.del"])
boxplot(ER[con.l.,"ratio.exclu.inclusion.tumoral"]~ER[con.l.,"min.del"])
dev.off()

