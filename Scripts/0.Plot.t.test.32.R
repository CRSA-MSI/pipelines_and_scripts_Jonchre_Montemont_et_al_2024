##########################################
####      Length/ Distance/pvalue    #####
##########################################
library("vioplot")
library("stringr")
library("eulerr")

EX <- read.delim("3.Exon.aberrant.txt", sep="\t", header=T)
EX$exon.id <- as.vector(EX$exon.id)
ER <- read.delim("RNASeq.Full.Output.txt", sep="\t", header=T)

colnames(ER)[c(17,26,18,27)] <- paste("Average",colnames(ER)[c(17,26,18,27)],sep=".") 
# Define minimal frequency
freq <- read.delim("Frequence.global.10.2019.txt", sep="\t", header=T)
freq$ratio2 <- 0
colnames(freq)[c(2:5)] <- c("freq.T.N.High.I","freq.T.N.low.I","freq.T.N.High.S","freq.T.N.low.S")

# Filter frequency below 5%
freq[which(freq$freq.T.N.High.I<5),"freq.T.N.High.I"] <- 0
freq[which(freq$freq.T.N.low.I.I<5),"freq.T.N.low.I"] <- 0
freq[which(freq$freq.T.N.High.S<3),"freq.T.N.High.S"] <- 0
freq[which(freq$freq.T.N.low.S<3),"freq.T.N.low.S"] <- 0

ER <- merge(ER, freq, by="exon.id")

######### EXON
ex <- read.table("0.exonsnb.bed" ,sep="\t", header=F)
colnames(ex) <- c("chr","Start","Stop","strand", "gene", "nb.exons","nb.total.exons")

ex$exon.id <- data.frame(str_split_fixed(as.character(ex$gene), " ", 4))[,4]
ex$exon.id <- gsub(";","",ex$exon.id)
ex$exon.id <- paste(ex$exon.id,ex$nb.exons, sep=".")
ex$length.exon <- ex$Stop-ex$Start+1

ex$exon.id. <- data.frame(str_split_fixed(as.character(ex$gene), " ", 4))[,4]
ex$exon.id. <- gsub(";","",ex$exon.id.)


ex <- ex[order(as.numeric(ex$nb.exons)),] 
ex <- ex[order((ex$exon.id.)),] 

#Compute Exon length
tp. <- tp <- NULL

for (i in 1:nrow(ex)){
  tp <- ex[i,"Stop"]-ex[i+1,"Start"]
  tp. <- c(tp., tp)
}

tp.A <- tp.
tp.B <- c("NA",tp.A[])

#Compute Intron length 
ex$intro.n.less.1 <-ex$intro.n.plus.1 <-NA
ex$intro.n.plus.1[which(ex$strand=="-")] <- tp.A[which(ex$strand=="-")]
ex$intro.n.less.1[which(ex$strand=="-")] <- tp.B[which(ex$strand=="-")]

tp. <- tp <- NULL

for (i in 1:nrow(ex)){
  tp <- ex[i+1,"Start"]-ex[i,"Stop"]
  tp. <- c(tp., tp)
}

tp.A. <- tp.
tp.B. <- c("NA",tp.A.)
ex$intro.n.plus.1[which(ex$strand=="+")] <- tp.A.[which(ex$strand=="+")]
ex$intro.n.less.1[which(ex$strand=="+")] <- tp.B.[which(ex$strand=="+")]
ex$intro.n.plus.1[which(ex$nb.exons==ex$nb.total.exons)] <-"NA"
ex$intro.n.less.1[which(ex$nb.exons==1)] <-"NA"

ex$rank <- ex$nb.exon/ex$nb.total.exons*100
ex$frameshift <- "pb" 
ex[which(ex$length.exon%%3==0),"frameshift"] <- "frameshift"
ex[which(ex$length.exon%%3!=0),"frameshift"] <- "Inframe"

################################
##### Without duplicated exon.id
ER.. <- ER[-which(duplicated(ER$exon.id)),]
ER.. <- merge(ER..,ex, by="exon.id", all.x=T)
ER.. <- merge(ER.., EX, by="exon.id", all.x=T)
#Check HSP
ER..[ER..$exon.id=="NM_005590.3.5","freq.T.N.High.I"]

ms <- read.delim("MS.Candidats.Slicing.simple.version.txt", sep="\t")
ms$exon.id <- paste(ms$refseq,ms$nb.exon, sep=".")
#MS.Candidats.Slicing.txt
ER.. <- merge(ER.., ms, by="exon.id",all.x=T, all.y=F)

colnames(ER..)[ which(colnames(ER..)=="type")] <- "orientation"
ER.. <-  ER..[-which(is.na(ER..$bases)),]
  
##### MSI / NORMAL
msiN <- ER..[which(ER..$padjust.t.test.MSI.Normal<=0.05 & ER..$freq.T.N.High.I!=0),]
colnames(msiN)[which(colnames(msiN)=="type.y")] <- "exon.type"
length(which(msiN$Average.MSI.ratio.ES.IS<msiN$Average.Normal.msi.ratio.ES.IS))
#277
length(which(msiN$Average.MSI.ratio.ES.IS>msiN$Average.Normal.msi.ratio.ES.IS))
#1700
msiN$skipping <- "pb"
msiN[which(msiN$Average.MSI.ratio.ES.IS<msiN$Average.Normal.msi.ratio.ES.IS & msiN$orientation=="-"),"skipping"] <- "DOWN.MSI.acc"
msiN[which(msiN$Average.MSI.ratio.ES.IS<msiN$Average.Normal.msi.ratio.ES.IS & msiN$orientation=="+"),"skipping"] <- "DOWN.MSI.don"
msiN[which(msiN$Average.MSI.ratio.ES.IS>msiN$Average.Normal.msi.ratio.ES.IS & msiN$orientation=="-"),"skipping"] <- "UP.MSI.acc"
msiN[which(msiN$Average.MSI.ratio.ES.IS>msiN$Average.Normal.msi.ratio.ES.IS & msiN$orientation=="+"),"skipping"] <- "UP.MSI.don"

tp <- intersect(as.vector(msiN[which(msiN$skipping=="DOWN.MSI.acc"),"exon.id"]),as.vector(msiN[which(msiN$skipping=="DOWN.MSI.don"),"exon.id"]))

for( i in tp) {
  msiN[which(msiN$exon.id==i),"skipping"] <- "DOWN.MSI.don.acc"
}

tp <- intersect(as.vector(msiN[which(msiN$skipping=="UP.MSI.acc"),"exon.id"]),as.vector(msiN[which(msiN$skipping=="UP.MSI.don"),"exon.id"]))

for( i in tp){
  msiN[which(msiN$exon.id==i),"skipping"] <- "UP.MSI.don.acc"
}
    
msiN <-    msiN[-which(duplicated(msiN$exon.id)),]
tab <- table(msiN$exon.type,msiN$skipping)
tab <- tab[,c(4,1,5,2,6,3)]

par(oma=c(3,1,1,0))
svg("4.Exon.Skipping.MSI.Normal.svg")
par(mfrow = c(1,1))
barplot(tab, las=2, ylab="number of events", ylim=c(0,600))
dev.off()

##### MSS / NORMAL
mssN <- ER..[which(ER..$padjust.t.test.MSS.Normal<=0.05 & ER..$freq.T.N.High.S!=0 ),]
colnames(mssN)[which(colnames(mssN)=="type.y")] <- "exon.type"
length(which(mssN$Average.MSS.ratio.ES.IS<mssN$Average.Normal.mss.ratio.ES.IS))
#39
length(which(mssN$Average.MSS.ratio.ES.IS>mssN$Average.Normal.mss.ratio.ES.IS))
#222
mssN$skipping <- "UP.MSS"
mssN[which(mssN$Average.MSS.ratio.ES.IS<mssN$Average.Normal.mss.ratio.ES.IS & mssN$orientation=="-"),"skipping"] <- "DOWN.MSS.acc"
mssN[which(mssN$Average.MSS.ratio.ES.IS<mssN$Average.Normal.mss.ratio.ES.IS & mssN$orientation=="+"),"skipping"] <- "DOWN.MSS.don"
mssN[which(mssN$Average.MSS.ratio.ES.IS>mssN$Average.Normal.mss.ratio.ES.IS & mssN$orientation=="-"),"skipping"] <- "UP.MSS.acc"
mssN[which(mssN$Average.MSS.ratio.ES.IS>mssN$Average.Normal.mss.ratio.ES.IS & mssN$orientation=="+"),"skipping"] <- "UP.MSS.don"

tp <- intersect(as.vector(mssN[which(mssN$skipping=="DOWN.MSS.acc"),"exon.id"]),as.vector(mssN[which(mssN$skipping=="DOWN.MSS.don"),"exon.id"]))

for( i in tp){
  mssN[which(mssN$exon.id==i),"skipping"] <- "DOWN.MSS.don.acc"
}

tp <- intersect(as.vector(mssN[which(mssN$skipping=="UP.MSS.acc"),"exon.id"]),as.vector(mssN[which(mssN$skipping=="UP.MSS.don"),"exon.id"]))

for( i in tp){
  mssN[which(mssN$exon.id==i),"skipping"] <- "UP.MSS.don.acc"
}

mssN <-    mssN[-which(duplicated(mssN$exon.id)),]  
tab <- table(mssN$exon.type,mssN$skipping)
tab <- tab[,c(4,1,5,2,6,3)]

par(oma=c(3,1,1,0))
svg("4.Exon.Skipping.MSS.Normal.svg")
barplot(tab, las=2, ylab="number of events", ylim=c(0,600))
dev.off()

##### MSI.MSS / NORMAL
msi.sN <- ER..[which(ER..$padjust.t.test.MSS.Normal<=0.05 & ER..$padjust.t.test.MSI.Normal<=0.05  & ER..$freq.T.N.High.I!=0),]
colnames(msi.sN)[which(colnames(msi.sN)=="type.y")] <- "exon.type"

length(which(msi.sN$Average.MSS.ratio.ES.IS<msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS<msi.sN$Average.Normal.msi.ratio.ES.IS))
#63
length(which(msi.sN$Average.MSS.ratio.ES.IS>msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS>msi.sN$Average.Normal.msi.ratio.ES.IS))
#152
msi.sN$skipping <- "opposit"
msi.sN[which(msi.sN$Average.MSS.ratio.ES.IS<msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS<msi.sN$Average.Normal.msi.ratio.ES.IS & msi.sN$orientation=="-"),"skipping"] <- "DOWN.MSI.S.acc"
msi.sN[which(msi.sN$Average.MSS.ratio.ES.IS<msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS<msi.sN$Average.Normal.msi.ratio.ES.IS & msi.sN$orientation=="+"),"skipping"] <- "DOWN.MSI.S.don"
msi.sN[which(msi.sN$Average.MSS.ratio.ES.IS>msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS>msi.sN$Average.Normal.msi.ratio.ES.IS & msi.sN$orientation=="-"),"skipping"] <- "UP.MSI.S.acc"
msi.sN[which(msi.sN$Average.MSS.ratio.ES.IS>msi.sN$Average.Normal.mss.ratio.ES.IS & msi.sN$Average.MSI.ratio.ES.IS>msi.sN$Average.Normal.msi.ratio.ES.IS & msi.sN$orientation=="+"),"skipping"] <- "UP.MSI.S.don"

tp <- intersect(as.vector(msi.sN[which(msi.sN$skipping=="DOWN.MSI.S.acc"),"exon.id"]),as.vector(msi.sN[which(msi.sN$skipping=="DOWN.MSI.S.don"),"exon.id"]))

for( i in tp){
  msi.sN[which(msi.sN$exon.id==i),"skipping"] <- "DOWN.MSS.don.acc"
}

tp <- intersect(as.vector(msi.sN[which(msi.sN$skipping=="UP.MSI.S.acc"),"exon.id"]),as.vector(msi.sN[which(msi.sN$skipping=="UP.MSI.S.don"),"exon.id"]))

for( i in tp){
  msi.sN[which(msi.sN$exon.id==i),"skipping"] <- "UP.MSS.don.acc"
}

msi.sN <-    msi.sN[-which(duplicated(msi.sN$exon.id)),]
tab <- table(msi.sN$exon.type,msi.sN$skipping)
tab <- tab[,c(5,1,6,2,7,3)]
par(oma=c(3,1,1,0))
svg("4.Exon.Skipping.MSI.MSS.Normal.svg")
barplot(tab, las=2, ylab="number of events", ylim=c(0,600))
dev.off()
  
##### MSI / MSS
msimss <- ER..[which(ER..$padjust.t.test.MSI.MSS<=0.05 & ER..$freq.T.N.High.I!=0),]
colnames(msimss)[which(colnames(msimss)=="type.y")] <- "exon.type"

length(which(msimss$Average.MSI.ratio.ES.IS<msimss$Average.Normal.msi.ratio.ES.IS))
#104
length(which(msimss$Average.MSI.ratio.ES.IS>msimss$Average.Normal.msi.ratio.ES.IS))
#519
msimss$skipping <- "pb"
msimss[which(msimss$Average.MSI.ratio.ES.IS<msimss$Average.MSS.ratio.ES.IS & msimss$orientation=="-"),"skipping"] <- "DOWN.MSI.acc"
msimss[which(msimss$Average.MSI.ratio.ES.IS<msimss$Average.MSS.ratio.ES.IS & msimss$orientation=="+"),"skipping"] <- "DOWN.MSI.don"
msimss[which(msimss$Average.MSI.ratio.ES.IS>msimss$Average.MSS.ratio.ES.IS & msimss$orientation=="-"),"skipping"] <- "UP.MSI.acc"
msimss[which(msimss$Average.MSI.ratio.ES.IS>msimss$Average.MSS.ratio.ES.IS & msimss$orientation=="+"),"skipping"] <- "UP.MSI.don"

tp <- intersect(as.vector(msimss[which(msimss$skipping=="DOWN.MSI.acc"),"exon.id"]),as.vector(msimss[which(msimss$skipping=="DOWN.MSI.don"),"exon.id"]))

for( i in tp){
  msimss[which(msimss$exon.id==i),"skipping"] <- "DOWN.MSI.don.acc"
}

tp <- intersect(as.vector(msimss[which(msimss$skipping=="UP.MSI.acc"),"exon.id"]),as.vector(msimss[which(msimss$skipping=="UP.MSI.don"),"exon.id"]))

for( i in tp){
  msimss[which(msimss$exon.id==i),"skipping"] <- "UP.MSI.don.acc"
}

msimss <-    msimss[-which(duplicated(msimss$exon.id)),]
tab <- table(msimss$exon.type,msimss$skipping)
tab <- tab[,c(4,1,5,2,6,3)]

par(oma=c(3,1,1,0))
svg("4.Exon.Skipping.MSI.MSS.svg")
barplot(tab, las=2, ylab="number of events", ylim=c(0,600))
dev.off()

#Pie Chart
table(msiN$exon.type,msiN$skipping)[,c(4,6)]
table(msimss$exon.type,msimss$skipping)[,c(4,6)]
table(msiN$exon.type,msiN$skipping)[,c(1,3,4,6)]
table(msimss$exon.type,msimss$skipping)[,c(1,3,4,6)]

#remove duplicates in msiN and msimss 
rownames(msiN) <- msiN$exon.id
rownames(msimss) <- msimss$exon.id
msimss <-  msimss[setdiff (msimss$exon.id, msiN[which(msiN$skipping=="UP.MSI.acc"|msiN$skipping=="UP.MSI.don.acc"),"exon.id"]),]
ab.alt <- apply(cbind(table(msiN$exon.type,msiN$skipping)[,c(4,6)],table(msimss$exon.type,msimss$skipping)[,c(4,6)]),1, sum)

svg("4.PieChart.svg")
par(mfrow=c(1,2))
pie(c(41219,5293), labels=c("Constitutif", "Alternatif"), col=c("#ffaa00ff","#ffeeaaff"), main="Exons Type in Healthy Colorectal Tissue", cex.main=0.8)
pie(ab.alt, labels=c("Constitutif", "Alternatif"), col=c("#ffaa00ff","#ffeeaaff"), main="", cex.main=0.8)

dev.off()

chisq.test(matrix(c(41219,5293,293,692),2,2, byrow=TRUE), correct=TRUE) 
#p-value < 2.2e-16 

write.table(msiN,"1.msiN.output.txt", sep="\t", row.names=F)
write.table(mssN,"1.mssN.output.txt", sep="\t", row.names=F)
write.table(msimss,"1.msimss.output.txt",  sep="\t", row.names=F)

#########     msiN UP Acc    #################
svg("Freq.figure.MSI.NORMAL.UP.MSI.acc.svg")
par(mfrow=c(2,2))
hist(main="UP.MSI.ACC.constit",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc" & msiN$exon.type=="aberrant"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,60))
plot(main="UP.MSI.ACC.constit",density(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc" & msiN$exon.type=="aberrant"),"freq.T.N.High.I"], bw=7), add=T, las=1, ylim=c(0,0.035))
hist(main="UP.MSI.ACC.Alter",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc" & msiN$exon.type=="altern"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,60))
plot(main="UP.MSI.ACC.Alter",density(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc" & msiN$exon.type=="altern"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
dev.off()

svg("Freq.figure.MSI.NORMAL.UP.MSI.don.acc.svg")
par(mfrow=c(2,2))
hist(main="UP.MSI.ACC.don.constit",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don.acc" & msiN$exon.type=="aberrant"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50))
plot(main="UP.MSI.don.acc.constit",density(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don.acc" & msiN$exon.type=="aberrant"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
hist(main="UP.MSI.don.acc.Alter",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don.acc" & msiN$exon.type=="altern"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50))
plot(main="UP.MSI.don.acc.Alter",density(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don.acc" & msiN$exon.type=="altern"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
dev.off()

svg("Freq.figure.MSI.MSS.UP.MSI.acc.svg")
par(mfrow=c(2,2))
hist(main="UP.MSI.ACC.constit",msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.acc" & msimss$exon.type=="aberrant"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,17))
plot(main="UP.MSI.ACC.constit",density(msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.acc" & msimss$exon.type=="aberrant"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
hist(main="UP.MSI.ACC.Alter",msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.acc" & msimss$exon.type=="altern"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,17))
plot(main="UP.MSI.ACC.aberrant",density(msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.acc" & msimss$exon.type=="altern"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
dev.off()

svg("Freq.figure.MSI.MSS.UP.MSI.don.acc.svg")
par(mfrow=c(2,2))
hist(main="UP.MSI.don.acc.constit",msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.don.acc" & msimss$exon.type=="aberrant"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,17))
plot(main="UP.MSI.don.acc.constit",density(msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.don.acc" & msimss$exon.type=="aberrant"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
hist(main="UP.MSI.don.acc.Alter",msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.don.acc" & msimss$exon.type=="altern"),"freq.T.N.High.I"], breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,17))
  plot(main="UP.MSI.don.acc.aberrant",density(msimss[which(msimss$freq.T.N.High.I>1 & msimss$skipping=="UP.MSI.don.acc" & msimss$exon.type=="altern"),"freq.T.N.High.I"], bw=7), add=T, las=2, ylim=c(0,0.035))
dev.off()

############ FIGURE ALL
ab <- c(msiN[which(msiN$freq.T.N.High.I>1 & ( msiN$skipping=="UP.MSI.don.acc" | msiN$skipping=="UP.MSI.acc") & msiN$exon.type=="aberrant"),"freq.T.N.High.I"],
msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="UP.MSI.don.acc"| msimss$skipping=="UP.MSI.acc") & msimss$exon.type=="aberrant"),"freq.T.N.High.I"])
alt <- c(msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="UP.MSI.don.acc" | msimss$skipping=="UP.MSI.acc") & msimss$exon.type=="altern"),"freq.T.N.High.I"],
msiN[which(msiN$freq.T.N.High.I>1 & (msiN$skipping=="UP.MSI.don.acc"| msiN$skipping=="UP.MSI.acc") & msiN$exon.type=="altern"),"freq.T.N.High.I"])

ab.d <- c(msiN[which(msiN$freq.T.N.High.I>1 & (msiN$skipping=="DOWN.MSI.don.acc" | msiN$skipping=="DOWN.MSI.acc") & msiN$exon.type=="aberrant"),"freq.T.N.High.I"],
msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="DOWN.MSI.don.acc"| msimss$skipping=="DOWN.MSI.acc") & msimss$exon.type=="aberrant"),"freq.T.N.High.I"])
alt.d <- c(msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="DOWN.MSI.don.acc" | msimss$skipping=="DOWN.MSI.acc") & msimss$exon.type=="altern"),"freq.T.N.High.I"],
msiN[which(msiN$freq.T.N.High.I>1 & ( msiN$skipping=="DOWN.MSI.don.acc"| msiN$skipping=="DOWN.MSI.acc") & msiN$exon.type=="altern"),"freq.T.N.High.I"])

svg("Freq.figure.MSI.full.B.svg")
par(mfrow=c(2,2))
hist(main="UP.MSI.don.acc.constit",ab, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
hist(main="UP.MSI.don.acc.Alter",alt, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
hist(main="UP.MSI.don.acc.constit",ab.d, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
hist(main="UP.MSI.don.acc.Alter",alt.d, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
dev.off()

ab <- c(msiN[which(msiN$freq.T.N.High.I>1 & ( msiN$skipping=="UP.MSI.don.acc" | msiN$skipping=="UP.MSI.acc" | msiN$skipping=="DOWN.MSI.don.acc" | msiN$skipping=="DOWN.MSI.acc") & msiN$exon.type=="aberrant"),"freq.T.N.High.I"],
msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="UP.MSI.don.acc"| msimss$skipping=="UP.MSI.acc" | msimss$skipping=="DOWN.MSI.don.acc"| msimss$skipping=="DOWN.MSI.acc") & msimss$exon.type=="aberrant"),"freq.T.N.High.I"])
alt <- c(msimss[which(msimss$freq.T.N.High.I>1 & (msimss$skipping=="UP.MSI.don.acc" | msimss$skipping=="UP.MSI.acc" | msimss$skipping=="DOWN.MSI.don.acc" | msimss$skipping=="DOWN.MSI.acc") & msimss$exon.type=="altern"),"freq.T.N.High.I"],
msiN[which(msiN$freq.T.N.High.I>1 & (msiN$skipping=="UP.MSI.don.acc"| msiN$skipping=="UP.MSI.acc"| msiN$skipping=="DOWN.MSI.don.acc"| msiN$skipping=="DOWN.MSI.acc") & msiN$exon.type=="altern"),"freq.T.N.High.I"])
 
svg("Freq.figure.MSI.full.A.svg")
par(mfrow=c(1,2))
hist(main="UP.MSI.don.acc.constit",ab, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
hist(main="UP.MSI.don.acc.Alter",alt, breaks=50, density=F, las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events", ylim=c(0,50),xlim=c(0,100))
dev.off()

t.test(ab, alt)
# p-value < 2.2e-16

par(mfrow=c(2,2))
boxplot(main="UP.MSI.ACC",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc"),"freq.T.N.High.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events")
boxplot(main="Down.MSI.ACC",msiN[which(msiN$freq.T.N.low.I>1 & msiN$skipping=="DOWN.MSI.acc"),"freq.T.N.low.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events")
boxplot(main="UP.MSI.Don",msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don"),"freq.T.N.High.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events")
boxplot(main="Down.MSI.don",msiN[which(msiN$freq.T.N.low.I>1 & msiN$skipping=="DOWN.MSI.don"),"freq.T.N.low.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events")
    
par(mfrow=c(1,1))
plot(2,2,xlim=c(0,5), ylim=c(0,100) )
vioplot(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.acc"),"freq.T.N.High.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events",ylim=c(0,100),col="#ffcc00ff", at=1, add=T)
vioplot(msiN[which(msiN$freq.T.N.low.I>1 & msiN$skipping=="DOWN.MSI.acc"),"freq.T.N.low.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events",ylim=c(0,100),col="#ffcc00ff", at=2, add=T)
vioplot(msiN[which(msiN$freq.T.N.High.I>1 & msiN$skipping=="UP.MSI.don"),"freq.T.N.High.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events",ylim=c(0,100),col="#ffcc00ff", at=3, add=T)
vioplot(msiN[which(msiN$freq.T.N.low.I>1 & msiN$skipping=="DOWN.MSI.don"),"freq.T.N.low.I"], las=2, xlab="Frequency (Percentage of tumors)", ylab="number of skipping events",ylim=c(0,100),col="#ffcc00ff", at=4, add=T)
 
NIu <- unique(c(as.vector(msiN[which(msiN$skipping=="UP.MSI.acc"),"id"]), as.vector(msiN[which(msiN$skipping=="UP.MSI.don.acc"),"id"])))
length(NIu)
#839

NId <- unique(c(as.vector(msiN[which(msiN$skipping=="DOWN.MSI.acc"),"id"]), as.vector(msiN[which(msiN$skipping=="DOWN.MSI.don.acc"),"id"])))
length(NId)
#142

NI <- unique(c(NIu,NId))
length(NI)
#981

SSu <- unique(c(as.vector(msimss[which(msimss$skipping=="UP.MSI.acc"),"id"]), as.vector(msimss[which(msimss$skipping=="UP.MSI.don.acc"),"id"])))
length(SSu)
#276
SSd <- unique(c(as.vector(msimss[which(msimss$skipping=="DOWN.MSI.acc"),"id"]), as.vector(msimss[which(msimss$skipping=="DOWN.MSI.don.acc"),"id"])))
length(SSd)
#20
SS <- unique(c(SSu,SSd))
length(SS)
#296   

length(intersect(NIu, NId))
#0
length(intersect(NIu, SSu))
#130
length(intersect(NIu, SSd))
#1    
length(intersect(NId, SSu))
#2
length(intersect(NId, SSd)) 
#5
length(intersect(SSu, SSd))
#0

length(setdiff(NIu, NId))
#839
length(setdiff(NIu, SSu))
#709
length(setdiff(NIu, SSd))
#838    
length(setdiff(NId, SSu))
#140
length(setdiff(NId, SSd)) 
#137
length(setdiff(SSu, SSd))
#276

length(setdiff( NId, NIu))
#839
length(setdiff( SSu, NIu))
#709
length(setdiff(SSd, NIu))
#838    
length(setdiff(SSu, NId))
#140
length(setdiff(SSd, NId)) 
#137
length(setdiff(SSd,SSu))
#276

NId <- as.data.frame(NId, as.is=T)
NId$NId<- as.character(NId$NId)
NId$to <- T

NIu<- as.data.frame(NIu, as.is=T)
NIu$NIu<- as.character(NIu$NIu)
NIu$ti <- T

SSd <- as.data.frame(SSd, as.is=T)
SSd$SSd<- as.character(SSd$SSd)
SSd$tu <- T

SSu<- as.data.frame(SSu, as.is=T)
SSu$SSu<- as.character(SSu$SSu)
SSu$ta <- T

tab <- merge(NId, SSd, by.x="NId", by.y="SSd", all=T)
tab <- merge(tab, NIu, by.x="NId", by.y="NIu", all=T)
tab <- merge(tab, SSu, by.x="NId", by.y="SSu", all=T)

colnames(tab) <- c("id","MSI/N.down", "MSS/N.down", "MSI/N.Up", "MSS/N.Up")

for (i in 2:5){
  tab[which(is.na(tab[,i])),i] <- F
}

tab <- tab[,-1]

logi(tab[,2])

set.seed(5)
svg("Venn.dia2.svg")
plot(euler(tab[, 1:4], shape = "ellipse"), quantities = TRUE)
dev.off()

plot(venn(tab[1:100, 1:2]))

length(intersect(SS, NI))
#138
length(setdiff(SS, NI))
#158
length(setdiff(NI, SS))
#843

tp <- read.delim("Venn.dia1.txt")
set.seed(1)
svg("Venn.dia1.svg")
plot(euler(tp[, 1:2], shape = "ellipse"), quantities = TRUE)
dev.off()

