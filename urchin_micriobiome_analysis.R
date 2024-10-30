# Author: Mia M. Bengtsson
# Contact: mia.bengtsson@uni-greifswald.de
# Date: October 30th 2024

#load data
seq<-read.csv("/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/urchin_seqtab.csv", header=T, row.names = 1)
meta<-read.csv("/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/urchin_map.csv", header=T, row.names = 1)
tax<-read.csv("/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/urchin_tax.csv", header=T, row.names = 1)

#sequencing depth
meta$libsize<-rowSums(seq)

hist(meta$libsize, breaks=20)
min(meta$libsize)

#cut out samples with low sequencing depth
meta1<-meta[meta$libsize>10000,] 
seq1<-seq[meta$libsize>10000,] 

#cut out seqs with zero abundance
seltax<-which(colSums(seq1)>0)
seq2<-seq1[,seltax]
tax1<-tax[seltax,]
#and cut out eukaryotes
seq3<-seq2[,tax1$level1!="Eukaryota"]
tax2<-tax1[tax1$level1!="Eukaryota",]

#cut out samples which were dominated by plastid sequences and now have low sequencing depth
meta1$libsize2<-rowSums(seq3)
meta2<-meta1[meta1$libsize2>10000,] 
seq4<-seq3[meta1$libsize2>10000,]

#observed sequence variants (richness)
meta2$nseq<-apply(seq4>0, 1, sum)

boxplot(meta2$dna~meta2$treatment)

library(vegan)

#calculate rarified richness and evenness
meta2$Rrichness<-rarefy(seq4, sample = min(meta2$libsize2))
meta2$evenness<-(diversity(seq4))/log(meta2$Rrichness)

#plot richness and evenness
boxplot(meta2$Rrichness~meta2$treatment , vertical=T , main="Richness", ylab="no. of seqs (rarefied)")
boxplot(meta2$evenness~meta2$treatment , vertical=T , main="Evenness", ylab="Pielous evenness")

#anova richness and evenness
hist(meta2$Rrichness) #check for normality
hist(log(meta2$Rrichness))
summary(aov(log(meta2$Rrichness)~meta2$treatment))
TukeyHSD(aov(log(meta2$Rrichness)~meta2$treatment))

hist(meta2$evenness)
#hist(log(meta2$evenness))
summary(aov(log(meta2$evenness)~meta2$treatment))
TukeyHSD(aov(log(meta2$evenness)~meta2$treatment))

#hellinger-transform dataset
seq4n<-decostand(seq4, method = "hellinger")

mds1<-metaMDS(seq4n)

plot(mds1, display="sites", type = "none")
points(mds1, col=meta2$alga, pch=21, cex=3)
text(mds1, cex=0.5)

#statistical test PERMANOVA
adonis2(seq4n~treatment, meta2) 
adonis2(seq4n[7:36,]~treatment*tank, meta2[7:36,]) 
adonis2(seq4n~alga, meta2) 

adist<-vegdist(seq4n, "bray")

#pairwise permanova post-hoc test
amod<-betadisper(adist, meta2$treatment)
permutest(amod)
plot(amod)
boxplot(amod)
mod.hsd<-TukeyHSD(amod)
mod.hsd
plot(mod.hsd)

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pmod<-pairwise.adonis2(seq4n~treatment, data = meta2)
pmod2<-pairwise.adonis2(seq4n[7:36,]~treatment*tank, meta2[7:36,]) 

p.adjust(p=c(pmod$before_vs_fucus[1,5], pmod$before_vs_control[1,5], pmod$before_vs_palmaria[1,5],pmod$before_vs_saccharina[1,5],pmod$fucus_vs_control[1,5],pmod$fucus_vs_palmaria[1,5], pmod$fucus_vs_saccharina[1,5], pmod$control_vs_palmaria[1,5], pmod$control_vs_saccharina[1,5], pmod$palmaria_vs_saccharina[1,5]),"bonferroni")

plot(anosim(seq4n, grouping = meta2$treatment))

#barplots

# aggregate seqs to class level (approximately)

agg<-aggregate(t(seq4), list(phylum=tax2$level3), sum)

ix<-as.character(agg[,1])
ix<-c("Unclassified",ix[2:length(ix)])
phylum<-data.frame(t(agg[,-1]))

names(phylum)<-ix
rownames(phylum)<-rownames(meta2)
#get rid of low abundance groups

Other<-phylum[,colSums(phylum)<2000]
phylum1<-phylum[,colSums(phylum)>2000]
phylum1$Other<-rowSums(Other)

phylum_norm<-data.frame(apply(phylum1,FUN="/",MARGIN=2, rowSums(phylum1)))
summary(rowSums(phylum_norm))

#sort according to sites
#phylum_norm$ids<-meta1$tnitrate
#phylum_norm<-phylum_norm[ order(phylum_norm$ids), ]
#phylum_norm<-phylum_norm[,1:ncol(phylum_norm)-1]

#plot
quartz(7,14)
par(mfrow=c(1, 2), mar=c(2,3,2,2))
rb<-c("#AF6BE7","#6CD6D1","#7EE043","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#C48783","#CDCCD8","#9D7CAF","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF","#AF6BE7","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF","#AF6BE7","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF")
#picked from http://tools.medialab.sciences-po.fr/iwanthue/
# terrible, but at least contrasting
rb1<-rb[1:23]
barplot(as.matrix(t(phylum_norm)), col=rb1, cex.names = 0.5)
#legend
x1<-1:10
x2<-11:19
plot(1:19, type="n", axes=F, xlab="", ylab="", xlim=c(1, 20))

points(rep(1, length(x1)),x1, pch=15, col=rb[x1], cex=2)
text(rep(6, length(x1)),x1, names(phylum_norm)[x1])
points(rep(1, length(x2)),x2, pch=15, col=rb[x2], cex=2)
text(rep(6, length(x2)),x2, names(phylum_norm)[x2])

# aggregate seqs to genus level (approximately)

agg<-aggregate(t(seq4), list(phylum=tax2$level6), sum)

ix<-as.character(agg[,1])
ix<-c("Unclassified",ix[2:length(ix)])
phylum<-data.frame(t(agg[,-1]))

names(phylum)<-ix
rownames(phylum)<-rownames(meta2)
#get rid of low abundance groups

Other<-phylum[,colSums(phylum)<2000]
phylum1<-phylum[,colSums(phylum)>2000]
phylum1$Other<-rowSums(Other)

phylum_norm<-data.frame(apply(phylum1,FUN="/",MARGIN=2, rowSums(phylum1)))
summary(rowSums(phylum_norm))

#sort according to sites
#phylum_norm$ids<-meta1$tnitrate
#phylum_norm<-phylum_norm[ order(phylum_norm$ids), ]
#phylum_norm<-phylum_norm[,1:ncol(phylum_norm)-1]

#plot
quartz(7,14)
par(mfrow=c(1, 2), mar=c(2,3,2,2))
rb<-c("#AF6BE7","#6CD6D1","#7EE043","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#C48783","#CDCCD8","#9D7CAF","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF","#AF6BE7","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF","#AF6BE7","#6CD6D1","#7EE043","#E15F71","#DCA9D6","#C6DB7F","#E85530","#67A0D3", "#659976","#68DE93","#5DA53E","#7D9098","#7684E6", "#C0DDB6","#E84A93","#D071BC","#C48783","#E04BD6","#CDCCD8","#9D7CAF")
#picked from http://tools.medialab.sciences-po.fr/iwanthue/
# terrible, but at least contrasting
rb1<-rb[1:23]
barplot(as.matrix(t(phylum_norm)), col=rb1, cex.names = 0.5)
#legend
x1<-1:10
x2<-11:19
plot(1:19, type="n", axes=F, xlab="", ylab="", xlim=c(1, 20))

points(rep(1, length(x1)),x1, pch=15, col=rb[x1], cex=2)
text(rep(6, length(x1)),x1, names(phylum_norm)[x1])
points(rep(1, length(x2)),x2, pch=15, col=rb[x2], cex=2)
text(rep(6, length(x2)),x2, names(phylum_norm)[x2])

#how much percentage is otu 1 out of gammaproteobacteria?

rowSums(agg[18,2:37])
sum(seq1$X1)

sum(seq1$X1)/rowSums(agg[18,2:37])


#differential abundance analysis
library(phyloseq)
library(DESeq2)



meta1$before<-0
meta1$before[meta1$time=="before"]<-1
meta1$after<-0
meta1$after[meta1$time=="after"]<-1
meta1$fucus<-0
meta1$fucus[meta1$treatment=="fucus"]<-1
meta1$palmaria<-0
meta1$palmaria[meta1$treatment=="palmaria"]<-1
meta1$control<-0
meta1$control[meta1$treatment=="control"]<-1
meta1$saccharina<-0
meta1$saccharina[meta1$treatment=="saccharina"]<-1

#converting treatments to factors due to issues downstream with DESeq2
meta1$before<-as.factor(meta1$before)
meta1$after<-as.factor(meta1$after)
meta1$fucus<-as.factor(meta1$fucus)
meta1$palmaria<-as.factor(meta1$palmaria)
meta1$control<-as.factor(meta1$control)
meta1$saccharina<-as.factor(meta1$saccharina)

meta1$treatment<-as.factor(meta1$treatment)

#making phyloseq objects

cbind(rownames(tax),colnames(seq1)) #seem to match
rownames(tax)<-colnames(seq1)



phy_urchin<-phyloseq(otu_table(as.matrix(seq1), taxa_are_rows = F), tax_table(as.matrix(tax)), sample_data(meta1))
#using the whole 16S count data including plastids to meet requirements of the DESeq model
phy_urchin

#deseq2 test for all treatments together

test_urchin<- phyloseq_to_deseq2(phy_urchin, ~ treatment)

test_urchin = DESeq(test_urchin, test="Wald", fitType="parametric")
alpha<-0.01

before_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "before","control"))
sigtab_before = before_res[which(before_res$padj < alpha), ]
sigtab_before = cbind(as(sigtab_before, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_before), ], "matrix"))

fucus_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "fucus","control"))
sigtab_fucus = fucus_res[which(fucus_res$padj < alpha), ]
sigtab_fucus = cbind(as(sigtab_fucus, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_fucus), ], "matrix"))

saccharina_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "saccharina","control"))
sigtab_saccharina = saccharina_res[which(saccharina_res$padj < alpha), ]
sigtab_saccharina = cbind(as(sigtab_saccharina, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_saccharina), ], "matrix"))

palmaria_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "palmaria","control"))
sigtab_palmaria = palmaria_res[which(palmaria_res$padj < alpha), ]
sigtab_palmaria = cbind(as(sigtab_palmaria, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_palmaria), ], "matrix"))

#write.csv(sigtab_before, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_before.csv")
#write.csv(sigtab_fucus, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_fucus.csv")
#write.csv(sigtab_palmaria, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_palmaria.csv")
#write.csv(sigtab_saccharina, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_saccharina.csv")

fucus_saccharina_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "fucus","saccharina"))
sigtab_fucus_saccharina = fucus_saccharina_res[which(fucus_saccharina_res$padj < alpha), ]
sigtab_fucus_saccharina = cbind(as(sigtab_fucus_saccharina, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_fucus_saccharina), ], "matrix"))

fucus_palmaria_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "fucus","palmaria"))
sigtab_fucus_palmaria = fucus_palmaria_res[which(fucus_palmaria_res$padj < alpha), ]
sigtab_fucus_palmaria = cbind(as(sigtab_fucus_palmaria, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_fucus_palmaria), ], "matrix"))

saccharina_palmaria_res = results(test_urchin, cooksCutoff = FALSE, alpha=alpha, contrast = c("treatment", "saccharina","palmaria"))
sigtab_saccharina_palmaria = saccharina_palmaria_res[which(saccharina_palmaria_res$padj < alpha), ]
sigtab_saccharina_palmaria = cbind(as(sigtab_saccharina_palmaria, "data.frame"), as(tax_table(phy_urchin)[rownames(sigtab_saccharina_palmaria), ], "matrix"))

#write.csv(sigtab_fucus_saccharina, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_fucus_saccharina.csv")
#write.csv(sigtab_fucus_palmaria, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_fucus_palmaria.csv")
#write.csv(sigtab_saccharina_palmaria, "/Users/Mia/Dropbox/jobb/Greifswald/Sea_urchins/R-files/sigtab_saccharina_palmaria.csv")


#make a table with all the differentially abundant ASVs
sigASV_all<-unique(c(rownames(sigtab_fucus), rownames(sigtab_palmaria), rownames(sigtab_saccharina), rownames(sigtab_before)))
sigtax<-tax[sigASV_all,]
#remove plastids
sigtax2<-sigtax[sigtax$level2!="Plastid",]
#calculate relative abundance in each treatment
sigtax2$per_total<-0
sigtax2$per_control<-0
sigtax2$per_field<-0
sigtax2$per_fucus<-0
sigtax2$per_saccharina<-0
sigtax2$per_palmaria<-0

#custom plotting order
sigtab_before$phylum<-13
sigtab_before$phylum[sigtab_before$level3=="Gammaproteobacteria"]<-2
sigtab_before$phylum[sigtab_before$level3=="Alphaproteobacteria"]<-1
sigtab_before$phylum[sigtab_before$level3=="Epsilonproteobacteria"]<-4
sigtab_before$phylum[sigtab_before$level3=="Deltaproteobacteria"]<-3
sigtab_before$phylum[sigtab_before$level2=="Bacteroidetes"]<-5
sigtab_before$phylum[sigtab_before$level2=="Verrucomicrobia"]<-6
sigtab_before$phylum[sigtab_before$level2=="Planctomycetes"]<-7
sigtab_before$phylum[sigtab_before$level2=="Firmicutes"]<-8
sigtab_before$phylum[sigtab_before$level2=="Actinobacteria"]<-9
sigtab_before$phylum[sigtab_before$level2=="Lentisphaerae"]<-10
sigtab_before$phylum[sigtab_before$level2=="Spirochaetes"]<-11
sigtab_before$phylum[sigtab_before$level2=="Plastid"]<-12

sigtab_fucus$phylum<-13
sigtab_fucus$phylum[sigtab_fucus$level3=="Gammaproteobacteria"]<-2
sigtab_fucus$phylum[sigtab_fucus$level3=="Alphaproteobacteria"]<-1
sigtab_fucus$phylum[sigtab_fucus$level3=="Epsilonproteobacteria"]<-4
sigtab_fucus$phylum[sigtab_fucus$level3=="Deltaproteobacteria"]<-3
sigtab_fucus$phylum[sigtab_fucus$level2=="Bacteroidetes"]<-5
sigtab_fucus$phylum[sigtab_fucus$level2=="Verrucomicrobia"]<-6
sigtab_fucus$phylum[sigtab_fucus$level2=="Planctomycetes"]<-7
sigtab_fucus$phylum[sigtab_fucus$level2=="Firmicutes"]<-8
sigtab_fucus$phylum[sigtab_fucus$level2=="Actinobacteria"]<-9
sigtab_fucus$phylum[sigtab_fucus$level2=="Lentisphaerae"]<-10
sigtab_fucus$phylum[sigtab_fucus$level2=="Spirochaetes"]<-11
sigtab_fucus$phylum[sigtab_fucus$level2=="Plastid"]<-12

sigtab_saccharina$phylum<-13
sigtab_saccharina$phylum[sigtab_saccharina$level3=="Gammaproteobacteria"]<-2
sigtab_saccharina$phylum[sigtab_saccharina$level3=="Alphaproteobacteria"]<-1
sigtab_saccharina$phylum[sigtab_saccharina$level3=="Epsilonproteobacteria"]<-4
sigtab_saccharina$phylum[sigtab_saccharina$level3=="Deltaproteobacteria"]<-3
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Bacteroidetes"]<-5
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Verrucomicrobia"]<-6
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Planctomycetes"]<-7
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Firmicutes"]<-8
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Actinobacteria"]<-9
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Lentisphaerae"]<-10
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Spirochaetes"]<-11
sigtab_saccharina$phylum[sigtab_saccharina$level2=="Plastid"]<-12

sigtab_palmaria$phylum<-13
sigtab_palmaria$phylum[sigtab_palmaria$level3=="Gammaproteobacteria"]<-2
sigtab_palmaria$phylum[sigtab_palmaria$level3=="Alphaproteobacteria"]<-1
sigtab_palmaria$phylum[sigtab_palmaria$level3=="Epsilonproteobacteria"]<-4
sigtab_palmaria$phylum[sigtab_palmaria$level3=="Deltaproteobacteria"]<-3
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Bacteroidetes"]<-5
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Verrucomicrobia"]<-6
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Planctomycetes"]<-7
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Firmicutes"]<-8
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Actinobacteria"]<-9
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Lentisphaerae"]<-10
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Spirochaetes"]<-11
sigtab_palmaria$phylum[sigtab_palmaria$level2=="Plastid"]<-12

library(ggplot2)

quartz()
par(mfrow=c(5,1), mar=c(1,5,0.5,2), cex=1.1)

radius<-sqrt(sigtab_before$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="", xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_before$phylum[sigtab_before$log2FoldChange<0]), sigtab_before$log2FoldChange[sigtab_before$log2FoldChange<0], circles=radius[sigtab_before$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("khaki1", 0.6), add=T)
symbols( (sigtab_before$phylum[sigtab_before$log2FoldChange>0]), sigtab_before$log2FoldChange[sigtab_before$log2FoldChange>0], circles=radius[sigtab_before$log2FoldChange>0], inches=F , xlim = c(0,13  ), bg=alpha("aquamarine3", 0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49, adj=0,"Field")
text(-0.3, -35, adj=0,"Control")

radius<-sqrt(sigtab_fucus$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="", xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_fucus$phylum[sigtab_fucus$log2FoldChange<0]), sigtab_fucus$log2FoldChange[sigtab_fucus$log2FoldChange<0], circles=radius[sigtab_fucus$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("khaki1", 0.6), add=T)
symbols( (sigtab_fucus$phylum[sigtab_fucus$log2FoldChange>0]), sigtab_fucus$log2FoldChange[sigtab_fucus$log2FoldChange>0], circles=radius[sigtab_fucus$log2FoldChange>0], inches=F , xlim = c(0,13), bg=alpha("thistle", 0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49, adj=0,"Fucus")
text(-0.3, -35, adj=0,"Control")

radius<-sqrt(sigtab_saccharina$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="",  xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_saccharina$phylum[sigtab_saccharina$log2FoldChange<0]), sigtab_saccharina$log2FoldChange[sigtab_saccharina$log2FoldChange<0], circles=radius[sigtab_saccharina$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13 ), bg=alpha("khaki1", 0.6), add=T)
symbols( (sigtab_saccharina$phylum[sigtab_saccharina$log2FoldChange>0]), sigtab_saccharina$log2FoldChange[sigtab_saccharina$log2FoldChange>0], circles=radius[sigtab_saccharina$log2FoldChange>0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("steelblue",0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49,adj=0, "Saccharina")
text(-0.3, -35,adj=0, "Control")

radius<-sqrt(sigtab_palmaria$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="",  xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols((sigtab_palmaria$phylum[sigtab_palmaria$log2FoldChange<0]), sigtab_palmaria$log2FoldChange[sigtab_palmaria$log2FoldChange<0], circles=radius[sigtab_palmaria$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13 ), bg=alpha("khaki1", 0.6), add=T)
symbols((sigtab_palmaria$phylum[sigtab_palmaria$log2FoldChange>0]), sigtab_palmaria$log2FoldChange[sigtab_palmaria$log2FoldChange>0], circles=radius[sigtab_palmaria$log2FoldChange>0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("tomato", 0.6), ylim = c(-40,40), xaxt="n", add=T)
#text( jitter(sigtab_palmaria$phylum[sigtab_palmaria$log2FoldChange<0]), sigtab_palmaria$log2FoldChange[sigtab_palmaria$log2FoldChange<0], labels=rownames(sigtab_palmaria), cex=0.5)
#text( jitter(sigtab_palmaria$phylum[sigtab_palmaria$log2FoldChange>0]), sigtab_palmaria$log2FoldChange[sigtab_palmaria$log2FoldChange>0], labels=rownames(sigtab_palmaria), cex=0.5)


text(-0.3, 49, adj=0,"Palmaria")
text(-0.3, -35,adj=0, "Control")

plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "", xlab="", axes=F, xlim=c(0,13), ylim=c(-40,50))
axis(1, at=1:13, pos=50,cex=0.9, labels = c("Alphaproteobacteria","Gammaproteobacteria","Deltaproteobacteria","Epsilonproteobacteria","Bacteroidetes","Verrucomicrobia","Planctomycetes","Firmicutes","Actinobacteria","Lentisphaerae", "Spirochaetes", "Plastid", "Other"), las=2)


#plot for algal comparisons
sigtab_fucus_saccharina$phylum<-13
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level3=="Gammaproteobacteria"]<-2
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level3=="Alphaproteobacteria"]<-1
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level3=="Epsilonproteobacteria"]<-4
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level3=="Deltaproteobacteria"]<-3
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Bacteroidetes"]<-5
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Verrucomicrobia"]<-6
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Planctomycetes"]<-7
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Firmicutes"]<-8
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Actinobacteria"]<-9
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Lentisphaerae"]<-10
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Spirochaetes"]<-11
sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$level2=="Plastid"]<-12

sigtab_fucus_palmaria$phylum<-13
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level3=="Gammaproteobacteria"]<-2
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level3=="Alphaproteobacteria"]<-1
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level3=="Epsilonproteobacteria"]<-4
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level3=="Deltaproteobacteria"]<-3
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Bacteroidetes"]<-5
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Verrucomicrobia"]<-6
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Planctomycetes"]<-7
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Firmicutes"]<-8
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Actinobacteria"]<-9
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Lentisphaerae"]<-10
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Spirochaetes"]<-11
sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$level2=="Plastid"]<-12

sigtab_saccharina_palmaria$phylum<-13
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level3=="Gammaproteobacteria"]<-2
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level3=="Alphaproteobacteria"]<-1
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level3=="Epsilonproteobacteria"]<-4
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level3=="Deltaproteobacteria"]<-3
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Bacteroidetes"]<-5
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Verrucomicrobia"]<-6
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Planctomycetes"]<-7
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Firmicutes"]<-8
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Actinobacteria"]<-9
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Lentisphaerae"]<-10
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Spirochaetes"]<-11
sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$level2=="Plastid"]<-12

quartz()
par(mfrow=c(4,1), mar=c(1,5,0.5,2), cex=1.1)

radius<-sqrt(sigtab_fucus_saccharina$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="", xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$log2FoldChange<0]), sigtab_fucus_saccharina$log2FoldChange[sigtab_fucus_saccharina$log2FoldChange<0], circles=radius[sigtab_fucus_saccharina$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("steelblue", 0.6), add=T)
symbols( (sigtab_fucus_saccharina$phylum[sigtab_fucus_saccharina$log2FoldChange>0]), sigtab_fucus_saccharina$log2FoldChange[sigtab_fucus_saccharina$log2FoldChange>0], circles=radius[sigtab_fucus_saccharina$log2FoldChange>0], inches=F , xlim = c(0,13  ), bg=alpha("thistle", 0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49, adj=0,"Fucus")
text(-0.3, -35, adj=0,"Saccharina")

radius<-sqrt(sigtab_fucus_palmaria$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="", xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$log2FoldChange<0]), sigtab_fucus_palmaria$log2FoldChange[sigtab_fucus_palmaria$log2FoldChange<0], circles=radius[sigtab_fucus_palmaria$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("tomato", 0.6), add=T)
symbols( (sigtab_fucus_palmaria$phylum[sigtab_fucus_palmaria$log2FoldChange>0]), sigtab_fucus_palmaria$log2FoldChange[sigtab_fucus_palmaria$log2FoldChange>0], circles=radius[sigtab_fucus_palmaria$log2FoldChange>0], inches=F , xlim = c(0,13), bg=alpha("thistle", 0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49, adj=0,"Fucus")
text(-0.3, -35, adj=0,"Palmaria")

radius<-sqrt(sigtab_saccharina_palmaria$baseMean/pi)/30
plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "log2-Fold Change", xlab="",  xaxt="n", xlim=c(0,13), ylim=c(-40,55))
abline(a=0, b=0, lty=2)
symbols( (sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$log2FoldChange<0]), sigtab_saccharina_palmaria$log2FoldChange[sigtab_saccharina_palmaria$log2FoldChange<0], circles=radius[sigtab_saccharina_palmaria$log2FoldChange<0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13 ), bg=alpha("tomato", 0.6), add=T)
symbols( (sigtab_saccharina_palmaria$phylum[sigtab_saccharina_palmaria$log2FoldChange>0]), sigtab_saccharina_palmaria$log2FoldChange[sigtab_saccharina_palmaria$log2FoldChange>0], circles=radius[sigtab_saccharina_palmaria$log2FoldChange>0], inches=F , ylab = "log2-Fold Change", xlab="", xlim = c(0,13  ), bg=alpha("steelblue",0.6), ylim = c(-40,40),xaxt="n", add = T)
text(-0.3, 49,adj=0, "Saccharina")
text(-0.3, -35,adj=0, "Palmaria")

plot(x=1:13, y=seq(-40, 40, length.out=13), type="n",ylab = "", xlab="", axes=F, xlim=c(0,13), ylim=c(-40,50))
axis(1, at=1:13, pos=50,cex=0.9, labels = c("Alphaproteobacteria","Gammaproteobacteria","Deltaproteobacteria","Epsilonproteobacteria","Bacteroidetes","Verrucomicrobia","Planctomycetes","Firmicutes","Actinobacteria","Lentisphaerae", "Spirochaetes", "Plastid", "Other"), las=2)


#core ASV analysis

#core community identification (frequency-abundance)

freq16<-apply(seq6>0, 2, sum)
abund16<-colSums(seq6)
plot(abund16, freq16, log = "x")
tax2[which(freq16==nrow(seq6)),]
tax2[which(freq16==nrow(seq6)-1),]
tax2[which(freq16==nrow(seq6)-2),]
tax2[which(freq16==nrow(seq6)-3),]


#how much percentage are core ASVs out of all?

sum(seq6$X1)/sum(seq6)
sum(seq6$X9)/sum(seq6)
sum(seq6$X4)/sum(seq6)
sum(seq6$X7)/sum(seq6)
sum(seq6$X8)/sum(seq6)

# how many of the reads are eukaryotic?

sum(seq[tax$level1=="Eukaryota"])/sum(seq)
