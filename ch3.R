library(RSQLite)

con=dbConnect(dbDriver("SQLite"), "SNPsmall")

animids = dbGetQuery(con, "select distinct animal from snps")
animids = as.vector(animids$animal)

hold = dbGetQuery(con, paste("select * from snps where animal = '", animids[1], "'", sep=""))

snpids = as.vector(dbGetQuery(con, "select distinct name from snpmap")[,1])

dbDisconnect(con)


first_snp = dbGetQuery(con, paste("select * from snps where snp='",snpids[1], "'",sep=""))
first_snp$allele1 = factor(first_snp$allele1)
first_snp$allele2 = factor(first_snp$allele2)

snp = data.frame(first_snp, genotype=factor(paste(first_snp$allele1, first_snp$allele2, sep=""),
                                            levels=c("AA","AB", "BB")))
plot(snp$x, snp$y, col=snp$genotype, pch=as.numeric(snp$genotype), 
     xlab="x", ylab="y",
     main=snp$snp[1], cex.main=0.9)
legend("bottomleft",paste(levels(snp$genotype), " (",summary(snp$genotype),")",sep=""),
  col= 1:length(levels(snp$genotype)),
  pch= 1:length(levels(snp$genotype)),
  cex=0.7)

alleles=factor(c(as.character(snp$allele1), as.character(snp$allele2)), levels=c("A","B"))
summary(alleles) / sum(summary(alleles)) * 100

obs = summary(factor(snp$genotype))
hwal = summary(factor(c(as.character(snp$allele1), as.character(snp$allele2))))
hwal = hwal / sum(hwal)

exp = c(hwal[1]^2, 2 * hwal[1] * hwal[2], hwal[2]^2) * sum(obs)
names(exp) = c("AA","AB","BB")
xtot = sum((abs(obs-exp) - c(0.5, 1, 0.5))^2 / exp)
pval = 1 - pchisq(xtot, 1)
print(pval)

sumslides = matrix(NA, 83, 4)
rownames(sumslides) = animids
colnames(sumslides) = c ("-/-", "A/A", "A/B", "B/B")
numgeno = matrix(9, 54977, 83)
for (i in 1:83){
  hold = dbGetQuery(con, paste("select * from snps where animal='",animids[i],"'",sep=""))
  hold = data.frame(hold, genotype = factor(paste(hold$allele1,hold$allele2,sep=""), levels = c("--", "AA", "AB", "BB")))
  hold = hold[order(hold$snp),]
  sumslides[i, ] = summary(hold$genotype)
  temp = hold$genotype
  levels(temp) = c(9,0,1,2)
  numgeno[,i] = as.numeric(as.character(temp))
  # change to 9 genotypes under GC score cutoff
  numgeno[which(hold$gcscore<0.6),i] = 9
}

rownames(numgeno) = hold$snp
colnames(numgeno) = animids

samplehetero=sumslides[,3] / (sumslides[,2] + sumslides[,3] + sumslides[,4])
up = mean(samplehetero) + 3 * sd(samplehetero)
down = mean(samplehetero) - 3 * sd(samplehetero)
hsout = length(which(samplehetero > up))
hsout = hsout + length(which(samplehetero < down))

plot(sort(samplehetero),1:83,col="blue",cex.main=0.9,
  cex.axis=0.8,cex.lab=0.8,
  ylab="sample",xlab="heterozygosity",
  main=paste("Sample heterozygosity\nmean:",
  round(mean(samplehetero),3)," sd:",
  round(sd(samplehetero),3)),
  sub=paste("mean: black line ",3,
  "SD: red line number of outliers:",hsout),
  cex.sub=0.8)
abline(v=mean(samplehetero))
abline(v=mean(samplehetero)-3*sd(samplehetero),col="red")
abline(v=mean(samplehetero)+3*sd(samplehetero),col="red")

animcor = cor(numgeno)
library("gplots")
hmcol = greenred(256)
heatmap(animcor,col = hmcol,symm = T,labRow = " ",labCol = " ", trace="none")

genotypes=read.table("SNPxSample.txt",
  header=T,sep="\t",na.strings = "9",
  colClasses = "factor")
dim(genotypes)

for (i in 1:length(genotypes[1,]))
  levels(genotypes[,i])=c("AA","AB","BB",NA)

indexsnp = apply(genotypes, 1, function(x) length(which(is.na(x) == T)))
indexsnp = which(indexsnp == length(genotypes[1,]))
indexsample = apply(genotypes, 2, function(x) length(which(is.na(x) == T)))
indexsample = which(indexsample == length(genotypes[,1]))

genotypes = genotypes[-indexsnp,]

weight=rnorm(83,mean=50,sd=10)

plot(density(weight),col="blue",main="Density plot of weights")
abline(v=mean(weight),col="red")
lines(density(rnorm(83000,mean=50,sd=10)),col="green",lty=2)

singlesnp = function(trait, snp){
  if (length(levels(snp)) > 1) lm(trait~snp)
  else NA
}

results = apply(genotypes, 1, function(x) singlesnp(weight, factor(t(x))))
pvalfunc = function(model){
  if(class(model)=="lm") anova(model)[[5]][1]
  else NA
}
pvals = lapply(results, function(x) pvalfunc(x))
names(results) = row.names(genotypes)
pvals = data.frame(snp = row.names(genotypes), pvalue=unlist(pvals))

index=sort(pvals$pvalue,index.return=T)[[2]][1:5]

estimates=NULL
for (i in 1:5){
  estimates = rbind(estimates, coefficients(summary(results[[index[i]]])))
}
estimates = cbind(rep(c("AA mean","AB dev","BB dev"), 5), estimates, rep(names(results)[index], each=3))
estimates = data.frame(estimates)
names(estimates) = c("genotype","effect","stderror","t-value","p-value","snp")
for (i in 2:5){
  estimates[,i] = signif(as.numeric(as.character(estimates[,i])), 2)
}
print(estimates)

map=dbGetQuery(con, paste("select * from snpmap",sep=""))
merged=merge(pvals,map,by.x=1,by.y=1)

plot(merged$position[which(merged$chromosome==1)],
  -log(merged$pvalue[which(merged$chromosome==1)]),
  xlab="map position",ylab="-log odds",
  col="blue",pch=20,main="Chromosome 1")
abline(h=-log(0.01),col="red")

length(which(pvals$pvalue<0.01))
length(which(pvals$pvalue<0.01/length(pvals$pvalue)))
sort(pvals$pvalue)[1:5]
