sires = read.table("~/Documents/workspace/pagdur/chapter2/siredata.txt", header=T, sep="\t", skip=3)
prog = read.table("~/Documents/workspace/pagdur/chapter2/progdata.txt", header=T, sep="\t", skip=3)

prog = prog[-which(prog$weight=="-"),]
prog$weight = as.numeric(as.character(prog$weight))
prog = prog[-which(prog$weight > 400),]

index=grep("m", names(prog))
missing = numeric()
for (i in 1:length(index))
  missing=c(missing, which(prog[,index[i]] == "-"))
print(missing)

missingU = unique(missing)
prog = prog[-missingU,]

for(i in 1:length(index))
  prog[,index[i]] = factor(prog[,index[i]])

indexms = grep("m", names(sires))
indexms = matrix(indexms, length(indexms)/2,2,byrow=T)

indexm = grep("m", names(prog))
indexm = matrix(indexm, length(indexms)/2,2,byrow=T)

library(made4)

compatible=matrix(NA,length(sires$id), length(indexms[,1]))

for(j in 1 : dim(indexms)[1]){
  for(i in 1 : length(sires$id)){
    indexs=which(prog$sire==sires$id[i])
    sirealleles=sires[i,indexms[j,]]
    sirealleles=c(as.character(sirealleles[,1]),
    as.character(sirealleles[,2]))
    
    hold=prog[indexs,indexm[j,]]
    hold=factor(c(as.character(hold[,1]),
    as.character(hold[,2])))
    hold=sort(summary(hold),decreasing=T)
    topalleles=names(hold)[1:2]
    
    compatible[i,j]=length(comparelists(sirealleles,topalleles)$Set.Dif)
    
    if(i==1 & j==1){
      cat("allele counts in offspring\n")
      print(hold)
      cat("most common alleles in offspring\n")
      print(topalleles)
      cat("sire alleles\n")
      print(sirealleles)
    }
  }
}

compatible = matrix(NA, length(prog$id), length(indexms[,1]))
for (j in 1:length(indexms[,1])) {
  for (i in 1:length(sires$id)) {
    indexs = which(prog$sire==sires$id[i])
    
    sirealleles = sires[i,indexms[j,]]
    sirealleles = c(as.character(sirealleles[,1]), as.character(sirealleles[,2]))
    
    for (k in 1:length(indexs)){
      hold = prog[indexs[k],indexm[j,]]
      topalleles = c(as.character(hold[,1]), as.character(hold[,2]))
      compatible[indexs[k],j] = length(comparelists(sirealleles, topalleles)$intersect)
    }
  }
}
compatible = data.frame(compatible)
for (i in 1:length(compatible[1,]))
  compatible[,i] = factor(as.character(compatible[,i]))
cat("\nSummary of alleles in sires and offspring\n")
summary(compatible)

index = which(compatible[,1] == 0)
prog = prog[-index,]

alleles = summary(factor(c(as.character(prog$m11), as.character(prog$m12))))
alleles = alleles/sum(alleles)

hold = data.frame(m11 = as.character(prog$m11), m12 = as.character(prog$m12))
hold[,1] = as.character(hold[,1])
hold[,2] = as.character(hold[,2])

sorted=character()

for (i in 1 : length(hold[,1]))
  sorted = rbind(sorted,sort(as.character(hold[i,])))
genotypes = paste(as.character(sorted[,1]), as.character(sorted[,2]),sep="_")
genotypes = summary(factor(genotypes))
genotypes = genotypes / sum(genotypes)

allgeno = NULL
for(i in 1 : length(indexm[,1])){
  hold = data.frame(prog[indexm[i,]])
  hold[,1] = as.character(hold[,1])
  hold[,2] = as.character(hold[,2])
  sorted = character()
  for (j in 1 : length(hold[,1])){
    sorted = rbind(sorted, sort(as.character(hold[j,])))
  }
  genotypes = paste(as.character(sorted[,1]), as.character(sorted[,2]), sep="_")
  allgeno = cbind(allgeno,genotypes)
}
colnames(allgeno) = c("M1", "M2", "M3", "M4", "M5")
markers = data.frame(prog[,1:4], allgeno)

model1 = lm(weight~M1, data=markers)
model2 = lm(weight~sex+M1, data=markers)
model3 = lm(weight~sex+M1+sire, data=markers)

library(nlme)
linear=coefficients(lm(weight~sex+sire+M5, data=markers))[c(1,2,12:28)]
random=fixef(lme(weight~sex+M5,random=~1|sire,data=markers))
linear=data.frame(effect=names(linear),fixed=linear)
random=data.frame(effect=names(random),random=random)
comparison=merge(linear,random,by="effect")
comparison=data.frame(comparison,difference=comparison$fixed-comparison$random)
print(comparison)