setwd("/Users/siakhnin/Documents/workspace/pagdur/")

geno = readRDS("chapter4/genotypes.rds")

freqA = numeric(10000)
freqB = numeric(10000)

for (i in 1:10000){
  hold = 0
  for (j in 1:2000)
    hold = hold + geno[i, j]
  freqB[i] = hold/4000
  freqA[i] = 1-freqB[i]
}


freqA = numeric(10000)
freqB = numeric(10000)
for (i in 1:10000){
  freqB[i] = sum(geno[i,]) / 4000
  freqA[i] = 1 - freqB[i]
}

freqB = rowSums(geno) / 4000
freqA = 1 - freqB

pheno =rnorm(2000,mean=100,sd=10)
pvals = numeric(10000)
for (i in 1:10000)
  pvals[i] = coef(summary(lm(pheno~geno[i,])))[2, 4]


pvals = apply(geno,1, function(snp) coef(summary(lm(pheno~snp)))[2,4])

M = geno - 1
M2 = M*geno

gwas = readRDS("chapter4/gwasData.rds")
pheno = read.table("chapter4/phenotypes.txt",header=T,sep="\t")$Pheno
map = read.table("chapter4/map.txt",header=T,sep="\t")
dim(gwas)

effect = numeric(10000) # effect sizes (coefficients)
pval = numeric(10000) # p-values
for (i in 1:10000){
  res = coef(summary(lm(pheno~gwas[i,])))[2,c(1,4)]
  effect[i] = res[1]
  pval[i] = res[2]
}
effect[1:4]
pval[1:4]

which(pval == min(pval))
which(pval == max(pval))

par(mfrow = c(2,1))

plot(gwas[8577,],pheno,
xlab = "genotypes", ylab="phenotypes",
main = paste("effect size:",round(effect[8577],2)))
mod = lm(pheno~gwas[8577,])
abline(mod,lwd=2,col="blue") # add regression line

# least significant SNP
plot(gwas[7296,],pheno,
xlab = "genotypes",ylab="phenotypes",
main = paste("effect size:",round(effect[7296],2)))
mod = lm(pheno~gwas[7296,])
abline(mod,lwd=2,col="blue")

y = pheno
x = gwas[1,]
mod = lm(y~x)
summary(mod)

b_hat = sum((x - mean(x)) * (y-mean(y))) / sum((x - mean(x))^2)
a_hat = mean(y) - b_hat*mean(x)
y_hat = a_hat + b_hat * x
y_hat[1:4]

TSS = sum((y - mean(y))^2)
SSR = sum((y_hat - mean(y))^2)
SSE = sum((y_hat - y)^2)

R_sq = SSR/TSS
Fval = SSR / (SSE / (length(y) - 2))
Pval = 1 - pf(Fval, 1, (length(y) - 2))
c(R_sq, Fval, Pval)

X = matrix(0, 2000, 2)
X[,1] = 1
interceptM = numeric(10000)
effectM = numeric(10000)

for (i in 1:10000){
  X[,2] = gwas[i,]
  XtX = t(X) %*% X
  lhs = solve(XtX)
  rhs = t(X) %*% y
  sol = lhs %*% rhs
  interceptM[i] = sol[1, 1]
  effectM[i] = sol[2, 1]
}

# snpBLUP
h2=0.5 # heritability
y=pheno
p=rowMeans(gwas)/2 # frequency of second allele
# lambda (equation 4.7)
d=2*sum(p*(1-p))
ve=var(y)*(1-h2) # residual variance
va=var(y)*h2 # additive variance
lambda = d * (ve/va)
 # equivalent for lambda:
   #lambda=(1-h2)/(h2/d)
   # snpBLUP (equation 4.6)
X=t(gwas-(p*2)) # freq. adjusted X matrix
XtX=t(X) %*% X # X'X
diag(XtX)=diag(XtX)+lambda # X'X+lambda
ones=rep(1,length(y)) # 1n
oto=t(ones)%*%ones # 1n'1n
otX=t(ones)%*%X # 1n'X
Xto=t(X)%*%ones # X'1n
 # build full matrix:
lhs=rbind(cbind(oto,otX),cbind(Xto,XtX))
oty=t(ones)%*%y # 1n'y
Xty=t(X)%*%y # X'y
rhs=rbind(oty,Xty)
effectBLUP=solve(lhs)%*%rhs # SNP solutions
