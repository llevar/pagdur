library(RSQLite)

con=dbConnect(dbDriver("SQLite"), dbname = "SNPsmall")

animids = dbGetQuery(con, "select distinct animal from snps")
animids = as.vector(animids$animal)

hold = dbGetQuery(con, paste("select * from snps where animal = '", animids[1], "'", sep=""))

snpids = as.vector(dbGetQuery(con, "select distinct name from snpmap")[,1])

first_snp = dbGetQuery(con, paste("select * from snps where snp='",snpids[1], "'",sep=""))
first_snp$allele1 = factor(first_snp$allele1)
first_snp$allele2 = factor(first_snp$allele2)
dbDisconnect(con)

