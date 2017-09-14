#!/bin/sh

#NOTE: run this script from within the directory containing the QC'ed files
#NOTE: this file must be run after QC.sh as it uses the trimmed_pruned files
echo "Have you run QC.sh? Respond Y or N, followed by [ENTER]: "
read QCdone
#NOTE: there must be an individual list with FID and IID columns (no column labels needed)
echo "Do you have an individual list with two columns (FID and IID) in this folder? Respond Y or N, followed by [ENTER]: "
read indlist

if [ $QCdone == "Y" ] && [ $indlist == "Y" ]
then
	echo "good to go"
else
	echo "You need to prepare the files for this script to run properly"
fi

module load plink
module load r

echo "Enter name of study population (e.g. WSC, MrOS, APOE), followed by [ENTER]: "
read study

echo "Enter number of genotype/chip/array pseudocohorts, followed by [ENTER]: "
read cohortnum

# gather genotype/chip/array pseudocohort(s) names into an array called "list"
for i in {1..$cohortnum}
do
	echo "Enter name(s) of genotype/chip/array pseudocohort(s), separated by spaces, followed by [ENTER]: "
	read cohortnames
	list=($cohortnames)
done

# determine phenotype name (as designated in the ${cohortname}_pheno.txt file)
#	can generate the pheno file, using the fam2pheno.R script: https://www.dropbox.com/s/g8e5zzwkvdc8ny0/fam2pheno.R?dl=0
echo "Enter BINARY phenotype name as it appears in the {cohortname}_pheno.txt file, followed by [ENTER]: "
read pheno

echo "Enter CONTINUOUS/QUANTITATIVE phenotype name as it appears in the {cohortname}_pheno.txt file, followed by [ENTER]: "
read quant

alpha=5e-8
echo "Enter significance threshold (alpha), followed by [ENTER]: "
read alpha

workDir=$PWD

# loop over cohort names in "list" to submit jobs
for k in "${list[@]}"
do
	# make directory to receive files for SNPtest output
	mkdir -p $workDir/${k}_gabel
	cp $workDir/${k}_subset/${k}.* ${k}_gabel
	cp $workDir/${k}_subset/${k}_pheno.txt ${k}_gabel
	
cat <<EOF >pheno01forGABEL.R
fam <- read.table("$workDir/${k}_gabel/${k}.fam",header=F)
fam <- transform(fam, V6 = ifelse(V6 == -9 , "NA", V6))
fam <- transform(fam, V6 = ifelse(V6 == 0 , "NA", V6))
fam <- transform(fam, V6 = ifelse(V6 == 1 , 0, V6))
fam <- transform(fam, V6 = ifelse(V6 == 2 , 1, V6))
nonames <- unname(fam)
write.table(fam,file="$workDir/${k}_gabel/${k}.fam",row.names=F,col.names=F,quote=F,sep=" ")
pheno <- read.table("$workDir/${k}_gabel/${k}_pheno.txt",header=T)
pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == -9 , "NA", ${pheno}))
pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 0 , "NA", ${pheno}))
pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 1 , 0, ${pheno}))
pheno <- transform(pheno, ${pheno} = ifelse(${pheno} == 2 , 1, ${pheno}))
write.table(pheno,file="$workDir/${k}_gabel/${k}_pheno.txt",row.names=F,quote=F,sep=" ")
q()
EOF
	
	R CMD BATCH pheno01forGABEL.R
	
	plink --bfile $workDir/${k}_gabel/${k} --pheno $workDir/${k}_gabel/${k}_pheno.txt --pheno-name ${pheno} --1 --transpose --recode --out $workDir/${k}_gabel/${k}_gabel
	plink --bfile $workDir/${k}_gabel/${k} --pheno $workDir/${k}_gabel/${k}_pheno.txt --pheno-name ${quant} --transpose --recode --out $workDir/${k}_gabel/${k}_gabel_qtl
	
	cat <<EOF >$workDir/${k}_gabel/${k}_gabel.R
#GenABEL package information: http://www.genabel.org/
install.packages("GenABEL", repos='http://cran.us.r-project.org')
library("GenABEL")
pheno<-read.table("$workDir/${k}_gabel/${k}_pheno.txt",header=T)
colnames(pheno)[2] <- "id"
pheno<-pheno[,c(2:dim(pheno)[2])]
binfam<-read.table("$workDir/${k}_gabel/${k}_gabel.tfam",header=F)
binsub <- binfam[,c(2,5)]
names(binsub)<-c("id","sex")
contfam<-read.table("$workDir/${k}_gabel/${k}_gabel_qtl.tfam",header=F)
contsub <- contfam[,c(2,5)]
names(contsub)<-c("id","sex")
#generate praw files
#-first change sex designations from M-1, F-2, unkwn-0 to M-1, F-0, unkwn-NA
binpraw <- transform(binsub, sex = ifelse(sex == 0, "NA", sex))
binpraw <- transform(binpraw, sex = ifelse(sex == 2, 0, sex))
contpraw <- transform(contsub, sex = ifelse(sex == 0, "NA", sex))
contpraw <- transform(contpraw, sex = ifelse(sex == 2, 0, sex))
#-next merge in pheno data
binpraw<-merge.data.frame(binpraw, pheno, by = "id", sort = F)
contpraw<-merge.data.frame(contpraw, pheno, by = "id", sort = F)
#--change missing data (-9) to "NA"
binpraw <- transform(binpraw, ${pheno} = ifelse(${pheno} == -9 , "NA", ${pheno}))
contpraw <- transform(contpraw, ${pheno} = ifelse(${pheno} == -9 , "NA", ${pheno}))
#--put age into 3rd column
bin_idx <- grep("age", names(binpraw), ignore.case=T)
cont_idx <- grep("age", names(contpraw), ignore.case=T)
if (length(bin_idx)==0) {
	print("No age column in ${pheno}_pheno.txt file")
	} else {
	binpraw <- binpraw[, c(bin_idx, (1:ncol(binpraw))[-bin_idx])]
	binpraw <- binpraw[,c(2,3,1,4:ncol(binpraw))]
	colnames(binpraw)[3]<-"age"
}
if (length(cont_idx)==0) {
	print("No age column in ${pheno}_pheno.txt file")
	} else {
	contpraw <- contpraw[, c(cont_idx, (1:ncol(contpraw))[-cont_idx])]
	contpraw <- contpraw[,c(2,3,1,4:ncol(contpraw))]
	colnames(contpraw)[3]<-"age"
}
write.table(binpraw,file="$workDir/${k}_gabel/${k}_gabel.praw",quote=F,row.names=F)
write.table(contpraw,file="$workDir/${k}_gabel/${k}_gabel_qtl.praw",quote=F,row.names=F)

#first analyze the CONTINUOUS/QUANTITATIVE phenotype
convert.snp.tped(tped = "$workDir/${k}_gabel/${k}_gabel_qtl.tped", tfam = "$workDir/${k}_gabel/${k}_gabel_qtl.tfam", out = "$workDir/${k}_gabel/${k}_gabel_qtl.raw", strand = "u")
g.dat <- load.gwaa.data(phen = "$workDir/${k}_gabel/${k}_gabel_qtl.praw", gen = "$workDir/${k}_gabel/${k}_gabel_qtl.raw", force = T)
slotNames(g.dat)
slotNames(g.dat@gtdata)
colnames(g.dat@phdata)
sample.size <- g.dat@gtdata@nids
snps.total <- g.dat@gtdata@nsnps
print(c(sample.size, snps.total))
summary(g.dat@phdata[,"${quant}"])
options(bitmapType='cairo')
png("$workDir/${k}_gabel/${k}_${quant}_summary.png",height=1000,width=1000)
hist(g.dat@phdata[,"${quant}"], main="Quantitative phenotype data summary", xlab = "${quant}", freq = F,breaks=20, col="gray")
rug(g.dat@phdata[,"${quant}"])
dev.off()
test.snp <- scan.glm('${quant} ~ CRSNP', family = gaussian(), data = g.dat)
names(test.snp)
alpha <- ${alpha}
attach(test.snp)
length(snpnames)==length(P1df)
snpnames[which(P1df<alpha)]
P1df[which(P1df<alpha)]
test.qt <- qtscore(${quant}, data = g.dat, trait = "gaussian")
test.qt@lambda
descriptives.scan(test.qt)
row.names(results(test.qt))[results(test.qt)[,"Pc1df"] < alpha]
results(test.qt)[which(results(test.qt)[,"P1df"] < alpha),"P1df"]
results(test.qt)[which(results(test.qt)[,"Pc1df"] < alpha),"Pc1df"]
obs <- sort(results(test.qt)[,"P1df"]) 
ept <- ppoints(obs)
options(bitmapType='cairo')
png("$workDir/${k}_gabel/${k}_${quant}_QQ.png",height=1000,width=1000)
plot(-log10(ept), -log10(obs), main = "${k} QQ plot, ${quant}", xlab="Expected -log10(pvalue)", ylab="Observed -log10(pvalue)")
abline(0, 1, col = "red")
abline(h = 8, lty = 2)
dev.off()
png("$workDir/${k}_gabel/${k}_${quant}_GWAS.png",height=1000,width=1000)
plot(test.qt, col = "black")
dev.off()

test.qt.sex <- qtscore(${quant} ~ sex, data = g.dat, trait = "gaussian")
row.names(results(test.qt.sex))[which(results(test.qt)[,"P1df"] < alpha)]
summary(lm(${quant} ~ sex, data = g.dat))

#now doing the process on the BINARY phenotype
convert.snp.tped(tped = "$workDir/${k}_gabel/${k}_gabel.tped", tfam = "$workDir/${k}_gabel/${k}_gabel.tfam", out = "$workDir/${k}_gabel/${k}_gabel.raw", strand = "u")
b.dat <- load.gwaa.data(phen = "$workDir/${k}_gabel/${k}_gabel.praw", gen = "$workDir/${k}_gabel/${k}_gabel.raw", force = T)
slotNames(b.dat)
slotNames(b.dat@gtdata)
colnames(b.dat@phdata)
b.dat@gtdata@nids
case.size <- length(which(b.dat@phdata[,"${pheno}"] == 1))
control.size <- length(which(b.dat@phdata[,"${pheno}"] == 0))
case.size 
control.size 
snpsb.total <- b.dat@gtdata@nsnps
testb.snp <- scan.glm('${pheno} ~ CRSNP', family = binomial(), data = b.dat)
names(testb.snp)
attach(testb.snp)
length(snpnames)==length(P1df)
snpnames[which(P1df<alpha)]
P1df[which(P1df<alpha)]
testb.qt <- qtscore(${pheno}, data = b.dat, trait = "binomial")
slotNames(testb.qt)
descriptives.scan(testb.qt)
row.names(results(testb.qt))[results(testb.qt)[,"P1df"] < alpha]
results(testb.qt)[results(testb.qt)[,"P1df"] < 1e-4,"P1df"]
results(testb.qt)[results(testb.qt)[,"Pc1df"] < 1e-4,"Pc1df"]
gkin <- ibs(g.dat, weight = "freq")
gkin[1:10,1:10]
cps.full <- cmdscale(as.dist(.5 - gkin), eig = T, k = 10)
names(cps.full)
attach(cps.full)
cps <- points
png("$workDir/${k}_gabel/${k}_${pheno}_PC1xPC2.png",height=1000,width=1000)
plot(cps[,1], cps[,2], pch = 20)
#Eliminated the next 2 lines because they allow plotting different points by designated population "popn"
###plot(cps[,1], cps[,2], pch = g.dat@phdata$popn)
###legend("topright", c("TSI","MEX", "CEU"), pch = c(1,2,3))
dev.off()
colnames(cps)<-c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10') 
gpc.dat <- g.dat
gpc.dat@phdata<-cbind(g.dat@phdata, cps)
test.pc.a <- scan.glm('${pheno} ~ CRSNP + C1 + C2 + C3 + C4 + C5', family=gaussian(), data = gpc.dat)
attach(test.pc.a)
length(snpnames)==length(P1df)
snpnames[which(P1df<alpha)]
P1df[which(P1df<alpha)]
test.pc.b <- qtscore(${pheno} ~  C1 + C2 + C3 + C4 + C5, data = gpc.dat, trait = "gaussian")
test.pc.b@lambda
png("$workDir/${k}_gabel/${k}_${pheno}_MDSscree.png",height=1000,width=1000)
plot(eig[1:10]/sum(eig), axes = F, type = "b", xlab = "Components",  ylim = c(0,0.05), ylab = "Proportion of Variations", main = "MDS analysis scree plot") 
axis(1, 1:10)
axis(2)
dev.off()

png("$workDir/${k}_gabel/${k}_${pheno}_MDScum.png",height=1000,width=1000)
plot(cumsum(eig[1:10])/sum(eig), axes = F, type = "b", ylim = c(0,0.2), xlab = "Components", ylab = "Proportion of Variations", main = "MDS analysis cumulative plot") 
axis(1, 1:10)
axis(2)
dev.off()

row.names(results(test.qt))[results(test.qt)[,"Pc1df"] < alpha]
results(test.qt)[which(results(test.qt)[,"Pc1df"] < alpha),"Pc1df"]
test.qt@lambda
obs <- sort(results(test.qt)[,"chi2.1df"]) 
ept <- sort(qchisq(ppoints(obs), df = 1))

png("$workDir/${k}_gabel/${k}_${pheno}_preGC.png",height=1000,width=1000)
plot(ept, obs, main = "Genomic control (lambda = slope of the dashed line)", xlab="Expected chisq, 1df", ylab="Observed chisq, 1df")
abline(0, 1, col = "red")
abline(0, test.qt@lambda[1], lty = 2)
dev.off()

median(results(test.qt)[,"chi2.1df"])/0.456
obs <- sort(results(test.qt)[,"Pc1df"])
ept <- ppoints(obs)

png("$workDir/${k}_gabel/${k}_${pheno}_postGC.png",height=1000,width=1000)
plot(-log10(ept), -log10(obs), main = "GWAS QQ plot adj. via Genomic Control", xlab="Expected -log10(pvalue)", ylab="Observed -log10(pvalue)")
abline(0, 1, col = "red")
abline(h = 8, lty = 2)
dev.off()

adj.gkin = gkin
diag(adj.gkin) = hom(g.dat)[,"Var"]
test.eg <- egscore(${pheno}, data = g.dat, kin = adj.gkin, naxes = 2)
descriptives.scan(test.eg)
snp.eg <- row.names(results(test.eg))[results(test.eg)[,"P1df"] < alpha]
pvalue.eg <- results(test.eg)[which(results(test.eg)[,"P1df"] < alpha),"P1df"]
lambda.eg <- test.eg@lambda
snp.eg 
pvalue.eg
lambda.eg
for (k in 1:10){
	test.tmp <- egscore(${pheno}, data = g.dat, kin = adj.gkin, naxes = k)
	print(test.tmp@lambda["estimate"])
}
obs <- sort(results(test.eg)[,"Pc1df"])
ept <- ppoints(obs)

png("$workDir/${k}_gabel/${k}_${pheno}_QQ_EIGENSTRATadj.png",height=1000,width=1000)
plot(-log10(ept), -log10(obs), main = "GWAS QQ plot adj. w/ EIGENSTRAT", xlab="Expected -log10(pvalue)", ylab="Observed -log10(pvalue)")
abline(0, 1, col = "red")
abline(h = 8, lty = 2)
dev.off()

png("$workDir/${k}_gabel/${k}_${pheno}_Manhattanadj.png",height=1000,width=1000)
plot(test.qt, col = "black")
add.plot(test.eg, col = "gray", pch = 3)
legend("topright", c("Original plot","After correction w/ EIGENSTRAT"), pch = c(1,3))
dev.off()

q()
EOF

	R CMD BATCH $workDir/${k}_gabel/${k}_gabel.R
done