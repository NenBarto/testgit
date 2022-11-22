library(vcfR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(gtools)
library(Biostrings)
library(dplyr)
library(plyr)
library(dartR)

homedir="/Users/Nenad"
projectdir=paste0(homedir,"/projects/dingo")
#projectdir=paste0(homedir,"/projects/Anthony")
scriptsdir=paste0(projectdir,"/scripts")
tabledir=paste0(projectdir,"/tables")

system(paste0("mkdir -p ",tabledir))
figdir=paste0(projectdir,"/figures")
system(paste0("mkdir -p ",figdir))
popfile<-paste0(projectdir,"/annotation/popmap2")
populationsdir<-paste0(projectdir,"/results/populations/")

annotationdir<-"../annotation"

#rm(list=ls())
#vcfFile<-paste0("../vcfs/dingo_200_dog_wolf_default_722.3.vcf.all.gz")
vcfFile<-paste0("../vcfs/dingo_200_dog_wolf_default_unique_722.3.vcf.gz")
vcf<-read.vcfR(vcfFile)
head(vcf)

queryMETA(vcf)
gt <- extract.gt(vcf, IDtoRowNames = F)
fixed <- getFIX(vcf)
snps <- cbind(fixed[,1:5], gt)

head(snps)[,1:10]

snps.1 <- as.data.frame(as.matrix(snps))
snps.1$CHROM <- as.character(as.factor(snps.1$CHROM))
snps.1$POS <- as.character(as.factor(snps.1$POS))
snps.1$ID <- as.character(as.factor(snps.1$ID))
snps.1$identifier <- with(snps.1, paste0(CHROM, POS, ID))

rDataFiles<-list.files("../Robjects/",full.names=T)
rDataFiles<-rDataFiles[!grepl("canids8c",rDataFiles)]
for(rDataFile in rDataFiles[1]){
  fileName<-gsub("\\..*","",basename(rDataFile))
  gl<-gl.load(rDataFile)
  
  gl2plink(gl,plink_path="/opt/miniconda3/envs/stacks/bin/plink",outfile = "gl_vcf")
  
  #execute vcftools 
  outVCF<-paste0("../admixVCFs/",fileName,".vcf",)
  vcftoolsLine<-paste("vcftools --gzvcf",vcfFile,"--snps",SNPfile,"--keep",sampleFile,"--recode --recode-INFO-all --out",outVCF)
  system(vcftoolsLine)
}





