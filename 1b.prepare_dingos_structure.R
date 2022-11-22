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

plinkdir=paste0(projectdir,"/plink")
system(paste0("mkdir -p ",plinkdir))

figdir=paste0(projectdir,"/figures")
system(paste0("mkdir -p ",figdir))
popfile<-paste0(projectdir,"/annotation/popmap2")
populationsdir<-paste0(projectdir,"/results/populations/")

annotationdir<-"../annotation"

#rm(list=ls())
#vcfFile<-paste0("../vcfs/dingo_200_dog_wolf_default_722.3.vcf.all.gz")
#vcfFile<-paste0("../vcfs/dingo_200_dog_wolf_default_unique_722.3.vcf.gz")
#vcf<-read.vcfR(vcfFile)
#head(vcf)

#queryMETA(vcf)
#gt <- extract.gt(vcf, IDtoRowNames = F)
#fixed <- getFIX(vcf)
#snps <- cbind(fixed[,1:5], gt)
#
#head(snps)[,1:10]
#
#snps.1 <- as.data.frame(as.matrix(snps))
#snps.1$CHROM <- as.character(as.factor(snps.1$CHROM))
#snps.1$POS <- as.character(as.factor(snps.1$POS))
#snps.1$ID <- as.character(as.factor(snps.1$ID))
#snps.1$identifier <- with(snps.1, paste0(CHROM, POS, ID))

rDataFiles<-list.files("../Robjects/",full.names=T,pattern="dingos7k_nm")
rDataFiles<-rDataFiles[!grepl("canids8c",rDataFiles)]

e <- new.env(parent=as.environment("package:dartR"))
gl2plink2 <- function(x,
                     plink_path = getwd(),
                     bed_file = FALSE,
                     outfile = "gl_plink",
                     outpath = tempdir(),
                     chr_format = "character",
                     pos_cM = "0",
                     ID_dad = "0",
                     ID_mom = "0",
                     sex_code = "unknown",
                     phen_value = "0",
                     verbose = NULL) {
  local=e
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  outfilespec <- file.path(outpath, outfile)
  
  snp_temp <- as.data.frame(cbind(as.character(x$chromosome),x$position))
  colnames(snp_temp) <- c("chrom","snp_pos")
  
  
  if (chr_format == "numeric") {
    snp_temp$chrom <- as.numeric(snp_temp$chrom)
  }
  if (chr_format == "character") {
    snp_temp$chrom <- as.character(snp_temp$chrom)
  }
  
  snp_temp$snp_pos <- as.numeric(snp_temp$snp_pos)
  
  # Convert any NA values to 0
  snp_temp[is.na(snp_temp$snp_pos), "snp_pos"] <- 0
  # Convert any NA values to 0
  snp_temp[snp_temp$snp_pos == 0, "chrom"] <- 0
  
  # Chromosome code
  snp_chr <- snp_temp$chrom
  # Variant identifier
  var_id <- locNames(x)
  
  # Base-pair coordinate
  pos_bp <- snp_temp$snp_pos
  
  gl_map <- cbind(snp_chr, var_id, pos_cM, pos_bp)
  
  write.table(
    gl_map,
    file = paste0(outfilespec, ".map"),
    quote = F,
    row.names = F,
    col.names = F
  )
  
  ########## .fam (PLINK sample information file)
  
  sample.id_temp <- indNames(x)
  sample.id_temp <-
    gsub(" ", replacement = "_", sample.id_temp)
  
  # Family ID ('FID')
  FID <- as.character(x$pop)
  FID <- gsub(" ","_",FID)
  # Within-family ID ('IID'; cannot be '0')
  IID <- sample.id_temp
  
  ID_dad <- gsub(" ","_",ID_dad)
  ID_mom <- gsub(" ","_",ID_mom)
  
  # Sex code ('1' = male, '2' = female, '0' = unknown)
  if (length(sex_code) > 1) {
    sex_code <- as.character(sex_code)
    sex_code[startsWith(sex_code, "f") |
               startsWith(sex_code, "F")] <- "2"
    sex_code[startsWith(sex_code, "m") |
               startsWith(sex_code, "M")] <- "1"
    sex_code[startsWith(sex_code, "u") |
               startsWith(sex_code, "U")] <- "0"
    sex_code[nchar(sex_code) == 0] <- "0"
  }
  
  gl_fam <-
    cbind(FID, IID, ID_dad, ID_mom, sex_code, phen_value)
  
  x_mat <- as.matrix(x[, ])
  homs1 <- paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
  hets <- x@loc.all
  homs2 <- paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
  xx <- matrix(NA, ncol = ncol(x_mat), nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)) {
    for (ii in 1:ncol(x_mat)) {
      inp <- x_mat[i, ii]
      if (!is.na(inp)) {
        if (inp == 0)
          xx[i, ii] <- homs1[ii]
        else if (inp == 1)
          xx[i, ii] <- hets[ii]
        else if (inp == 2)
          xx[i, ii] <- homs2[ii]
      } else{
        xx[i, ii] = "0/0"
      }
    }
  }
  xx <- gsub("/", " ", xx)
  xx <- apply(xx, 1, paste, collapse = " ")
  xx <- cbind(gl_fam, xx)
  
  write.table(
    xx,
    file = paste0(outfilespec, ".ped"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if (bed_file) {
    prefix.in_temp <- outfilespec
    prefix.out_temp <- outfilespec
    
    allele_tmp <- gsub("/"," ", x$loc.all)
    allele_tmp <- strsplit(allele_tmp,split = " ")
    allele_tmp <- Reduce(rbind,allele_tmp)[,1]
    allele_tmp <- cbind(locNames(x), allele_tmp)
    write.table(allele_tmp,
                file = file.path(outpath,"mylist.txt"),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
    )
    
    make_plink <-
      function(plink.path,
               prefix.in = prefix.in_temp,
               prefix.out = prefix.out_temp,
               extra.options = "--dog") {
        bedfile.out <- paste0(prefix.out, ".bed")
        system_verbose(
          paste(
            plink.path,
            "--file",
            prefix.in,
            "--allow-no-sex",
            "--allow-extra-chr",
            # paste("--reference-allele",file.path(tempdir(),'mylist.txt')),
            "--out",
            prefix.out,
            extra.options
          )
        )
        bedfile.out
      }
    
    
    system_verbose = function(...) {
      report = system(..., intern = T)
      message(
        paste0(
          "\n\n----------Output of function start:\n\n",
          paste(report, collapse = "\n"),
          "\n\n----------Output of function finished...\n\n"
        )
      )
    }
    
    make_plink(plink.path = paste0(plink_path, "/plink"))
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(NULL)
  requireNamespace("dartR", quietly = TRUE)
}

environment(gl2plink2) <- asNamespace('dartR')

gl2vcf2<-function (x, plink_path = getwd(), outfile = "gl_vcf", outpath = tempdir(), 
          snp_pos = "0", snp_chr = "0", chr_format = "character", pos_cM = "0", 
          ID_dad = "0", ID_mom = "0", sex_code = "unknown", phen_value = "0", 
          verbose = NULL) 
{
  verbose <- gl.check.verbosity(verbose)
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", verbosity = verbose)
  datatype <- utils.check.datatype(x, verbose = verbose)
  gl2plink(x = x, outfile = "gl_plink_temp", outpath = outpath, 
           chr_format = chr_format, pos_cM = pos_cM, ID_dad = ID_dad, 
           ID_mom = ID_mom, sex_code = sex_code, phen_value = phen_value, 
           verbose = NULL)
  prefix.in_temp <- paste0(outpath, "/gl_plink_temp")
  prefix.out_temp <- file.path(outpath, outfile)
  allele_tmp <- gsub("/", " ", x$loc.all)
  allele_tmp <- strsplit(allele_tmp, split = " ")
  allele_tmp <- Reduce(rbind, allele_tmp)[, 2]
  allele_tmp <- cbind(locNames(x), allele_tmp)
  write.table(allele_tmp, file = file.path(tempdir(), "mylist.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  make_plink <- function(plink.path, prefix.in = prefix.in_temp, 
                         prefix.out = prefix.out_temp, autosome.only = FALSE, 
                         extra.options = "") {
    bedfile.out <- paste0(prefix.out, ".bed")
    system_verbose(paste(plink.path, "--file", prefix.in, 
                         "--recode", "vcf", if (autosome.only) 
                           "--autosome"
                         else "", "--dog --allow-no-sex", paste("--reference-allele", 
                                                          file.path(tempdir(), "mylist.txt")), "--out", 
                         prefix.out, extra.options))
  }
  system_verbose = function(...) {
    report = system(..., intern = T)
    message(paste0("\n\n----------Output of function start:\n\n", 
                   paste(report, collapse = "\n"), "\n\n----------Output of function finished...\n\n"))
  }
  make_plink(plink.path = paste0(plink_path, "/plink"), extra.options = "--aec")
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  invisible(NULL)
  requireNamespace("dartR", quietly = TRUE)
}

environment(gl2vcf2) <- asNamespace('dartR')


for(rDataFile in rDataFiles){
  fileName<-gsub("\\..*","",basename(rDataFile))
  gl<-gl.load(rDataFile)
  samples<-gl$ind.names
  VCFs<-gl$loc.names
  
  test <- gl.filter.callrate(gl)
  posColumn<-colnames(test$other$loc.metrics)[grepl("ChromPos",colnames(test$other$loc.metrics))][1]
  chrColumn<-colnames(test$other$loc.metrics)[grepl("Chrom_",colnames(test$other$loc.metrics))][2]
  
  test$position <- test$other$loc.metrics[,posColumn]
  test$chromosome <- test$other$loc.metrics[,chrColumn]
  test$chromosome<-factor(gsub(".*chromosome_","",test$chromosome))
  test$chromosome<-factor(gsub("_.*","",test$chromosome))
  
  test<-test[,test$chromosome %in% 1:38]
  
  anno<-data.frame(id=test$ind.names,pop=pop(test))
  anno$index<-as.numeric(factor(anno$pop))
  
  noSexData<-data.frame(id1=anno$id,id2=anno$id)
  
  
  #gl2faststructure(gl,outpath = ".",outfile=paste0(fileName,".str"))
  gl2plink2(test,bed_file = T,chr_format = "integer",plink_path = "/opt/miniconda3/envs/stacks/bin/",outpath = plinkdir,outfile=paste0(fileName))
  gl2vcf2(test,plink_path = "/opt/miniconda3/envs/stacks/bin/",outpath = plinkdir,outfile=paste0(fileName),snp_pos = test$position,snp_chr = test$chromosome)
  
  structPopFile<-paste0("../faststructure/",fileName,".populations.str.txt")
  write.table(anno,structPopFile,quote=F,sep="\t",row.names=F,col.names = F)
  write.table(anno,paste0(plinkdir,"/",fileName,".annotation.txt"),quote=F,sep="\t",row.names=F,col.names = F)
  write.table(noSexData,paste0(plinkdir,"/",fileName,".nosex"),quote=F,sep="\t",row.names=F,col.names = F)
}




