library(admixtools)
library(tidyverse)

#Run data for dingos + dogs + wolves
prefix="../plink_dingos7k_with_wolves/dingos7k_with_wolves"
my_f2_dir<-"../plink_dingo7k_ww/"

#first empty the out folder
system(paste0("rm -rf ",my_f2_dir))

#define populations
mypops<-c("alpine","desert","dog","Mallee","Wolves")

#extract plink file
extract_f2(prefix,my_f2_dir,pops = mypops,format = "plink",maxmiss=1)
f2_blocks = f2_from_precomp(my_f2_dir)
count_snps(f2_blocks)
dim(f2_blocks)

pop1<-"Wolves"
pop2<-"dog"
pop3<-"Mallee"
pop4<-"desert"
pop5<-"alpine"


#f2 measures the amount of genetic drift that separates two populations
f2_blocks[,,1]                # f2-statistics of the 1st SNP block
apply(f2_blocks, 1:2, mean)   # average across all blocks
f2_blocks[pop1, pop2, ] 

#in our case f2 is largest between wolves and Mallee, followed by
#wolves and desert, wolves and alpine, wolves and dog
#this is to be expected

#f3 measures 
f3(f2_blocks) 


#f4 measures the amount of drift that is shared between two population pairs
#if a pair (dingo+wolves) forms a clade relative to other pair (desert+alpine) then f4 
#should be 0

f4(f2_blocks, pop1, pop2, pop4, pop5)
#dog and wolves are not a clade to desert+alpine since the 
#correlation in allele frequencies between pairs is significantly
#different from 0

f4(f2_blocks, pop1, pop2, pop3, pop4)
f4(f2_blocks, pop1, pop2, pop3, pop5)
#dog and wolves are a clade to desert+Mallee, but not a clade to alpine and mallee

#all combinations are here
f4(f2_blocks)

#f3 tests whether population is admixed
#if negative, first is admixed between a second and third population
#in our cases, dingoes should not be admixed with dogs

f3res<-f3(f2_blocks)
write.table(f3res,"f3res_wolves.txt",sep="\t")

fst(my_f2_dir)

left = c('alpine','desert','Mallee')
right = c('dog','Wolves')
target = 'alpine'
pops = c(left, right, target)
results = qpadm(f2_blocks, left, right, target)
write.table(results,"qpadm_wolves.txt",sep="\t")

#create graph
#ADeMDoW=origin of alpine,desert,mallee,dog,wolves
graph<-data.frame(
  to  =c("ADeMDoW", "ADeMDoW","ADeMDo","ADeMDo","ADeM",  "ADeM","ADe",   "ADe",   "dog"),
  from=c("Wolves",  "ADeMDo", "dog",   "ADeM",  "Mallee","ADe", "desert","alpine","alpine"))


my_graph<-tibble(graph)
qpg_results = qpgraph(f2_blocks, my_graph)
plot_graph(qpg_results$edges,fix = T)



#Run data for dingos + dogs

prefix="../plink_dingos7k_nm/dingos7k_nm"
my_f2_dir<-"../plink_dingo7k_nmOut/"

#first empty the out folder
system(paste0("rm -rf ",my_f2_dir))

#define populations
mypops<-c("alpine","desert","dog","Mallee")

#extract plink file
extract_f2(prefix,my_f2_dir,pops = mypops,format = "plink",maxmiss=1)
f2_blocks = f2_from_precomp(my_f2_dir)
count_snps(f2_blocks)
dim(f2_blocks)

pop1<-"alpine"
pop2<-"desert"
pop3<-"dog"
pop4<-"Mallee"


#f2 measures the amount of genetic drift that separates two populations
f2_blocks[,,1]                # f2-statistics of the 1st SNP block
apply(f2_blocks, 1:2, mean)   # average across all blocks
f2_blocks[pop1, pop2, ] 

#in our case f2 is largest between wolves and Mallee, followed by
#wolves and desert, wolves and alpine, wolves and dog
#this is to be expected


#f4 measures the amount of drift that is shared between two population pairs
#if a pair (dingo+wolves) forms a clade relative to other pair (desert+alpine) then f4 
#should be 0

f4(f2_blocks, pop1, pop2, pop4, pop5)
#dog and wolves are not a clade to desert+alpine since the 
#correlation in allele frequencies between pairs is significantly
#different from 0

f4(f2_blocks, pop1, pop2, pop3, pop4)
f4(f2_blocks, pop1, pop2, pop3, pop5)
#dog and wolves are a clade to desert+Mallee, but not a clade to alpine and mallee

#all combinations are here
f4(f2_blocks)

#f3 tests whether population is admixed
#if negative, first is admixed between a second and third population
#in our cases, dingoes should not be admixed with dogs
f3res<-f3(f2_blocks) 
write.table(f3res,"f3res_onlyDingos.txt",sep="\t")


graph<-data.frame(
  to  =c("ADeMDo","ADeMDo","ADeM",  "ADeM","ADe",   "ADe",   "dog"),
  from=c("dog",   "ADeM",  "Mallee","ADe", "desert","alpine","alpine"))

my_graph<-tibble(graph)
qpg_results = qpgraph(f2_blocks, my_graph)
plot_graph(qpg_results$edges,fix = T)




