# gSEM

We performed gSEM in iPSYCH on 5 PGC summary statistics. 
Following instructions of [their github]() with a specific nod to their main-README note to use an `export` command to avoid paralelization issues. 

## step 1: munge summary statistics

```r
library(dplyr)
require(GenomicSEM)

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
# downloading from https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v/file/805005013708
hm3 <- "/home/jasmijn/MDD_heterogeneity/people/jolien/results/2025_02_14_gSEM_5disPGC/data/w_hm3.snplist"

#name the traits 
trait.names <- c("ADHD","AUT","BPD","MDD","SCZ")

#list the sample sizes. All but PTSD have SNP-specific sum of effective sample sizes so only its
#sample size is listed here
N=c(NA,NA,NA,NA,NA)


#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01



#can rename the headers of my files to what they're expecting with this list in the column.names within munge()
names <- list(N="N_col")

setwd("xxx/munge")

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter,column.names=names)
```

## step2: LDSC to get covariance matrix

```r
library(dplyr)
require(GenomicSEM)

traits <- c(ADHD, AUT, BPD, MDD, SCZ)
trait.names <- c("ADHD","AUT","BPD","MDD","SCZ")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev <- c(.5,.5,.5,.5)

prevsource <- "published"
population.prev <- c(0.03235,0.01965,0.0158,0.12285,0.01745)

#logfile name
ldsc.log <- paste0(indir,"/../log/",prevsource,"_ldsc.log")

## run ldsc - defaults of the pacakge.

#the folder of LD scores (see prev. step for source)
refbase <- "xxx/eur_w_ld_chr/"

ld <- refbase

#the folder of LD weights [typically the same as folder of LD scores]
wld <- refbase

#run LDSC
setwd("xxx/ldsc_out/")

LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#optional command to save the output as a .RData file for later use
ldsc_out <- paste0(indir,"/../ldsc_out/",prevsource,"_LDSCoutput.RData")
save(LDSCoutput,file=ldsc_out) 



```

## step 3: prep sumstats formats..

```r
library(dplyr)
require(GenomicSEM)
library(dplyr)
require(GenomicSEM)


# indir <- "/home/jasmijn/MDD_heterogeneity/people/jolien/results/2025_02_14_gSEM_5disPGC/ldsc_out"

trait.names <- c("ADHD","AUT","BPD","MDD","SCZ")

sumsdir <- "/xxx/data"

ADHD <- paste0(sumsdir,"/noB-PGC.ADHD.GWAS.MTAGheader.txt")
AUT <- paste0(sumsdir,"/noB-PGC.AUT.GWAS.MTAGheader.txt")
BPD <- paste0(sumsdir,"/noB-PGC.BPD.GWAS.MTAGheader.txt")
MDD <- paste0(sumsdir,"/noB-PGC.MDD.GWAS.MTAGheader_NeffAdj.txt")
SCZ <- paste0(sumsdir,"/noB-PGC.SCZ.GWAS.MTAGheader.txt")


#create vector of the summary statistics files
files <- c(ADHD,AUT,BPD,MDD,SCZ)

N=c(NA,NA,NA,NA,NA)

#using 1000G reference LD
ref= "xxx/data/reference.1000G.maf.0.005.txt"


info.filter=.9
maf.filter=0.05 # diffrent than default 0.01!


se.logit=c(T,T,T,T,T)
linprob=c(F,F,F,F,F)

#step 3: prep
PSYCH_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,N=N,info.filter=info.filter,se.logit=se.logit,linprob=linprob,maf.filter=maf.filter) #,column.names=names

sumsout <- "/home/jasmijn/MDD_heterogeneity/people/jolien/results/2025_02_14_gSEM_5disPGC/outputsums/psych_sumstats.RData"
save(PSYCH_sumstats, file= sumsout)




```

## step 4: common-factor GWAS

!!! IMPORTANT !!!
to avoid issues, make sure to run this on clusters before doing interactive R (after srun) or to include it in submission scripts. otherwise it takes unneccessary much resources and time!

`export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`


```r
library(dplyr)
require(GenomicSEM)

#load the PSYCH_sumstats from prev. step
sumsout <- "xxx/psych_sumstats.RData"
load(sumsout)


	prevsource <- "published"

	indir <- "/home/jasmijn/MDD_heterogeneity/people/jolien/results/2025_02_14_gSEM_5disPGC/munge"
	ldsc_out <- paste0(indir,"/../ldsc_out/",prevsource,"_LDSCoutput.RData")
	load(ldsc_out)
	#step 4: getting the common factor GWAS out
	#run the multivariate GWAS using parallel processing
	
	setwd("/home/jasmijn/MDD_heterogeneity/people/jolien/results/2025_02_14_gSEM_5disPGC/outputsums/")
	#IMPORTANT: cores requires a number, default NULL doesn't work.prob because parallel default = TRUE

#run the multivariate GWAS using parallel processing

PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats, estimation="DWLS", cores = 99, toler=FALSE,SNPSE=FALSE,parallel=TRUE,GC="standard",MPI=FALSE,TWAS=FALSE,smooth_check=FALSE)

outfile <- paste0("xxx/outputsums/",prevsource,"_PSYCH_factor.RData")

save(PSYCH_factor, file=outfile)

# # output psych_factor to a file -- this is the meta-analyzed gwas output
outfile <- paste0("xxx/outputsums/",prevsource,"_PSYCH_factor.txt")
write.table(PSYCH_factor, outfile, quote=FALSE, col.names=TRUE, row.names=FALSE)


# ## if you want the sample size for this meta-analysis for followup analyses: 
# #restrict to MAF of 40% and 10%
PSYCH_factor2<-subset(PSYCH_factor, PSYCH_factor$MAF <= .4 & PSYCH_factor$MAF >= .1)
#calculate expected sample size (N_hat)
N_hat<-mean(1/((2*PSYCH_factor2$MAF*(1-PSYCH_factor2$MAF))*PSYCH_factor2$se_c^2))

#output that N_hat somewhere for later reference.
noutfile <- paste0("xxx/outputsums/",prevsource,"_PSYCH_factor-expectedN.txt")
write(N_hat, file=noutfile)
print(N_hat)


```


## step 5: obtain factor loadings

```r
library(dplyr)
require(GenomicSEM)

#load the PSYCH_sumstats thingy
sumsout <- "xxx/psych_sumstats.RData"
load(sumsout)



prevsource <- "published"
#prevsource <- "register"

indir <- "xxx/munge"
ldsc_out <- paste0(indir,"/../ldsc_out/",prevsource,"_LDSCoutput.RData")
load(ldsc_out)

# step 3: specify an estiamte a structural equation model
CommonFactor_DWLS <- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

CommonFactor_DWLS



```