# incremental MTAG (in UKB)


First, we need to rank the phenotypes by rG with the focal phenotype (LifetimeMDD)
To do this we use [LDSC](https://github.com/bulik/ldsc/). 

## step1: munge sumstats
```bash

#location of ldsc installation
ldscbase="xxx/ldsc"

#note to change --signed-sumstats based on the headers of sumstats / make sure to match sumstats with this. the 0 indiciates binary phenotype

outfile= #full basename of the output to be generated.

${ldscbase}/munge_sumstats.py \
--sumstats ${sumdir}/${pheno}.sumstats \
--out ${outfile} \
--snp snpid \
--a1 a1 \
--a2 a2 \
--p pval \
--signed-sumstats beta,0 \
--N-col n

```

## step 2: obtain rG's with focal phenotype

Loop over all phenotypes to be assessed, and save them pairwise going 
focalsumstats - sumstatsA
focalsumstats - sumstatsB
focalsumstats - sumstatsC

their file names are, per pair, comma-separated and fed to LDSC (`--rg` flag)

Reference files used as linked via the [LDSC github](https://github.com/bulik/ldsc/). 
```bash

ldscbase="/project/genomics/jolien.rietkerk/CoordinatedInteraction/bin/ldsc"

${ldscbase}/ldsc.py \
--rg ${sumstatsbase}/${pheno1}.sumstats.gz,${sumstatsbase}/${pheno2}.sumstats.gz \
--ref-ld-chr ${refbase}/20211027_refld_H2rG/ \
--w-ld-chr ${refbase}/20211027_refld_H2rG/ \
--out ${outdir}/${pheno1}_${pheno2}_LDSC_rG

```

## step 3: rank the summary statistics on rG with focal
Merging the output of LDSC into a file that can be read by R more easily: 

```bash
#loop over the pairs of sumstats and obtain only the first lines as useful for R.
tail -n 5 ${outfile}_LDSC_rG.log | head -n 2 > ${outdir}/${pheno1}_${pheno2}_formatted.txt

```

In R, obtain the output information on rG of the pairs of summary statistics and rank the phenotypes by abs(rG) with the focal phenotype

```r
disorderA <- "LifetimeMDD"
phenolist_full= #list of phenotypes that were paired

data <- NULL
for (disorderB in phenolist_full) {
	current_file <- paste0(wdir,"/input/",disorderA,"_",disorderB,"_formatted.txt")

	current_data <- read.table(current_file, header=TRUE)

	entry <- data.frame("p1"=disorderA, "p2"=disorderB,
	                    "rg"=current_data$rg,"se"=current_data$se,
	                    "z"=current_data$z,"p"=current_data$p,
	                    "h2_obs"=current_data$h2_obs, 
	                    "h2_obs_se"=current_data$h2_obs_se,
	                    "h2_int"=current_data$h2_int,
	                    "h2_int_se"=current_data$h2_int_se, 
	                    "gcov_int"=current_data$gcov_int,
	                    "gcov_int_se"=current_data$gcov_int_se)
	data <- rbind(data, entry)

	# pairs <- rbind(pairs, data.frame(pair = paste0(disorderA,"_",disorderB), stringsAsFactors = FALSE))

	}
}

## ranking abs(rG)

rankdata <- data.frame("p1"=character(),"p2"=character(),
                       "rg"=numeric(),"p"=numeric(),"rank"=numeric())
focus <- "LifetimeMDD"

current_data <- data
#rank using arrange() from dplyr
current_data <- arrange(current_data, desc(rg))
#create a rank column for plotting
current_data$rank <- seq(1,nrow(current_data))

rankdata <- rbind(rankdata,current_data)


```

To obtain incremental MTAG, for each increment, we need a file to reference to that lists the phenotypes included in that increment, to be matched with their summary statistics later. So, incrementalMTAG on 5 disorders, would require 4 files of meta-analyses, like so: 
focalpheno - phenoA
focalpheno - phenoA - phenoB
focalpheno - phenoA - phenoB - phenoC
focalpheno - phenoA - phenoB - phenoC - phenoD
focalpheno - phenoA - phenoB - phenoC - phenoD - phenoE

each of these will be named 'increment 1 - 4' of this incremental MTAG analysis, and the increment nr will be looped over lateron. 

```r
#### output increments files
outdir <- paste0(wdir,"/output")

focus <- "LifetimeMDD"

for (i in seq(1,nrow(current_data))) {
new_increment <- current_data[1:i,]$p2
print(new_increment)
current_increments <- c(focus, new_increment)

#write the increment to file  -- increment nr = i
outfile <- paste0(outdir,"/",focus,".increment",i,".txt")
write.table(current_increments, outfile, col.names=FALSE, row.names=FALSE, 
            quote=FALSE)
}
```

# step 4: perform MTAG on each of the increments.

For easy upscaling make sure that sumstats are in the same directory and have names that are easily referenced like `${indir}/${pheno}.sumstats`

IMPORTANT: The incrementalMTAG was also performed as part of 10xCV PRS, and so, only done on 90% of total individuals. Make sure to obtain input GWASs for 90% of the individuals in your dataset!

As such, we looped over 10 folds and the 14 increments to perform incrementalMTAG. 


To create the sumstats flag to go into MTAG, I created a small R script that takes the output text files from above (for each increment) and outputs a comma-separated line with the full filedirectories


```r
#scriptname: conversion_variable_increments.R

library(dplyr)

############ outputting the orders to use in MTAG command..
## requirements for MTAG flag:
# --sumstats [{File1},{File2}...]
# Specify the list of summary statistics files to
# perform multitrait analysis. Multiple files paths must
# be separated by ",". Please read the documentation to
# find the up-to-date set of acceptable file formats. A
# general guideline is that any files you pass into MTAG
# should also be parsable by ldsc and you should take
# the additional step of specifying the names of the
# main columns below to avoid reading errors.

# conversiontable will provide: incrementfilename and fold
# for the concept list, each focus will have 28 increment filenames
# and then go over each of them 10 times. 


args <- commandArgs(trailingOnly = TRUE)

indir <- args[1]
focus <- args[2]
increment <- args[3]
increments_file <- args[4]
fold <- args[5]


increments <- read.table(increments_file, header=FALSE, stringsAsFactors=FALSE)

flaginfo <- c()
for (incr in increments) {
  
  current_file <- paste0(indir,"/training",fold,".",incr,".mtag.txt")
  
  flaginfo <- append(flaginfo, current_file)
}

cat(flaginfo, sep=",") #output to std out to use in bash after converting. 


q()
```

Then use that to perform the final MTAG runs:

```bash
sumstats_flag=$(Rscript ${code}/conversion_variable_increments.R ${indir} ${focus} ${increment} ${incrementfile} ${fold} )

echo ${sumstats_flag}

## to run MTAG, activate that environment
conda activate mtag

## run MTAG
mtag_dir="xxx/mtag"

${mtag_dir}/mtag.py \
--sumstats ${sumstats_flag} \
--n_min 0 \
--std_betas \
--cores 1 \
--stream_stdout \
--verbose \
--force \
--use_beta_se \
--se_name betase \
--out ${outdir}/${focus}.increment${increment}.fold${fold} 

```
