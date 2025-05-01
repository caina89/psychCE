# psychCE

This repository documents the code used to implement the Coordinated Epsitasis (CE) framework as part of the investigations documented in [this preprint]().

environments used to run these analyses: for CE tests, GWAS and PRSice runs: `CE-env.yml`, for gSEM analyses `gSEM-env.yml`, for ldsc [this](https://github.com/bulik/ldsc/blob/master/environment.yml) `environment.yml` file.

these were loaded into submission scripts using 
`source ~/.bashrc`
`conda activate <ENVNAME>`

## GWAS
GWAS that lie at the basis of the PRS that are used as input for pc-CE and wg-CE tests were obtained using [plink2](https://www.cog-genomics.org/plink/2.0/)

```bash
bfile= #genotype files, plink format, excluding extension
phenofile= #binary phenotype file, header: FID IID PHENO
covar= #covariate file, header FID IID cov1 cov2 ...
outfile= #full basename of output file to be created
memory_use= #on cluster submission to match with specifications from job submission to avoid out-of-memory job failure

## plink
${plinkdir}/plink2 \
--bfile "$bfile" \
--pheno "${phenofile}" \
--1 \
--glm hide-covar \
--covar-variance-standardize \
--covar "${covar}" \
--out "${outfile}" \
--memory ${memory_use} \
--threads 1 

```

## GWAS meta-analyses
We performed multiple meta-analyses using both METAL and MTAG, the latter also in an incremental manner. 

### METAL
The [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) commanline tool was used to perform inverse-variance weigthed meta-analyses on PGC and iPSYCH summary statistics. The code used to run this is listed in the file `meta-analysis_METAL.md` in this Github.

### MTAG
The [MTAG] package was used to perform rG weighted meta-analysis across psychiatric disorders. Code used to run this on iPSYCH and PGC summary statistics: 

```bash
sumstats_flag= #comma-separated list of input summary statistics full filenames -- MTAG will run on all input sumstas as focals once, listing them as trait_1, trait_2 output, matching input list order.

outfile=${outdir}/MTAG_5disPGC #full basename of output sumstats

#location of MTAG installation
mtag_dir="xxx/mtag"

${mtag_dir}/mtag.py \
--sumstats ${sumstats_flag} \
--n_min 0 \
--std_betas \
--cores 1 \
--stream_stdout \
--verbose \
--force \
--snp_name ID \
--a1_name A1 \
--a2_name A2 \
--eaf_name EAF \
--z_name Z_STAT \
--n_name N_col \
--chr_name CHR \
--bpos_name POS \
--out ${outfile}

```

### incremental MTAG
To perform incremental MTAG, first phenotypes are ranked by their rG with the focal phenotype (LifetimeMDD in our case) and then incrementally meta-analyses are performed on those ranked phenotypes.
Code description in `incremental-MTAG.md`. 


## genomic Structural Equation Moddeling (gSEM)
to investigate whether p-factor is heterogeneous. 
is at the level of input sumstats, GWAS on common P factor was obtained, as well as factor loadings of the model. For this, code is described in `meta-analyses_gSEM.md`

## Genetic Liability Scores
We used two different methods to estimate genetic liability to an outcome. 
1. Polygenic Risk Scores
2. Pearson-Aitkens Family Genetic Risk Scores. 

### Polygenic Risk Scores
PRS were obtained from GWASs or meta-analyzed GWASs. We used PRSice-2 to obtain them. 
For those that were not munged before, there is a risk that there were still some duplicate SNPs and other minor issues. 
PRSice-2 checks for these, and the outputs a file that lists those SNPs to be excluded. We utilize this functionality by running PRSice-2 once, and if it creates the output file, re-run PRSice-2 with an `--exclude` flag added like so: 

```bash
Rscript ${packagebase}/PRSice/PRSice.R \
--prsice ${packagebase}/PRSice/PRSice_linux \
--base ${basefile} \
--a1 A1 \
--a2 REF \
--pvalue P \
--stat OR \
--snp ID \
--binary-target T \
--or \
--pheno ${phenobase} \
--cov ${covar} \
--base-maf MAF:0.05 \
--base-info INFO:0.9 \
--target-list ${source}/genolist_set${split}_${half}.txt \
--keep ${keepbase}/FID_IID_test_${fold}.txt \
--thread 1 \
--memory 30gb \
--out ${outdir}/prsice.${symp}.fold${fold}.${half}

## this detects duplicate snps, check if there are any created and then rerun the PRSice command with the --extract .valid file.
if test -f ${outdir}/prsice.${symp}.fold${fold}.${half}.valid 
then
Rscript ${packagebase}/PRSice/PRSice.R \
--prsice ${packagebase}/PRSice/PRSice_linux \
--base ${basefile} \
--a1 A1 \
--a2 REF \
--pvalue P \
--stat OR \
--snp ID \
--binary-target T \
--or \
--pheno ${phenobase} \
--cov ${covar} \
--base-maf MAF:0.05 \
--base-info INFO:0.9 \
--target-list ${source}/genolist_set${split}_${half}.txt \
--keep ${keepbase}/FID_IID_test_${fold}.txt \
--thread 1 \
--memory 30gb \
--out ${outdir}/prsice.${symp}.fold${fold}.${half} \
--extract ${outdir}/prsice.${symp}.fold${fold}.${half}.valid 
fi

```

### PA-FGRS
PA-FGRS were obtained following paper and [instructions listed here. ](https://github.com/BioPsyk/PAFGRS)
Code to obtain them in iPSYCH in our application is provided in `PAFGRS.R`


## Coordinated Epistasis implementation (and mundlak correction)


There are four different implementations of the CE framework described in that work: 

1. As published: Even-Odd partition CE (EO-CE)
2. As published: per-chromosome partition CE (pc-CE)
3. Whole genome CE (wg-CE)
4. family risk score-based CE (fgrs-CE)

For implementation 1, 2 and 3, polygenic risk scores are used as pathway proxies within the CE framework. 
xxx files list how those PRS were obtained in the UK Biobank
xxx files list how those PRs were obatined in iPSYCH

Mundlak correction for 10-fold cross-validated whole-genome PRS and per-chromosome PRS are listed in xxx and xxx respectively.

xxx lists how implementation 4 obtained [PA-FGRS]() scores in the Danish Register

Final CE tests all come down to a loglikelihood ratio test between two regression models, as this github describes in detail. 
Code used for these final tests are listed xx and xxx. 


-- also note logit/probit/lm 

## Simulations
Code to run the two simulations presented in the pre-print are xxx and xxx. 


