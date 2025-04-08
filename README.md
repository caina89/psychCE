# psychCE

This repository documents the code used to implement the Coordinated Epsitasis (CE) framework as part of the investigations documented in [this preprint]().

environments used to run these analyses: for CE tests, GWAS and PRSice runs: xxx.yml
and for gSEM analyses xxx.yml

these were loaded into submission scripts using 
source ~/.bashrc
conda activate (ENVNAME)

## GWAS
GWAS that lie at the basis of the PRS that are used as input for pc-CE and wg-CE tests were obtained using plink2 {insert link}

```bash
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
to boost power, multiple meta-analyses were used. For the first set of MTAG analyses, with MTAG.All, sumstats were obtained from this paper (cite).
However, incremental MTAG analyses using those input phenotypes in that meta-analyses was run for this project specifically and is described in xxx file. 

furthermore, five psychiatric disorders were meta-analyzed using mtag as well, to compare with gSEM method and the code for this is listed xxx.

## genomic Structural Equation Moddeling (gSEM)
to investigate whether p-factor is heterogeneous. 
is at the level of input sumstats, GWAS on common P factor was obtained, code listed xxxx. 
factor loadings obtained using xxx. 

## Genetic Liability Scores
### Polygenic Risk Scores
were obtained from GWAS (see above), to deal with genotype quality issues in iPSYCH, upscaled to check if it first error'd out and then re-run with the adjusted file. 

if QC already done, like analyses that include munging with LDSC (such as gSEM), or where genotype quality control already excluded these variants, a single PRSice command was used:
```bash
Rscript ${prsicedir}/PRSice.R \
--prsice ${prsicedir}/PRSice_linux \
--base ${sumstats} \
--cov  ${covar} \
--extract ${QCfile} \
--target ${bfile} \
--pheno ${phenofile} \
--score "avg" \
--binary-target T \
--beta \
--a1 ${a1} \
--a2 ${a2} \
--pvalue ${p} \
--snp ${rsid} \
--chr ${chr} \
--stat ${beta} \
--seed 8 \
--bar-levels 1 \
--fastscore \
--num-auto 22 \
--print-snp \
--thread 1 \
--out ${outdir}/${meta}_5disPGC_${target}_${pheno}_${prev}_wholegenome

```

all other analyses would include an if statement to check whether QC was needed and would use PRSice to do that QC. 

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
were obtained by morten krebs following their paper and instructions (cite).
full code is listed in the xxx file. 


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


