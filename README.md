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

#### Mundlak correction of PRS for CE
In the UKB, we obtained 10xCV PRS. We found that this introduces a bias in the PRS that needs to be corrected. The correction we implement is the Mundak method, which is described in detail in the Supplementary Text and Methods of the paper. 

The code used to perform mundlak correction on the PRS data is listed in `mundlak-correction.md`. 

### PA-FGRS
PA-FGRS were obtained following paper and [instructions listed here. ](https://github.com/BioPsyk/PAFGRS)
Code to obtain them in iPSYCH in our application is provided in `PAFGRS.R`


## Coordinated Epistasis implementation (and mundlak correction)
We implement CE framework in different ways, but the code to run in essence is the same for all instances. What changes per analysis is the input. 
It is always a single dataframe, with header `FID IID PHENO cov1 cov2...`
and then the PRS are added to that, either individual PRS, or per-chromosome PRS, or PA-FGRS. 

And then the one thing that changes is the definition of the formulas given to the regression models. 
Depending different nr's of PC's included it changes the `cov1 cov2 ... etc`
and changes also as follows: 
- EO-CE: `PHENO ~ cov1.. covX + PRSeven + PRSodd + PRSeven * PRSodd`
- pc-CE: `PHENO ~ cov1 .. covX + PRS1 + PRS2 ... + PRS22 + PRS1 * PRS2 + PRS1 * PRS3 .. + PRS21 * PRS 22` which can be obtained with small R code: 
```r
prs_times <- c()
pairs=c()

for (first in prs[1:22]) {
    for (second in prs[1:22]) {
        if (first != second & !paste0(second,":",first) %in% pairs ) {
                    firstchr <- strsplit(first,"_")[[1]][2]
        secondchr <- strsplit(second,"_")[[1]][2]

        prs_times <- append(prs_times, paste0(first," * ",second))

        pairs <- append(pairs,paste0(first,":", second))

        }
    }
}


nr_times <- length(prs_times) # to use for mean gamma. 
prs_times <- paste0(prs_times, collapse=" + ")

formula_H1 <- paste0("PHENO ~ " , covs_form, " + ", prs_plus, " + ", prs_times)

```
Also note that in UKB pc-CE 10xCV included Mundlak correction, adding `+ panelmean` to the formula.
- wg-CE and PA-FGRS CE are the same: `PHENO ~ cov1 ... covX + wgPRS + wgPRS * wgPRS `
- wg-CE * PA-FGRS CE: `PHENO ~ cov1 .. covX + wgPRS + PAFGRS + wgPRS * PAFGRS`.

formula_H0 always excludes all the `*` terms, formula_H1 always is the full formula as I list above. 

Obtaining the Standard Error's for the different models does change. For pc-CE we have a bit of code to obtain that. Note here the `nr_times` variable obtained from obtaining the formula above
```r
    #calculate meangammaSE -- original code from Andy, its different from std.err().. 
	Sig <- vcov(model_H1) ## Returns the variance-covariance matrix of the main parameters of a fitted model object.
	sub <- seq(nrow(Sig) - nr_times + 1 , nrow(Sig)) #rows corresponding to interaction effects between chr pairs #data has FID and IID, so ncoldata +1 if it did not, but -1 if it does. 
	Sig <- Sig[ sub, sub ] #only part of the covariate matrix we need
	sig2 <- mean( Sig ) ### this is the SE of the overall gamma
```

For all other interaction tests that don't generate more than one interaction term tested, the SE is just obtained from the regression model output. 

```r
    #output stderr
    all_err <- as.data.frame(summary(model_H1)$coefficients[,"Std. Error"])
    all_err$variables <- rownames(all_err)
    colnames(all_err) <- c("SE","variables")
    all_err <- dplyr::filter(all_err, grepl(paste0(":"),all_err$variables))
	SE <- all_err$SE

	stats[paste0(disorder,"_",model), "gammaSE"] <- SE

```

Then, for the three regression models: linear, and logistic regression with probit and logit link functions, the following is run on the input data: 

```r
## ---- create output matrix -------------

stats <- as.data.frame(matrix(nrow=3,ncol=4,NA))
colnames(stats) <- c("pheno_model","meangamma","meangammaSE","loglikpval")
rownames(stats) <- c(paste0(pheno,"_logit"),paste0(pheno,"_probit"),paste0(pheno,"_lm"))
# adding the rownames as a column as well for correct output formatting, the rownames are for indexing.
stats$pheno_model <- c(paste0(pheno,"_logit"),paste0(pheno,"_probit"),paste0(pheno,"_lm"))


# ---- training the models ---
## testing bmi, so only lm. 
for (model in c("logit","probit","lm")) {
	if (model == "logit") {
		fam <- binomial(link="logit")
	} else if (model == "probit") {
		fam <- binomial(link="probit")
	} else if (model == "lm") {
		fam <- gaussian()
	}
	
	model_H1 <- glm(formula_H1,data=data,family=fam)
	model_H0 <- glm(formula_H0,data=data,family=fam)


	# the statistics
	## log likelyhood ratio test
	#overallPval <- logLikPVal(model_H0,model_H1)
	overallPval <- anova(model_H1,model_H0,test="LRT")$"Pr(>Chi)"[2]
	
	#write pvalue to output matrix
	stats[paste0(pheno,"_",model), "loglikpval"] <- overallPval

	## calculating the mean gamma and its Standard Error.
	allgamma <- as.data.frame(summary(model_H1)$coefficients[,"Estimate"])
	allgamma$variables <- rownames(allgamma)

	allgamma <- filter(allgamma, grepl(paste0(":"),allgamma$variables))
	colnames(allgamma) <- c("Estimate","varaibles")
	meangamma <- mean(allgamma$"Estimate") ## this is the mean gamma

    #!!! OBTAIN SE HERE, SEE ABOVE!!!

	#write meangamma and SE to output matrix
	stats[paste0(pheno,"_",model), "meangamma"] <- meangamma
	stats[paste0(pheno,"_",model), "meangammaSE"] <- sqrt(sig2)
	
	# document the specific model output -- for H1
	model_out <- as.data.frame(summary(model_H1)$coefficients)
	colnames(model_out) <- c("Estimate","Std_error","z_value","Pvalue(z)")
	model_out$variables <- rownames(model_out)

	modeloutfile <- paste0(outdir,"/",PRSname,".base-",base,".target-",target,".",pheno,".chr.",model,"_full_H1_output.txt")
	write.table(model_out, modeloutfile, quote=FALSE, col.names=TRUE, row.names=FALSE)

}


```



## Simulations
All code to replicate our simulations is given in the folder `simulations`
