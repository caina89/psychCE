# METAL combining iPSYCH and PGC summary statistics

Note to make sure that dataformat matches that for METAL. 
Mainly, if sumstats have OR, make sure to perform analysis on -log10(OR)

For submission upscaling, first create a 'lookup table' to reference the looped variables

```bash

wdir="xxx"
src=${wdir}/src

disorder_list=$(echo "SCZ" "BPD" "AUT" "ADHD" "ANO" "MDD" )
slurmvariable=1



for cohort in "2012" "2015i"
do
for pheno in ${disorder_list}
do
echo ${cohort} ${pheno} ${slurmvariable} >> ${src}/METAL_conversion.txt
slurmvariable=$((slurmvariable+1))
done
done

```

Next, use that in a submission bash script and submit. 

```bash
#!/bin/bash
#SBATCH --output=xxx.out
#SBATCH --error=xxx.err
#SBATCH -J metal
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10000 #in MB
#SBATCH --ntasks=1
#SBATCH --account xxx
#SBATCH --chdir xxx/output 


wdir="xxx"
source="${wdir}/src"

outdir="${wdir}/output"

src="xxx/src"

conversiontable=${source}/METAL_conversion.txt
cohort=$(cat ${conversiontable} | awk -v x="${SLURM_ARRAY_TASK_ID}" '{if ( $3 == x ) { print $1 } }')
pheno=$(cat ${conversiontable} | awk -v x="${SLURM_ARRAY_TASK_ID}" '{if ( $3 == x ) { print $2 } }')


sumsdir="xxx/data"

#these sumstats were formatted and header was edited to work with both METAL and MTAG
PGC_sumstats="${sumsdir}/PGC.${pheno}.GWAS.MTAGheader.txt"

## different location for iPSYCH sumstats, don't have beta column, made it.
ipsychdir="xxx/beta"
iPSYCH_sumstats=${ipsychdir}/iPSYCH${cohort}.${pheno}.GWAS.METALheader.txt

outfile="${pheno}_PGC-${cohort}_METAL_ .tbl" #IMPORTANT: note the space between the basename and extension here, if not added, METAL will fail silently

#location of METAL installation
metal_dir="xxx/generic-metal"


if [[ ${pheno} == "MDD" ]]; then
    FRQ="MAF"
elif [[ ${pheno} == "SCZ" ]]; then
    FRQ="MAF"
else 
    FRQ="EAF"
fi


# # ipsych doesnt have neff column, so calculated one for the whole thing. loading here

ipsych_Neff=$(cat ${wdir}/data/Neff_case_control_EurUnrel_perdisorder.txt | awk -v x=${cohort} -v y=${pheno} '{if ($1 == x && $2 == y ) { print $6 } }')

# echo ${outfile}

###### make the metal execution script
echo -e "
SCHEME STDERR \n\
CUSTOMVARIABLE TotalSampleSize \n\
LABEL TotalSampleSize as N_col \n\
\n\
MARKER ID \n\
ALLELE A1 A2 \n\
EFFECT B \n\
PVALUE P \n\
STDERR SE \n\
WEIGHT N_col \n\
PROCESS ${PGC_sumstats} \n\
\n\
MARKER ID \n\
ALLELE A1 A2 \n\
EFFECT BETA \n\
PVALUE P \n\
STDERR LOG.OR._SE \n\
WEIGHTLABEL DONTUSECOLUMN \n\
DEFAULTWEIGHT ${ipsych_Neff} \n\
PROCESS ${iPSYCH_sumstats} \n\
\n\
OUTFILE ${outfile} \n\
ANALYZE \n\
\n\
QUIT
" > ${src}/2024_05_23_METAL_execution.sh

#### execute
${metal_dir}/metal < ${src}/2024_05_23_METAL_execution.sh
####

exit 0
```
