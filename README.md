# psychCE

## Coordinated Epistasis implementation
This repository documents the code used to implement the Coordinated Epsitasis (CE) framework as part of the investigations documented in [this preprint]().
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

## Simulations
Code to run the two simulations presented in the pre-print are xxx and xxx. 


## Meta-analyses
MTAG and METAL here

### incremental MTAG
rG estimation to order the phenotypes and then MTAG details for running the incremental MTAG on lifetimeMDD outcome with MTAG.All phenotypes.

## genomic Structural Equation Modelling to obtain P-factor PRS into wg-CE

