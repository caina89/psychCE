rm(list=ls())


args <- commandArgs(trailingOnly = T)

print(args)


pheno_data=""
stam_data=""
design_data=""
kinship_data="" #kin_mean_diff_wide.RData -- would need code on how you made this

pheno_mame <- ""
h2=0.4
date_pheno <- FALSE # to indicate that the outcome was 0/1 not dates.


library(data.table)
library(PAFGRS)
library(Matrix)
library(survival)

pheno <- fread(pheno_data)
stam <- fread(stam_data)
design <- fread(design_data)
kin <- fread(kinship_data)


# Population sample:
kontrol2015 <- design[ kontrol2015I==1, pid ]
kontrol_children <- unique( unlist( dt_r2[ id1 %in% kontrol2015 & child==T, rel_ids] ) )
kontrol_parents <- unique( unlist( dt_r2[ id1 %in% kontrol2015 & parent==T, rel_ids] ) )
kontrol_gr_parents <- unique( unlist( dt_r2[ id1 %in% kontrol_parents & parent==T, rel_ids] ) )
kontrol_gr_gr_parents <- unique( unlist( dt_r2[ id1 %in% kontrol_gr_parents & parent==T, rel_ids] ) )
kontrol_gr_gr_gr_parents <- unique( unlist( dt_r2[ id1 %in% kontrol_gr_gr_parents & parent==T, rel_ids] ) )

pop_sample_1 <- c( kontrol2015, kontrol_children, kontrol_parents, kontrol_gr_parents, kontrol_gr_gr_parents, kontrol_gr_gr_gr_parents)
pop_sample_2 <- unique(unlist(c( kontrol2015, unique( dt_r2[id1 %in% kontrol2015, rel_ids ] ) )))

## Build data
dx_col <- which(names(pheno) %in% c("pid", pheno_name))
pheno <- pheno[,..dx_col]
data <- merge( stam[ ,c(1:5)], pheno, by="pid", all.x=T )

# Compute birth date
data[,birth_date := as.Date(fdato,format="%d/%m/%Y")]


dx_col <- which(names(data)==pheno_name)
data$dx <- data[[dx_col]]

kin <- fread("~/iPSYCH_rels/data/kinship_2.0.csv")

load("cc_pheno_dx.Rdata") ## what's in this?
design <- fread(design_data)


K = kin[,.(i=id1,j=id2,x=r_mean_diff/2)]
out[is.na(k_p)|k_p==0,k_p:=out[k_p>0,min(k_p,na.rm=T)]]
pheno = out[,.(id=pid,aff=dx_ind,thr=qnorm(1-k_pop),w=k_p/k_pop)]
pheno[!is.na(w)&!w==0&aff==1,.N,w]
pheno[aff==1,.N,w]
pheno[!is.na(w)&!w==0&aff==1,w:=1]
pheno[aff==1,.N,w]
probands <- design$pid
rm(design,kin,par_mat)


PAFGRS_est <- PAFGRS::FGRS_wrapper(proband_ids = probands,
                            K = K,
                            pheno = pheno[!is.na(aff),.(id,aff,thr,w)],
                            h2=h2, method="PAFGRS")

  
  
