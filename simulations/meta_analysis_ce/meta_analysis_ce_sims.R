rm( list=ls() )

N <- 1e3
M <- 1e3
P <- 2
rg <- +.7
re <- -.2

r2out <- matrix( NA, 1000, 3 )
pvout <- array( NA, dim=c( 1000, 20, 20, 3 ) )

for( zzz in 1:1000 ){
  print(zzz)
  
  ###1 two phenotypes, purely additive, with some rg 
  ##LIES: its P=4, not 2 here.
  #(TODO: bonus point: with local rg)
  h2s   <- seq( .4, .9, len=P )
  Sig_g <- diag(sqrt(h2s  )) %*% ( diag(P)*(1-rg) + rg ) %*% diag(sqrt(h2s  ))
  Sig_e <- diag(sqrt(1-h2s)) %*% ( diag(P)*(1-re) + re ) %*% diag(sqrt(1-h2s))
  G <- matrix(rnorm(N*M),N,M)
  B <- matrix(rnorm(M*P),M,P) %*% chol( Sig_g/M ) #### can be done with mvnorm
  Y <- G %*% B + matrix(rnorm(N*P),N,P) %*% chol( Sig_e ) # note: chol is the matrix version of sqrt
  
  Gtst <- matrix(rnorm(N*M),N,M)
  Ytst <- Gtst %*% B + matrix(rnorm(N*P),N,P) %*% chol( Sig_e )
  
  #2 run mtag on them (I think I can do this in R, I just use the equation from Turley NatGen, right? no need for their software, right?)
  #2a: run GWASes
  gwas <- array( NA, dim=c(4,M,P), dimnames=list( c('beta','sd','t','p'), 1:M, 1:P ) )
  for( j in 1:M )
    for( p in 1:P )
      gwas[,j,p] <- summary( lm( Y[,p] ~ G[,j] ) )$coef[2,]
  
  #2b mtag the gwas
  mtag_fxn <- function(betahats,Sigma,Omega){
    omt <- as.numeric(Omega[1,])
    omtt <- Omega[1,1]
    num   <- (omt/omtt) %*% solve( Omega - (omt %o% omt)/omtt + Sigma ) %*% t(betahats)
    denom <- (omt/omtt) %*% solve( Omega - (omt %o% omt)/omtt + Sigma ) %*% (omt/omtt)
    as.numeric(num)/as.numeric(denom)
  }
  b_mtag <- mtag_fxn( betahats=gwas[1,,], Sigma=cov( gwas[1,,] - B ), Omega=Sig_g/M ) # empirical variance = something andy did
  b_mtag2<- mtag_fxn( betahats=gwas[1,,], Sigma=Sig_g/M + Sig_e/N, Omega=Sig_g/M ) #theoretical variance, = method
  
  #3: test CE
  #3a: build PGS
  prs <- matrix( NA, N, 20 )
  prsmtag <- matrix( NA, N, 20 )
  prsmtag2<- matrix( NA, N, 20 )
  for( chr in 1:20 ){
    snps_per_chr <- M / 20
    snps_chr <- 1:snps_per_chr + (chr-1) * snps_per_chr
    prs     [,chr] <- Gtst[,snps_chr] %*% gwas[1,snps_chr,1]
    prsmtag [,chr] <- Gtst[,snps_chr] %*% b_mtag[snps_chr]
    prsmtag2[,chr] <- Gtst[,snps_chr] %*% b_mtag2[snps_chr]
  }
  
  #3b: then run CE tests on the focal trait
  run_ce <- function(y,prs){
    res <- array( NA, dim=c(4,20,20), dimnames=list( c('beta','sd','t','p'), 1:20, 1:20 ) )
    for( i in 1:20 )
      for( j in setdiff( 1:20, i ) )
        res[,i,j] <- summary( lm( y ~ prs[,i] * prs[,j] ) )$coef[4,]
    res
  }
  
  pvout[zzz,,,1] <-  run_ce( Ytst[,1], prs )['p',,] #standard PRS (basic)
  pvout[zzz,,,2] <-  run_ce( Ytst[,1], prsmtag )['p',,] #mtag PRS. (focal 1)
  pvout[zzz,,,3] <-  run_ce( Ytst[,1], prsmtag2 )['p',,] #mtag PRS2 (focal 2)
  
  # par( mfrow=c(1,3) )
  # hist( pvout[,,,1] )
  # hist( pvout[,,,2] )
  # hist( pvout[,,,3] )
  # 
  
  #### this part is just for checking MTAG
  r2out[zzz,] <-  c(
    summary( lm( Ytst[,1] ~ prs ))$r.squared,
    summary( lm( Ytst[,1] ~ prsmtag ))$r.squared,
    summary( lm( Ytst[,1] ~ prsmtag2 ))$r.squared
  )
  ## correlation phenotype and three PRS
  
  # lims <- range( r2out, na.rm=T )
  # plot( r2out[,2], r2out[,1], xlim=lims, ylim=lims )
  # points( r2out[,2], r2out[,3], col=2, pch=16 )
  # abline( a=0, b=1 )
}

save(r2out, file="./1000iterations_r2out.RData")

save(pvout, file="./1000iterations_pvout_pval3.RData")



