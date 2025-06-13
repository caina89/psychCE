rm( list=ls() )

N <- 1e4
M <- 1e2
P <- 4
h2 <- .3
rg <- .7
re <- .7

pdf( 'ExtendedDataFigure2a.pdf', w=12, h=5 )
par( mfcol=c(2,6), mar=c(4,4,1,1) )
for( zzz in 1:6 ){
  
  # simulate multiple traits
  Sig_g <- h2     * ( diag(P)*(1-rg) + rg )
  Sig_e <- (1-h2) * ( diag(P)*(1-re) + re )
  G <- matrix(rnorm(N*M),N,M)
  B <- matrix(rnorm(M*P),M,P) %*% chol( Sig_g/M )
  Y <- G %*% B + matrix(rnorm(N*P),N,P) %*% chol( Sig_e )
  
  # run GWAS
  gwas <- array( NA, dim=c(4,M,P), dimnames=list( c('beta','sd','t','p'), 1:M, 1:P ) )
  for( j in 1:M )
    for( p in 1:P )
      gwas[,j,p] <- summary( lm( Y[,p] ~ G[,j] ) )$coef[2,]
  
  # run mtag
  mtag_fxn <- function(betahats,Sigma,Omega){
    omt <- as.numeric(Omega[1,])
    omtt <- Omega[1,1]
    num   <- (omt/omtt) %*% solve( Omega - (omt %o% omt)/omtt + Sigma ) %*% t(betahats)
    denom <- (omt/omtt) %*% solve( Omega - (omt %o% omt)/omtt + Sigma ) %*% (omt/omtt)
    as.numeric(num)/as.numeric(denom)
  }
  b_mtag <- mtag_fxn( betahats=gwas[1,,], Sigma=cov( gwas[1,,] - B ), Omega=Sig_g/M ) 
  
  # build PGS in test data
  Gtst <- matrix(rnorm(N*M),N,M)
  Ytst <- Gtst %*% B + matrix(rnorm(N*P),N,P) %*% chol( Sig_e )
  
  prs     <- matrix( NA, N, 20 )
  prsmtag <- matrix( NA, N, 20 )
  for( chr in 1:20 ){
    snps_chr <- 1:(M/20) + (chr-1) * (M/20)
    prs     [,chr] <- Gtst[,snps_chr] %*% gwas[1,snps_chr,1]
    prsmtag [,chr] <- Gtst[,snps_chr] %*% b_mtag[snps_chr]
  }
  
  # run CE test
  run_ce <- function(y,prs){
    res <- array( NA, dim=c(4,20,20), dimnames=list( c('beta','sd','t','p'), 1:20, 1:20 ) )
    for( i in 1:20 )
      for( j in setdiff( 1:20, i ) )
        res[,i,j] <- summary( lm( y ~ prs[,i] * prs[,j] ) )$coef[4,]
    res
  }
  
  hist( run_ce( Ytst[,1], prs     )['p',,], breaks=11, main='', ylab='Frequency', xlab='pc-CE p-values (no MTAG)' )
  hist( run_ce( Ytst[,1], prsmtag )['p',,], breaks=11, main='', ylab='Frequency', xlab='pc-CE p-values (with MTAG)' )
  
}
dev.off()