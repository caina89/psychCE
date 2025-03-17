rm( list=ls() )

N <- 1e4
M <- 1e2

prev <- .2
h2s <- c(.4,.4)
rg <- 0
re <- 0
Sig_g <- diag(sqrt(h2s  )) %*% ( diag(2)*(1-rg) + rg ) %*% diag(sqrt(h2s  ))
Sig_e <- diag(sqrt(1-h2s)) %*% ( diag(2)*(1-re) + re ) %*% diag(sqrt(1-h2s))

niter <- 1e3
EOtest  <- array( NA, dim=c(4,4,niter) )
Dirtest <- array( NA, dim=c(3,4,niter) )
EOtest_oracle  <- array( NA, dim=c(4,4,niter) )
Dirtest_oracle <- array( NA, dim=c(4,4,niter) )
for( it in 1:niter ){
  
  # simulate two subtype liabilities
  G <- matrix(rnorm(N*M),N,M)
  B <- matrix(rnorm(M*2),M,2) %*% chol( Sig_g/M )
  Y <- G %*% B + matrix(rnorm(N*2),N,2) %*% chol( Sig_e ) 
  
  # define disease as people exceeding either liability limit
  maxLiabs <- apply( Y, 1, max )
  cc <- as.numeric( maxLiabs > quantile(maxLiabs,1-prev) )
  
  # run GWAS on training data
  gwas <- sapply( 1:M, function(j) summary( glm( cc ~ G[,j], subset=1:(N/2), family='binomial' ) )$coef[2,] )
  
  # build PRS on test data
  PS <- G[1:(N/2)+N/2,           ] %*% gwas[1,]
  PSe<- G[1:(N/2)+N/2,1:(M/2)    ] %*% gwas[1,1:(M/2)    ]
  PSo<- G[1:(N/2)+N/2,1:(M/2)+M/2] %*% gwas[1,1:(M/2)+M/2]
  
  cctst <- cc[1:(N/2)+N/2]
  EOtest [,,it] <- summary( glm( cctst ~ PSe * PSo   , family='binomial' ) )$coef
  Dirtest[,,it] <- summary( glm( cctst ~ PS + I(PS^2), family='binomial' ) )$coef
  
  #### now repeat, using subtypes
  # run GWAS on training data
  cc1  <- as.numeric( Y[,1] > quantile(maxLiabs,1-prev) )
  cc2  <- as.numeric( Y[,2] > quantile(maxLiabs,1-prev) )
  gwas1 <- sapply( 1:M, function(j) summary( glm( cc1 ~ G[,j], subset=1:(N/2), family='binomial' ) )$coef[2,] )
  gwas2 <- sapply( 1:M, function(j) summary( glm( cc2 ~ G[,j], subset=1:(N/2), family='binomial' ) )$coef[2,] )
  
  # build PRS on test data
  PS1 <- G[1:(N/2)+N/2,           ] %*% gwas1[1,]
  PS2 <- G[1:(N/2)+N/2,           ] %*% gwas2[1,]
  PSe1<- G[1:(N/2)+N/2,1:(M/2)    ] %*% gwas1[1,1:(M/2)    ]
  PSo2<- G[1:(N/2)+N/2,1:(M/2)+M/2] %*% gwas2[1,1:(M/2)+M/2]
  
  EOtest_oracle [,,it] <- summary( glm( cctst ~ PSe1 * PSo2   , family='binomial' ) )$coef
  Dirtest_oracle[,,it] <- summary( glm( cctst ~ PS1 * PS2, family='binomial' ) )$coef
  
  par( mfrow=c(2,2) )
  myhist <- function( x1, x2, xlab, breaks ){
    if( missing( breaks ) ){
      xlim <- c(-1,1)*max(abs(c(x1,x2)),na.rm=T)
      breaks <- seq( xlim[1], xlim[2], length=21 )
    } else {
      xlim <- range( breaks )
    }
    y1 <- hist( x1, breaks=breaks, plot=F )
    y2 <- hist( x2, breaks=breaks, plot=F )
    plot( xlim, range( c( y1$density, y2$density ) ), type='n', xlab=xlab, main='', ylab='Density' )
    lines( y1$mids, y1$density, lwd=2, col=1 )
    lines( y2$mids, y2$density, lwd=2, col=2 )
  }
  myhist( EOtest [4,4,], EOtest_oracle [4,4,], seq( 0, 1, len=21 ), xlab='EO test p-value' )
  myhist( Dirtest[3,4,], Dirtest_oracle[4,4,], seq( 0, 1, len=21 ), xlab='Direct test p-value' )
    legend( 'topright', fill=1:2, leg=c( 'Disease PRS', 'Subtype PRS' ), bty='n' )
  myhist( EOtest [4,1,], EOtest_oracle [4,1,], xlab='EO test Effect Size' ); abline( v=0, col=3 )
  myhist( Dirtest[3,1,], Dirtest_oracle[4,1,], xlab='Direct test Effect Size' );  abline( v=0, col=3 )
  
}
