rm( list=ls() )
set.seed(34)

if( file.exists( 'ExtendedDataFigure1.Rdata' ) ){
  load( 'ExtendedDataFigure1.Rdata' )
} else {
  niter <- 100
  N     <- 1e4
  S     <- 100
  h2    <- .3
  lams  <- seq(-1,1,len=9)
  out <- array( dim=c(2,4,length(lams),niter) )
}

for( it in 1:niter )
  for( lam.i in seq( along=lams ) )
{
  if( any( !is.na( out[1,,lam.i,it] ) ) ) next
  print( it )
  G     <- scale(matrix( rnorm(2*N*S), 2*N, S )) ## scale not needed but makes G_j^T G_j exactly 1
  betas <- sqrt(h2/(2*S)) * matrix( rnorm(S*2), S, 2 )
  P <- G %*% betas
  y <- P[,1] + P[,2] + lams[lam.i] * P[,1] * P[,2] + rnorm(2*N)*sqrt(1 - h2 - lams[lam.i]^2 * h2^2 / 4 )

  betahat <- t(G[1:N,]) %*% y[1:N]

  P1   <- G[N+1:N,1:(S/2)    ] %*% betahat[1:(S/2)    ]
  P2   <- G[N+1:N,1:(S/2)+S/2] %*% betahat[1:(S/2)+S/2]
  P    <- P1 + P2
  ytst <- y[N+1:N]

  P1  <- scale(P1) * sqrt( h2/4 )
  P2  <- scale(P2) * sqrt( h2/4 )
  P   <- scale(P)  * sqrt( h2/4 )

  out[1,,lam.i,it]  <- summary( lm( ytst ~ P1  + P2 + P1:P2 ) )$coef['P1:P2',]
  out[2,,lam.i,it]  <- summary( lm( ytst ~ P + I(P^2) ) )$coef['I(P^2)',]

  rm( G )
  save.image( 'ce_sims.Rdata' )
}

lamarray <- Zarray <- array( dim=c(2,length(lams),niter) )
for( lam.i in seq( along=lams ) )
  lamarray[,lam.i,]<- lam.i
for( j in 1:2 )
  Zarray  [j,,]<- j

pdf( 'ExtendedDataFigure1.pdf', w=10, h=5 )
par( mfrow=c(1,2), mar=c(4,4,1,1) )

boxplot( out[,1,,] ~ Zarray + lamarray, col=rep(c(2,4),3), axes=F, xlab='Gamma', ylab='Gamma-hat', ylim=c(-1.5,1.5) )
axis(1,at=.5+seq(1,length(lams)*2,2),labels=lams)
axis(2)
legend( 'topleft', fill=c(2,4,3), leg=c( 'Even/Odd', 'Whole Genome', 'Truth' ), bty='n' )
for( lam in lams )
  lines( which(lams==lam)*2 + c(-1,0), rep(lam,2), col=3, lwd=3 )
abline( h = 0, col='grey' )

boxplot( -log10(out[,4,,]+1e-20) ~ Zarray + lamarray, col=rep(c(2,4),3), axes=F, xlab='Gamma', ylab='CE p-value (-log10)' )
axis(1,at=.5+seq(1,length(lams)*2,2),labels=lams)
axis(2)

dev.off()
