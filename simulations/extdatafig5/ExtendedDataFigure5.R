rm( list=ls() )
set.seed(34)
source( 'fxns.R' )

if( file.exists( 'ExtendedDataFigure5.Rdata' ) ){
  load( 'ExtendedDataFigure5.Rdata' )
} else {
  niter <- 100
  N     <- 1e4
  S     <- 100
  h2    <- .3
  prev  <- .2
  rho_gs<- seq(0,1,len=5)
  
  out <- array( dim=c(3,4,4,length(rho_gs),niter) )
  Ptrues    <- list()
  liabtrues <- list()
}

for( it in 1:niter )
  for( lam.i in seq( along=rho_gs ) )
    for( k in 1:3 )
try({
  if( any( !is.na( out[k,1,,lam.i,it] ) ) ) next
  print( it )
  G     <- scale(matrix( rnorm(2*N*S), 2*N, S )) 
  betas <- sqrt(h2/(2*S)) * matrix( rnorm(S*2), S, 2 ) %*% mat.sqrt( matrix( c( 1, rho_gs[lam.i], rho_gs[lam.i], 1 ), 2, 2 ) )
  Ptrue <- G %*% betas

  #eps   <- matrix( rnorm(2*N*2), 2*N, 2 ) %*% mat.sqrt( (1-h2) * matrix( c( 1, rho_e, rho_e, 1 ), 2, 2 ) )
  eps   <- matrix( rnorm(2*N*2), 2*N, 2 ) * sqrt( 1-h2 )
  liabs <- Ptrue + eps
  y    <-  threshfxn(liabs,k)
  if( it == 1 & lam.i == 2 ){
    Ptrues[[k]] <- Ptrue[1:1e4,] 
    liabtrues[[k]] <- liabs[1:1e4,] 
  }
  
  ytst <- y[N+1:N]
  betahat <- t(G[1:N,]) %*% y[1:N]

  P    <- scale( G[N+1:N,] %*% betahat )
  P1  <- scale(Ptrue[N+1:N,1])
  P2  <- scale(Ptrue[N+1:N,2])
  P12 <- scale(rowSums(Ptrue)[N+1:N])

  out[k,1,,lam.i,it]  <- summary( glm( ytst ~ P  + I(P^2)    , family = binomial(link = "probit") ) )$coef['I(P^2)',]
  out[k,2,,lam.i,it]  <- summary( glm( ytst ~ P1 + P2 + P1:P2, family = binomial(link = "probit") ) )$coef['P1:P2',]
  out[k,3,,lam.i,it]  <- summary( glm( ytst ~ P1 + I(P1^2)   , family = binomial(link = "probit") ) )$coef['I(P1^2)',]
  out[k,4,,lam.i,it]  <- summary( glm( ytst ~ P12+ I(P12^2)  , family = binomial(link = "probit") ) )$coef['I(P12^2)',]

  save.image( file='ExtendedDataFigure5.Rdata' )
  myplot('ExtendedDataFigure5.pdf', Ptrues )
})