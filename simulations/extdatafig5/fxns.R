## from phenix
mat.sqrt <- function(X,inv=FALSE){
  if( length(X) == 1 )
    if( inv ){
      return(1/sqrt(X))
    } else {
      return(sqrt(X))
    }
  eig <- eigen(X,symmetric=TRUE)
  Q   <- eig$vec
  
  negs  <- which( eig$val < 0 )
  if( length( negs ) > 0 ){
    warning( paste( 'There are', length( negs ), 'negative eigenvalues; they are being set to 1e-8' ) )
    eig$val[ negs ] <- 1e-8
  }
  
  if( ! inv ){
    Lam <- diag( sqrt(eig$val) )
  } else {
    Lam <- diag( 1/sqrt(eig$val) )
  }
  
  return( Q %*% Lam %*% t(Q) )
}

myplot <- function(pdfname,Ptrues, ylim=c(-.11,.11)){
  #rho_gs <- rho_gs * 100 
  nx <- dim(out)[2]
  pdf( pdfname, w=8, h=7 )
  par( mfrow=c(3,3), mar=c(4,4,1,1) )
  for( k in 1:3 ){
    P <- Ptrues[[k]]
    plot( c(-2,2), c(-2,2), type='n', xlab='Pathway 1', ylab='Pathway 2' )
    points( P[,1], P[,2], pch=16, col=c('black','purple')[threshfxn(P,k)+1] )
    
    plot( range(rho_gs), ylim, type='n', xlab='Genetic Correlation', ylab='Gamma-hat', axes=F )
    axis(1,at=rho_gs,labels=round(rho_gs,2))
    axis(2)
    abline( h = 0, col='grey' )
    for( i in 1:nx )
      mylines( rho_gs, out[k,i,1,,], col=i )
    
    ys <- (out[k,,4,,]<.05)
    plot( range(rho_gs), c(0,1), type='n', xlab='Genetic Correlation', ylab='CE power (p<.05)', axes=F )
    abline(h=.05,col='grey')
    for( i in 1:nx )
      lines( rho_gs, rowMeans( ys[i,,], na.rm=T ), col=i, lwd=2 )
    axis(1,at=rho_gs,labels=round(rho_gs,2))
    axis(2)
    if( k == 1 )
      legend( bty='n', 'topright', leg=c( 'Whole Genome', 'True P1*P2', 'True P1*P1', 'True wg-CE' ), fill=1:nx, cex=1.2 )
  }
  dev.off()
}

threshfxn <- function(liabs,k){
  if( k == 1 ){
    liabs <- apply( liabs, 1, sum )
  } else if( k == 2 ){
    liabs <- apply( liabs, 1, min )
  } else {
    liabs <- apply( liabs, 1, max )
  }
  as.numeric( liabs >  quantile(liabs,1-prev) )
}

mylines <- function( x, y, col, ... ){
  mu <- apply( y, 1, mean, na.rm=T )
  se <- apply( y, 1, function(z) sd( z, na.rm=T ) )
  points( x, mu, col=col, pch=16, cex=1.3, ... )
  lines ( x, mu, col=col, pch=16, lwd=1.7, ... )
  arrows( x, mu+se, x, mu-se, code=3, angle=90, length=0.03, col=col, lwd=1.4 )
} 