## Functions
pij = function( i, px, N, s ) 
{
  pxprime = px * ( s*(px+1) + 1 )/( 2*s*px + 1 )
  return( dbinom( i, 2*N, pxprime ) )
} 
#comment


probability.matrix = function( N, s, ng )
{
  mat = matrix(NA, nrow = 2*N+1, ncol = 2*N+1)
  
  ## Assigning probability values.
  for( i in 1:nrow(mat) )
    for( j in 1:ncol(mat) )
    {
      px = (j-1)/2/N
      mat[i,j] = dbinom( i-1, 2*N, px * ( s*(px+1) + 1 )/( 2*s*px + 1 ) )
    }
  
  ## Matrix multiplication ~ res = mat^ng
  res = mat
  if( ng >= 2 )
    for( g in 2:ng )
      res = res %*% mat
    
  return(res)
}

     
impl.sampling.fn = function( fn.opt )
{
  if( is.character(fn.opt) )
  {
    if( fn.opt == "dbinom" )
      fn = function(N,i,S,k) dbinom( i, 2*S, k/2/N )
    else if( fn.opt == "dhyper" )
      fn = function(N,i,S,k) dhyper( i, k, 2*N-k, 2*S )
    else
      stop("unknown sampling distribution.")
  } else if( is.function(fn.opt) )
  {
    fn = fn.opt
  } else
    stop( "sampling option is not correct." )
  
  return( fn )
}
      

impl.likelihood = function( p0, i1, S1, i2, S2, N, ng, s, 
                            ps1.fn = "dbinom", ps2.fn = "dbinom" )
{
  ## Computing the probability matrix.
  mat = probability.matrix(N, s, ng)
  
  ## Sampling function
  fn1 = impl.sampling.fn( ps1.fn )
  fn2 = impl.sampling.fn( ps2.fn )
  
  ## Computing the likelihood
  L = 0
  for( k1 in 0:(2*N) ) 
  {
    # probability of sampling at t1
    ps1 = p0[k1+1] * fn1(N,i1,S1,k1)
    
    # probability of sampling at t2
    ps2 = 0
    for( k2 in 0:(2*N) )
      ps2 = ps2 + mat[k2+1,k1+1] * fn2(N,i2,S2,k2)
    
    L = L + ps1 * ps2
  }   
  
  return(L)
}

## Main function
likelihood = function( data, N, s, 
                       p0 = NULL, ps1.fn = "dbinom", ps2.fn = "dbinom" )
{
  if( is.null(p0) )
    p0 = rep(1/(2*N+1),2*N+1)
  
  p2 = 1
  
  for( k in 1:(nrow(data)-1) )
    p2 = p2 * impl.likelihood( p0,  
                               data[k,2],   # i1
                               data[k,3],   # S1
                               data[k+1,2], # i2
                               data[k+1,3], # S2
                               N, 
                               data[k+1,1]-data[k,1], # t2-t1
                               s, 
                               ps1.fn, ps2.fn )
  
  return(p2)
}


## Generating pseudo data
Ne = 10
S1 = Ne
S2 = Ne

t1 = 0
t2 = 2
    
i1 = 10
i2 = 10
    
s = 0
    
data = matrix(NA,ncol = 3, nrow = 2)
data[1,] = c(t1,i1,S1)
data[2,] = c(t2,i2,S2)
    
likelihood(data,Ne,s)
# [1] 0.00390656



