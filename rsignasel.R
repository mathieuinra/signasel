## Functions
pij = function( i, px, N, s ) 
{
  pxprime = px * ( s*(px+1) + 1 )/( 2*s*px + 1 )
  return( dbinom( i, 2*N, pxprime ) )
}


probability.matrix = function( N, s, ng )
{
  mat = matrix(NA, nrow = 2*N+1, ncol = 2*N+1)
  
  ## Assigning probability values.
  for( i in 1:nrow(mat) )
    for( j in 1:ncol(mat) )
      mat[i,j] = pij( i-1, 1.*(j-1)/2/N, N, s )
  
  ## Matrix multiplication ~ res = mat^(ng-1)
  res = mat
  if( ng >= 2 )
  {
    for( g in 2:ng )
      res = res %*% mat
  }
    
  return(res)
}


impl.likelihood = function( p0, i1, S1, i2, S2, N, ng, s, ps1.replacement = T, ps2.replacement = T )
{
  mat = probability.matrix(N, s, ng)
  
  L = 0.
  
  for( k1 in 0:(2*N) ) 
  {
    # probability of sampling at t1
    if(ps1.replacement)
    {
      ps1 = p0[k1+1] * dbinom( i1, 2*S1, k1/2/N ) # with replacement
    } else
      ps1 = p0[k1+1] * dhyper( i1, k1, 2*N-k1, 2*S1 ) # without replacement
      
      
    # probability of sampling at t2
    ps2 = 0.
    for( k2 in 0:(2*N) )
    {
      if(ps2.replacement)
      {
        tmp = dbinom( i2, 2*S2, k2/2/N ) 
      } else
        tmp = dhyper( i2, k2, 2*N-k2, 2*S2 )
        
      ps2 = ps2 + mat[k2+1,k1+1] * tmp
     }
    
    L = L + ps1 * ps2
  }   
  
  return(L)
}

## Main function
likelihood = function( data, N, s, p0 = NULL )
{
  if( is.null(p0) )
    p0 = rep(1/(2*N+1),2*N+1)
  
  p2 = 1.
  
  for( k in 1:(nrow(data)-1) )
    p2 = p2 * impl.likelihood( p0,  
                               data[k,2],   # i1
                               data[k,3],   # S1
                               data[k+1,2], # i2
                               data[k+1,3], # S2
                               N, 
                               data[k+1,1]-data[k,1], # t2-t1
                               s )
  
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
    
likelihood(data,N,s)
# print 0.002341363


