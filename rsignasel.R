## Functions
pij = function( i, px, N, s ) 
{
  pxprime = 1. * px * ( s*(px+1.) + 1. )/(2.*s*px + 1.)
  return( dbinom( i, 2*N, pxprime ) )
}


probability_matrix = function( N, s, ng )
{
  mat = matrix(NA, nrow = 2*N+1, ncol = 2*N+1)
  
  for( i in 1:nrow(mat) )
    for( j in 1:ncol(mat) )
      mat[i,j] = pij( i-1, 1.*(j-1)/2/N, N, s )
    
    res = mat
    
    for( g in 2:ng )
      res = res %*% mat
    
    return(res)
}


r.__likelihood = function( i1, S1, i2, S2, N, ng, s )
{
  mat = probability_matrix(N, s, ng)
  
  prob_i1 = 1./(2*N+1-i1)
  
  # Declaring miscellaneous variables.
  L = 0.
  
  for( i in i1:(2*N) ) 
  {
    # probability of sampling at t1
    ps1 = dhyper( i1, i, 2*N-i, 2*S1 )
    
    # probability of sampling at t2
    ps2 = 0.
    for( k in i2:(2*N) )
      ps2 = ps2 + mat[k+1,i+1] * dhyper( i2, k, 2*N-k, 2*S2 )
    
    L = L + ps1 * ps2
  }
  
  return(L * prob_i1)
}

## Main function
r.likelihood = function( data, N, s )
{
  p2 = 1.
  
  for( k in 1:(dim(data)[1]-1) )
    p2 = p2 * r.__likelihood( data[k,2],   # i1
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

t1 = 1
t2 = 10

i1.true = 5
s = 0.1

i.tmp = i1.true

for( i in 1:(t2-t1) )
  i.tmp = rbinom(1, 2*Ne, i.tmp/2/Ne)

i1 = rhyper(1, i1.true, 2*Ne-i1.true, 2*S1)
i2 = rhyper(1, i.tmp, 2*Ne-i.tmp, 2*S2)

data = matrix(c(t1,t2,i1,i2,S1,S2),ncol=3)

like = sapply( 10:20, function(x) r.likelihood(data,x,0) )
plot(like)





