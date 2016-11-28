#include <cmath>
#include <gsl/gsl_randist.h>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

using matrix_type = vector< vector< double > >; // col then row!

auto make_matrix( int nrow, int ncol )
    -> matrix_type
{
    return matrix_type( ncol, vector<double>(nrow) );
}


matrix_type mult_matrix( matrix_type const& A, matrix_type const& B )
{
  auto res = make_matrix( A[0].size(), B.size() ); 
  
  for( auto i = size_t{0}, end = res.size(); i < end; ++i ) // for each col
    for( auto j = size_t{0}, endj = res[0].size(); j < endj; ++j ) // for each row
    {
      auto& tmp = res[i][j];
      
      for( auto k = 0ul, endk = res[0].size(); k < endk; ++k )
        tmp += A[k][j] * B[i][k];
    }
  
  return res;
}


auto probability_vector( double px, int N, double s) 
{
    /***********************************************************
    Gives the vector of probabilities of all possible allele
    numbers (0 to 2N) in a population of size N, when the frequency
    at previous generation was px and selection coefficient is s,
    according to a Wrigt-Fisher model. This vector will serve as a
    column of the recursion matrix, so f(i,j) is element a(j,i) of
    the matrix.
    
    px: allele frequency at generations t.
    N: population sizes at generations (t+1)
    s: coefficient of selection
    ***********************************************************/
    
    
    // Fitnesses of the 3 genotypes A1A1, A1A2, A2A2 (A1 is the selected allele).    
    auto w11 = 1. + 2 * s;
    auto w12 = 1. + s;
    auto w22 = 1.;
    
    auto a = w11 - w12;
    auto b = w12 - w22;
    auto c = a - b;
    
    // Allele frequency after selection at generations t.
    auto pxprime = (px * (px * a + w12)) / (2 * px * b + w22 + px * px * c);
    
    // Probabilities of all allele numbers at generations (t+1)
    auto nn2 = 2*N+1;
    auto col = vector<double>(nn2);
    
    for( auto i = size_t{0}; i < nn2 ; ++i )
        col[i] = gsl_ran_binomial_pdf( i, pxprime, nn2 );
    
    return col;
}


auto probability_matrix( int N, double s, int ng )
{
    /***********************************************************
    Create the Wright-Fisher recursion matrix over ng generations
    s: coefficient of selection
    N: population size (Ne), for generations 1 to ng-1
    ng: number of generations (iterations)
    ***********************************************************/
    
    // Making the (2N+1)*(2N+1) recursion matrix from one generation to the
    // next (pop size is always N)
    auto res = make_matrix(2*N+1, 2*N+1);
    
    auto k = 0;
    for( auto& vec: res )
        vec = probability_vector( double{k++}/2/N, N, s );
    
    // Iterating the matrix over generations
    for( auto g = 1; g < ng; ++g )
        res = mult_matrix( res, mat );
        
    return mat;
}



auto likelihood_impl( int i1, int S1, int i2, int S2, int N, int ng, double s )
    -> double
{
    /***********************************************************
    Calculating the Likelihood, i.e. simply multiplies the recursion
    matrix by the vector of initial frequencies.
    
    p0 is the supposedly known allele frequency in the pop at
    generation 1, from which the (known) sample (i1,S1) was drawn,
    and from which the iteration starts. It is the responsability
    of the calling function to choose p0. 
    
    p0 must be passed as a matrix of the form [[prob,value]]
    recall slatkin:
    L=Pr(p1|N)Pr(i1|p1,S1)Pr(i2|i1,p1,S1,S2,N)
    ***********************************************************/
    
    // Generating the recursion matrix.
    auto mat = MatProb(N, s, ng);
    
    // Slatking stuffs.
    probp0 <- matp0[,1] // ??
    valp0 <- matp0[,2]  // ??
    ## normalization
    probp0 <- probp0 / sum(probp0)
    stopifnot(sum(probp0) == 1)
        
    // Declaring miscellaneous variables.
    auto L = 0.;
    
    for (kp0 in seq_along(probp0)) 
    {
        auto p0 = valp0[kp0]
        
        // probability of sampling at t1
        auto ps1 = gsl_ran_binom_pdf( i1, 2*S1, p0 );
        
        // probability pop at t1 (N) -> pop at t2 (N) (vector)
        // v1 <- rep(0, 2*N+1)
        // i0 <- round(p0 * 2 * N)
        // v1[i0+1] <- 1 ## i+1 here because i is in 0:(2n) but the first
        //              ## element of a vector is 1 in R
        // v2 <- mat %*% v1
        auto v2 = mat[size_t(p0*2*N)];
        
        // probability of sampling at t2
        auto ps2 = 0.;
        for( auto k = 0, end = 2*N+1; k < end; ++k )
            ps2 += gsl_ran_binom_pdf( i2, 2*S2, v2[k] );
        
        L += probp0[kp0] * ps1 * ps2
    }
    
    return L;
}



auto likelihood( IntegerMatrix const& data, double s ) 
    -> double
{
    /***********************************************************
    Calculates the Likelihood for multiple time samples.
    Data is assumed to have the format of a matrix with nrow = 
    number of samples and for each sample the columns: g, i, S, N:
    g1 i1 S1 N1
    g2 i2 S2 N2
    g3 i3 S3 N3
    etc.
    with g: generation, i: number of allele copies, S: sample size, 
    N: (effective) population size.
    We simply call %likelihood from g1 to g2 (ng=g2-g1), then from
    g2 to g3, etc. and multiply the resulting probabilities. s is the
    same for all samples and is maximized over all generations.
    ***********************************************************/
    
    // Getting the data.
    auto g = data.column(0);
    auto i = data.column(1);
    auto S = data.column(3);
    auto n = data.column(4);
    
    auto p2 = 1.;
    
    for(auto k = 0, end = data.rows()-1; k < end; ++k )
        p2 *= likelihood_impl( i[k], S[k], i[k+1], S[k+1], N[k+1], g[k+1]-g[k], s );
  }
  
  return p2;
}


################################################################
WFMaxiLike2S <- function(data, maxs) {
################################################################
  ## Finds the maximum likelihood on s, returns L(smax), smax
  ## Uses the 'optimize' function of R
  ## WARNING: this function does not find the max in some
  ## cases. Must be checked
  ##
  ## IMPORTANT: maxs gives min and max possible values for s 
  ## must be realistic (??)  
################################################################
  res <- optimize(WFLike2S, interval=c(-maxs, maxs), 
                  data,
                  maximum = TRUE)
  warn <- 0
  ## Issue warning message if maxs is reached
  if ((res$maximum - maxs)**2 < 1e-6) {
    print(paste("Warning: Maximum value smax = ", maxs,
                  " was reached."))
    warn <- 1
    ## it would be possible to launch something here (exhaustive search?)
    ## however i doubt it improves, because likelihood is continuous
    ## and has a single maximum (this must be checked!) hence we are
    ## at the extreme s value. Still it is BAD (unrealistic?)
    ## return value 
  }
  return(c(res$objective, res$maximum, warn))
}
################################################################


auto checkparam( IntegerMatrix const& data )
    -> void
{
    /***********************************************************
    Data is assumed to have the format of a matrix with nrow = 
    number of samples and for each sample the columns: 
    g1 i1 S1 N1
    g2 i2 S2 N2
    g3 i3 S3 N3
    etc.
    with g: generation, i: number of allele copies, S: sample 
    size, N: (effective) population size
    ***********************************************************/
    
    // Getting the data.
    auto g = data.column(0);
    auto i = data.column(1);
    auto S = data.column(3);
    auto n = data.column(4);
    
    // Checking the data.
    for( auto k = 0, end = data.rows(); k < end; ++k )
    {
        if( g[k] <= 0 || S[k] <= 0 || n[k] <= 0 || 
            i[k] < 0 || S[k] > n[k] || i[k] > 2*S[k] )
            throw runtime_error( "unvalid parameters" );
    }
}


auto signseltest( IntegerMatrix const& data, maxs = 1 )
    -> NumericMatrix
{
    /***************************************************************
    Compute the test statistics to detect selection. 
    Data is a matrix with as much rows as the number of samples. For each sample (row) the columns are:
    g1 i1 S1 N1
    g2 i2 S2 N2
    g3 i3 S3 N3
    etc.
    with g: generation, i: number of allele copies, S: sample size, N: (effective) population size.
    ***************************************************************/
  
    // Checking the parameter values.
    checkparam( data );
    
    // Computing the likelihood of the null hypothesis (s=0).
    auto L0 = WFLike2S( data, 0 );
    
    // Computing the maximum likelihood (?): what is maximised?
    auto x = WFMaxiLike2S(data, maxs);
    auto Lmax = x[1];
    auto smax = x[2];
    auto warn = x[3];
    
    // Computing the likelihood ratio.
    auto LRT = -2 * log(L0 / Lmax);
    
    // Computing p-value assuming LRT follows a Chi-square low with 1 df.
    auto pvalue = -log10(1 - pchisq(LRT, 1)); /// TODO: pchisq
    auto res = NumericMatrix(1,5);
    res[0] = L0, res[1] = Lmax, res[2] = smax, res[3] = LRT, res[4] = pvalue, res[5] = warn;
    colnames(res) = CharacterVector::create("L0", "Lmax", "smax", "LRT", "-log10pvalue", "warn");
    // rownames(res) <- ""
    
    return res
}
//////////////////////////////////////////////////////////////////


