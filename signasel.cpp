#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/****************************************************************
    Matrix functions 
****************************************************************/
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
      
      for( auto k = size_t{0}, endk = res[0].size(); k < endk; ++k )
        tmp += A[k][j] * B[i][k];
    }
  
  return res;
}


/****************************************************************
    Simulation functions
****************************************************************/

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
    N: population size at generations (t+1)
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
        col[i] = gsl_ran_binomial_pdf( i, pxprime, 2*N );
    
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
    auto mat = make_matrix(2*N+1, 2*N+1);
    
    auto k = 0;
    for( auto& vec: mat )
        vec = probability_vector( double(k++)/2/N, N, s );
    
    auto res = mat;
    
    // Iterating the matrix over generations
    for( auto g = 1; g < ng; ++g )
        res = mult_matrix( res, mat );
        
    return res;
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
    auto mat = probability_matrix(N, s, ng);
    
    // Defining the prior of i1* (true value of i1), two vectors of:
    // - values of i1*.
    // - its probability.
    auto val_i1  = vector<int>(N+1-i1);
    iota( val_i1.begin(), val_i1.end(), i1 );
    
    auto prob_i1 = vector<double>( val_i1.size(), 1./val_i1.size() ); 
        
    // Declaring miscellaneous variables.
    auto L = 0.;
    
    for( auto i = size_t{0}, end = prob_i1.size(); i < end; ++i )
    {
        auto p0 =  double(val_i1[i])/2/N;
        
        // probability of sampling at t1
        auto ps1 = gsl_ran_binomial_pdf( i1, p0, 2*S1 );
        
        // probability pop at t1 (N) -> pop at t2 (N) (vector)
        auto v2 = mat[val_i1[i]];
        
        // probability of sampling at t2
        auto ps2 = 0.;
        for( auto k = i2, end = 2*N+1; k < end; ++k )
            ps2 += gsl_ran_binomial_pdf( i2, v2[k], 2*S2 );
        
        L += prob_i1[i] * ps1 * ps2;
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
    auto S = data.column(2);
    auto N = data.column(3);
    
    auto p2 = 1.;
    
    for(auto k = 0, end = data.rows()-1; k < end; ++k )
        p2 *= likelihood_impl( i[k], S[k], i[k+1], S[k+1], N[k+1], g[k+1]-g[k], s );
  
    return p2;
}


auto maximised_likelihood_s( IntegerMatrix const& data ) 
{
    /**********************************************************
    Finds the maximum likelihood on s. Returns L(smax) and smax.
    
    Using an easy grid computation for maximum searching.
    **********************************************************/
    
    // Declaring variables.
    auto step = 0.01;
    auto smax = 0.;
    auto Lmax = likelihood( data, 0 );
    
    // Searching the maximum likelihood of s.
    for( auto i = size_t{1}, end = size_t(1/step)+1; i < end; ++i )
    {
        auto Ltmp = likelihood(data, i*step);
        if( Lmax < Ltmp )
        {
            smax = i*step;
            Lmax = Ltmp;
        }
    }
    
    return make_pair(smax, Lmax);
}


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
    auto S = data.column(2);
    auto n = data.column(3);
    
    // Checking the data.
    for( auto k = 0, end = data.rows(); k < end; ++k )
    {
        if( g[k] <= 0 || S[k] <= 0 || n[k] <= 0 || 
            i[k] < 0 || S[k] > n[k] || i[k] > 2*S[k] )
            throw runtime_error( "unvalid parameters" );
    }
}


auto signaseltest( IntegerMatrix const& data )
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
    auto L0 = likelihood( data, 0 );
    
    // Computing the maximum likelihood value of s.
    auto x = maximised_likelihood_s( data );
    auto smax = x.first;
    auto Lmax = x.second;
    
    // Computing the likelihood ratio.
    auto LRT = -2 * log(L0 / Lmax);
    
    // Computing p-value assuming LRT follows a Chi-square low with 1 df.
    auto pvalue = -log10( 1 - gsl_cdf_chisq_P(LRT, 1) ); 
    auto res = NumericMatrix(1,5);
    res[0] = L0, res[1] = Lmax, res[2] = smax, res[3] = LRT, res[4] = pvalue;
    colnames(res) = CharacterVector::create("L0", "Lmax", "smax", "LRT", "-log10pvalue");
    // rownames(res) <- ""
    
    return res;
}
//////////////////////////////////////////////////////////////////


