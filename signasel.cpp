#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Rcpp.h>

#include "max_f.hpp"

using namespace Rcpp;
using namespace std;


namespace impl
{ // namespace impl

    /****************************************************************
        Matrix
    ****************************************************************/
    
    template < typename T >
    class matrix
    {
    private:
      
      std::unique_ptr< T[] > _data;
      size_t _n_row;
      size_t _n_col;
      
    public:
      
      matrix( size_t n_row, size_t n_col ):
      _data( std::make_unique<T[]>(n_col*n_row) ),
      _n_row( n_row ),
      _n_col( n_col )
      {
        /* No further construction here. */
      }
      
      matrix( matrix const& other ):
      _data( std::make_unique<T[]>(other._n_col*other._n_row) ),
      _n_row( other._n_row ),
      _n_col( other._n_col )
      {
        std::copy( &other._data[0],
                   &other._data[_n_col*_n_row],
                               &_data[0]
        );
      }
      
      auto operator*=( matrix<T> const& other )
        -> matrix<T>&
      {
        auto res_line = std::make_unique<T[]>(_n_col);
        
        for( auto i = size_t{0}; i < _n_row; ++i )
        {
          for( auto j = size_t{0}; j < _n_col; ++j )
          {
            auto& tmp = res_line[j] = 0;
            for( auto k = size_t{0}; k < _n_col; ++k )
              tmp += _data[i*_n_col+k] * other._data[k*_n_col+j];
          }
          
          auto it = &res_line[0];
          for( auto k = i*_n_col, end = (i+1)*_n_col; k < end; ++k )
            _data[k] = *it++;
        }
        
        return *this;
      }
      
      auto operator()( size_t i, size_t j ) const -> T { return _data[i*_n_col+j]; }
      auto operator()( size_t i, size_t j ) -> T& { return _data[i*_n_col+j]; }
      
      auto nrow() const { return _n_row; }
      auto ncol() const { return _n_col; }
    };
    


    /****************************************************************
        Simulation functions
    ****************************************************************/
    inline
    auto pij( size_t i, double px, int N, double s ) 
      -> double
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
        // auto w11 = 1. + 2 * s;
        // auto w12 = 1. + s;
        // auto w22 = 1.;
        // 
        // auto a = w11 - w12; 
        // auto b = w12 - w22;
        // auto c = a - b;
        // 
        // // Allele frequency after selection at generations t.
        // auto pxprime = (px * (px * a + w12)) / (px * (2 * b + px * c) + w22);
        
        auto pxprime = px * ( s*(px+1.) + 1. )/(2.*s*px + 1.);

        return gsl_ran_binomial_pdf( i, pxprime, 2*N );
    }


    auto probability_matrix( size_t N, double s, int ng )
    {
        /***********************************************************
        Create the Wright-Fisher recursion matrix over ng generations
        s: coefficient of selection
        N: population size (Ne), for generations 1 to ng-1
        ng: number of generations (iterations)
        ***********************************************************/

        // Making the (2N+1)*(2N+1) recursion matrix from one generation to the
        // next (pop size is always N)
        auto mat = matrix<double>(2*N+1, 2*N+1);

        for( auto i = size_t{0}, end = mat.nrow(); i < end; ++i )
          for( auto j = size_t{0}, jend = mat.ncol(); j < jend; ++j )
            mat(i,j) = pij( i, double(j)/2/N, N, s );

        auto res = mat;

        // Iterating the matrix over generations
        for( auto g = 1; g < ng; ++g )
            res *= mat;

        return res;
    }

    
    auto likelihood( size_t i1, size_t S1, size_t i2, size_t S2, 
                     size_t N, size_t ng, double s )
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

        // Generating the recursion matrix
        /// @dev: critically slow function.
        auto mat = probability_matrix(N, s, ng);

        // Defining the prior of i1* (true value of i1) as a uniform prior.
        auto prob_i1 = 1./(N+1-i1);

        // Declaring miscellaneous variables.
        auto L = 0.;

        for( auto i = i1, end = 2*N+1; i < end; ++i )
        {
            auto p0 =  double(i)/2/N;

            // probability of sampling at t1
            auto ps1 = gsl_ran_binomial_pdf( i1, p0, 2*S1 );

            // // probability pop at t1 (N) -> pop at t2 (N) (vector)
            // auto& v2 = mat[i];

            // probability of sampling at t2
            auto ps2 = 0.;
            for( auto k = i2, end = 2*N+1; k < end; ++k )
                ps2 += gsl_ran_binomial_pdf( i2, mat(k,i), 2*S2 );

            L += prob_i1 * ps1 * ps2;
        }

        return L;
    }

} // namespace impl



auto likelihood( IntegerMatrix const& data, int N, double s ) 
    -> double
{
    /***********************************************************
    Calculates the Likelihood for multiple time samples.
    Data is assumed to have the format of a matrix with nrow = 
    number of samples and for each sample the columns: g, i, S:
    g1 i1 S1
    g2 i2 S2
    g3 i3 S3
    etc.
    with g: generation, i: number of allele copies, S: sample size.
    
    We simply call %likelihood from g1 to g2 (ng=g2-g1), then from
    g2 to g3, etc. and multiply the resulting probabilities. s is the
    same for all samples and is maximized over all generations.
    ***********************************************************/
    
    auto p2 = 1.;
    
    for( auto k = 0, end = data.rows()-1; k < end; ++k )
        p2 *= impl::likelihood( data(k,1),   // i1
                                data(k,2),   // S1
                                data(k+1,1), // i2
                                data(k+1,2), // S2
                                N, 
                                data(k+1,0)-data(k,0), // t2-t1
                                s );
  
    return p2;
}


auto max_like_s( IntegerMatrix const& data, int N ) 
{
    /**********************************************************
    Finds the maximum likelihood on s. Returns L(smax) and smax.
    
    Using an easy grid computation for maximum searching.
    **********************************************************/
    using namespace std::placeholders;
    return max_f( bind(likelihood, data, N, _1), 0., 1. );
}


auto max_like_Ne( IntegerMatrix const& data, 
                  double s,
                  size_t min = 2, size_t max = 100
                  )
{
    using namespace std::placeholders;
    return max_f( bind(likelihood, data, _1, s), min, max );
}


auto checkparam( IntegerMatrix const& data )
    -> void
{
    /***********************************************************
    Data is assumed to have the format of a matrix with nrow = 
    number of samples and for each sample the columns: 
    g1 i1 S1
    g2 i2 S2
    g3 i3 S3
    etc.
    with g: generation, i: number of allele copies, S: sample 
    size.
    ***********************************************************/
       
    // Checking the data.
    for( auto k = 0, end = data.rows(); k < end; ++k )
    {
        if( data(k,0) <= 0 || data(k,2) <= 0 || /*n[k] <= 0 ||*/ 
            data(k,1) < 0 || /*S[k] > n[k] ||*/ data(k,1) > 2*data(k,2) )
            throw runtime_error( "unvalid parameters" );
    }
}


//[[Rcpp::export]]
NumericMatrix estimate_s( IntegerMatrix const& data, int N )
{
    /***************************************************************
    Compute the test statistics to detect selection. 
    Data is a matrix with as much rows as the number of samples. For each sample (row) the columns are:
    g1 i1 S1
    g2 i2 S2
    g3 i3 S3
    etc.
    with g: generation, i: number of allele copies, S: sample size, N: (effective) population size.
    ***************************************************************/
  
    // Checking the parameter values.
    checkparam( data );
    
    // Computing the likelihood of the null hypothesis (s=0).
    auto L0 = likelihood( data, N, 0 );
    
    // Computing the maximum likelihood value of s.
    auto x = max_like_s( data, N );
    auto smax = x.value;
    auto Lmax = x.likelihood;
    
    // Computing the likelihood ratio.
    auto LRT = -2 * log(L0 / Lmax);
    
    // Computing p-value assuming LRT follows a Chi-square low with 1 df.
    auto pvalue = -log10( 1 - gsl_cdf_chisq_P(LRT, 1) ); 
    auto res = NumericMatrix(1,5);
    res[0] = L0, res[1] = Lmax, res[2] = smax, res[3] = LRT, res[4] = pvalue;
    colnames(res) = CharacterVector::create("L0", "Lmax", "smax", "LRT", "-log10pvalue");
    rownames(res) = CharacterVector::create( "" );
    
    return res;
}


//[[Rcpp::export]]
NumericMatrix estimate_Ne( IntegerMatrix const& data )
{
    /***************************************************************
    Compute the test statistics to detect selection. 
    Data is a matrix with as much rows as the number of samples. For each sample (row) the columns are:
    g1 i1 S1
    g2 i2 S2
    g3 i3 S3
    etc.
    with g: generation, i: number of allele copies, S: sample size, N: (effective) population size.
    ***************************************************************/
  
    // Checking the parameter values.
    checkparam( data );
 
    // Computing the likelihood of the null hypothesis (s=0).
    auto x = max_like_Ne( data, 0 );
    
    auto res = NumericMatrix(1,2);
    res[0] = x.value, res[1] = x.likelihood;
    colnames(res) = CharacterVector::create("Ne", "Lmax" );
    rownames(res) = CharacterVector::create( "" );
    
    return res;
}


