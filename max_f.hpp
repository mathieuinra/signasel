#pragma once 

#include <cmath>
#include <functional>

template < typename T >
struct like
{
    T value;
    double likelihood;
};
    

template < typename F, typename T >
auto max_f( F f, T min, T max, int max_counter = 1000 )
  -> like<T>
{
    // Declaring variables.
    auto gr = (std::sqrt(5) + 1) / 2;
    auto tol = 1e-5;
    auto a = min;
    auto b = max;
    
    auto c = b - (b - a) / gr;
    auto d = a + (b - a) / gr;
    
    auto counter = -1;
    
    while( std::abs(c - d) > tol && ++counter < max_counter )
    {
      if( f(c) > f(d) )
        b = d;
      else
        a = c;
      
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
    }
        
    return {(b + a) / 2, f((b + a) / 2)};
}


template < typename F >
auto max_f( F f, size_t min, size_t max )
  -> like<size_t>
{
    auto val = f(min);
    auto x = min;
    
    while( ++min < max )
    {
      auto tmp = f(min);
      if( tmp > val )
      {
        val = tmp;
        x = min;
      }
    }
    
    return {x, val};
    
}




