#ifndef VECT_H
#define VECT_H

// This file defines some templated mathematical functions on std::vectors.

// C++ includes
#include <iostream>
#include <vector>
#include <assert.h>

// Print newline-separated vector contents.
template<class T>
void print(const std::vector<T> & u)
{
  for (const auto & val : u)
    std::cout << val << std::endl;
}

// Computes the dot product x.y, returning the result as a scalar.
template<class T>
T dot(const std::vector<T> & x,
      const std::vector<T> & y)
{
  assert(x.size() == y.size());
  T dot_prod(0);
  for (unsigned int i=0; i<x.size(); ++i)
    dot_prod += x[i] * y[i];
  return dot_prod;
}

// Compute the l2-norm of vector r.
template<class T>
T norm(const std::vector<T> & r)
{
  return sqrt(dot(r,r));
}

// Scale a vector using operators, y = a*u
template<class T>
std::vector<T> operator * (const T & a,
                           const std::vector<T> & u)
{
  const unsigned int n = u.size();
  std::vector<T> ret(n);

  for (unsigned int i=0; i<n; ++i)
    ret[i] = a * u[i];

  return ret;
}

// Subtract in-place, u -= v
template<class T>
void operator -= (std::vector<T> & u,
                  const std::vector<T> & v)
{
  const unsigned int n = u.size();
  for (unsigned int i=0; i<n; ++i)
    u[i] -= v[i];
}

// Subtract in-place, u -= a * v
template<class T>
void subtract_scaled(std::vector<T> & u,
                     const T & a,
                     const std::vector<T> & v)
{
  const unsigned int n = u.size();
  for (unsigned int i=0; i<n; ++i)
    u[i] -= a * v[i];
}

// Subtract and return by value
template<class T>
std::vector<T> operator - (const std::vector<T> & u,
                           const std::vector<T> & v)
{
  const unsigned int n = u.size();
  std::vector<T> ret(n);

  for (unsigned int i=0; i<n; ++i)
    ret[i] = u[i] - v[i];

  return ret;
}

#endif
