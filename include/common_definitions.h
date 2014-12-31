#ifndef __common_definitions__
#define __common_definitions__

#include <string>
#include <vector>
#include <algorithm> // std::partition

// For mpfr_class.  It may be possible to forward declare this, but
// I'm not sure how.
#include "gmpfrxx.h"

typedef double Real;

// Returns a string which contains the 32 decimal digits of x in the form:
// -9.4489927222288222340758013830322e-01L
// The user can optionally turn off the appending of a trailing L by passing
// 'false' for the second argument.
std::string fix_string(const mpfr_class& x,
                       bool append_L = true);

// iota fills up a range with the entries (value, value+1, value+2, ...)
// and is a handy way to initialize an array.
template <typename ForwardIter, typename T>
void iota (ForwardIter first, ForwardIter last, T value)
{
  while (first != last)
    {
      *first = value++;
      ++first;
    }
}


/**
 * Simple (non-recursive) factrial function implementation.  Using big
 * ints so we hopefully don't have to worry about integer overflow for
 * large n...
 */
inline
mpz_class factorial(mpz_class n)
{
  mpz_class factorial_n = 1;

  if (n==0)
    return factorial_n;

  for (unsigned int i=1; i<n; i++)
    factorial_n *= i+1;

  return factorial_n;
}



/**
 * A simple class representing an arbitrary-precision point in
 * three-dimensional space.  Templated on the underlying type, which
 * in this application is typically either mpfr_class or mpq_class.
 */
template <class T>
class Point
{
public:
  Point(const T & x = 0,
        const T & y = 0,
        const T & z = 0)
  {
    _coords[0] = x;
    _coords[1] = y;
    _coords[2] = z;
  }

  const T & operator()(unsigned i) const
  {
    return _coords[i];
  }

  T & operator()(unsigned i)
  {
    return _coords[i];
  }

private:
  T _coords[3];
};



/**
 * An indirect unary comparison object to be used with the
 * indirect_partition function.  An obvious generalization would be to
 * also template on a custom unary comparison object.
 */
template <class RandomAccessIterator>
struct indirect_unary_is_negative
{
  // ctor
  indirect_unary_is_negative(RandomAccessIterator r) :
    _iterator(r)
  {}

  // unary comparison operator
  bool operator()(size_t lhs)
  {
    return _iterator[lhs] < 0;
  }

private:
  RandomAccessIterator _iterator;
};



/**
 * A binary less-than comparison object.  An obvious generalization
 * would be to also template on a custom comparison object while doing
 * the indirect comparison.
 */
template <class RandomAccessIterator>
struct indirect_less_than
{
  // ctor
  indirect_less_than(RandomAccessIterator r) :
    _iterator(r)
  {}

  // comparison operator - calls the user's comparison function on
  // v[lhs] and v[rhs]
  bool operator()(size_t lhs, size_t rhs)
    {
      // Note: operator[] is defined for random access iterators!
      return _iterator[lhs] < _iterator[rhs];
    }

private:
  RandomAccessIterator _iterator;
};



/**
 * A binary less-than comparison object which compares absolute values instead of values.
 */
template <class RandomAccessIterator>
struct indirect_abs_less_than
{
  // ctor
  indirect_abs_less_than(RandomAccessIterator r) :
    _iterator(r)
  {}

  // comparison operator - calls the user's comparison function on
  // v[lhs] and v[rhs]
  bool operator()(size_t lhs, size_t rhs)
    {
      // Note: operator[] is defined for random access iterators!
      return abs(_iterator[lhs]) < abs(_iterator[rhs]);
    }

private:
  RandomAccessIterator _iterator;
};



/**
 * A binary greater-than comparison object which compares absolute values.
 */
template <class RandomAccessIterator>
struct indirect_abs_greater_than
{
  // ctor
  indirect_abs_greater_than(RandomAccessIterator r) :
    _iterator(r)
  {}

  // comparison operator - calls the user's comparison function on
  // v[lhs] and v[rhs]
  bool operator()(size_t lhs, size_t rhs)
    {
      // Note: operator[] is defined for random access iterators!
      return abs(_iterator[lhs]) > abs(_iterator[rhs]);
    }

private:
  RandomAccessIterator _iterator;
};



/**
 * Calls std::partition() using an indirect comparison operator so
 * that we can partition multiple vectors at the same time.
 */
template <class RandomAccessIterator>
void indirect_partition(RandomAccessIterator beg,
                        RandomAccessIterator end,
                        std::vector<size_t>& indices)
{
  // Initialize the indices array
  indices.resize(std::distance(beg, end));
  iota(indices.begin(), indices.end(), 0);

  // Construct indirect comparison object
  indirect_unary_is_negative<RandomAccessIterator> comp(beg);

  // Call std::partition on the indices with the indirect comparison object
  std::partition(indices.begin(), indices.end(), comp);
}



/**
 * Similar to std::partition but calls std::sort indirectly.
 */
template <class RandomAccessIterator>
void indirect_sort(RandomAccessIterator beg,
                   RandomAccessIterator end,
                   std::vector<size_t>& indices)
{
  // Initialize the indices array
  indices.resize(std::distance(beg, end));
  iota(indices.begin(), indices.end(), 0);

  // Construct comparator object.  There are a couple available to
  // choose from...
  indirect_less_than<RandomAccessIterator> comp(beg);
  // indirect_abs_less_than<RandomAccessIterator> comp(beg);
  // indirect_abs_greater_than<RandomAccessIterator> comp(beg);

  // Sort the indices, based on the data
  std::sort(indices.begin(), indices.end(), comp);
}



// Interleave the entries of a vector (a,b,c,d,e) such that the result is (a,e,b,d,c)
template <typename T>
void interleave(std::vector<T> & v)
{
  std::vector<T> result(v.size());

  unsigned n = v.size();
  for (unsigned i=0, c=0; i<n; i+=2, ++c)
    {
      result[i] = v[c];
      if (i+1 < n)
        result[i+1] = v[n-1-c];
    }

  // Copy result into the input vector
  v = result;
}

#endif // __common_definitions__
