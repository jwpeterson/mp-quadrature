#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <stdlib.h> // exit
#include <stdexcept>
#include <assert.h>

/**
 * Templated matrix class implementing operations needed for
 * solving quadrature rules.  Templated so that we can use
 * the 'mpfr_class' data type to get multi-precision matrices.
 */
template <class T>
class Matrix
{
public:
  /**
   * Constructor
   */
  Matrix(unsigned n_rows=0, unsigned n_cols=0);

  /**
   * Destructor.  Don't forget to make virtual if you ever derive from this.
   */
  ~Matrix();

  /**
   * Resize the Matrix so that it is n_rows * n_cols
   */
  void resize(unsigned n_rows, unsigned n_cols);

  /**
   * Clears (zeros) any existing values and sets the number of rows/cols to 0.
   */
  void clear();

  /**
   * Returns the number of rows (cols) in the Matrix
   */
  unsigned n_rows() const { return _n_rows; }
  unsigned n_cols() const { return _n_cols; }

  /**
   * Returns the \p (i,j) element of the matrix as a writeable or const reference.
   */
  T& operator() (unsigned i, unsigned j);
  const T& operator() (unsigned i, unsigned j) const;

  /**
   * Print the matrix, by default to \p std::cout.
   */
  void print(std::ostream& os = std::cout) const;

  /**
   * Same print function as above, but allows you to send the matrix
   * directly to the stream.  Note the strange "additional" template
   * declaration... this is necessary for friend functions for some
   * reason, otherwise you get a compiler warning that you're declaring
   * a non-template function.  Also note, this does have to be a 'friend'
   * function, otherwise you get a different 'must take exactly one argument'
   * error.
   */
  template <class U>
  friend std::ostream& operator<< (std::ostream& os, const Matrix<U>& m);

  /**
   * Replaces this matrix by the pair of triangular matrices L and U,
   * such that A=LU, and solves the linear system Ax=b for x.  Note: this
   * method overwrites the matrix, so if you need to use it again for
   * some reason after lu_solve, be sure to make a copy of it first!
   */
  void lu_solve(std::vector<T>& x, const std::vector<T>& b);

  /**
   * Only for symmetric positive definite matrices. Factors the
   * matrix into LL^T and then does back substitution. This overwrites
   * the matrix.
   */
  void cholesky_solve(std::vector<T>& x, const std::vector<T>& b);

  /**
   * Matrix-vector product
   * Returns y = (*this) * x
   */
  std::vector<T> operator * (const std::vector<T> & x);

  /**
   * Transposed Matrix-vector product
   * Returns y = (*this)^T * x
   */
  std::vector<T> matvec_transpose (const std::vector<T> & x);

  /**
   * Compute a measure of the matrix asymmetry given by:
   * |a(i,j) - a(j,i)|
   *  Useful for debugging.
   */
  T asymmetry() const;

private:
  /**
   * Ingredients of the lu_solve() function.
   */
  void _lu_decompose();
  void _lu_back_substitute(std::vector<T>& x, const std::vector<T>& b) const;

  /**
   * Ingredients of the cholesky_solve() function.
   */
  void _cholesky_decompose();
  void _cholesky_back_substitute(std::vector<T>& x, const std::vector<T>& b) const;

  /**
   * Private data members.
   */
  unsigned _n_rows;
  unsigned _n_cols;
  std::vector<T> _val;
  std::vector<unsigned> _pivots;
};





template <class T>
Matrix<T>::Matrix(unsigned n_rows, unsigned n_cols)
  : _n_rows(n_rows),
    _n_cols(n_cols)
{
  _val.resize(_n_rows * _n_cols);
}



template <class T>
Matrix<T>::~Matrix() = default;



template<class T>
void Matrix<T>::resize(unsigned n_rows, unsigned n_cols)
{
  _n_rows = n_rows;
  _n_cols = n_cols;
  _val.resize(_n_rows * _n_cols);
}



template<class T>
void Matrix<T>::clear()
{
  _n_rows = 0;
  _n_cols = 0;
  _val.clear();
}



template <class T>
T & Matrix<T>::operator() (unsigned i, unsigned j)
{
  return _val[_n_cols*i + j];
}



template <class T>
const T & Matrix<T>::operator() (unsigned i, unsigned j) const
{
  return _val[_n_cols*i + j];
}




template <class T>
void Matrix<T>::print(std::ostream& os) const
{
  // A convenient reference to *this
  const Matrix<T>& A = *this;

  // save the initial format flags
  std::ios_base::fmtflags os_flags = os.flags();

  // Print out the entries
  for (unsigned i=0; i<_n_rows; i++)
    {
      // There seems to be something wrong with writing mpfr_class
      // values and std::setw?  It puts a weird gap between the digits
      // and the exponent...  It also seems to ignore std::fixed and
      // just does a dynamic mix of scientific and fixed-point formatting.
      // At least it does seem to follow the precision flag.
      for (unsigned j=0; j<_n_cols; j++)
        {
          // Don't print values smaller than epsilon=1.e-30
          T val = A(i,j);

          if (abs(val) < 1.e-30)
            val = 0.;

          os << std::setprecision(6)
             << std::showpos
             << val << " ";
        }
      os << std::endl;
    }

  // restore the original format flags
  os.flags(os_flags);
}



// This is the friend function, but note that we don't use 'friend' here
template <class U>
std::ostream& operator<< (std::ostream& os, const Matrix<U>& m)
{
 m.print(os);
 return os;
}






template <class T>
void Matrix<T>::lu_solve(std::vector<T>& x, const std::vector<T>& b)
{
  if (_n_rows != _n_cols)
    {
      std::cerr << "Error, cannot LU solve non-square matrix!" << std::endl;
      exit(1);
    }

  this->_lu_decompose();
  this->_lu_back_substitute(x, b);
}



template <class T>
void Matrix<T>::cholesky_solve(std::vector<T>& x, const std::vector<T>& b)
{
  if (_n_rows != _n_cols)
    {
      std::cerr << "Error, cannot Cholesky solve non-square matrix!" << std::endl;
      exit(1);
    }

  this->_cholesky_decompose();
  this->_cholesky_back_substitute(x, b);
}




template <class T>
void Matrix<T>::_lu_decompose()
{
  // For doing partial pivoting.  We also need these values when doing
  // back substitution!
  _pivots.resize(_n_rows);

  // A convenient reference to *this
  Matrix<T>& A = *this;

  for (unsigned i=0; i<_n_rows; ++i)
    {
      // Find the pivot row by searching down the i'th column
      _pivots[i] = i;

      // Note: we can't use std::abs here since it's not defined for
      // the mpfr_class type...  We also can't catch the return as a
      // "Real" since there is no automatic conversion to that type,
      // and it would lose precision.
      T max = abs( A(i,i) );
      for (unsigned j=i+1; j<_n_rows; ++j)
        {
          T candidate_max = abs( A(j,i) );
          if (max < candidate_max)
            {
              max = candidate_max;
              _pivots[i] = j;
            }
        } // end for(j)

      // If the max was found in a different row, interchange rows.
      // Here we interchange the *entire* row, in Gaussian elimination
      // you would only interchange the subrows A(i,j) and A(p(i),j), for j>i
      if (_pivots[i] != i)
        {
          for (unsigned j=0; j<_n_rows; ++j)
            std::swap( A(i,j), A(_pivots[i], j) );
        }

      // If the max abs entry found is zero, the matrix is singular
      if (A(i,i) == T(0))
        throw std::runtime_error("Matrix is singular!");

      // Scale upper triangle entries of row i by the diagonal entry
      // Note: don't scale the diagonal entry itself!
      const T diag_inv = T(1) / A(i,i);
      for (unsigned j=i+1; j<_n_rows; ++j)
        A(i,j) *= diag_inv;

      // Update the remaining sub-matrix A[i+1:_n_rows][i+1:_n_rows]
      // by subtracting off (the diagonal-scaled)
      // upper-triangular part of row i, scaled by the
      // i'th column entry of each row.  In terms of
      // row operations, this is:
      // for each r > i
      //   SubRow(r) = SubRow(r) - A(r,i)*SubRow(i)
      //
      // If we were scaling the i'th column as well, like
      // in Gaussian elimination, this would 'zero' the
      // entry in the i'th column.
      for (unsigned row=i+1; row<_n_rows; ++row)
        for (unsigned col=i+1; col<_n_rows; ++col)
          A(row,col) -= A(row,i) * A(i,col);

    } // end for(i)
}




template <class T>
void Matrix<T>::_lu_back_substitute(std::vector<T>& x, const std::vector<T>& b) const
{
  x.clear();
  x.resize (_n_rows);

  // A convenient reference to *this
  const Matrix<T>& A = *this;

  // Temporary vector storage.  We use this instead of
  // modifying the RHS.
  std::vector<T> z = b;

  // Lower-triangular "top to bottom" solve step, taking into account pivots
  for (unsigned i=0; i<_n_rows; ++i)
    {
      // Swap
      if (_pivots[i] != i)
        std::swap( z[i], z[_pivots[i]] );

      x[i] = z[i];

      for (unsigned j=0; j<i; ++j)
        x[i] -= A(i,j)*x[j];

      x[i] /= A(i,i);
    }

  // Upper-triangular "bottom to top" solve step.
  // Note integer looping, since the loop index goes negative here.
  const unsigned int last_row = _n_rows-1;

  for (int i=last_row; i>=0; --i)
    {
      for (int j=i+1; j<static_cast<int>(_n_rows); ++j)
        x[i] -= A(i,j)*x[j];
    }
}



template <class T>
void Matrix<T>::_cholesky_decompose()
{
  // A convenient reference to *this
  Matrix<T> & A = *this;

  for (unsigned int i=0; i<_n_rows; ++i)
    for (unsigned int j=i; j<_n_cols; ++j)
      {
        for (unsigned int k=0; k<i; ++k)
          A(i,j) -= A(i,k) * A(j,k);

        if (i == j)
          {
            if (A(i,i) <= T(0))
              throw std::runtime_error("Matrix is not SPD!");

            A(i,i) = sqrt(A(i,j));
          }
        else
          A(j,i) = A(i,j) / A(i,i);
      }
}



template <class T>
void Matrix<T>::_cholesky_back_substitute(std::vector<T>& x, const std::vector<T>& b) const
{
  // A convenient reference to *this
  const Matrix<T> & A = *this;

  // Now compute the solution to Ax =b using the factorization.
  x.resize(_n_rows);

  // Solve for Ly=b
  for (unsigned int i=0; i<_n_cols; ++i)
    {
      T temp = b[i];

      for (unsigned int k=0; k<i; ++k)
        temp -= A(i,k) * x[k];

      x[i] = temp / A(i,i);
    }

  // Solve for L^T x = y
  for (unsigned int i=0; i<_n_cols; ++i)
    {
      const unsigned int ib = (_n_cols - 1) - i;

      for (unsigned int k=(ib+1); k<_n_cols; ++k)
        x[ib] -= A(k,ib) * x[k];

      x[ib] /= A(ib,ib);
    }
}



template <class T>
std::vector<T>
Matrix<T>::operator * (const std::vector<T> & x)
{
  std::vector<T> y(_n_rows);
  for (unsigned int i=0; i<_n_rows; ++i)
    for (unsigned int j=0; j<_n_cols; ++j)
      y[i] += (*this)(i,j) * x[j];

  return y;
}



template <class T>
std::vector<T> Matrix<T>::matvec_transpose (const std::vector<T> & x)
{
  std::vector<T> y(_n_cols);
  for (unsigned int i=0; i<_n_cols; ++i)
    for (unsigned int j=0; j<_n_rows; ++j)
      y[i] += (*this)(j,i) * x[j];

  return y;
}



template <class T>
T Matrix<T>::asymmetry() const
{
  // Only for square matrices
  assert(_n_rows == _n_cols);

  T asymmetry(0);

  for (unsigned int i=0; i<_n_rows; ++i)
    for (unsigned int j=i+1; j<_n_cols; ++j)
      asymmetry += abs((*this)(i,j) - (*this)(j,i));

  return asymmetry;
}

#endif // __matrix_h__
