#include <sstream>
#include <stdexcept>

#include "exceptions/matrix_exceptions.h"
#include "matrix.h"

#include "utils/comparison.h"
#include "utils/matrix_funcs.h"

using namespace supmath;

template<Numerical T>
supmath::Matrix<T>::Matrix() : _row_size(0), _col_size(0)
{}

template<Numerical T>
supmath::Matrix<T>::Matrix(const size_t row_size, const size_t col_size) :
    _row_size(row_size), _col_size(col_size), _elements(row_size, std::vector<T>(col_size, 0))
{}

template<Numerical T>
supmath::Matrix<T>::Matrix(const std::vector<std::vector<T>> &elements)
{
  size_t row_size = static_cast<size_t>(elements.size());
  if ( elements.size() <= 0 ) {
    throw InvalidMatrixDimensionException("The row of the Matrix should not be empty.");
  }

  size_t col_size = static_cast<size_t>(elements[0].size());
  for ( const auto &row : elements ) {
    if ( row.size() != col_size ) {
      throw InvalidMatrixDimensionException("All rows must have the same number of columns.");
    }
  }

  _col_size = col_size;
  _row_size = row_size;
  _elements = elements;
}

template<Numerical T>
template<Numerical U>
supmath::Matrix<T>::Matrix(const Matrix<U> &other) :
    _row_size(other.getRowSize()), _col_size(other.getColSize()),
    _elements(other.getRowSize(), std::vector<T>(other.getColSize()))
{

  for ( size_t row = 0; row < other.getRowSize(); row++ ) {
    for ( size_t col = 0; col < other.getColSize(); col++ ) {
      _elements[row][col] = static_cast<T>(other[row][col]);
    }
  }
}

template<Numerical T>
std::vector<T> &supmath::Matrix<T>::operator[](const size_t index)
{
  if ( index >= _row_size ) {
    throw std::out_of_range(
        "The matrix has " + std::to_string(_row_size) + " rows, but tried to access row at index "
        + std::to_string(index)
    );
  }
  return _elements[index];
}

template<Numerical T>
const std::vector<T> &supmath::Matrix<T>::operator[](const size_t index) const
{
  if ( index >= _row_size ) {
    throw std::out_of_range(
        "The matrix has " + std::to_string(_row_size) + " rows, but tried to access row at index "
        + std::to_string(index)
    );
  }
  return _elements[index];
}

template<Numerical T>
template<Numerical U>
auto supmath::Matrix<T>::operator+(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>
{
  if ( _row_size != matrix.getRowSize() || _col_size != matrix.getColSize() ) {
    throw MatrixMismatchException(
        "A " + std::to_string(_row_size) + "x" + std::to_string(_col_size) + " matrix is trying to add with a "
        + std::to_string(matrix.getRowSize()) + "x" + std::to_string(matrix.getColSize()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(_row_size, _col_size);

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < _col_size; col++ ) {
      result[row][col] = _elements[row][col] + matrix[row][col];
    }
  }

  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Matrix<T>::operator-(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>
{
  if ( _row_size != matrix.getRowSize() || _col_size != matrix.getColSize() ) {
    throw MatrixMismatchException(
        "A " + std::to_string(_row_size) + "x" + std::to_string(_col_size) + " matrix is trying to subtract with a "
        + std::to_string(matrix.getRowSize()) + "x" + std::to_string(matrix.getColSize()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(_row_size, _col_size);

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < _col_size; col++ ) {
      result[row][col] = _elements[row][col] - matrix[row][col];
    }
  }

  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Matrix<T>::operator*(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>
{
  if ( _col_size != matrix.getRowSize() ) {
    throw MatrixMismatchException(
        "The first Matrix has a column size of " + std::to_string(_col_size)
        + ", and the second Matrix has a row size of " + std::to_string(matrix.getRowSize())
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(_row_size, matrix.getColSize());

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < matrix.getColSize(); col++ ) {
      for ( size_t factor_row = 0; factor_row < matrix.getRowSize(); factor_row++ ) {
        result[row][col] += _elements[row][factor_row] * matrix[factor_row][col];
      }
    }
  }

  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Matrix<T>::operator*(const U factor) const -> Matrix<std::common_type_t<T, U>>
{
  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(_row_size, _col_size);

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < _col_size; col++ ) {
      result[row][col] = factor * _elements[row][col];
    }
  }
  int    a = 24;
  double b = 3;
  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Matrix<T>::elementWiseProduct(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>
{
  if ( _row_size != matrix.getRowSize() || _col_size != matrix.getColSize() ) {
    throw MatrixMismatchException(
        "A " + std::to_string(_row_size) + "x" + std::to_string(_col_size)
        + " matrix is trying to element-wise product with a " + std::to_string(matrix.getRowSize()) + "x"
        + std::to_string(matrix.getColSize()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(_row_size, _col_size);

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < _col_size; col++ ) {
      result[row][col] = _elements[row][col] * matrix[row][col];
    }
  }

  return result;
}

template<Numerical T>
Matrix<T> supmath::Matrix<T>::transpose() const
{
  Matrix<T> transposed(_col_size, _row_size);

  for ( size_t row = 0; row < _col_size; row++ ) {
    for ( size_t col = 0; col < _row_size; col++ ) {
      transposed[row][col] = _elements[col][row];
    }
  }

  return transposed;
}

template<Numerical T>
double supmath::Matrix<T>::determinant() const
{
  if ( _row_size != _col_size ) {
    throw MatrixNonSquareException("The matrix is not a square matrix.");
  }

  if ( _row_size == 1 ) {
    return static_cast<double>(_elements[0][0]);
  }

  if ( _row_size == 2 ) {
    return static_cast<double>(_elements[0][0] * _elements[1][1] - _elements[0][1] * _elements[1][0]);
  }

  if ( _row_size == 3 ) {
    return static_cast<double>(
        _elements[0][0] * (_elements[1][1] * _elements[2][2] - _elements[1][2] * _elements[2][1])
        - _elements[0][1] * (_elements[1][0] * _elements[2][2] - _elements[1][2] * _elements[2][0])
        + _elements[0][2] * (_elements[1][0] * _elements[2][1] - _elements[1][1] * _elements[2][0])
    );
  }

  try {
    auto [eliminated_matrix, nu_swaps] = matrix_funcs::gaussianElimination(*this);

    double det                         = 1;
    for ( size_t pivot = 0; pivot < _row_size; pivot++ ) {
      det *= eliminated_matrix[pivot][pivot];
    }

    return std::pow(-1, nu_swaps) * det;
  }
  catch ( MatrixSingularException & ) {
    return 0;
  }
}

template<Numerical T>
Matrix<double> supmath::Matrix<T>::inverse() const
{
  if ( _row_size != _col_size ) {
    throw MatrixNonSquareException("The matrix is not a square matrix.");
  }

  const double det = determinant();

  if ( det == 0 ) {
    throw MatrixNonInvertibleException("The matrix is singular and hence is not invertible.");
  }

  Matrix<double> inversed(_row_size, _col_size);

  if ( _row_size == 1 ) {
    inversed[0][0] = 1.0 / static_cast<double>(_elements[0][0]);
    return inversed;
  }

  if ( _row_size == 2 ) {
    inversed[0][0] = _elements[1][1] / det;
    inversed[0][1] = -_elements[0][1] / det;
    inversed[1][0] = -_elements[1][0] / det;
    inversed[1][1] = _elements[0][0] / det;
    return inversed;
  }

  return matrix_funcs::gaussianEliminationInverse(*this);
}

template<Numerical T>
void supmath::Matrix<T>::swap(const size_t chosen_index, const size_t swapped_index, MatrixOrder order)
{
  if ( order == MatrixOrder::Row ) {
    if ( chosen_index >= _row_size || swapped_index >= _row_size ) {
      throw std::out_of_range(
          "The matrix has " + std::to_string(_row_size) + " rows, but tried to swap row at index "
          + std::to_string(chosen_index) + " and " + std::to_string(swapped_index)
      );
    }
    _elements[chosen_index].swap(_elements[swapped_index]);
  }
  else {
    if ( chosen_index >= _col_size || swapped_index >= _col_size ) {
      throw std::out_of_range(
          "The matrix has " + std::to_string(_col_size) + " columns, but tried to swap column at index "
          + std::to_string(chosen_index) + " and " + std::to_string(swapped_index)
      );
    }

    for ( size_t row = 0; row < _row_size; row++ ) {
      std::swap(_elements[row][chosen_index], _elements[row][swapped_index]);
    }
  }
}

template<Numerical T>
template<Numerical U>
bool supmath::Matrix<T>::operator==(const Matrix<U> &matrix) const
{
  if ( _row_size != matrix.getRowSize() || _col_size != matrix.getColSize() ) {
    return false;
  }

  double epsilon = std::is_same<T, float>::value ? 1e-5 : 1e-9;

  for ( size_t row = 0; row < _row_size; row++ ) {
    for ( size_t col = 0; col < _col_size; col++ ) {
      if ( !comparison::numeric_almost_equal(matrix[row][col], _elements[row][col], epsilon) ) {
        return false;
      }
    }
  }

  return true;
}

template<Numerical T>
supmath::Matrix<T>::operator std::string() const
{
  std::ostringstream oss;

  oss << "Matrix([";
  for ( size_t row = 0; row < getRowSize(); ++row ) {
    oss << "[";
    for ( size_t col = 0; col < getColSize(); ++col ) {
      oss << _elements[row][col];
      if ( col != getColSize() - 1 ) {
        oss << ", ";
      }
    }
    oss << (row == getRowSize() - 1 ? "]" : "], ");
  }
  oss << "])";

  return oss.str();
  ;
}

template<Numerical T, Numerical U>
auto supmath::operator*(const T factor, const Matrix<U> &matrix) -> Matrix<std::common_type_t<T, U>>
{
  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(matrix.getRowSize(), matrix.getColSize());

  for ( size_t row = 0; row < matrix.getRowSize(); row++ ) {
    for ( size_t col = 0; col < matrix.getColSize(); col++ ) {
      result[row][col] = factor * matrix[row][col];
    }
  }

  return result;
}

template<Numerical T>
std::ostream &supmath::operator<<(std::ostream &os, const Matrix<T> &matrix)
{
  os << static_cast<std::string>(matrix);

  return os;
}

// Explicit Instantiation
template class supmath::Matrix<short>;
template class supmath::Matrix<int>;
template class supmath::Matrix<long>;
template class supmath::Matrix<long long>;
template class supmath::Matrix<float>;
template class supmath::Matrix<double>;

template std::ostream &supmath::operator<<(std::ostream &, const Matrix<short> &);
template std::ostream &supmath::operator<<(std::ostream &, const Matrix<int> &);
template std::ostream &supmath::operator<<(std::ostream &, const Matrix<long> &);
template std::ostream &supmath::operator<<(std::ostream &, const Matrix<long long> &);
template std::ostream &supmath::operator<<(std::ostream &, const Matrix<float> &);
template std::ostream &supmath::operator<<(std::ostream &, const Matrix<double> &);

INSTANTIATION_MATRIX_CONSTRUCTOR(short, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, float)

INSTANTIATION_MATRIX_FUNCTIONS(short, short)
INSTANTIATION_MATRIX_FUNCTIONS(short, int)
INSTANTIATION_MATRIX_FUNCTIONS(short, long)
INSTANTIATION_MATRIX_FUNCTIONS(short, long long)
INSTANTIATION_MATRIX_FUNCTIONS(short, float)
INSTANTIATION_MATRIX_FUNCTIONS(short, double)
INSTANTIATION_MATRIX_FUNCTIONS(int, short)
INSTANTIATION_MATRIX_FUNCTIONS(int, int)
INSTANTIATION_MATRIX_FUNCTIONS(int, long)
INSTANTIATION_MATRIX_FUNCTIONS(int, long long)
INSTANTIATION_MATRIX_FUNCTIONS(int, float)
INSTANTIATION_MATRIX_FUNCTIONS(int, double)
INSTANTIATION_MATRIX_FUNCTIONS(long, short)
INSTANTIATION_MATRIX_FUNCTIONS(long, int)
INSTANTIATION_MATRIX_FUNCTIONS(long, long)
INSTANTIATION_MATRIX_FUNCTIONS(long, long long)
INSTANTIATION_MATRIX_FUNCTIONS(long, float)
INSTANTIATION_MATRIX_FUNCTIONS(long, double)
INSTANTIATION_MATRIX_FUNCTIONS(long long, short)
INSTANTIATION_MATRIX_FUNCTIONS(long long, int)
INSTANTIATION_MATRIX_FUNCTIONS(long long, long)
INSTANTIATION_MATRIX_FUNCTIONS(long long, long long)
INSTANTIATION_MATRIX_FUNCTIONS(long long, float)
INSTANTIATION_MATRIX_FUNCTIONS(long long, double)
INSTANTIATION_MATRIX_FUNCTIONS(float, short)
INSTANTIATION_MATRIX_FUNCTIONS(float, int)
INSTANTIATION_MATRIX_FUNCTIONS(float, long)
INSTANTIATION_MATRIX_FUNCTIONS(float, long long)
INSTANTIATION_MATRIX_FUNCTIONS(float, float)
INSTANTIATION_MATRIX_FUNCTIONS(float, double)
INSTANTIATION_MATRIX_FUNCTIONS(double, short)
INSTANTIATION_MATRIX_FUNCTIONS(double, int)
INSTANTIATION_MATRIX_FUNCTIONS(double, long)
INSTANTIATION_MATRIX_FUNCTIONS(double, long long)
INSTANTIATION_MATRIX_FUNCTIONS(double, float)
INSTANTIATION_MATRIX_FUNCTIONS(double, double)
