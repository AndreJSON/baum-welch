#ifndef MATRIX_H
#define MATRIX_H
#include "assert.h"
#include <iostream>

using std::size_t;
using std::cout;
using std::cin;
using std::endl;

template <class T>
class Matrix {
 private:
  // Row and columns, size = rows * cols
  size_t _rows, _cols, _size;
  // elements of matrix
  T* _mat;

  // check so that matrices has the same order
  bool same_order(const Matrix<T>&) const;
  void check_order(const Matrix<T>&) const;

  // check so indices are within the matrix
  void check_bounds(size_t, size_t) const;
  void check_bounds(size_t) const;

 public:
  // default 1x1 matrix
  Matrix();
  
  // rows x cols matrix
  Matrix(size_t rows, size_t cols);
  // copy constructor
  Matrix(const Matrix<T>&);
  ~Matrix();

  // get # rows and columns. length is for vectors
  size_t rows() const;
  size_t cols() const;
  size_t length() const;

  // matrix multiplication operators
  Matrix<T> operator*(const Matrix<T>&) const;

  // safe get functions (slow, using checkBounds)
  // will throw exceptions and debug info iff when out of bounds
  T& operator()(size_t, size_t);
  T& operator()(size_t);
  T operator()(size_t, size_t) const;
  T operator()(size_t) const;

  // assignment op
  Matrix<T>& operator=(const Matrix<T>&);

  /** 
  * print the matrix
  * rows cols (print row by row)
  * Example: if we have matrix
  * 1 2 3
  * 4 5 6
  * 7 8 9
  * output will be
  * 3 3 1 2 3 4 5 6 7 8 9
  */
  void print() const;
};

template <class T>
Matrix<T> readMatrix();

#include "matrix.hpp"

/**
  PRIVATE FUNCTIONS
*/
template <class T>
bool Matrix<T>::same_order(const Matrix<T>& other) const {
  return _rows == other._rows && _cols == other._cols;
}

template <class T>
void Matrix<T>::check_order(const Matrix<T>& other) const {
  assert(same_order(other));
}

template <class T>
void Matrix<T>::check_bounds(size_t r, size_t c) const {
  assert(r >= 0 && r < _rows && c >= 0 && c < _cols);
}

template <class T>
void Matrix<T>::check_bounds(size_t i) const {
  assert(i >= 0 && i < _size);
}

/**
  CONSTUCTORS
*/
template <class T>
Matrix<T>::Matrix()
    : _rows(1), _cols(1), _size(1) {
  _mat = new T[1]{0};
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols)
    : _rows(rows), _cols(cols), _size(rows * cols) {
  _mat = new T[_size]; //Could cause errors without {0}
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& other)
    : _rows(other._rows), _cols(other._cols), _size(other._rows * other._cols) {
  _mat = new T[_rows * _cols];
  for (int i = 0; i < _size; ++i) {
    _mat[i] = other._mat[i];
  }
}

/**
  DESTRUCTOR
*/
template <class T>
Matrix<T>::~Matrix() {
  delete[] _mat;
}

/**
  PUBLIC FUNCTIONS
*/

template <class T>
size_t Matrix<T>::length() const {
  return _size;
}

template <class T>
size_t Matrix<T>::rows() const {
  return _rows;
}

template <class T>
size_t Matrix<T>::cols() const {
  return _cols;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const {
  assert(_cols == rhs._rows);
  Matrix<T> res(_rows, rhs._cols);
  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < rhs._cols; ++j) {
      for (int k = 0; k < _cols; k++) {
        res(i, j) += operator()(i, k) * rhs(k, j);
      }
    }
  }
  return res;
}

template <class T>
T& Matrix<T>::operator()(size_t r, size_t c) {
  check_bounds(r, c);
  return _mat[r * _cols + c];
}

template <class T>
T& Matrix<T>::operator()(size_t i) {
  check_bounds(i);
  return _mat[i];
}

template <class T>
T Matrix<T>::operator()(size_t r, size_t c) const {
  check_bounds(r, c);
  return _mat[r * _cols + c];
}

template <class T>
T Matrix<T>::operator()(size_t i) const {
  check_bounds(i);
  return _mat[i];
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
  if (this == &other) return *this;
  
  _rows = other._rows;
  _cols = other._cols;
  if( _size != other._size) {
     _size = other._size;
     delete[] _mat;
     _mat = new T[_size];
  }

  for (int i = 0; i < _size; ++i) {
    _mat[i] = other._mat[i];
  }
  return *this;
}

template <class T>
void Matrix<T>::print() const {
  cout << _rows << " " << _cols;
  for (int i = 0; i < _size; ++i) {
    cout << " " << _mat[i];
  }
  cout << std::endl;
}

/**
  NON CLASS FUNCTIONS
*/

template <class T>
Matrix<T> readMatrix() {
  int rows, cols;
  cin >> rows >> cols;
  Matrix<T> res(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    cin >> res(i);
  }
  return res;
}

#endif
