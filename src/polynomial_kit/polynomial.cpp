/**
 * Polynomial math library implementation.
 *
 * Giorgio Manca <giorgio.manca.97@gmail.com>
 *
 * June 30, 2023
 */

/**
 * Copyright 2024 dotX Automation s.r.l.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <polynomial_kit/polynomial.hpp>

namespace PolynomialKit
{
  template <typename T>
  Polynomial<T>::Polynomial() {
    zero_ = T();
    reserve(1);
    size_ = 1;
    poly_(0,0) = zero_;
  }

  template <typename T>
  Polynomial<T>::Polynomial(const T& value) {
    zero_ = T();
    reserve(1);
    size_ = 1;
    poly_(0,0) = value;
  }

  template <typename T>
  Polynomial<T>::Polynomial(const MatrixX<T>& matrix) {
    if(matrix.rows() > 2 && matrix.cols() > 2) {
      std::invalid_argument("Argument is not a vector.");
    }
    zero_ = T();
    if(matrix.rows() >= 1 && matrix.cols() >= 1) {
      if(matrix.rows() > 1) {
        reserve(matrix.rows());
        size_ = matrix.rows();
        poly_ = MatrixX<T>::Map(matrix.data(), 1, matrix.rows());
      } else {
        reserve(matrix.cols());
        size_ = matrix.cols();
        poly_ = MatrixX<T>::Map(matrix.data(), 1, matrix.cols());
      }
    }
  }

  template <typename T>
  Polynomial<T>::Polynomial(const Polynomial<T>& other) {
    zero_ = T();
    reserve(other.size());
    size_ = capacity();
    poly_ = other.poly_(0, seq(0, other.degree()));
  }

  template <typename T>
  Polynomial<T>::~Polynomial() {}

  template <typename T>
  inline unsigned int Polynomial<T>::degree() const {
    return size()-1;
  }

  template <typename T>
  inline unsigned int Polynomial<T>::size() const {
    return size_;
  }

  template <typename T>
  inline unsigned int Polynomial<T>::capacity() const {
    return poly_.cols();
  }

  template <typename T>
  void Polynomial<T>::reserve(unsigned int capacity, bool force) {
    if(capacity == 0) {
      std::invalid_argument("capacity cannot be zero");
    }
    if(force || capacity > this->capacity()) {
      poly_.conservativeResize(1, capacity);
      if(this->capacity() < this->size()) {
        size_ = this->capacity();
      }
    }
  }

  template <typename T>
  void Polynomial<T>::regrade(unsigned int degree) {
    if(degree > this->degree()) {
      reserve(degree+1);
    } else {
      this->size_ = degree + 1;
    }
  }

  template <typename T>
  void Polynomial<T>::reset(unsigned int degree, bool clean) {
    if(degree < this->size()){
      this->size_ = degree + 1;
    }
    if(degree == 0) {
      poly_(0, 0) = zero_;
    }
    if(clean) {
      this->clean();
    }
  }

  template <typename T>
  void Polynomial<T>::clean() {
    for(unsigned int i = this->size(); i < capacity(); i++) {
      poly_(0, i) = zero_;
    }
  }

  template <typename T>
  void Polynomial<T>::trim() {
    if(capacity() > this->size()) {
      poly_.conservativeResize(1, this->size());
    }
  }

  template <typename T>
  void Polynomial<T>::set(unsigned int degree, const T& coeff) {
    if(degree > this->degree() && coeff != zero_){
      reserve(degree + 1);
      for(unsigned int i = this->size(); i < degree; i++) {
        this->poly_(0,i) = zero_;
      }
      size_ = degree+1;
    } else if (degree == this->degree() && coeff == zero_) {
      size_--;
    }
    poly_(0, degree) = coeff;
  }

  template <typename T>
  T Polynomial<T>::get(unsigned int degree) const {
    if(degree >= this->size()) {
      std::invalid_argument("degree out of bound.");
    }
    return poly_(0, degree);
  }

  template <typename T>
  T Polynomial<T>::coeff(unsigned int degree) const {
    if(degree < this->size()){
      return poly_(0, degree);
    } else {
      return zero_;
    }
  }

  template <typename T>
  T Polynomial<T>::eval(const T& point) const {
    T x = point;
    T y = poly_(0, 0);
    for(unsigned int i = 1; i < this->size(); i++) {
      y += poly_(0, i) * x;
      x *= point;
    }
    return y;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::diff(unsigned int degree) {
    Polynomial res;
    if(degree == 0) {
      res = *this;
    } if(degree <= this->degree()) {
      res.reserve(this->size() - degree);
      res.size_ = res.capacity();
      res.poly_(0, seq(0, res.size()-1)) = this->poly_(0, seq(degree, this->degree()));
      for(unsigned int i = 0; i < res.size(); i++) {
        unsigned int j = i + degree;
        unsigned int m = j;
        while(true) {
          if(j == i) {
            break;
          } else {
            res.poly_(0, i) *= m;
            m = (--j);
          }
        }
      }
    }
    return res;
  }

  template <typename T>
  VectorX<T> Polynomial<T>::row_vector(unsigned int mindegree) {
    return row_matrix(mindegree);
  }

  template <typename T>
  VectorX<T> Polynomial<T>::col_vector(unsigned int mindegree) {
    return col_matrix(mindegree);
  }

  template <typename T>
  MatrixX<T> Polynomial<T>::row_matrix(unsigned int mindegree) {
    MatrixX<T> res = MatrixX<T>(1, std::max(this->size(), mindegree+1));
    res(0, seq(0, degree())) = poly_(0, seq(0, degree()));
    return res;
  }

  template <typename T>
  MatrixX<T> Polynomial<T>::col_matrix(unsigned int mindegree) {
    MatrixX<T> res = MatrixX<T>(std::max(this->size(), mindegree+1), 1);
    res(seq(0, degree()), 0) = poly_(0, seq(0, degree())).transpose();
    return res;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>& other) {
    this->reserve(other.size());
    this->size_ = other.size();
    this->poly_(0, seq(0, other.degree())) = other.poly_(0, seq(0, other.degree()));
    return *this;
  }

    template <typename T>
  Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& other) {
    this->reserve(other.size());
    if(other.size() > this->size()) {
      this->poly_(0, seq(0, this->degree())) += other.poly_(0, seq(0, this->degree()));
      this->poly_(0, seq(this->size(), other.degree())) = other.poly_(0, seq(this->size(), other.degree()));
      this->size_ = other.size();
    } else {
      this->poly_(0, seq(0, other.degree())) += other.poly_(0, seq(0, other.degree()));
    }
    while(this->size() > 1 && this->poly_(0, this->degree()) == zero_) {
      this->size_--;
    }
    return *this;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& other) {
    this->reserve(other.size());
    if(other.size() > this->size()) {
      this->poly_(0, seq(0, this->degree())) -= other.poly_(0, seq(0, this->degree()));
      this->poly_(0, seq(this->size(), other.degree())) = -other.poly_(0, seq(this->size(), other.degree()));
      this->size_ = other.size();
    } else {
      this->poly_(0, seq(0, other.degree())) -= other.poly_(0, seq(0, other.degree()));
    }
    while(this->size() > 1 && this->poly_(0, this->degree()) == zero_) {
      this->size_--;
    }
    return *this;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& other) {
    unsigned int size = this->degree() + other.degree() + 1;
    MatrixX<T> res = MatrixX<T>::Zero(1, size);
    for(unsigned int i = 0; i < this->size(); i++) {
      for(unsigned int j = 0; j < other.size(); j++) {
        res(0, i+j) += this->poly_(0, i) * other.poly_(0, j);
      }
    }
    this->reserve(size);
    this->size_ = size;
    this->poly_(0, seq(0, this->degree())) = res;
    return *this;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator^=(unsigned int p) {
    if(p == 0) {
      *this = Polynomial<T>(std::pow(zero_, zero_));
    } else if (p > 1) {
      Polynomial<T> square(*this);
      bool first = true;
      while(true) {
        if(p % 2 == 1) {
          std::cout << 1 << std::endl;
          if(first) {
            *this = square;
            first = false;
          } else {
            *this *= square;
          }
        } else {
          std::cout << 0 << std::endl;
        }
        p /= 2;
        if(p == 0) {
          break;
        }
        square *= square;
      }
    }
    return *this;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator<<=(unsigned int s) {
    if(s > 0) {
      this->reserve(this->size() + s);
      for(unsigned int i = this->degree()+s; i >= s; i--) {
        this->poly_(0, i) = this->poly_(0, i-s);
      }
      for(unsigned int i = 0; i < s; i++) {
        this->poly_(0, i) = zero_;
      }
      this->size_ += s;
    }
    return *this;
  }

  template <typename T>
  Polynomial<T>& Polynomial<T>::operator>>=(unsigned int s) {
    if(s > 0) {
      for(unsigned int i = s; i < this->size(); i++) {
        this->poly_(0, i-s) = this->poly_(0, i);
      }
      this->size_ -= s;
    }
    return *this;
  }

  template <typename T>
  bool Polynomial<T>::operator==(const Polynomial<T>& other) const {
    if(this->degree() == other.degree()) {
      return this->poly_(0, seq(0, this->degree())) == other.poly_(0, seq(0, other.degree()));
    } else {
      return false;
    }
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator+() const {
    Polynomial<T> res(*this);
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator-() const {
    Polynomial<T> res(*this);
    res.poly_ = -res.poly_;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator+(const Polynomial<T>& other) const {
    Polynomial<T> res(*this);
    res += other;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator-(const Polynomial<T>& other) const {
    Polynomial<T> res(*this);
    res -= other;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator*(const Polynomial<T>& other) const {
    Polynomial<T> res(*this);
    res *= other;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator^(unsigned int p) const {
    Polynomial<T> res(*this);
    res ^= p;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator<<(unsigned int s) const {
    Polynomial<T> res(*this);
    res <<= s;
    return res;
  }

  template <typename T>
  Polynomial<T> Polynomial<T>::operator>>(unsigned int s) const {
    Polynomial<T> res(*this);
    res >>= s;
    return res;
  }

  template <typename T>
  T Polynomial<T>::operator()(const T& point) const {
    return this->eval(point);
  }

  template class Polynomial<double>;
  template class Polynomial<float>;
  template class Polynomial<long>;
  template class Polynomial<int>;
  template class Polynomial<short>;

  template class Polynomial<std::complex<double>>;
  template class Polynomial<std::complex<float>>;
  template class Polynomial<std::complex<long>>;
  template class Polynomial<std::complex<int>>;
  template class Polynomial<std::complex<short>>;
}
