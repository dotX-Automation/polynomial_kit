/**
 * Polynomial math library.
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

#ifndef POLYNOMIAL_KIT__POLYNOMIAL_HPP_
#define POLYNOMIAL_KIT__POLYNOMIAL_HPP_

#include "visibility_control.h"

#include <cmath>
#include <iostream>

#include <Eigen/Core>

using namespace Eigen;


namespace PolynomialKit
{
  template <typename T>
  class POLYNOMIAL_KIT_PUBLIC Polynomial {
    public:
      Polynomial();
      Polynomial(const T& value);
      Polynomial(const MatrixX<T>& matrix);
      Polynomial(const Polynomial<T>& matrix);
      ~Polynomial();

      inline unsigned int degree() const;
      inline unsigned int size() const;
      inline unsigned int capacity() const;
      void reserve(unsigned int capacity, bool force = false);
      void regrade(unsigned int degree);
      void reset(unsigned int degree = 0, bool clean = false);
      void clean();
      void trim();

      void set(unsigned int degree, const T& coeff);
      T get(unsigned int degree) const;
      T coeff(unsigned int degree) const;
      T eval(const T& value) const;

      Polynomial diff(unsigned int degree);

      VectorX<T> row_vector(unsigned int mindegree = 0);
      VectorX<T> col_vector(unsigned int mindegree = 0);
      MatrixX<T> row_matrix(unsigned int mindegree = 0);
      MatrixX<T> col_matrix(unsigned int mindegree = 0);

      Polynomial& operator=(const Polynomial& other);

      Polynomial& operator+=(const Polynomial& other);
      Polynomial& operator-=(const Polynomial& other);
      Polynomial& operator*=(const Polynomial& other);
      Polynomial& operator^=(unsigned int p);
      Polynomial& operator<<=(unsigned int s);
      Polynomial& operator>>=(unsigned int s);

      bool operator==(const Polynomial& other) const;

      Polynomial operator+() const;
      Polynomial operator-() const; 

      Polynomial operator+(const Polynomial& other) const;
      Polynomial operator-(const Polynomial& other) const;
      Polynomial operator*(const Polynomial& other) const;
      Polynomial operator^(unsigned int p) const;
      Polynomial operator<<(unsigned int s) const;
      Polynomial operator>>(unsigned int s) const;

      T operator()(const T& point) const;

      template <typename U>
      explicit operator Polynomial<U>() const {
        Polynomial<U> res;
        res.capacity_ = this->capacity_;
        res.size_ = this->size_;
        res.poly_ = this->poly_.template cast<U>();
        return res;
      }

      friend std::ostream &operator<<(std::ostream &output, const Polynomial<T> &poly) {
        unsigned int i = poly.degree();
        while(true) {
          output << poly.poly_(0, i);
          if(i > 0) {
            output << " ";
            i--;
          } else {
            break;
          }
        }
        return output;
      }

    private:
      T zero_ = T();
      MatrixX<T> poly_ = MatrixX<T>(1,1);
      unsigned int size_ = 1;
  };

  using Polynomiald = Polynomial<double>;
  using Polynomialf = Polynomial<float>;
  using Polynomiall = Polynomial<long>;
  using Polynomiali = Polynomial<int>;
  using Polynomials = Polynomial<short>;

  using Polynomialcd = Polynomial<std::complex<double>>;
  using Polynomialcf = Polynomial<std::complex<float>>;
  using Polynomialcl = Polynomial<std::complex<long>>;
  using Polynomialci = Polynomial<std::complex<int>>;
  using Polynomialcs = Polynomial<std::complex<short>>;
}

#endif  //POLYNOMIAL_KIT__POLYNOMIAL_HPP_
