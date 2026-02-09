#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>

template <class T>
struct Vec3 {
  T data[3];

  constexpr Vec3() : data{T{}, T{}, T{}} {}
  constexpr Vec3(T x, T y, T z) : data{x, y, z} {}

  template <class U>
  constexpr explicit Vec3(const Vec3<U>& other)
      : data{static_cast<T>(other.x()), static_cast<T>(other.y()), static_cast<T>(other.z())} {}

  constexpr T& x() { return data[0]; }
  constexpr const T& x() const { return data[0]; }
  constexpr T& y() { return data[1]; }
  constexpr const T& y() const { return data[1]; }
  constexpr T& z() { return data[2]; }
  constexpr const T& z() const { return data[2]; }

  constexpr T& operator[](std::size_t i) {
    assert(i < 3);
    return data[i];
  }

  constexpr const T& operator[](std::size_t i) const {
    assert(i < 3);
    return data[i];
  }

  constexpr void setZero() { data[0] = data[1] = data[2] = T{}; }

  constexpr Vec3 operator+() const { return *this; }
  constexpr Vec3 operator-() const { return Vec3(-data[0], -data[1], -data[2]); }

  constexpr Vec3& operator+=(const Vec3& rhs) {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    return *this;
  }

  constexpr Vec3& operator-=(const Vec3& rhs) {
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    return *this;
  }

  constexpr Vec3& operator*=(const T& scalar) {
    data[0] *= scalar;
    data[1] *= scalar;
    data[2] *= scalar;
    return *this;
  }

  constexpr Vec3& operator/=(const T& scalar) {
    data[0] /= scalar;
    data[1] /= scalar;
    data[2] /= scalar;
    return *this;
  }

  constexpr T dot(const Vec3& rhs) const {
    return data[0] * rhs.data[0] + data[1] * rhs.data[1] + data[2] * rhs.data[2];
  }

  constexpr Vec3 cross(const Vec3& rhs) const {
    return Vec3(data[1] * rhs.data[2] - data[2] * rhs.data[1],
                data[2] * rhs.data[0] - data[0] * rhs.data[2],
                data[0] * rhs.data[1] - data[1] * rhs.data[0]);
  }

  constexpr T squaredNorm() const { return dot(*this); }

  T norm() const { return static_cast<T>(std::sqrt(squaredNorm())); }
};

template <class T>
constexpr Vec3<T> operator+(Vec3<T> lhs, const Vec3<T>& rhs) {
  lhs += rhs;
  return lhs;
}

template <class T>
constexpr Vec3<T> operator-(Vec3<T> lhs, const Vec3<T>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <class T>
constexpr Vec3<T> operator*(Vec3<T> lhs, const T& scalar) {
  lhs *= scalar;
  return lhs;
}

template <class T>
constexpr Vec3<T> operator*(const T& scalar, Vec3<T> rhs) {
  rhs *= scalar;
  return rhs;
}

template <class T>
constexpr Vec3<T> operator/(Vec3<T> lhs, const T& scalar) {
  lhs /= scalar;
  return lhs;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
  os << v.x() << ' ' << v.y() << ' ' << v.z();
  return os;
}
