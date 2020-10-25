#pragma once

namespace imgalg {
using Real = float;

constexpr static int iround(Real const f) {
  // only works for positive f
  return static_cast<int>(f + 0.5f);
}

constexpr static int iround(int const i) { return i; }

template <typename T> static constexpr int div_round(T const a, T const b) {
  return static_cast<int>((a + b / 2) / b);
}
} // namespace imgalg

#ifdef USE_OPENCV
#include "ImageAlgo_opencv.h"
#else
#include "ImageAlgo_gil.h"
#endif