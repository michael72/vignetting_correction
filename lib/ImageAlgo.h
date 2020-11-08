#pragma once

#include <algorithm>
#include <ciso646>
#include <cmath>
#include <type_traits>

#include "ImageAlgoBase.h"

namespace imgalg {

constexpr static PixelT clamp(int const i) {
  return static_cast<PixelT>(((0xff - i) >> 31) | (i & ~(i >> 31)));
}

constexpr static PixelT clamp(Real const f) { return clamp(iround(f)); }

class ImageAlgo : public ImageAlgoBase {
public:
  using HistogramType = Real;

  template <typename T> static constexpr T square(T const x) { return x * x; }

  template <typename T> static Real sqrt(T sq) {
    return sqrtf(static_cast<Real>(sq));
  }

  static Real dist(Point const &p) { return sqrt(square(p.x) + square(p.y)); }

  template <typename T, typename View, typename Fun>
  static bool reduce(T &akku, View &view,
                     std::function<void(int)> const &row_fun, Fun const &fun) {
    auto const [cols, rows] = dimensions(view);

    using PixelIt = decltype(row_begin<Pixel>(view, 0));
    static auto constexpr bool_result =
        std::is_same<std::invoke_result_t<Fun, T &, PixelIt &, int>,
                     bool>::value;

    for (int row = 0; row < rows; ++row) {
      auto row_it = row_begin<Pixel>(view, row);
      row_fun(row);
      for (int col = 0; col < cols; ++col) {
        if constexpr (bool_result) {
          if (not fun(akku, row_it, col)) {
            return false;
          }
        } else {
          fun(akku, row_it, col);
        }
      }
    }

    return true;
  }
};
} // namespace imgalg
