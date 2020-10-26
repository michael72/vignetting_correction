#include "VignettingCorrection.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <ciso646>
#include <sys/stat.h>
#include <iostream>

#include "GaussianBlur.h"

namespace vgncorr {

using namespace imgalg;

static auto constexpr DebugPrint = false;
static auto constexpr MeasureTime = true;

constexpr float log2f(float val) {
  union {
    float val;
    int32_t x;
  } u = {val};
  auto const log_2i = (((u.x >> 23) & 255) - 128);
  u.x = (u.x & ~(255 << 23)) + (127 << 23);
  return log_2i + ((-0.34484843f) * u.val + 2.02466578f) * u.val - 0.67487759f;
}

constexpr int log2i(int const val) {
  if (val == 1) {
    return 0;
  }
  return 1 + log2i(val / 2);
}

/// Polynomial to calculate the vignetting correction with.
/// This is of the form:
/// g_a,b,c(r) = 1 + ar^2 + br^4 + cr^6.
class Poly : public imgalg::ImageAlgo {
 public:
  static auto constexpr num_coefficients = 3;
  using CoeffType = Real;
  using Coefficients = std::array<CoeffType, num_coefficients>;

  Poly(Coefficients const &coefficients, Point const &mid_point)
      : coeffs_(coefficients),
        mid_point_(mid_point),
        d2_(square(dist(mid_point))) {}

  void set_row(int const row) const;
  /// Check if the polynomial is increasing. I.e. the derivative is positive for
  /// the current coefficients.
  bool is_increasing() const;
  Real calc_at(int const col) const;

 private:
  Coefficients const coeffs_;
  Point const mid_point_;
  Real const d2_;
  Real mutable sq_row_dist_{};
};

void Poly::set_row(int const row) const {
  sq_row_dist_ = static_cast<Real>(square(mid_point_.y - row));
}

Real Poly::calc_at(int const col) const {
  auto const r2 = (sq_row_dist_ + square(mid_point_.x - col)) / d2_;
  auto const r4 = r2 * r2;
  auto const r6 = r4 * r2;
  auto const g = 1 + coeffs_[0] * r2 + coeffs_[1] * r4 + coeffs_[2] * r6;
  return g;
}

/**
The derivate of the given polynomial
g =  g_a,b,c(r) = 1 + ar^2 + br^4 + cr^6 is:
g'_a,b,c(r) = 2*ar + 4*br^3 + 6*cr^5 has to be > 0,
so we get - replacing q = r^2 and dividing by 2r:

a + 2*b*q + 3*c*q^2 > 0

defining:

d = 3 * a * c
b2 = b * b

q_plus = (-2*b + sqrt(4*b^2 - 12*ac)) / (6*c)  # d = 3ac
q_minus = (-2*b - sqrt(4*b^2 - 12*ac)) / (6*c)
=>
q_plus = (-b + sqrt(b - d)) / (3*c)
q_minus = (-b - sqrt(b - d)) / (3*c)

Any of the Conditions C1 - C9 must hold true for the derivative to be positive:

[Horizontal positive line]
    C1 = (a > 0 and b == c == 0);
[Increasing line with non-positive root]
    C2 = (a >= 0 and b > 0 and c = 0);
[Decreasing line with root >= 1]
    C3 = (c = 0 and b < 0 and -a <= 2b);
[Convex parab. without roots]
    C4 = (c > 0 and b2 < d);
[Convex parab., only one non-positive root]
    C5 = (c > 0 and b2 = d and b >= 0);
[Convex parab., only one root >= 1]
    C6 = (c > 0 and b2 = d and -b >= 3c);
[Convex parab., non-positive highest root]
    C7 = (c > 0 and b2 > d and q_plus <= 0);
[Convex parab., lowest root >= 1]
    C8 = (c > 0 and b2 > d and q_minus >= 1);
[Concave parab., lowest root <= 0; highest root >= 1]
C9 = (c < 0 and b2 > d and q_plus >= 1 and q_minus <= 0);
*/
bool Poly::is_increasing() const {
  auto const [a, b, c] = coeffs_;

  if ((a > 0 and b == 0 and c == 0) or  // C1
      (a >= 0 and b > 0 and c == 0)) {  // C2
    return false;
  }
  auto const d = 3 * a * c;
  auto const b2 = b * b;
  if ((c > 0 and b2 < d) or                   // C4
      (c > 0 and b2 == d and b >= 0) or       // C5
      (c > 0 and b2 == d and -b >= 3 * c)) {  // C6
    return true;
  }
  if (c == 0) {
    return b < 0 and -a <= 2 * b;  // shortcut + C3
  }
  auto const e = b2 - d;
  if (e <= 0) {
    return false;
  }
  // b2 > d ensured (with e > 0)
  auto const sqrt_e = sqrt(e);
  auto const c_3 = 3 * c;
  auto const q_plus = (-b + sqrt_e) / c_3;
  if (c > 0 and q_plus <= 0) {  // C7
    return true;
  }
  auto const q_minus = (-b - sqrt_e) / c_3;
  if (c > 0 and q_minus >= 1) {  // C8
    return true;
  }
  if (c < 0 and q_plus >= 1 and q_minus <= 0) {  // C9
    return true;
  }
  return false;
}

VignettingCorrection::VignettingCorrection(ImgOrig const &input_image, Factors const& factors)
    : factors_(factors),
	input_image_orig_(const_view(input_image)),
      input_image_(
          GaussianBlur::blur(scaled_down_gray(input_image_orig_, factors_.scale), factors_.blur)) {}

VignettingCorrection::~VignettingCorrection() {}

template <int SmoothRadius>
void VignettingCorrection::_smooth_histogram(
    HistogramType (&histogram)[HistogramSize + 1]) {
  HistogramType histo_orig[HistogramSize + 2 * SmoothRadius];

  // add mirrored values at the edges
  for (int i = 0; i < SmoothRadius; ++i) {
    // 0 -> 4, ... 3 -> 1
    histo_orig[i] = histogram[SmoothRadius - i];
    // 260 -> 254, 261 -> 253 ...
    histo_orig[HistogramSize + SmoothRadius + i] =
        histogram[HistogramSize - 2 - i];
  }
  memcpy(histo_orig + SmoothRadius, histogram,
         HistogramSize * sizeof(HistogramType));

  auto constexpr factor_sum = square(SmoothRadius + 1);
  //  smooth the histogram
  for (int i = 0; i < HistogramSize; ++i) {
    HistogramType sum = 0;
    auto *orig = &histo_orig[i];
    for (int j = 0; j < SmoothRadius; ++j) {
      sum += (j + 1) * orig[j];
    }
    sum += (SmoothRadius + 1) * orig[SmoothRadius];
    for (int j = 0; j < SmoothRadius; ++j) {
      sum += (SmoothRadius - j) * orig[SmoothRadius + j + 1];
    }
    histogram[i] = sum / factor_sum;
  }
}

Real VignettingCorrection::_calc_H(Poly const &poly) const {
  static auto constexpr fact = (Depth - 1) / static_cast<float>(log2i(Depth)) /
                               (MaxAllowedBrightness / 255.f);

  struct CalcH {
    CalcH() { memset(histogram, 0, sizeof(histogram)); }

    HistogramType histogram[MaxAllowedBrightness + 1];
    bool ok{true};
  } c{};

  auto const calc_pixel = [&poly](CalcH &c, auto const row_it,
                                  int const col) -> bool {
    auto const g = poly.calc_at(col);

    auto const log_img = fact * log2f(1 + g * row_it[col]);
	if (auto const k_lower =
		static_cast<unsigned>(log_img); k_lower < MaxAllowedBrightness) {
      // add discrete histogram value for actual floating point
      // example: 86.1 is to 0.9 in 86 and to 0.1 in 87
      auto const k = static_cast<HistogramType>(log_img) - k_lower;
      auto *hist = &c.histogram[k_lower];
      hist[0] += 1 - k;
      hist[1] += k;
      return true;
    }
    return false;
  };
  if (auto const result = reduce<CalcH>(
          c, const_view(input_image_),
          [&poly](int const row) { poly.set_row(row); }, calc_pixel);
      result) {
    return _calc_entropy(c.histogram);
  } else {
    return std::numeric_limits<Real>::max();
  }
}

Real VignettingCorrection::_calc_entropy(
    HistogramType (&histogram)[MaxAllowedBrightness + 1]) {
  _smooth_histogram<HistogramSmoothFactor>(histogram);
  float sum = 0;
  for (int i = 0; i < MaxAllowedBrightness; ++i) {
    sum += histogram[i];
  }
  Real H = 0;
  for (int i = 0; i < MaxAllowedBrightness; ++i) {
    if (auto const pk = histogram[i] / sum; pk) {
      // pk is < 1 (as divided by sum) => log < 0
      H -= pk * log2f(pk);
    }
  }
  return H;
}

Poly VignettingCorrection::_calc_best_poly() const {
  auto const mid = _center_of_mass();
  Poly::Coefficients best_coefficients = {0, 0, 0}, current_best = {0, 0, 0};
  auto H_min = _calc_H(Poly(best_coefficients, mid));
  if constexpr (DebugPrint) {
    std::cout << "Hmin " << H_min << " @ (0, 0, 0) " << std::endl;
  }

  // check the new coefficients are minimizing the H value
  auto const chk_H = [&](int const coeff_idx, int const sign, float const delta) -> bool {
    auto coeff{best_coefficients};
    coeff[coeff_idx] += sign * delta;
    auto const poly = Poly{coeff, mid};
    auto found_min = false;

    if (poly.is_increasing()) {
      if (auto const H = _calc_H(poly); H < H_min) {
        current_best = coeff;
        H_min = H;
        if constexpr (DebugPrint) {
          std::cout << "Hmin " << H_min << " @ (" << coeff[0] << ", " << coeff[1] << ", "
                    << coeff[2] << ")" << std::endl;
        }
        found_min = true;
      }
    }
    return found_min;
  };

  auto walkback = 0;
  float delta = DeltaStart;
  auto const find_dir = [&]() {
    auto found_dir = 0;
    for (int idx : {0, 1, 2}) {
      for (int sign : {-1, 1}) {
        if (auto const dir = (idx + 1) * sign;
            dir != walkback and chk_H(idx, sign, delta)) {
          // dir != walkback => do not walk the same way back again
          found_dir = dir;
        }
      }
    }
    walkback = -found_dir;
    return found_dir;
  };
  while (delta * DeltaMinDivider >= 1) {
    auto found_dir = find_dir();
    if (found_dir) {	  
	  best_coefficients = current_best;
	} else {
      delta /= 2;
    }
  }

  return Poly{best_coefficients, {mid.x * factors_.scale, mid.y * factors_.scale}};
}

Point VignettingCorrection::_center_of_mass() const {
  struct Params {
    float weight_x{0.f};
    float weight_y{0.f};
    float sum{0.f};
    int row;
  } p{};

  reduce<Params>(
      p, const_view(input_image_), [&p](int const row) { p.row = row; },
      [](Params &p, auto const &row_it, int const col) {
        auto const d = row_it[col];
        p.sum += d;
        p.weight_y += (p.row + 1) * d;
        p.weight_x += (col + 1) * d;
      });
  auto const s2 = p.sum / 2;
  return {div_round(p.weight_x, p.sum), div_round(p.weight_y, p.sum)};
}

ImgOrig VignettingCorrection::correct() {
  [[maybe_unused]] auto const start = std::chrono::high_resolution_clock::now();
  auto const dim = dimensions(input_image_orig_);
  auto const[cols, rows] = dim;
  auto const best_poly = _calc_best_poly();

  if constexpr (MeasureTime) {
    auto const end = std::chrono::high_resolution_clock::now();

    auto const duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "took: " << duration << " ms" << std::endl;
  }
  ImgOrig result(create_img(dim));

  auto v = view(result);
  const auto nc = num_channels(input_image_orig_);

  for (auto row = 0; row < rows; ++row) {
    auto row_it = row_begin<PixelOut>(input_image_orig_, row);
    auto output_it = row_begin<PixelOut>(v, row);
    best_poly.set_row(row);
    for (auto col = 0; col < cols; ++col) {
      auto const g = best_poly.calc_at(col);
      auto &pix_out = output_it[col];
      auto const &pix_in = row_it[col];
      for (int i = 0; i < nc; ++i) {
		pix_out[i] = clamp(g * pix_in[i]);
      }
    }
  }

  return result;
}

Factors VignettingCorrection::default_factors(imgalg::ImgOrig const &img) {
  auto const [cols, rows] = dimensions(img);
  auto factors = Factors{16, 19};
  if (cols < 2000) {
    factors.blur = 11;
    if (cols < 1000) {
      factors.scale = 8;
      factors.blur = (cols <= 500) ? 5 : 7;
    }
  }
  return factors;
}

}  // namespace vgncorr

int main(int argc, char *argv[]) {
  using namespace vgncorr;
  if (argc < 2) {
    std::cerr << "Please provide a filename\n";
    return -1;
  }
  std::string const path = argv[1];
 
  struct stat buf;
  if (stat(path.c_str(), &buf) != 0) {
	std::cerr << "File not found!\n";
	return -1;
  }

  ImgOrig orig;
  using namespace imgalg;

  ImageAlgo::load_image(orig, path);

  auto const factors = VignettingCorrection::default_factors(orig);
  vgncorr::VignettingCorrection corr(orig, factors);
  auto const out = corr.correct();
  auto out_path = path;
  out_path = out_path.replace(path.find("."), 1, "_corr.");
 
  ImageAlgo::save_image(out, out_path);

  return 0;
}
