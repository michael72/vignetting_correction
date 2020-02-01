#include "VignettingCorrection.h"

#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

#include <algorithm>
#include <array>
#include <ciso646>
#include <iostream>

namespace vgncorr {

static auto constexpr DebugPrint = false;
static auto constexpr MeasureTime = false;

auto constexpr _square = [](auto const x) -> auto { return x * x; };

template <typename T> static Real _sqrt(T sq) {
  return sqrtf(static_cast<Real>(sq));
}
static Real _dist(Point const &p) { return _sqrt(_square(p.x) + _square(p.y)); }

#if CALCULATE_EXACT == 0
// Helper for simple logarithmic calculation
static constexpr int _calc_log2i(int const num) {
  auto res = -1;
  for (int i = num; i; i >>= 1) {
    ++res;
  }
  if (res > 0 and num & (1 << (res - 1))) {
    // round up
    ++res;
  }
  return res;
}
// pre-calculated lookup table
struct LogTable {
  auto static constexpr Size = 1024;
  uchar table[Size];
  LogTable() {
    for (int i = 0; i < Size; ++i) {
      table[i] = _calc_log2i(i) & 0xffu;
    }
  }
};

static int _log2i(int const num) {
  static const LogTable t{};
  if (num < LogTable::Size) {
    return t.table[num];
  }
  return _calc_log2i(num);
}
#endif // ! CALCULATE_EXACT

/// Polynomial to calculate the vignetting correction with.
/// This is of the form:
/// g_a,b,c(r) = 1 + ar^2 + br^4 + cr^6.
class Poly {
public:
  static auto constexpr num_coefficients = 3;
  using CoeffType = Real;
  using Coefficients = std::array<CoeffType, num_coefficients>;

  Poly(Coefficients const &coefficients, Point const &mid_point)
      : coeffs_(coefficients), mid_point_(mid_point),
        d2_(_square(_dist(mid_point))) {}

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
  sq_row_dist_ = static_cast<Real>(_square(mid_point_.y - row));
}

Real Poly::calc_at(int const col) const {
  auto const r2 = (sq_row_dist_ + _square(mid_point_.x - col)) / d2_;
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

  if ((a > 0 and b == 0 and c == 0) or // C1
      (a >= 0 and b > 0 and c == 0)) { // C2
    return false;
  }
  auto const d = 3 * a * c;
  auto const b2 = b * b;
  if ((c > 0 and b2 < d) or                  // C4
      (c > 0 and b2 == d and b >= 0) or      // C5
      (c > 0 and b2 == d and -b >= 3 * c)) { // C6
    return true;
  }
  if (c == 0) {
    return b < 0 and -a <= 2 * b; // shortcut + C3
  }
  auto const e = b2 - d;
  if (e <= 0) {
    return false;
  }
  // b2 > d ensured (with e > 0)
  auto const sqrt_e = _sqrt(e);
  auto const c_3 = 3 * c;
  auto const q_plus = (-b + sqrt_e) / c_3;
  if (c > 0 and q_plus <= 0) { // C7
    return true;
  }
  auto const q_minus = (-b - sqrt_e) / c_3;
  if (c > 0 and q_minus >= 1) { // C8
    return true;
  }
  if (c < 0 and q_plus >= 1 and q_minus <= 0) { // C9
    return true;
  }
  return false;
}

VignettingCorrection::VignettingCorrection(cv::Mat const &input_image)
    : input_image_orig_(input_image) {
  cv::resize(input_image, input_image_, cv::Size(), 1. / ScaleFactor,
             1. / ScaleFactor);
}

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

  auto constexpr factor_sum = _square(SmoothRadius + 1);
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
  auto const rows = input_image_.rows;
  auto const cols = input_image_.cols;

  static auto const fact =
      (Depth - 1) / log2f(Depth) / (MaxAllowedBrightness / 256.f);

  HistogramType histogram[MaxAllowedBrightness + 1];
  memset(histogram, 0, sizeof(histogram));

  for (int row = 0; row < rows; ++row) {
    auto const *row_buffer = input_image_.ptr<uchar>(row);
    poly.set_row(row);
    for (int col = 0; col < cols; ++col) {
      auto const g = poly.calc_at(col);

#if CALCULATE_EXACT
      auto const log_img = fact * log2f(1 + g * row_buffer[col]);
      auto k_lower = static_cast<unsigned>(log_img);
      auto k_upper = k_lower + 1;
      if (k_lower >= MaxAllowedBrightness) {
        return std::numeric_limits<float>::max();
      }
      // add discrete histogram value for actual floating point
      // example: 86.1 is to 0.9 in 86 and to 0.1 in 87
      histogram[k_lower] += (k_upper - static_cast<HistogramType>(log_img));
      histogram[k_upper] += (log_img - static_cast<HistogramType>(k_lower));
#else
      auto const f = static_cast<int>(1 + g * row_buffer[col]);
      auto const log_img = fact * _log2i(f);
      auto const k = static_cast<unsigned>(log_img + 0.5f);
      if (k >= MaxAllowedBrightness) {
        return std::numeric_limits<float>::max();
      }
      ++histogram[k];
#endif
    }
  }

  return _calc_entropy(histogram);
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
  float delta = DeltaStart;
  auto const mid =
      _center_of_mass(); // Point{input_image_.cols / 2, input_image_.rows / 2};
  Poly::Coefficients best_coefficients = {0, 0, 0}, current_best = {0, 0, 0};
  auto H_min = _calc_H(Poly(best_coefficients, mid));

  auto const chk_H = [&](int const coeff_idx, int const sign) -> bool {
    auto coeff{best_coefficients};
    coeff[coeff_idx] += sign * delta;
    auto const poly = Poly{coeff, mid};
    auto found_min = false;

    if (poly.is_increasing()) {
      auto const H = _calc_H(poly);
      if (H < H_min) {
        current_best = coeff;
        H_min = H;
        if constexpr (DebugPrint) {
          std::cout << "amin " << coeff[0] << ", bmin " << coeff[1] << ", cmin "
                    << coeff[2] << ", Hmin " << H_min << std::endl;
        }
        found_min = true;
      }
    }
    return found_min;
  };

  auto last_dir = 0;
  while (delta * DeltaMinDivider > 1) {
    auto found_dir = 0;
    for (int idx : {0, 1, 2}) {
      for (int sign : {-1, 1}) {
        auto dir = (idx + 1) * sign;
        if (-dir != last_dir and chk_H(idx, sign)) {
          // -dir != last_dir => do not walk the same way back again
          found_dir = dir;
        }
      }
    }
    last_dir = found_dir;
    if (found_dir) {
      best_coefficients = current_best;
    } else {
      delta = delta / 2.f;
    }
  }

  return Poly{best_coefficients, {mid.x * ScaleFactor, mid.y * ScaleFactor}};
}

Point VignettingCorrection::_center_of_mass() const {
  cv::Mat img(input_image_.size(), CV_8UC1);
  cv::GaussianBlur(input_image_, img, cv::Size(GaussianSize, GaussianSize), 0);

  auto weight_x = 0.f;
  auto weight_y = 0.f;
  auto sum = 0.f;
  unsigned const rows = img.rows;
  unsigned const cols = img.cols;
  for (auto row = 0u; row < rows; ++row) {
    auto const *data = img.ptr<uchar>(row);
    for (auto col = 0u; col < cols; ++col) {
      unsigned const d = data[col];
      sum += d;
      weight_y += (row + 1) * d;
      weight_x += (col + 1) * d;
    }
  }
  auto const s2 = sum / 2;
  return {static_cast<int>((weight_x + s2) / sum),
          static_cast<int>((weight_y + s2) / sum)};
}

cv::Mat VignettingCorrection::correct() {

  [[maybe_unused]] auto const start = std::chrono::high_resolution_clock::now();
  auto const best_poly = _calc_best_poly();

  if constexpr (MeasureTime) {
    auto const end = std::chrono::high_resolution_clock::now();

    auto const duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "took: " << duration << " ms" << std::endl;
  }
  cv::Mat result(input_image_orig_.size(), CV_8UC1);

  for (auto row = 0; row < result.rows; ++row) {
    auto const *const row_data = input_image_orig_.ptr<uchar>(row);
    auto *const output_row = result.ptr<uchar>(row);
    best_poly.set_row(row);
    for (auto col = 0; col < result.cols; ++col) {
      auto const g = best_poly.calc_at(col);
      // TODO apply to color image
      output_row[col] = std::clamp(
          static_cast<int>(g * row_data[col] / MaxBrightnessFactor + 0.5f), 0,
          255);
    }
  }

  return result;
}

} // namespace vgncorr

int main(int argc, char *argv[]) {
  using namespace cv;
  if (argc < 2) {
	std::cerr << "Please provide a filename\n";
	return -1;
  }
  std::string const path = argv[1];
  Mat input_image = cv::imread(path, IMREAD_UNCHANGED);
  cv::Mat gray(input_image.size(), CV_8UC1);
  cv::cvtColor(input_image, gray, cv::COLOR_BGR2GRAY);
  imshow("raw", input_image);

  vgncorr::VignettingCorrection corr(gray);
  auto const out = corr.correct();

  imshow("corrected", out);

  waitKey(0);

  return 0;
}
