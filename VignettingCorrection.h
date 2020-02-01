#pragma once

#include "opencv2/core.hpp"

#define CALCULATE_EXACT 0

namespace vgncorr {
using Point = cv::Point2i;
using Real = float;
#if CALCULATE_EXACT
using HistogramType = Real;
static constexpr auto exact = true;
#else
using HistogramType = int;
static constexpr auto exact = false;
#endif

class Poly;

class VignettingCorrection {
public:
  VignettingCorrection(cv::Mat const &input_image);
  ~VignettingCorrection();

  cv::Mat correct();
  static auto constexpr Depth = 256;
  static auto constexpr MaxBrightnessFactor = 1.1f;
  static_assert(MaxBrightnessFactor > 1.f);
  static auto constexpr MaxAllowedBrightness =
      static_cast<int>(256 * MaxBrightnessFactor + 0.5f);
  static auto constexpr HistogramSize = MaxAllowedBrightness;
  static auto constexpr ScaleFactor = 4;
  static auto constexpr GaussianSize = 80 / ScaleFactor + 1;
  static_assert((GaussianSize & 1) == 1); // Gaussian must be odd sized

  static auto constexpr DeltaStart = 8.f;
  static auto constexpr DeltaMinDivider = 256; // smallest delta is 1 / DeltaMinDivider

  static auto constexpr HistogramSmoothFactor = 8;

private:
  Real _calc_H(Poly const &poly) const;
  static Real _calc_entropy(HistogramType(&histogram)[MaxAllowedBrightness + 1]);

  Point _center_of_mass() const;
  Poly _calc_best_poly() const;
  template <int SmoothRadius = 4>
  static void _smooth_histogram(HistogramType (&histogram)[HistogramSize + 1]);

  cv::Mat const &input_image_orig_;
  // scaled and gray version of the input image
  cv::Mat input_image_;
};
} // namespace vgncorr
