#pragma once

#include "ImageAlgo.h"

namespace vgncorr {
using Real = float;
using HistogramType = Real;

class Poly;

struct Factors {
  int scale;
  int blur;
};

class VignettingCorrection : public imgalg::ImageAlgo {
 public:
  VignettingCorrection(imgalg::ImgOrig const &img, Factors const& factors);
  ~VignettingCorrection();

  static Factors default_factors(imgalg::ImgOrig const& img);

  imgalg::ImgOrig correct();
  static auto constexpr Depth = 256;
  static auto constexpr MaxAllowedBrightness = 255;
  static auto constexpr HistogramSize = MaxAllowedBrightness;
  
  static auto constexpr DeltaStart = 1.f;
  static auto constexpr DeltaMinDivider =
      256;  // smallest delta is 1 / DeltaMinDivider

  static auto constexpr HistogramSmoothFactor = 16;

 private:
  Real _calc_H(Poly const &poly) const;
  static Real _calc_entropy(
      HistogramType (&histogram)[MaxAllowedBrightness + 1]);

  imgalg::Point _center_of_mass() const;
  Poly _calc_best_poly() const;
  template <int SmoothRadius = 4>
  static void _smooth_histogram(HistogramType (&histogram)[HistogramSize + 1]);
  
  Factors factors_;
  imgalg::ImgViewOrig const input_image_orig_;
  // scaled and gray version of the input image
  imgalg::Img input_image_;
};
}  // namespace vgncorr
