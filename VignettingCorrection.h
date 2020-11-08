#pragma once

#include "ImageAlgo.h"

namespace vgncorr {
using Real = float;
using HistogramType = Real;

class Poly;

struct Config {
  // if scale or blur are 0 both values are determined automatically
  int scale = 0;
  int blur = 0;
  int histogram_smoothing = 16;
  float delta_start_divider = 1.f;
  int delta_max_precision = 1024;
};

class VignettingCorrection : public imgalg::ImageAlgo {
 public:
  VignettingCorrection(imgalg::ImgOrig const &img, Config const& config);
  ~VignettingCorrection();

  static Config default_factors(imgalg::ImgOrig const& img);

  imgalg::ImgOrig correct();
  static auto constexpr Depth = 256;
  static auto constexpr MaxAllowedBrightness = 255;
  static auto constexpr HistogramSize = MaxAllowedBrightness;

 private:
  Real _calc_H(Poly const &poly) const;
  Real _calc_entropy(
      HistogramType (&histogram)[MaxAllowedBrightness + 1]) const;

  imgalg::Point _center_of_mass() const;
  Poly _calc_best_poly() const;
  void _smooth_histogram(HistogramType (&histogram)[HistogramSize + 1]) const;
  
  Config config_;
  imgalg::ImgViewOrig const input_image_orig_;
  // scaled and gray version of the input image
  imgalg::Img input_image_;
};
}  // namespace vgncorr
