#pragma once
#include <algorithm>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

namespace imgalg {

class ImageAlgoBase {
public:
  using Point = cv::Point2i;
  using ImgOrig = cv::Mat;
  using ImgViewOrig = cv::Mat;
  using Img = ImgOrig;
  using Pixel = uint8_t;
  using PixelOrig = uint8_t;
  using PixelT = uint8_t;

  static inline cv::Mat const &const_view(cv::Mat const &mat) { return mat; }
  static inline cv::Mat &view(cv::Mat &mat) { return mat; }
  using Size = cv::Size;
  static Size dimensions(cv::Mat const &mat) { return mat.size(); }
  static Img create_img(Size const &size) { return cv::Mat{size, CV_8UC1}; }
  static Img scaled_down(ImgViewOrig const &input_image, int const SF) {
    Img result;
    cv::resize(input_image, result, cv::Size(), 1. / SF, 1. / SF,
               cv::INTER_LINEAR_EXACT);
    return result;
  }
  template <typename P, typename V> 
  static P *row_begin(V &img, int const row) {
    return img.ptr(row);
  }

  template <typename P, typename V>
  static P const *row_begin(V const &img, int const row) {
    return img.ptr(row);
  }

  static void load_image(ImgOrig &orig, cv::String const &path) {
    orig = cv::imread(path, cv::IMREAD_UNCHANGED);
  }

  static bool save_image(ImgOrig const &img, cv::String const &path) {
    return cv::imwrite(path, (cv::InputArray)img);
  }
};
} // namespace imgalg
