#pragma once
#include <algorithm>
#include <string>
#include <cctype>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)  // suppress warning to use wcstomb_s in
#pragma warning(disable : 4244)
#endif
#include <boost/gil.hpp>
#include <boost/gil/extension/io/bmp.hpp>
#include <boost/gil/extension/io/jpeg.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil/extension/numeric/convolve.hpp>
#include <boost/gil/extension/numeric/kernel.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

namespace imgalg {
using Pixel = boost::gil::gray8_pixel_t;
using PixelT = uint8_t;
using PixelOrig = boost::gil::rgb8_pixel_t;
using PixelOut = PixelOrig;
using Img = boost::gil::image<Pixel>;
using ImgOrig = boost::gil::image<PixelOrig>;
using ImgView = boost::gil::image<Pixel>::const_view_t;
using ImgViewOrig = boost::gil::image<PixelOrig>::const_view_t;
using Point = ImgView::point_t;

using Size = Point;

class ImageAlgoBase {
 public:
  template <typename V>
  static Size dimensions(V const &v) {
    return v.dimensions();
  }
  template <typename V> 
  static int num_channels(V const &v) {
    return boost::gil::num_channels<ImgViewOrig>::value;
  }

  static Size create_img(Size const &size) { return size; }

  template <typename P, typename V>
  static P *row_begin(V &img, int const row) {
    return img.row_begin(row);
  }

  template <typename P, typename V>
  static P const *row_begin(V const &img, int const row) {
    return img.row_begin(row);
  }

  static Img scaled_down_gray(ImgViewOrig const &input_image, int const SF) {
    auto const &orig_dim = input_image.dimensions();
    auto const dim = Point{orig_dim.x / SF, orig_dim.y / SF};
    Img result(dim);
    Img gray(orig_dim);
    auto v = view(result);
    copy_pixels(::boost::gil::color_converted_view<Pixel>(input_image),
                view(gray));
    auto const [cols, rows] = dim;
    auto static const SF2 = SF * SF;
    for (int row = 0; row < rows; ++row) {
      auto out = row_begin<Pixel>(v, row);
      for (int col = 0; col < cols; ++col) {
        auto sum = 0u;
        for (int y = 0; y < SF; ++y) {
          auto const &row_it = const_view(gray).row_begin(row * SF + y);
          for (int x = 0; x < SF; ++x) {
            sum += row_it[x + col * SF];
          }
        }
        out[col] = (sum + SF2 / 2 - 1) / SF2;
      }
    }
    return result;
  }

  template <typename Fun>
  static void img_action(std::string const &path, Fun const &fun) {
	auto const lower = [](auto const& s) {
		std::string data = s;
		std::transform(data.cbegin(), data.cend(), data.begin(),
		  [](unsigned char c) { return std::tolower(c); });
		return data;
	};
    auto const has_ending = [&](std::string const &ending) -> bool {
      return path.length() >= ending.length() and
             not lower(path).compare(path.length() - ending.length(), ending.length(),
                              ending);
    };
    if (has_ending(".png")) {
      fun(boost::gil::png_tag());
    } else if (has_ending(".jpg") or has_ending(".jpeg")) {
      fun(boost::gil::jpeg_tag());
    } else if (has_ending(".bmp")) {
      fun(boost::gil::bmp_tag());
    } else {
      throw std::runtime_error(std::string("Unrecognized file type in ") +
                               path);
    }
  }

  template <typename I>
  static void load_image(I &orig, std::string const &path) {
    img_action(path,
               [&](auto tag) { boost::gil::read_image(path, orig, tag); });
  }

  template <typename I>
  static void save_image(I const &img, std::string const &path) {
    img_action(path, [&](auto tag) {
      boost::gil::write_view(path, const_view(img), tag);
    });
  }
};

}  // namespace imgalg
