
#include "GaussianBlur.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace imgalg;
using namespace std;

constexpr static int iround(float const f) {
  // only works for positive f
  return static_cast<int>(f + 0.5f);
}

constexpr static int iround(int const i) { return i; }

constexpr static PixelT clamp(int const i) {
  return static_cast<PixelT>(((0xff - i) >> 31) | (i & ~(i >> 31)));
}

constexpr static PixelT clamp(float const f) { return clamp(iround(f)); }

// gaussian blur algorithm by Ivan Kutskir
// http://blog.ivank.net/fastest-gaussian-blur.html

static vector<int> boxesForGauss(float sigma,
                                 int n)  // standard deviation, number of boxes
{
  auto const wIdeal =
      sqrt((12 * sigma * sigma / n) + 1);  // Ideal averaging filter width
  int wl = floor(wIdeal);
  if (wl % 2 == 0) --wl;
  int wu = wl + 2;

  auto mIdeal =
      (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
  int m = iround(mIdeal);

  vector<int> sizes(n);
  for (auto i = 0; i < n; i++) {
    sizes[i] = i < m ? wl : wu;
  }
  return sizes;
}

static void boxBlurH(vector<PixelT> const& scl, vector<PixelT>& tcl, int w,
                     int h, int r) {
  float const iarr = 1.f / (2 * r + 1);
  for (auto i = 0; i < h; i++) {
    auto ti = i * w, li = ti, ri = ti + r;
    auto const fv = scl[ti], lv = scl[ti + w - 1];
    auto val = (r + 1) * fv;
    for (auto j = 0; j < r; ++j) {
      val += scl[ti + j];
    }
    for (auto j = 0; j <= r; ++j) {
      val += scl[ri++] - fv;
      tcl[ti++] = clamp(val * iarr);
    }
    for (auto j = r + 1; j < w - r; ++j) {
      val += scl[ri++] - scl[li++];
      tcl[ti++] = clamp(val * iarr);
    }
    for (auto j = w - r; j < w; ++j) {
      val += lv - scl[li++];
      tcl[ti++] = clamp(val * iarr);
    }
  }
}

static void boxBlurT(vector<PixelT> const& scl, vector<PixelT>& tcl, int w,
                     int h, int r) {
  float const iarr = 1.f / (2 * r + 1);
  for (auto i = 0; i < w; i++) {
    auto ti = i, li = ti, ri = ti + r * w;
    auto const fv = scl[ti], lv = scl[ti + w * (h - 1)];
    auto val = (r + 1) * fv;
    for (auto j = 0; j < r; ++j) {
      val += scl[ti + j * w];
    }
    for (auto j = 0; j <= r; ++j) {
      val += scl[ri] - fv;
      tcl[ti] = clamp(val * iarr);
      ri += w;
      ti += w;
    }
    for (auto j = r + 1; j < h - r; ++j) {
      val += scl[ri] - scl[li];
      tcl[ti] = clamp(val * iarr);
      li += w;
      ri += w;
      ti += w;
    }
    for (auto j = h - r; j < h; ++j) {
      val += lv - scl[li];
      tcl[ti] = clamp(val * iarr);
      li += w;
      ti += w;
    }
  }
}

static void boxBlur(vector<PixelT>& scl, vector<PixelT>& tcl, int w, int h,
                    int r) {
  tcl = scl;
  boxBlurH(tcl, scl, w, h, r);
  boxBlurT(scl, tcl, w, h, r);
}

vector<PixelT> gaussBlur(vector<PixelT> const& source, int w, int h, int radius,
                         int n) {
  auto scl = source;
  auto tcl = vector<PixelT>(source.size());

  auto bxs = boxesForGauss(radius, n);
  boxBlur(scl, tcl, w, h, (bxs[0] - 1) / 2);
  boxBlur(tcl, scl, w, h, (bxs[1] - 1) / 2);
  boxBlur(scl, tcl, w, h, (bxs[2] - 1) / 2);

  return tcl;
}

Img GaussianBlur::blur(Img const& img, int radius, int n) {
  auto pixels = std::vector<PixelT>();
  pixels.reserve(img.width() * img.height());

  reduce(
      pixels, const_view(img), [](int const row) {},
      [](auto& p, auto const row_it, int const col) {
        p.push_back(row_it[col]);
      });

  struct Param {
    std::vector<PixelT> blurred;
    int idx;
  };
  Param p = {gaussBlur(pixels, img.width(), img.height(), radius, n), 0};

  Img result = img;

  reduce(
      p, view(result), [](int const row) {},
      [](auto& p, auto& row_it, int const col) {
        const_cast<PixelT&>(row_it[col][0]) = p.blurred[p.idx++];
      });
  return result;
}