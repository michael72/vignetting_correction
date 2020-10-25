
#include "GaussianBlur.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace imgalg;
using namespace std;

// gaussian blur algorithm by Ivan Kutskir
// http://blog.ivank.net/fastest-gaussian-blur.html

static vector<int> boxesForGauss(int sigma,
                                 int n)  // standard deviation, number of boxes
{
  auto const wIdeal =
      sqrt((12 * sigma * sigma / n) + 1);  // Ideal averaging filter width
  int wl = static_cast<int>(wIdeal);
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

static void boxBlurH(PixelT const* scl, PixelT* tcl, int w,
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

static void boxBlurT(PixelT const* scl, PixelT* tcl, int w,
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

static void boxBlur(PixelT* scl, PixelT* tcl, int w, int h, int r) {
  memcpy(tcl, scl, w * h * sizeof(PixelT));
  boxBlurH(tcl, scl, w, h, r);
  boxBlurT(scl, tcl, w, h, r);
}

vector<PixelT> gaussBlur(vector<PixelT> const& source, int w, int h, int radius,
                         int n) {
  auto scl = source;
  auto tcl = vector<PixelT>(source.size());
  auto *pscl = &scl[0];
  auto *ptcl = &tcl[0];

  auto bxs = boxesForGauss(radius, n);
  boxBlur(pscl, ptcl, w, h, (bxs[0] - 1) / 2);
  boxBlur(ptcl, pscl, w, h, (bxs[1] - 1) / 2);
  boxBlur(pscl, ptcl, w, h, (bxs[2] - 1) / 2);

  return tcl;
}

Img GaussianBlur::blur(Img const& img, int radius, int n) {
  auto pixels = std::vector<PixelT>();
  auto const [cols, rows] = dimensions(img);
  auto const height = static_cast<int>(rows);
  auto const width = static_cast<int>(cols);

  pixels.reserve(width * height);

  reduce(
      pixels, const_view(img), [](int const row) {},
      [](auto& p, auto const row_it, int const col) {
        p.push_back(row_it[col]);
      });

  struct Param {
    std::vector<PixelT> blurred;
    int idx;
  };
  Param p = {gaussBlur(pixels, width, height, radius, n), 0};

  Img result = img;

  reduce(
      p, view(result), [](int const row) {},
      [](auto &p, auto &row_it, int const col) {
#ifdef USE_OPENCV
        row_it[col] = p.blurred[p.idx++];
#else
        const_cast<PixelT &>(row_it[col][0]) = p.blurred[p.idx++];
#endif
      });
  return result;
}
