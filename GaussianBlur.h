#pragma once

#include "ImageAlgo.h"

class GaussianBlur : public imgalg::ImageAlgo {
 public:
  static imgalg::Img blur(imgalg::Img const& src, int radius, int n = 3);
};
