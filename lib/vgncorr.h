#pragma once

#include <string>

namespace vgncorr 
{
struct Config 
{
  // if scale or blur are 0 both values are determined automatically
  int scale = 0;
  int blur = 0;
  int histogram_smoothing = 16;
  float delta_start_divider = 1.f;
  int delta_max_precision = 1024;
};

void correct(Config const& config, std::string const& filename_in, std::string const& filename_out);

}