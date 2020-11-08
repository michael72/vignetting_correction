#include "vgncorr.h"

#include "VignettingCorrection.h"

namespace vgncorr 
{
using namespace imgalg;

void correct(Config const& config, std::string const& filename_in, std::string const& filename_out)
{
  ImgOrig orig;
  using namespace imgalg;

  ImageAlgo::load_image(orig, filename_in);

  auto const factors = VignettingCorrection::default_factors(orig);
  vgncorr::VignettingCorrection corr(orig, factors);
  auto const out = corr.correct();
 
  ImageAlgo::save_image(out, filename_out);
}
}