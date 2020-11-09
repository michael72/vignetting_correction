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

  vgncorr::VignettingCorrection corr(orig, config);
  auto const out = corr.correct();
 
  ImageAlgo::save_image(out, filename_out);
}
}