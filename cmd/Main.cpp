#include "../lib/vgncorr.h"

#include <iostream>

int main(int argc, char *argv[]) {
  using namespace vgncorr;
  if (argc < 2) {
    std::cerr << "Please provide a filename\n";
    return -1;
  }
  std::string const path = argv[1];
 
  struct stat buf;
  if (stat(path.c_str(), &buf) != 0) {
  	std::cerr << "File not found!\n";
	  return -1;
  }
  auto out_path = path;
  out_path = out_path.replace(path.find("."), 1, "_corr.");

  vgncorr::Config config {};
  vgncorr::correct(config, path, out_path);

  return 0;
}
