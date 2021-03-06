# Vignetting Correction

Automatic Vignetting Correction algorithm in C++ after Laura Lopez-Fuentes's work
[Revisiting Image Vignetting Correction by Constrained Minimization of log-Intensity Entropy](https://www.researchgate.net/publication/300786398_Revisiting_Image_Vignetting_Correction_by_Constrained_Minimization_of_Log-Intensity_Entropy) which is based on the paper [Single-Image Vignetting Correction by Constrained Minimization of log-Intensity Entropy](https://www.semanticscholar.org/paper/Single-Image-Vignetting-Correction-by-Constrained-Torsten/e355fffc31fa0a7c5309bd2b90da84810e5ffb70)

See also
* https://github.com/HJCYFY/Vignetting-Correction
* https://github.com/dajuric/dot-devignetting

## Build the project

### Prerequisites

Install:
- [Python](https://python.org)
- [Conan](https://conan.io)
- [CMake](https://cmake.org)

On Ubuntu:
```shell
./configure.sh
./build.sh
./bin/vgncorr <img-filename>
```

On Windows:
```shell
build.bat
bin\vgncorr.exe <img-filename>
```

If you want to use opencv library instead of the default boost gil, then execute the corresponding `configure_opencv` and `build_opencv` files.

The resulting image is copied to the same folder.

