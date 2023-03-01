#pragma once

#include <Magick++.h>

namespace image_processing {
Magick::Image Negative(const Magick::Image &image);
Magick::Image GammaCorrection(const Magick::Image &image, float C, float y);
Magick::Image LogCorrection(const Magick::Image &image, float C);
} // namespace image_processing
