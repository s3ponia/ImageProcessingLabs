#pragma once

#include <Magick++.h>

namespace image_processing {
Magick::Image Mult(const Magick::Image &image, float C);
Magick::Image AddConstant(const Magick::Image &image, uint8_t C);
} // namespace image_processing