#include "image_processing.hpp"

#include <cmath>
#include <limits>

#include "utility.hpp"

namespace image_processing {
Magick::Image Mult(const Magick::Image &image, float C) {
  return ProcessGrayPixel(image, [&](auto gray_pixel) {
    return Magick::ColorGray{std::min(1., gray_pixel.shade() * C)};
  });
}
Magick::Image AddConstant(const Magick::Image &image, uint8_t C) {
  return ProcessGrayPixel(image, [&](auto gray_pixel) {
    return Magick::ColorGray{
        std::min(1., gray_pixel.shade() +
                         (float)C / std::numeric_limits<decltype(C)>::max())};
  });
}
} // namespace image_processing