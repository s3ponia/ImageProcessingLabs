#include "image_correction.hpp"

#include <cmath>

#include "utility.hpp"

namespace image_processing {
Magick::Image Negative(const Magick::Image &image) {
  return ProcessGrayPixel(image, [](Magick::ColorGray gray_pixel) {
    return Magick::ColorGray{1 - gray_pixel.shade()};
  });
}

Magick::Image GammaCorrection(const Magick::Image &image, float C, float y) {
  return ProcessGrayPixel(image, [&](Magick::ColorGray gray_pixel) {
    return Magick::ColorGray{std::min(1., C * std::pow(gray_pixel.shade(), y))};
  });
}

Magick::Image LogCorrection(const Magick::Image &image, float C) {
  return ProcessGrayPixel(image, [&](Magick::ColorGray gray_pixel) {
    return Magick::ColorGray{
        std::min(1., C * std::log(gray_pixel.shade() + 1))};
  });
}

} // namespace image_processing