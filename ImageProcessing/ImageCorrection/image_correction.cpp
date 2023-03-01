#include "image_correction.hpp"

#include <cmath>

namespace {
template <class T>
Magick::Image ProcessGrayPixel(const Magick::Image &image, T &&functor) {
  Magick::Image result = image;
  const auto geometry = result.size();

  result.modifyImage();
  result.type(Magick::GrayscaleType);
  auto pixels = result.getPixels(0, 0, geometry.width(), geometry.height());

  const auto get_off = [](const Magick::Geometry &geo, int x, int y) {
    return y * geo.width() + x;
  };

  for (std::size_t y = 0; y < geometry.height(); ++y) {
    for (std::size_t x = 0; x < geometry.width(); ++x) {
      auto &pixel = *(pixels + get_off(geometry, x, y));
      auto gray_pixel = Magick::ColorGray{pixel};
      pixel = functor(gray_pixel);
    }
  }

  result.syncPixels();

  return result;
}
} // namespace

namespace image_processing {
Magick::Image Negative(const Magick::Image &image) {
  return ProcessGrayPixel(image, [](Magick::ColorGray gray_pixel) {
      return Magick::ColorGray{1 - gray_pixel.shade()};
  });
}

Magick::Image GammaCorrection(const Magick::Image &image, float C, float y) {
  return ProcessGrayPixel(image, [&](Magick::ColorGray gray_pixel) {
      return Magick::ColorGray{std::min(1., C * std::pow(gray_pixel.shade(), y)) };
  });
}

Magick::Image LogCorrection(const Magick::Image &image, float C) {
  return ProcessGrayPixel(image, [&](Magick::ColorGray gray_pixel) {
      return Magick::ColorGray{std::min(1., C * std::log(gray_pixel.shade() + 1)) };
  });
}

} // namespace image_processing