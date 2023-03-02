#pragma once

#include <Magick++.h>

namespace image_processing {

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
} // namespace image_processing
