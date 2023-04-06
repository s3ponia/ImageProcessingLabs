#pragma once

#include <Magick++.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace image_processing {

template <class T, class... Images>
Magick::Image ProcessGrayPixel(T &&functor, const Magick::Image &image,
                               const Images &...additional_images) {
  static_assert((std::is_same_v<Images, Magick::Image> && ...));
  Magick::Image result = image;
  const auto geometry = result.size();

  assert((geometry == additional_images.size()) && ...);

  result.modifyImage();
  result.type(Magick::GrayscaleType);
  auto pixels = result.getPixels(0, 0, geometry.width(), geometry.height());

  std::tuple additional_pixels = {additional_images.getConstPixels(
      0, 0, geometry.width(), geometry.height())...};

  const auto get_off = [](const Magick::Geometry &geo, int x, int y) {
    return y * geo.width() + x;
  };

  for (std::size_t y = 0; y < geometry.height(); ++y) {
    for (std::size_t x = 0; x < geometry.width(); ++x) {
      auto &pixel = *(pixels + get_off(geometry, x, y));
      auto gray_pixel = Magick::ColorGray{pixel};
      auto add_pixels = std::apply(
          [&](auto... args) {
            return std::tuple{
                Magick::ColorGray{*(args + get_off(geometry, x, y))}...};
          },
          additional_pixels);

      std::apply(
          [&](auto... pixels) { pixel = functor(gray_pixel, pixels...); },
          add_pixels);
    }
  }

  result.syncPixels();

  return result;
}

template <class Container>
std::pair<double, double> FourierStep(const Container &data, int n,
                                      bool normalized) {
  double re = 0;
  double im = 0;

  auto N = std::size(data);

  for (size_t k = 0; k < N; ++k) {
    re += data[k] * std::cos(2 * M_PI * n * k / N);
    im += data[k] * std::sin(2 * M_PI * n * k / N);
  }

  if (normalized) {
    re /= N;
    im /= N;
  }

  return {re, im};
}

template <class Container> Container Fourier(const Container &data) {
  Container res(data);
  for (size_t i = 0; i < data.size(); ++i) {
    auto [re, im] = FourierStep(data, i, true);
    res[i] = std::sqrt(re * re + im * im);
  }
  return res;
}

template <class Container>
auto ComplexFourierWithoutSquare(const Container &data, bool normalize = true) {
  std::vector<std::complex<double>> res(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    auto [re, im] = FourierStep(data, i, normalize);
    res[i] = std::complex<double>{re, im};
  }
  return res;
}

template <class Container>
Container FourierWithoutSquare(const Container &data, bool normalize = true) {
  Container res(data);
  for (size_t i = 0; i < data.size(); ++i) {
    auto [re, im] = FourierStep(data, i, normalize);
    res[i] = re + im;
  }
  return res;
}

template <class Container> Container InverseFourier(const Container &data) {
  Container res(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    auto [re, im] = FourierStep(data, i, false);
    res[i] = re + im;
  }
  return res;
}

inline Magick::Image VecToImage(const std::vector<std::vector<double>> &vec2d) {
  Magick::Image result{Magick::Geometry(vec2d[0].size(), vec2d.size()),
                       Magick::ColorGray{1}};
  const auto geometry = result.size();

  result.modifyImage();
  result.type(Magick::GrayscaleType);
  auto pixels = result.getPixels(0, 0, geometry.width(), geometry.height());

  const auto get_off = [](const Magick::Geometry &geo, int x, int y) {
    return y * geo.width() + x;
  };

  double max = 0;
  double min = 99999999999;

  for (std::size_t y = 0; y < geometry.height(); ++y) {
    for (std::size_t x = 0; x < geometry.width(); ++x) {
      max = std::max(max, vec2d[y][x]);
      min = std::min(min, vec2d[y][x]);
    }
  }

  for (std::size_t y = 0; y < geometry.height(); ++y) {
    for (std::size_t x = 0; x < geometry.width(); ++x) {
      auto &pixel = *(pixels + get_off(geometry, x, y));
      pixel = Magick::ColorGray((vec2d[y][x] - min) / (max - min));
    }
  }

  result.syncPixels();

  return result;
}

inline std::vector<std::vector<double>>
ImageTo2dVec(const Magick::Image &image) {
  const auto geometry = image.size();
  std::vector<std::vector<double>> res(geometry.height(),
                                       std::vector<double>(geometry.width()));

  auto pixels = image.getConstPixels(0, 0, geometry.width(), geometry.height());

  const auto get_off = [](const Magick::Geometry &geo, int x, int y) {
    return y * geo.width() + x;
  };

  for (std::size_t y = 0; y < geometry.height(); ++y) {
    for (std::size_t x = 0; x < geometry.width(); ++x) {
      auto &pixel = *(pixels + get_off(geometry, x, y));
      res[y][x] = Magick::ColorGray{pixel}.shade() * 65536;
    }
  }

  return res;
}

template <class Container>
Container Fourier2DWithoutSquare(const Container &data) {
  Container res(data);
  for (size_t i = 0; i < data[0].size(); ++i) {
    std::vector<double> curCol(data.size());
    for (size_t x = 0; x < curCol.size(); ++x) {
      curCol[x] = data[x][i];
    }
    auto fft = FourierWithoutSquare(curCol);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[j][i] = fft[j];
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<double> curRow(data[0].size());
    for (size_t x = 0; x < curRow.size(); ++x) {
      curRow[x] = res[i][x];
    }
    auto fft = FourierWithoutSquare(curRow);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[i][j] = fft[j];
    }
  }

  return res;
}

template <class Container> Container Fourier2D(const Container &data) {
  Container res(data);
  for (size_t i = 0; i < data[0].size(); ++i) {
    std::vector<double> curCol(data.size());
    for (size_t x = 0; x < curCol.size(); ++x) {
      curCol[x] = data[x][i];
    }
    auto fft = FourierWithoutSquare(curCol);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[j][i] = fft[j];
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<double> curRow(data[0].size());
    for (size_t x = 0; x < curRow.size(); ++x) {
      curRow[x] = res[i][x];
    }
    auto fft = Fourier(curRow);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[i][j] = fft[j];
    }
  }

  return res;
}

template <class Container> Container InverseFourier2D(const Container &data) {
  Container res(data);
  for (size_t i = 0; i < data[0].size(); ++i) {
    std::vector<double> curCol(data.size());
    for (size_t x = 0; x < curCol.size(); ++x) {
      curCol[x] = data[x][i];
    }
    auto fft = FourierWithoutSquare(curCol, false);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[j][i] = fft[j];
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<double> curRow(data[0].size());
    for (size_t x = 0; x < curRow.size(); ++x) {
      curRow[x] = res[i][x];
    }
    auto fft = FourierWithoutSquare(curRow, false);
    for (size_t j = 0; j < fft.size(); ++j) {
      res[i][j] = fft[j];
    }
  }

  return res;
}

inline double HarmFunctionImpl(double A, double f, size_t i, double dt) {
  return A * std::sin(2 * M_PI * f * i * dt);
}

inline std::vector<double> HarmFunction(double A, double f, size_t size,
                                        double dt) {
  std::vector<double> res(size);

  for (size_t i = 0; i < res.size(); ++i) {
    res[i] = HarmFunctionImpl(A, f, i, dt);
  }

  return res;
}

inline std::pair<std::vector<double>, std::vector<double>>
SpecterFourier(const std::vector<double> &data, double dt) {
  auto borderF = 1. / (2. * dt);
  auto deltaF = borderF / (data.size() / 2.);
  std::vector<double> x(data.size() / 2);
  std::vector<double> y(data.size() / 2);

  for (size_t i = 0; i < x.size(); ++i) {
    x[i] = deltaF * i;
    y[i] = data[i];
  }

  return {std::move(x), std::move(y)};
}

inline std::vector<std::vector<double>>
DiagonalShift2d(std::vector<std::vector<double>> matrix, int rowShift,
                int colShift) {
  std::vector<std::vector<double>> res(matrix);

  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[0].size(); ++j) {
      auto posX = (i + matrix.size() - rowShift) % matrix.size();
      auto posY = (j + matrix[0].size() - colShift) % matrix.size();
      res[i][j] = matrix[posX][posY];
    }
  }

  return res;
}

inline std::vector<double> MultModel(const std::vector<double> &data1,
                                     const std::vector<double> &data2) {
  std::vector<double> res(data1.size());
  for (size_t i = 0; i < res.size(); ++i) {
    res[i] = data1[i] * data2[i];
  }
  return res;
}

inline std::vector<double> Convolution(const std::vector<double> &data1,
                                       const std::vector<double> &data2) {
  std::vector<double> res;
  const auto N = data2.size();
  const auto M = data1.size();

  for (size_t k = 0; k < N + M; ++k) {
    auto val = 0.;
    for (size_t m = 0; m < M; ++m) {
      if (k - m < 0 || k - m >= N) {
        continue;
      }
      val += data1[m] * data2[k - m];
    }
    res.push_back(val);
  }

  return res;
}

inline std::vector<double> InverseFilter(const std::vector<double> &main_data,
                                         std::vector<double> h_func,
                                         bool noised = false,
                                         double alpha = 0) {
  h_func.resize(std::max(h_func.size(), main_data.size()));

  auto complex_h_func = ComplexFourierWithoutSquare(h_func);
  auto complex_main = ComplexFourierWithoutSquare(main_data);

  if (!noised) {
    for (size_t i = 0; i < complex_main.size(); ++i) {
      complex_main[i] /= complex_h_func[i];
    }
  } else {
    auto div = complex_main;

    for (size_t i = 0; i < div.size(); ++i) {
      div[i] *= std::conj(complex_h_func[i]);
    }

    for (size_t i = 0; i < complex_main.size(); ++i) {
      auto divider = std::abs(complex_h_func[i]) * std::abs(complex_h_func[i]) +
                     alpha * alpha;
      complex_main[i] = div[i] / divider;
    }
  }

  std::vector<double> new_line(complex_main.size());

  for (size_t i = 0; i < new_line.size(); ++i) {
    new_line[i] = std::real(complex_main[i]) + std::imag(complex_main[i]);
  }

  return FourierWithoutSquare(new_line);
}

inline std::pair<std::vector<double>, std::vector<double>> CreateHeartBeat() {
  const auto A = 1;
  const auto f = 7;
  const auto dt = 0.005;
  const auto M = 200;
  const auto a = 30;
  const auto b = 1;

  std::vector<double> x(M);

  for (size_t i = 0; i < M; ++i) {
    x[i] = i;
  }

  std::vector<double> h1(x.size());

  for (size_t i = 0; i < M; ++i) {
    h1[i] = A * std::sin(2 * M_PI * f * dt * x[i]);
  }

  std::vector<double> h2(x.size());

  for (size_t i = 0; i < M; ++i) {
    h2[i] = b * std::exp(-a * dt * x[i]);
  }

  auto h = MultModel(h1, h2);
  const double max = *std::max_element(h.begin(), h.end());

  for (auto &el : h) {
    el = el / max * 120.;
  }

  for (auto &el : x) {
    el *= dt;
  }

  return {std::move(x), std::move(h)};
}

inline std::vector<std::vector<double>>
ImageRestore(std::vector<std::vector<double>> data, std::vector<double> h_func,
             bool noised = false, double alpha = 0) {
  std::vector<std::vector<double>> res(data);

  for (size_t i = 0; i < data.size(); ++i) {
    std::vector<double> curRow(data[0].size());
    for (size_t x = 0; x < curRow.size(); ++x) {
      curRow[x] = res[i][x];
    }

    auto inversed = InverseFilter(curRow, h_func, noised, alpha);
    for (size_t j = 0; j < curRow.size(); j++)
      res[i][j] = inversed[j];
  }

  return res;
}

inline std::vector<double> ReadDataFile(std::string path) {
  if (std::ifstream is{path, std::ios::binary | std::ios::ate}) {
    auto size = is.tellg();
    std::vector<float> res(size);
    is.seekg(0);
    if (is.read((char *)res.data(), size)) {
      return {res.begin(), res.end()};
    }
  }

  return {};
}

inline std::vector<std::vector<double>> FormImage(std::vector<double> data, int rows,
                                           int cols) {
  std::vector<std::vector<double>> res(rows, std::vector<double>(cols));

  int pos = 0;

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      res[i][j] = data[pos++];
    }
  }

  return res;
}
} // namespace image_processing
