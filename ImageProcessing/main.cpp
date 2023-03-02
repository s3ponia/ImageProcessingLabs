#define FMT_HEADER_ONLY

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include <Magick++.h>
#include <argparse/argparse.hpp>
#include <fmt/format.h>

#include "ImageCorrection/image_correction.hpp"
#include "ImageProcessing/image_processing.hpp"
#include "ImageResize/bilinear_interpolation.hpp"
#include "ImageResize/nearest_neighbourd.hpp"

void PrintImageGeometry(const Magick::Image &image) {
  Magick::Geometry geometry = image.size();

  std::cout << fmt::format("width: {}\nheight: {}", geometry.width(),
                           geometry.height())
            << std::endl;
}

Magick::Image readXCRFile(const std::filesystem::path &path,
                          const Magick::Geometry &geometry) {
  static const int PIXEL_SIZE = 2;
  Magick::Image image;
  image.size(geometry);
  image.magick("GRAY");
  image.depth(8);

  std::ifstream input_file(path, std::ios::binary);

  input_file.seekg(2048);

  const auto data_size = geometry.height() * geometry.width();

  auto data = std::make_unique<unsigned short[]>(data_size);

  input_file.read((char *)data.get(), data_size * PIXEL_SIZE);

  for (std::size_t i = 0; i < data_size; ++i) {
    auto *bytes = (unsigned char *)&data[i];
    std::swap(bytes[0], bytes[1]);
  }

  auto recalced_data = std::make_unique<unsigned char[]>(data_size);
  const auto min = *std::min_element(data.get(), data.get() + data_size);
  const auto max = *std::max_element(data.get(), data.get() + data_size);

  for (std::size_t i = 0; i < data_size; ++i) {
    recalced_data[i] =
        (unsigned char)(1.f * (data[i] - min) / (max - min) *
                        std::numeric_limits<unsigned char>::max());
  }

  Magick::Blob blob;
  blob.updateNoCopy(recalced_data.release(), data_size);

  image.read(blob);

  return image;
}

Magick::Image readFile(const std::filesystem::path &path,
                       std::optional<Magick::Geometry> geometry) {
  if (path.extension() == ".xcr") {
    if (geometry) {
      return readXCRFile(path, *geometry);
    } else {
      throw std::logic_error{"geometry must be specified for .xcr files"};
    }
  } else {
    return Magick::Image{path};
  }
}

template <class T> Magick::Geometry geometryFromVec(const T &vec) {
  return {vec[0], vec[1]};
}

int main(int argc, char **argv) {
  Magick::InitializeMagick(nullptr);

  argparse::ArgumentParser program("ImageProcessing");

  program.add_argument("image_path").help("path to the image to work with");

  program.add_argument("-t", "--type")
      .required()
      .help("specify the programm's behaviour: geo - prints image's geometry; "
            "shift - add constant to all pixels in image; "
            "mult - multiply by constant all pixels in image; "
            "rotate - rotate images on 90 degree; "
            "recalc - recalculate image pixels to range [0;255]; "
            "resize - resize image; "
            "gamma - apply gamma correction");

  program.add_argument("-s", "--size")
      .help("specify width and height: --size width height. Required *.xcr "
            "images")
      .nargs(2)
      .scan<'u', unsigned int>();

  program.add_argument("-ns", "--new-size")
      .help("specify new size for image (required for resize type)")
      .nargs(2)
      .scan<'u', unsigned int>();

  program.add_argument("-rt", "--resize-type")
      .default_value("nearest")
      .help("specify type of resizing: nearest or bilinear. Default: nearest");

  program.add_argument("-gt", "--gamma-type")
      .default_value("negative")
      .help("specify function for gamma correction: negative (default), gamma, "
            "logarithm");

  program.add_argument("-C")
      .help("specify parameter for gamma correction, multipling and shifting")
      .scan<'g', double>()
      .default_value(1.);

  program.add_argument("-y")
      .help("specify parameter for gamma correction")
      .scan<'g', double>()
      .default_value(2.2);

  program.add_argument("-o", "--output")
      .help("must be specified with image modifing functions");

  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  const auto size_vec = program.present<std::vector<unsigned int>>("--size");
  const auto geometry = size_vec
                            ? std::optional{geometryFromVec(size_vec.value())}
                            : std::nullopt;
  const auto image_path = program.get<std::string>("image_path");
  const auto runtime_type = program.get<std::string>("--type");
  auto image = readFile(image_path, geometry);

  if (runtime_type == "geo") {
    PrintImageGeometry(image);
  } else if (runtime_type == "shift") {
    image = image_processing::AddConstant(image,
                                          (uint8_t)program.get<double>("-C"));
  } else if (runtime_type == "mult") {
    image = image_processing::Mult(image, program.get<double>("-C"));
  } else if (runtime_type == "rotate") {
    image.rotate(270);
  } else if (runtime_type == "recalc") {
  } else if (runtime_type == "resize") {
    if (program.is_used("-ns")) {
      const auto res_type = program.get<std::string>("--resize-type");
      if (res_type == "nearest") {
        image = image_processing::ResizeNearestNeighbourd(
            image,
            geometryFromVec(program.get<std::vector<unsigned int>>("-ns")));
      } else if (res_type == "bilinear") {
        image = image_processing::ResizeBilinear(
            image,
            geometryFromVec(program.get<std::vector<unsigned int>>("-ns")));
      } else {
        std::cerr << "Unknown resize type" << std::endl;
        return 1;
      }
    } else {
      std::cerr << "--new-size: required for resize action" << std::endl;
      return 1;
    }
  } else if (runtime_type == "gamma") {
    const auto gamma_type = program.get<std::string>("--gamma-type");
    const auto c_param = program.get<double>("-C");
    const auto y_param = program.get<double>("-y");
    if (gamma_type == "negative") {
      image = image_processing::Negative(image);
    } else if (gamma_type == "gamma") {
      image = image_processing::GammaCorrection(image, c_param, y_param);
    } else if (gamma_type == "logarithm") {
      image = image_processing::LogCorrection(image, c_param);
    } else {
      std::cerr << "Unhandled --gt: " << runtime_type << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Unhandled --type: " << runtime_type << std::endl;
    return 1;
  }

  if (program.is_used("--output")) {
    image.write(program.get<std::string>("--output"));
  } else {
    std::cout << "--output: can be required for some actions" << std::endl;
  }

  return 0;
}
