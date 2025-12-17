//
// Created by Fuad Hasan on 10/19/25.
//

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "MeshConfig.h"

void MeshConfig::printParams() {
  printf(" Mesh Generation Parameters:\n");
  printf("  Input Wall File: %s\n", input_wall_file.c_str());
  printf("  Output Mesh File: %s\n", output_name.c_str());
  printf("  Gradation Rate: %.4f\n", gradation_rate);
  if (size_type == SizeType::Uniform) {
    printf("  Size Field Type: Uniform\n");
    printf("  Uniform Size: %.4f\n", uniform_size);
  } else if (size_type == SizeType::Gaussian) {
    printf("  Size Field Type: Gaussian\n");
    printf("  Gaussian Max Size: %.4f\n", gaussian_params.maxSize);
    printf("  Gaussian Min Size: %.4f\n", gaussian_params.minSize);
    printf("  Gaussian Sigma: %.4f\n", gaussian_params.sigma);
    if (!gaussian_params.userdefinedCenter) {
      printf("  Gaussian Center: calculate automatically\n");
    } else {
      printf("  Gaussian Center: (%.4f, %.4f, %.4f)\n",
             gaussian_params.center[0], gaussian_params.center[1],
             gaussian_params.center[2]);
    }
  }

  if (!userdefined_verts.empty()) {
    printf("  User Defined Vertices (%zu):\n", userdefined_verts.size());
    for (size_t i = 0; i < userdefined_verts.size(); ++i) {
      printf("   [%zu]: (%.4f, %.4f, %.4f)\n", i, userdefined_verts[i][0],
             userdefined_verts[i][1], userdefined_verts[i][2]);
    }
  } else {
    printf("  No User Defined Vertices.\n");
  }

  printf("\n");
}

MeshConfig::MeshConfig(const std::string &filename) {
  std::ifstream ifs(filename);
  if (!ifs)
    throw std::runtime_error("Cannot open config file: " + filename);

  std::unordered_map<std::string, std::string> kv;
  std::string line;
  size_t lineno = 0;
  while (std::getline(ifs, line)) {
    ++lineno;
    trim(line);
    if (line.empty())
      continue;
    if (line[0] == '#' || line[0] == ';')
      continue;

    auto pos = line.find('=');
    if (pos == std::string::npos) {
      throw std::runtime_error("Invalid line (no '=') at " + filename + ":" +
                               std::to_string(lineno));
    }
    std::string key = line.substr(0, pos);
    std::string val = line.substr(pos + 1);
    trim(key);
    trim(val);
    if (key.empty())
      throw std::runtime_error("Empty key at " + filename + ":" +
                               std::to_string(lineno));
    kv[toLower(key)] = val;
  }

  // required fields
  if (kv.count("input_wall") == 0 && kv.count("input_wall_file") == 0)
    throw std::runtime_error("Missing input_wall or input_wall_file in config");
  input_wall_file = getValue(kv, {"input_wall", "input_wall_file"});

  if (kv.count("input_model_file") != 0) {
    input_model_file = getValue(kv, {"input_model_file", "input_model"});
  }

  if (kv.count("output") == 0)
    throw std::runtime_error("Missing output in config");
  output_name = kv.at("output");

  if (kv.count("gradation_rate") == 0) {
    throw std::runtime_error("Missing gradation_rate\n");
  }
  gradation_rate = requireDouble(kv, "gradation_rate");
  if (gradation_rate < 0.0 || gradation_rate > 1.0) {
    throw std::runtime_error(
        "Invalid gradation_rate: %f, needs value from 0.0 to 1.0.\n");
  }

  // size field type
  std::string stype = "uniform";
  if (kv.count("size_field.type"))
    stype = toLower(kv.at("size_field.type"));
  else if (kv.count("size_field"))
    stype = toLower(kv.at("size_field"));
  if (stype == "uniform") {
    size_type = SizeType::Uniform;
    if (kv.count("size_field.uniform")) {
      uniform_size =
          parseDouble(kv.at("size_field.uniform"), "size_field.uniform");
      if (uniform_size <= 0.0)
        throw std::runtime_error("size_field.uniform must be > 0");
    }
  } else if (stype == "gaussian" || stype == "gauss") {
    size_type = SizeType::Gaussian;
    // required gaussian params
    if (!kv.count("size_field.gaussian.center"))
      throw std::runtime_error(
          "Missing size_field.gaussian.center for gaussian size field");
    if (kv.at("size_field.gaussian.center") == "calculate") {
      gaussian_params.userdefinedCenter = false;
    } else {
      gaussian_params.userdefinedCenter = true;
      parseVec3(kv.at("size_field.gaussian.center"), gaussian_params.center,
                "size_field.gaussian.center");
    }

    gaussian_params.maxSize = requireDouble(kv, "size_field.gaussian.max");
    gaussian_params.minSize = requireDouble(kv, "size_field.gaussian.min");
    gaussian_params.sigma = requireDouble(kv, "size_field.gaussian.sigma");

    if (gaussian_params.maxSize <= 0.0)
      throw std::runtime_error("size_field.gaussian.max must be > 0");
    if (gaussian_params.minSize <= 0.0)
      throw std::runtime_error("size_field.gaussian.min must be > 0");
    if (gaussian_params.minSize > gaussian_params.maxSize)
      throw std::runtime_error("size_field.gaussian.min must be <= max");
    if (gaussian_params.sigma <= 0.0)
      throw std::runtime_error("size_field.gaussian.sigma must be > 0");
  } else {
    throw std::runtime_error("Unknown size_field.type: " + stype);
  }

  // read user defined vertices
  // two input variable, userdefined_verts.num and userdefined_verts.coords
  // only read coords if num != 0
  if (kv.count("userdefined_verts.num") > 0) {
    int num_verts =
        static_cast<int>(requireDouble(kv, "userdefined_verts.num"));
    printf("Number of user defined vertices: %d\n", num_verts);
    if (num_verts < 0)
      throw std::runtime_error("userdefined_verts.num must be >= 0");
    if (num_verts > 0) {
      if (!kv.count("userdefined_verts.coords"))
        throw std::runtime_error(
            "Missing userdefined_verts.coords for user defined vertices");
      std::string coords_str = kv.at("userdefined_verts.coords");
      std::istringstream iss(coords_str);
      for (int i = 0; i < num_verts; ++i) {
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
          throw std::runtime_error(
              "Not enough coordinates in userdefined_verts.coords. Needs " +
              std::to_string(num_verts * 3) + " values. But found " +
              std::to_string(i * 3) + ".");
        }
        userdefined_verts.push_back({x, y, z});
      }
    }
  }
}

std::string MeshConfig::toLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

void MeshConfig::trim(std::string &s) {
  auto notspace = [](const int ch) { return !std::isspace(ch); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
  s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
}

std::string
MeshConfig::getValue(const std::unordered_map<std::string, std::string> &kv,
                     const std::vector<std::string> &keys) {
  for (auto &k : keys) {
    auto it = kv.find(k);
    if (it != kv.end())
      return it->second;
  }
  throw std::runtime_error("Missing keys in config");
}

double MeshConfig::parseDouble(const std::string &s, const std::string &name) {
  try {
    size_t idx = 0;
    double v = std::stod(s, &idx);
    if (idx != s.size()) {
      // allow trailing spaces
      std::string rest = s.substr(idx);
      for (char c : rest)
        if (!std::isspace(static_cast<unsigned char>(c)))
          throw std::runtime_error("Extra chars in " + name);
    }
    return v;
  } catch (...) {
    throw std::runtime_error("Invalid double for " + name + ": " + s);
  }
}

double MeshConfig::requireDouble(
    const std::unordered_map<std::string, std::string> &kv,
    const std::string &key) {
  auto it = kv.find(key);
  if (it == kv.end())
    throw std::runtime_error("Missing " + key);
  return parseDouble(it->second, key);
}

void MeshConfig::parseVec3(const std::string &s, double out[3],
                           const std::string &name) {
  // accept comma or space separated
  std::string temp = s;
  for (char &c : temp)
    if (c == ',')
      c = ' ';
  std::istringstream iss(temp);
  double a, b, c_val;
  if (!(iss >> a >> b >> c_val))
    throw std::runtime_error("Invalid vec3 for " + name + ": " + s);
  out[0] = a;
  out[1] = b;
  out[2] = c_val;
}