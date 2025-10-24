//
// Input file for mesh configuration parameters.
// Created by Fuad Hasan on 10/19/25.
//

#ifndef DG2MESHWBOUNDARY_MESHCONFIG_H
#define DG2MESHWBOUNDARY_MESHCONFIG_H

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

struct GaussSizeParams {
  double center[3]{0.0, 0.0, 0.0};
  bool userdefinedCenter = false;
  double maxSize = 0.0;
  double minSize = 0.0;
  double sigma = 0.0;
};

class MeshConfig {
public:
  enum class SizeType { Uniform, Gaussian };

  // Parsed fields
  std::string input_wall_file;
  std::string output_name;

  SizeType size_type = SizeType::Uniform;
  double uniform_size = 0.05;
  double gradation_rate = 0.2;

  std::vector<std::array<double, 3>> userdefined_verts;

  GaussSizeParams gaussian_params;

  // Load and parse config file, throws std::runtime_error on error
  explicit MeshConfig(const std::string &filename);

  void printParams();

private:
  static std::string toLower(std::string s);
  static void trim(std::string &s);
  static std::string
  getValue(const std::unordered_map<std::string, std::string> &kv,
           const std::vector<std::string> &keys);
  static double parseDouble(const std::string &s, const std::string &name);
  static double
  requireDouble(const std::unordered_map<std::string, std::string> &kv,
                const std::string &key);
  static void parseVec3(const std::string &s, double out[3],
                        const std::string &name);
};

#endif // DG2MESHWBOUNDARY_MESHCONFIG_H