// *************************** DG2MESHWBOUNDARY.CPP ************************ //
// Creates simmetrix mesh from a model using a wall coordinates file to      //
// keep boundary nodes on the wall defined.                                  //
// ************************************************************************* //

/* ************* TODOS *******************
 * Remove helper functions to a header file
 * Add support for given XGC mesh instead of wall coordinates file
 * Add support for given model
 * ************************************** */

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "MeshSim.h"
#include "MeshTypes.h"
#include "ModelTypes.h"
#include "SimCreateModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimPList.h"
#include "SimUtil.h"

#include "MeshConfig.h"

// Gaussian-based size function: size = minSize + (maxSize - minSize) *
// exp(-dist^2/(2*sigma^2))
double centerMaxSizeExpr(const double gpt[3], void *userdata);

void messageHandler(int type, const char *msg);
std::vector<std::array<double, 2>> readWall(const MeshConfig &config, bool &ccw,
                                            pProgress progress);
std::vector<std::array<double, 2>>
readWallCoordinatesFile(const std::string &filename);
std::vector<std::array<double, 2>>
readWallFromXGCNodeFile(const std::string &filename);
std::vector<std::array<double, 2>>
readWallFromSimMesh(const std::string &model_filename,
                    const std::string &mesh_filename, pProgress progress);
bool endsWith(const std::string &fullString,
              const std::string &ending);
void findCornerXptYpt(const std::vector<std::array<double, 2>> &wallCoords,
                      double corner[3], double xpt[3], double ypt[3]);
bool isCounterClockwise(const std::vector<std::array<double, 2>> &wallCoords);
std::vector<pGVertex>
create_vertices(const std::vector<std::array<double, 2>> &wallCoords,
                int nWallPoints, const pGImporter &importer, bool ccw);
std::string usage_string = "Usage: %s <input_file>\n";
std::string help_string =
    usage_string +
    "This program creates a Simmetrix mesh with corresponding model using a "
    "wall coordinates file to keep boundary nodes on the wall defined.\n"
    "The wall coordinates file should have the following format:\n"
    "First line: number of wall points (n)\n"
    "Next n lines: x y coordinates of each wall point\n"
    "The wall points have to be *sequential along the wall*.\n\n";

void check_mesh_wall(const std::vector<pGEdge> &edges, const pMesh *mesh);

int main(int argc, char *argv[]) {
  bool asking_for_help =
      ((argc == 2 && std::string(argv[1]) == "--help") ||
       (argc == 2 && std::string(argv[1]) == "-h") || (argc != 2));
  if (asking_for_help) {
    printf("%s", help_string.c_str());
    return 0;
  }

  MeshConfig meshConfig(argv[1]);
  meshConfig.printParams();
  std::string output_name = meshConfig.output_name;

  try {
    Sim_logOn("simmeshWboundaryNodes.log");
    SimModel_start(); // Call before Sim_readLicenseFile
    // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
    // pass in the location of a file containing your keys.  For a release
    // product, use Sim_registerKey()
    Sim_readLicenseFile(0);
    // Tessellation of GeomSim geometry requires Meshing to have started
    MS_init();
    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    // ************* Read Wall Coordinates *****************//
    bool ccw = true;
    std::vector<std::array<double, 2>> wallCoords = readWall(meshConfig, ccw, progress);
    const int nWallPoints = wallCoords.size();

    // ************* Import Model *****************//
    pGImporter importer = GImporter_new();
    const auto vertices =
        create_vertices(wallCoords, nWallPoints, importer, ccw);

    // create userdefined_vertices
    std::vector<pGVertex> userdefined_vertices(
        meshConfig.userdefined_verts.size());
    if (!meshConfig.userdefined_verts.empty()) {
      for (size_t i = 0; i < meshConfig.userdefined_verts.size(); ++i) {
        auto &v = meshConfig.userdefined_verts[i];
        userdefined_vertices[i] = GImporter_createVertex(importer, v.data());
      }
    }

    std::vector<pGEdge> edges(nWallPoints);
    for (size_t i = 0; i < nWallPoints; ++i) {
      pGVertex v_start = vertices[i];
      pGVertex v_end =
          vertices[(i + 1) % nWallPoints]; // Wrap around to first vertex
      double p_start[3], p_end[3];
      GV_point(v_start, p_start);
      GV_point(v_end, p_end);
      pCurve linearCurve = SCurve_createLine(p_start, p_end);
      edges[i] =
          GImporter_createEdge(importer, v_start, v_end, linearCurve, 0, 1, 1);
    }

    double lowerleft_corner[3], xpt[3], ypt[3];
    findCornerXptYpt(wallCoords, lowerleft_corner, xpt, ypt);

    pSurface planarSurface = SSurface_createPlane(lowerleft_corner, xpt, ypt);
    std::vector<int> facedirs(nWallPoints, 1);
    int loopIndex[1] = {0};
    pGFace face =
        GImporter_createFace(importer, nWallPoints, edges.data(),
                             facedirs.data(), 1, loopIndex, planarSurface, 1);

    pGModel model = GImporter_complete(importer);
    GImporter_delete(importer);

    // check that the input model is topologically valid before meshing
    pPList modelErrors = PList_new();
    if (!GM_isValid(model, 0, modelErrors)) {
      std::cerr << "Generated model not valid. Try running the program in "
                   "debug mode and check printed parameter and check the log."
                << std::endl;
      std::cerr << "Number of errors returned: " << PList_size(modelErrors)
                << std::endl;
      GM_release(model);
      return 1;
    }
    PList_delete(modelErrors);

    GM_setDisplayTolerance(model, 0.5);
    std::cout << "Number of vertices in model: " << GM_numVertices(model)
              << std::endl;
    std::cout << "Number of edges in model: " << GM_numEdges(model)
              << std::endl;
    std::cout << "Number of faces in model: " << GM_numFaces(model)
              << std::endl;
    std::cout << "Number of regions in model: " << GM_numRegions(model)
              << std::endl;
    GM_write(model, (output_name + ".smd").c_str(), 0, progress);

    // ************* Create Mesh *****************//
    pMesh mesh = M_new(0, model);
    pACase meshCase = MS_newMeshCase(model);
    pMeshNex meshNex = MeshNex_new(mesh);

    // specify mesh nodes and edges on wall
    for (int i = 0; i < nWallPoints; ++i) {
      double pos[3];
      GV_point(vertices[i], pos);
      MS_specifyVertex(mesh, pos, NULL, vertices[i], i);
    }
    // userdefined_verts on the face
    for (int i = 0; i < meshConfig.userdefined_verts.size(); ++i) {
      const auto &v = meshConfig.userdefined_verts[i];
      // had to classify the mesh vertices on geometric vertices
      // classifying on the face was not working (got stuck during meshing)
      MS_specifyVertex(mesh, v.data(), 0, userdefined_vertices[i],
                       nWallPoints + i); // tag after wall points
    }

    // done separately to specify vertices first
    // they will not work together
    for (int i = 0; i < nWallPoints; ++i) {
      const int vertTags[2]{i, (i + 1) % nWallPoints};
      MS_specifyEdge(mesh, vertTags, edges[i], i);
    }

    // center of the mesh
    void *size_exp = nullptr;
    if (meshConfig.size_type == MeshConfig::SizeType::Gaussian) {
      double center[3];
      if (!meshConfig.gaussian_params.userdefinedCenter) {
        meshConfig.gaussian_params.center[0] =
            0.5 * (lowerleft_corner[0] + xpt[0]);
        meshConfig.gaussian_params.center[1] =
            0.5 * (lowerleft_corner[1] + ypt[1]);
        meshConfig.gaussian_params.center[2] = 0.0;
      }

      size_exp = MS_registerSizeExprFunc("centerMaxSize", centerMaxSizeExpr,
                                         &meshConfig.gaussian_params);
      MS_setMeshSize(meshCase, face, 2, 0.0, "centerMaxSize($x,$y,$z)");
    } else if (meshConfig.size_type == MeshConfig::SizeType::Uniform) {
      MS_setMeshSize(meshCase, face, 2, meshConfig.uniform_size, NULL);
    }

    // set global mesh size
    // MS_setMeshSize(meshCase, face, 2, 0.5, NULL);
    MS_setGlobalSizeGradationRate(meshCase, meshConfig.gradation_rate);

    std::vector<double> edge_lengths(nWallPoints);
    for (int i = 0; i < nWallPoints; ++i) {
      // calculate edge size based on edge length
      const auto &edge = edges[i];
      const double edge_len = GE_length(edge);
      edge_lengths[i] = edge_len;

      //MS_setMeshSize(meshCase, edges[i], 1, edge_len * 5.0, NULL);
      //  fixme propagation could be useful but causing the program to stall
      //  MS_setMeshSizePropagation(meshCase, edges[i], 2, 1, 0.3, 2.0);
    }

    double max_edge_length = *std::max_element(edge_lengths.begin(), edge_lengths.end());
    std::cout << "[INFO] Maximum wall edge length: " << max_edge_length
              << std::endl;

    for (int i = 0; i < nWallPoints; ++i) {
      const auto &edge = edges[i];
      MS_setMeshSize(meshCase, edge, 1, max_edge_length, NULL);
    }

    // generate mesh
    pSurfaceMesher surfmesh = SurfaceMesher_new(meshCase, mesh);
    SurfaceMesher_setEnforceSpatialGradation(surfmesh, 1);
    SurfaceMesher_execute(surfmesh, progress);

    MS_unregisterSizeExprFunc(size_exp);

    // validity check
    check_mesh_wall(edges, &mesh);

    M_write(mesh, (output_name + ".sms").c_str(), 0, progress);
    printf("[INFO] Mesh written to file: %s.sms\n", output_name.c_str());
    std::string mesh_stat = "Mesh statistics:\n"
                            "Number of nodes: %d\n"
                            "Number of faces: %d\n";
    printf(mesh_stat.c_str(), M_numVertices(mesh), M_numFaces(mesh));

    MeshNex_setNodeDefault(meshNex, 1, 0);

    MeshNex_write(meshNex, (output_name + ".nex").c_str());
    printf("Nex file written to file: %s.nex\n", output_name.c_str());

    // ************* Clean Up ***************** //
    MeshNex_delete(meshNex);
    SurfaceMesher_delete(surfmesh);
    MS_deleteMeshCase(meshCase);
    M_release(mesh);
    GM_release(model);
    Progress_delete(progress);
    MS_exit();
    Sim_unregisterAllKeys();
    SimModel_stop();
    Sim_logOff();
  } catch (pSimInfo err) {
    std::cerr << "SimModSuite error caught:" << std::endl;
    std::cerr << "  Error code: " << SimInfo_code(err) << std::endl;
    std::cerr << "  Error string: " << SimInfo_toString(err) << std::endl;
    SimInfo_delete(err);
    return 1;
  } catch (...) {
    std::cerr << "Unhandled exception caught" << std::endl;
    return 1;
  }

  return 0;
}

void messageHandler(int type, const char *msg) {
  switch (type) {
  case Sim_InfoMsg:
    std::cout << "Info: " << msg << std::endl;
    break;
  case Sim_DebugMsg:
    std::cout << "Debug: " << msg << std::endl;
    break;
  case Sim_WarningMsg:
    std::cout << "Warning: " << msg << std::endl;
    break;
  case Sim_ErrorMsg:
    std::cout << "Error: " << msg << std::endl;
    break;
  }
}

std::vector<std::array<double, 2>> readWallFromSimMesh(const std::string &model_filename, const std::string &mesh_filename, pProgress progress) {
  const auto model = GM_load(model_filename.c_str(), NULL, progress);
  const auto mesh = M_load(mesh_filename.c_str(), model, progress);

  // find the leftmost mesh vertex to start searching for wall vertices
  VIter v_iter = M_vertexIter(mesh);
  pVertex leftmost_vertex = nullptr;
  pVertex current_vertex = nullptr;

  while (current_vertex = VIter_next(v_iter)) {
    if (!leftmost_vertex) {
      leftmost_vertex = current_vertex;
    } else {
      double current_coords[3];
      double leftmost_coords[3];
      V_coord(current_vertex, current_coords);
      V_coord(leftmost_vertex, leftmost_coords);
      if (current_coords[0] < leftmost_coords[0]) {
        leftmost_vertex = current_vertex;
      }
    }
  }
  VIter_delete(v_iter);

  // traverse the boundary edges to collect wall vertices
  std::vector<std::array<double, 2>> wall_points;
  current_vertex = leftmost_vertex;
  pEdge last_boundary_edge = nullptr;
  while (true) {
    pEdge current_boundary_edge = nullptr;
    const int n_adj_edges = V_numEdges(current_vertex);
    for ( int i = 0; i < n_adj_edges; ++i) {
      current_boundary_edge = V_edge(current_vertex, i);
      if (E_numFaces(current_boundary_edge) == 1 && current_boundary_edge != last_boundary_edge) {
        last_boundary_edge = current_boundary_edge;
        break;
      }
    }
    if (current_boundary_edge == nullptr) {
      double current_coords[3];
      V_coord(current_vertex, current_coords);
      printf("Error: For vertex at (%.6f, %.6f, %.6f), no boundary edge found but it has %d adjacent edges.\n",
             current_coords[0], current_coords[1], current_coords[2], n_adj_edges);
      std::abort();
    }

    // find the next vertex on the boundary edge
    current_vertex = E_otherVertex(current_boundary_edge, current_vertex);

    // store the vertex coordinates
    double coords[3];
    V_coord(current_vertex, coords);
    wall_points.push_back({coords[0], coords[1]});

    // complete the loop
    if (current_vertex == leftmost_vertex) {
      break; // completed the loop
    }
  }

  return wall_points;
}


std::vector<std::array<double, 2>> readWallFromXGCNodeFile(
    const std::string &filename) {
  printf("[Warning] Reading wall coordinates from XGC %s file. As of"
         "this writing, .node file does not contain wall nodes in order.", filename.c_str());
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  size_t boundary_point_count = 0;
  size_t total_point_count = 0;
  size_t dim = 0;
  size_t n_attributes = 0;
  size_t n_boundary_attrs = 0;

  // Read Header
  infile >> total_point_count >> dim >> n_attributes >> n_boundary_attrs;
  if ( total_point_count < 3 || dim != 2 || n_attributes != 0 ||
      n_boundary_attrs != 1) {
    std::string error_msg = "In file " + filename +
                            ", invalid header values: ";
    error_msg += "total_point_count=" + std::to_string(total_point_count) +
                 ", dim=" + std::to_string(dim) +
                 ", n_attributes=" + std::to_string(n_attributes) +
                 ", n_boundary_attrs=" + std::to_string(n_boundary_attrs) + ".";
    throw std::runtime_error(error_msg);
  }

  std::vector<std::array<double, 2>> boundary_points;
  // will strip the extra later
  boundary_points.reserve(total_point_count);

  double x, y;
  bool is_boundary;
  size_t node_id;

  for (size_t i = 0; i < total_point_count; ++i) {
    if (!(infile >> node_id >> x >> y >> is_boundary)) {
      throw std::runtime_error("Error reading point at line " +
                               std::to_string(i + 2));
    }

    if (is_boundary) {
      boundary_points.push_back({x, y});
      boundary_point_count++;
    }
  }

  // Resize to actual boundary point count
  boundary_points.resize(boundary_point_count);
  return boundary_points;
}


std::vector<std::array<double, 2>>
readWallCoordinatesFile(const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  size_t n;
  infile >> n;
  std::vector<std::array<double, 2>> points;
  points.reserve(n);

  double x, y;
  for (size_t i = 0; i < n; ++i) {
    if (!(infile >> x >> y)) {
      throw std::runtime_error("Error reading point at line " +
                               std::to_string(i + 2));
    }
    points.push_back({x, y});
  }
  return points;
}

void findCornerXptYpt(const std::vector<std::array<double, 2>> &wallCoords,
                      double corner[3], double xpt[3], double ypt[3]) {
  // corner is the lower left corner
  // xpt is the highest of x values with corner y
  // ypt is the highest of y values with corner x
  double min_x = wallCoords[0][0];
  double min_y = wallCoords[0][1];
  double max_x = wallCoords[0][0];
  double max_y = wallCoords[0][1];

  for (const auto &point : wallCoords) {
    if (point[0] < min_x) {
      min_x = point[0];
    }
    if (point[1] < min_y) {
      min_y = point[1];
    }
    if (point[0] > max_x) {
      max_x = point[0];
    }
    if (point[1] > max_y) {
      max_y = point[1];
    }
  }
  // tol is 5% of the range
  double tol_x = 0.05 * (max_x - min_x);
  double tol_y = 0.05 * (max_y - min_y);

  corner[0] = min_x - tol_x;
  corner[1] = min_y - tol_y;
  corner[2] = 0.0;

  xpt[0] = max_x + tol_x;
  xpt[1] = min_y - tol_y;
  xpt[2] = 0.0;

  ypt[0] = min_x - tol_x;
  ypt[1] = max_y + tol_y;
  ypt[2] = 0.0;

#ifndef NDEBUG
  printf("Corner: (%.6f, %.6f, %.6f)\n", corner[0], corner[1], corner[2]);
  printf("Xpt: (%.6f, %.6f, %.6f)\n", xpt[0], xpt[1], xpt[2]);
  printf("Ypt: (%.6f, %.6f, %.6f)\n", ypt[0], ypt[1], ypt[2]);
#endif
}

bool isCounterClockwise(const std::vector<std::array<double, 2>> &wallCoords) {
  // not robust, just checks the first three points
  const auto n_points = wallCoords.size();
  const double first_loc = double(1. / 3.) * double(n_points);
  const double second_loc = double(2. / 3.) * double(n_points);

  double p1[2] = {wallCoords[0][0], wallCoords[0][1]};
  double p2[2] = {wallCoords[first_loc][0], wallCoords[first_loc][1]};
  double p3[2] = {wallCoords[second_loc][0], wallCoords[second_loc][1]};

  double cross_prod_mod =
      (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);

  // if coss_prod is about 0, points are collinear
  // then move to the next set of points
  double tol = 1e-8 * (p1[0] + p1[1] + p2[0] + p2[1] + p3[0] +
                       p3[1]); // scale tol with size of vectors
  if (cross_prod_mod > tol) {
    return true;
  } else if (cross_prod_mod < -tol) {
    return false;
  } else {
    // todo : implement more robust check
    throw std::runtime_error(
        "This check is not robust enough for collinear points.");
  }
}

std::vector<pGVertex>
create_vertices(const std::vector<std::array<double, 2>> &wallCoords,
                const int nWallPoints, const pGImporter &importer,
                const bool ccw) {
  std::vector<pGVertex> vertices(nWallPoints);
  for (size_t i = 0; i < nWallPoints; ++i) {
    double coords[3] = {wallCoords[i][0], wallCoords[i][1], 0.0};
    const size_t index = ccw ? i : nWallPoints - (i + 1);
    vertices[index] = GImporter_createVertex(importer, coords);
  }
  return vertices;
}

double centerMaxSizeExpr(const double gpt[3], void *userdata) {
  auto *p = static_cast<const GaussSizeParams *>(userdata);
  if (!p || p->sigma <= 0.0 || p->maxSize <= 0.0 || p->minSize <= 0.0) {
    std::string error_msg = "Invalid GaussSizeParams: ";
    if (!p)
      error_msg += "nullptr userdata. ";
    if (p && p->sigma <= 0.0)
      error_msg += "sigma <= 0.0. ";
    if (p && p->maxSize <= 0.0)
      error_msg += "maxSize <= 0.0. ";
    if (p && p->minSize <= 0.0)
      error_msg += "minSize <= 0.0. ";
    throw std::runtime_error(error_msg);
  }
  double dx = gpt[0] - p->center[0];
  double dy = gpt[1] - p->center[1];
  double dz = gpt[2] - p->center[2];
  double dist2 = dx * dx + dy * dy + dz * dz;
  double exponent = -dist2 / (2.0 * p->sigma * p->sigma);
  double value = p->minSize + (p->maxSize - p->minSize) * std::exp(exponent);
  // ensure positive and at least minSize
  if (value < p->minSize)
    value = p->minSize;
  if (value <= 0.0)
    value = 1e-6;
  return value;
}

void check_mesh_wall(const std::vector<pGEdge> &edges, const pMesh *mesh) {
  EIter edge_iter = M_edgeIter(*mesh);
  int boundary_edge_count = 0;

  while (const pEdge &edge = EIter_next(edge_iter)) {
    int n_adj_faces = E_numFaces(edge);
    if (n_adj_faces == 1) {
      boundary_edge_count++;
    }
  }
  EIter_delete(edge_iter);

  if (boundary_edge_count != edges.size()) {
    std::string error_message = "Mesh boundary edge count (" +
                                std::to_string(boundary_edge_count) +
                                ") does not match expected wall edge count (" +
                                std::to_string(edges.size()) + ").";
    throw std::runtime_error(error_message);
  }
}
// Function to check if a string ends with another string
bool endsWith(const std::string& fullString,
              const std::string& ending)
{
  // Check if the ending string is longer than the full
  // string
  if (ending.size() > fullString.size())
    return false;

  // Compare the ending of the full string with the target
  // ending
  return fullString.compare(fullString.size()
                                - ending.size(),
                            ending.size(), ending)
         == 0;
}

std::vector<std::array<double, 2>> readWall(const MeshConfig& config, bool &ccw, pProgress progress) {
  std::string wallFile = config.input_wall_file;

  std::vector<std::array<double, 2>> wallCoords;
  const bool is_wall_file_xgc_node = endsWith(wallFile, ".node");
  const bool is_wall_file_sim_mesh = endsWith(wallFile, ".sms");
  if (is_wall_file_xgc_node) {
    wallCoords = readWallFromXGCNodeFile(wallFile);
  } else if (is_wall_file_sim_mesh) {
    std::string model_file = config.input_model_file;
    if (model_file.empty()) {
      throw std::runtime_error("When wall file is a Simmetrix mesh (.sms), "
                               "the input model file must be provided in the "
                               "configuration.");
    }
    wallCoords = readWallFromSimMesh(model_file, wallFile, progress);
  } else {
    wallCoords = readWallCoordinatesFile(wallFile);
  }

  ccw = isCounterClockwise(wallCoords);
  int nWallPoints = wallCoords.size();
  std::cout << "[INFO] Wall coordinates read from file: " << wallFile
      << " with " << nWallPoints << " points." << std::endl;
  std::cout << "[INFO] Wall points are in "
      << (ccw ? "*counter-clockwise*" : "*clockwise*") << " order."
      << std::endl;

#ifndef NDEBUG
  printf("Wall Points:\n");
  for (size_t i = 0; i < nWallPoints; ++i) {
    double coords[3] = {wallCoords[i][0], wallCoords[i][1], 0.0};
    printf("[INFO] Wall Point %zu: (%.6f, %.6f, %.6f)\n", i + 1, coords[0],
           coords[1], coords[2]);
  }
#endif

  return wallCoords;
}

