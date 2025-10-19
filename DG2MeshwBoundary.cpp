// *************************** DG2MESHWBOUNDARY.CPP ************************ //
// Creates simmetrix mesh from a model using a wall coordinates file to      //
// keep boundary nodes on the wall defined.                                  //
// ************************************************************************* //

#include <array>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "MeshSim.h"
#include "MeshTypes.h"
#include "ModelTypes.h"
#include "SimCreateModel.h"
#include "SimDisplay.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimPList.h"
#include "SimUtil.h"

// user data passed to the size expression
struct SizeFieldParams {
  double center[3]; // cx, cy, cz
  double maxSize;   // size at center (must be > 0)
  double minSize;   // minimum allowed size (must be > 0 and <= maxSize)
  double sigma;     // gaussian std-dev (controls falloff), > 0
};

// Gaussian-based size function: size = minSize + (maxSize - minSize) *
// exp(-dist^2/(2*sigma^2))
double centerMaxSizeExpr(const double gpt[3], void *userdata);

void messageHandler(int type, const char *msg);
std::vector<std::array<double, 2>>
readWallCoordinatesFile(const std::string &filename);
void findCornerXptYpt(const std::vector<std::array<double, 2>> &wallCoords,
                      double corner[3], double xpt[3], double ypt[3]);
bool isCounterClockwise(const std::vector<std::array<double, 2>> &wallCoords);
std::vector<pGVertex>
create_vertices(const std::vector<std::array<double, 2>> &wallCoords,
                int nWallPoints, const pGImporter &importer, bool ccw);
std::string usage_string = "Usage: %s <wallfile_name> <output_name>\n";
std::string help_string =
    "This program creates a Simmetrix mesh with corresponding model using a "
    "wall coordinates file to keep boundary nodes on the wall defined.\n"
    "The wall coordinates file should have the following format:\n"
    "First line: number of wall points (n)\n"
    "Next n lines: x y coordinates of each wall point\n"
    "The wall points have to be *sequential along the wall*.\n\n" +
    usage_string;

int main(int argc, char *argv[]) {
  bool asking_for_help = (argc == 2 && std::string(argv[1]) == "--help");
  if (asking_for_help) {
    printf(help_string.c_str(), argv[0]);
    return 0;
  }
  if (argc != 3) {
    printf(usage_string.c_str(), argv[0]);
    return 1;
  }

  std::string wallFile = argv[1];
  std::string output_name = argv[2];

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

    const auto wallCoords = readWallCoordinatesFile(wallFile);
    bool ccw = isCounterClockwise(wallCoords);
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

    // ************* Import Model *****************//
    pGImporter importer = GImporter_new();
    const auto vertices =
        create_vertices(wallCoords, nWallPoints, importer, ccw);

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

    // specify mesh nodes and edges on wall
    for (int i = 0; i < nWallPoints; ++i) {
      double pos[3];
      GV_point(vertices[i], pos);
      MS_specifyVertex(mesh, pos, NULL, vertices[i], i);
    }
    // done separately to specify vertices first
    // they will not work together
    for (int i = 0; i < nWallPoints; ++i) {
      const int vertTags[2]{i, (i + 1) % nWallPoints};
      MS_specifyEdge(mesh, vertTags, edges[i], i);
    }
    // center of the mesh
    double center[3] = {0.5 * (lowerleft_corner[0] + xpt[0]),
                        0.5 * (lowerleft_corner[1] + ypt[1]), 0.0};
    SizeFieldParams params = {
        {center[0], center[1], center[2]}, // center
        0.1,                               // maxSize at center
        0.005,                             // minSize
        0.05                               // sigma
    };
    void *size_exp =
        MS_registerSizeExprFunc("centerMaxSize", centerMaxSizeExpr, &params);
    MS_setMeshSize(meshCase, face, 2, 0.0, "centerMaxSize($x,$y,$z)");
    // set global mesh size
    pModelItem modelDomain = GM_domain(model);
    // MS_setMeshSize(meshCase, face, 2, 0.5, NULL);
    // MS_setGlobalSizeGradationRate(meshCase, 0.2);

    for (int i = 0; i < nWallPoints; ++i) {
      // MS_setMeshSize(meshCase, edges[i], 2, 1.0, NULL);
      //  fixme propagation could be useful but causing the program to stall
      //  MS_setMeshSizePropagation(meshCase, edges[i], 2, 1, 0.3, 2.0);
    }

    // generate mesh
    pSurfaceMesher surfmesh = SurfaceMesher_new(meshCase, mesh);
    SurfaceMesher_setEnforceSpatialGradation(surfmesh, 1);
    SurfaceMesher_execute(surfmesh, progress);

    // write mesh to file
    M_write(mesh, (output_name + ".sms").c_str(), 0, progress);
    printf("[INFO] Mesh written to file: %s.sms\n", output_name.c_str());

    // ************* Clean Up *****************//
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
  auto *p = static_cast<const SizeFieldParams *>(userdata);
  if (!p || p->sigma <= 0.0 || p->maxSize <= 0.0 || p->minSize <= 0.0) {
    std::string error_msg = "Invalid SizeFieldParams: ";
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