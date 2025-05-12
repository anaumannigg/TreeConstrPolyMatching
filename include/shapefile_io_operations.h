#ifndef _shapefile_io_operations_included_
#define _shapefile_io_operations_included_
#include "cgal_includes.h"
#include "polygon_wh.h"
#include "shapefil.h"
#include "solution.h"
#include "polygon_operations.h"


#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

//reads shapefile of single polygons without holes into a vector
bool ReadShapeFile(const std::string& filename, std::vector<Polygon_wh>& polys, bool ignore_invalid_geometries = false);

//writes polygons into shapefile
void writeToShapeFile(std::vector<Polygon_wh> polys, std::string path);

//writes polygons into shapefile labeled with their group index
void writeToShapeFile(std::vector<Polygon_wh> polys, std::vector<int> group_index, std::string path);

//writes polygons into shapefile labeled with their group index and match weight
void writeToShapeFile(std::vector<Polygon_wh> polys, std::vector<int> group_index, std::vector<double> match_weight, std::string path);

//writes polygons into wkt-file (.txt)
void writeToWKT(std::vector<Polygon_wh> polys, std::string path);

//writes Segmets into wkt-file (.txt)
void writeToWKT(std::vector<Segment> segments, std::string path);

//writes solution data into csv
void writeToCSV(Solution s, std::string path);

//write singular objectives of each connected component to csv
void writeToCSV(std::vector<double>& obj_of_set, std::string path);

//writes analysis data into csv
void writeToCSV(std::vector<int> set_sizes, std::vector<std::pair<int, int>> set_sizes_after_precomp, std::vector<double> execution_times, std::vector<bool> completed_exploration, std::string path);

//helper function to find the matches, in which two solutions differ
void findDifferingMatches(std::string fileA, std::string fileB);

//writes arrangement to shapefile (edges as well as faces in two separate files)
void writeToShapeFile(Arrangement arr, std::string path, std::vector<int> poly_face_ids = std::vector<int>());




#endif