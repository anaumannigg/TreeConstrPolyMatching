
#ifndef TCPOLYGONMATCHING_BINARY_IO_H
#define TCPOLYGONMATCHING_BINARY_IO_H

#include "cgal_includes.h"

void writePolysToBinaryFile(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2 , const std::string& file_path, size_t set_id);

std::pair<std::vector<Polygon_wh>,std::vector<Polygon_wh>> readPolysFromBinaryFile(const std::string& file_path, size_t set_id);


//CLASS PolygonFileReader for Line by Line reading to maybe improve performance (to be tested)

class BinaryPolygonFileReader {
public:
    explicit BinaryPolygonFileReader(const std::string& file_path) {
        ifs.open(file_path, std::ios::binary);
        if (!ifs) {
            throw std::runtime_error("Failed to open file for reading");
        }
    }

    ~BinaryPolygonFileReader() {
        if (ifs.is_open()) {
            ifs.close();
        }
    }

    bool readNextSet(size_t& set_id, std::vector<Polygon_wh>& first_polygon_set, std::vector<Polygon_wh>& second_polygon_set);

private:
    std::ifstream ifs;
};

//END CLASS

#endif //TCPOLYGONMATCHING_BINARY_IO_H
